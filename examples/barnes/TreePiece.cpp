#include <sstream>
using namespace std;

#include "TreePiece.h"
#include "SphShadow.h"
#include "SplitterGroup.h"
#include "ParamStorage.h"
#include "TipsyFile.h"
#include "Utilities.h"
#include "NodePayload.h"
#include "Physics.h"
#include "Visitor_def.h"
#include "Resumption.h"

#include "Orb3dLB.h"
#include "TaggedVector3D.h"

// readonly's
extern ParamStorage parameters;
extern CProxy_SplitterGroup splitterGroupProxy;
extern CProxy_TreePiece treePieceProxy;
extern CProxy_SphShadow sphShadowProxy;


extern MyTreeHandleType treeHandle;
extern GravityVisitor::TraversalHandleType gravityTraversalHandle;
extern DensitySphVisitor::TraversalHandleType sphDensityTraversalHandle;
extern BallSphVisitor::TraversalHandleType sphBallTraversalHandle;

TreePiece::TreePiece(CkMigrateMessage *m){
}

TreePiece::TreePiece() : 
  myNumParticles_(0),
  nExpectedParticles_(-1),
  iteration_(0),
#ifdef DEBUG_TRAVERSAL
  checker_(thisIndex, Key(0))
#endif
  gravityTraversalData_(), 
  sphTraversalData_(), 
  gravityOutstanding_(-1),
  sphOutstanding_(-1),
  ballSphOutstanding_(-1),
  nLeavesInitiated_(-1),
  nSphLeavesInitiated_(-1),
  nBallSphLeavesInitiated_(-1)
{
  usesAtSync = true;
  findOrb3dLb();
}

CProxy_TreePiece TreePiece::myProxy() const {
  return CProxy_TreePiece(thisProxy);
}

void TreePiece::load(const CkCallback &cb){
  loadTipsy();
  contribute(sizeof(BoundingBox), &myBox_, MultiphaseDistree::Reduce<BoundingBox>::reducer(), cb);
}

void TreePiece::loadTipsy(){
  string filename(parameters.inFile);
  Tipsy::TipsyReader r(filename);
  if(!r.status()) {
    CkPrintf("[%u] Fatal: Couldn't open tipsy file! (%s)\n", thisIndex, filename.c_str());
    CkExit();
    return;
  }

  Tipsy::header tipsyHeader = r.getHeader();
  int nTotalParticles = tipsyHeader.nbodies;
  int nTotalSPH = tipsyHeader.nsph;
  int nTotalDark = tipsyHeader.ndark;
  int nTotalStar = tipsyHeader.nstar;
  int dStartTime = tipsyHeader.time;
  int excess;
  unsigned int startParticle;

  unsigned int numLoaders = parameters.nPieces;

  int np = nTotalParticles / numLoaders;
  excess = nTotalParticles % numLoaders;
  startParticle = np * thisIndex;
  if(thisIndex < (unsigned int) excess) {
    np++;
    startParticle += thisIndex;
  }
  else {
    startParticle += excess;
  }

  // allocate an array for myParticles_
  int nStore = np;
  myParticles_.resize(nStore);

  if(!r.seekParticleNum(startParticle)) {
    CkAbort("Couldn't seek to my particles!");
    return;
  }

  Tipsy::gas_particle gp;
  Tipsy::dark_particle dp;
  Tipsy::star_particle sp;

  int iSPH = 0;
  int iStar = 0;
  myBox_.reset();
  myBox_.pe = 0.0;
  myBox_.ke = 0.0;

  for(unsigned int i = 0; i < myParticles_.size(); ++i) {
    if(i + startParticle < (unsigned int) tipsyHeader.nsph) {
      if(!r.getNextGasParticle(gp)) {
        CkAbort("failed to read gas particle!");
      }
      myParticles_[i].mass() = gp.mass;
      myParticles_[i].position() = gp.pos;
      myParticles_[i].velocity() = gp.vel;

      iSPH++;
    } else if(i + startParticle < (unsigned int) tipsyHeader.nsph
        + tipsyHeader.ndark) {
      if(!r.getNextDarkParticle(dp)) {
        CkAbort("failed to read dark particle!");
      }
      myParticles_[i].mass() = dp.mass;
      myParticles_[i].position() = dp.pos;
      myParticles_[i].velocity() = dp.vel;
    } else {
      if(!r.getNextStarParticle(sp)) {
        CkAbort("failed to read star particle!");
      }
      myParticles_[i].mass() = sp.mass;
      myParticles_[i].position() = sp.pos;
      myParticles_[i].velocity() = sp.vel;
      iStar++;
    }
    myParticles_[i].order = i+startParticle; 
    myParticles_[i].potential = 0.0;
    myBox_.grow(myParticles_[i].position());

    myBox_.mass += myParticles_[i].mass();
    myBox_.ke += myParticles_[i].mass()*myParticles_[i].velocity().lengthSquared();
    myBox_.pe = 0.0;
  }
  myBox_.ke /= 2.0;
  myBox_.numParticles = myParticles_.size();

}

void TreePiece::assignKeys(BoundingBox &universe, const CkCallback &cb){
  Key prepend;
  prepend = 1L;
  prepend <<= (TREE_KEY_BITS-1);

  Vector3D<Real> sz = universe.box.greater_corner - universe.box.lesser_corner;
  Vector3D<Real> prel;

  for(unsigned int i = 0; i < myParticles_.size(); i++){
    Particle *p = &(myParticles_[i]);
    // Obtain the integer grid points on which the particle falls in each dimension

    prel = p->position()-universe.box.lesser_corner;
    prel /= sz;
    prel *= (1.0*BOXES_PER_DIM);

    Key xint = (Key) prel.x;
    Key yint = (Key) prel.y;
    Key zint = (Key) prel.z;

    // Interleave bits
    Key mask = Key(0x1);
    Key k = Key(0x0);
    int shiftBy = 0;
    for(int j = 0; j < BITS_PER_DIM; j++){
      k |= ((zint & mask) <<  shiftBy);
      k |= ((yint & mask) << (shiftBy+1));
      k |= ((xint & mask) << (shiftBy+2));
      mask <<= 1;
      // minus 1 because mask itself has shifted
      // left by one position
      shiftBy += (NDIMS-1);
    }
    // Prepend the key with a '1' bit.
    k |= prepend;
    myParticles_[i].key = k;
  }

  myParticles_.quickSort();

  contribute(cb);
}

void TreePiece::count(CkVec<Key> &keys, const CkCallback &cb){
  // search for first particle with key GEQ check
  // in the range [start, finish)
  int start = 0;
  int finish = myParticles_.size();

  convertKeysToSplitters(keys);

  CkVec<int> counts(keys.size() / 2);
  CkAssert(counts.size() * 2 == keys.size());

  ostringstream oss;
  if(myParticles_.size() > 0){
    for(int i = 0; i < counts.size(); i++){
      Key from = keys[2*i];
      Key to = keys[2*i+1];

      int begin = Utilities::binary_search_ge(from, &myParticles_[0], start, finish);
      int end = Utilities::binary_search_ge(to, &myParticles_[0], begin, finish);
      counts[i] = end - begin;
      
      start = end;
    }
  }
  else{
    for(int i = 0; i < counts.size(); i++){
      counts[i] = 0;
    }
  }

  contribute(sizeof(int) * counts.size(), &counts[0], CkReduction::sum_int, cb);
}

void TreePiece::convertKeysToSplitters(CkVec<Key> &keys){
  for(int i = 0; i < keys.size(); i++){
    // don't try to convert largest possible key to a splitter
    // because it already is a splitter.
    if(keys[i] == (~Key(0))) continue;

    Key &k = keys[i];
    k = Utilities::toSplitter(k);
  }
}

void TreePiece::flush(const CkCallback &cb){
  int start = 0;
  int finish = myParticles_.size();

  const CkVec<Splitter> &splitters = splitterGroupProxy.ckLocalBranch()->getSplitters();

  nPieces_ = splitters.size();

  if(thisIndex < nPieces_){
    // this is how many particles i'm expecting to receive
    nExpectedParticles_ = splitters[thisIndex].nParticles;
    // XXX instead of sending tpRootKey, could send the depth of the 'from' key in splitter
    tpRootKey_ = splitters[thisIndex].tpRootKey;
  }
  else{
    nExpectedParticles_ = 0;
    tpRootKey_ = Key(0);
  }

  int nFlushedParticles = 0;

  // flush particles held by self to rightful owners
  if(myParticles_.size() > 0){
    for(int i = 0; i < nPieces_; i++){
      Key from = splitters[i].from;
      Key to = splitters[i].to;

      int begin = Utilities::binary_search_ge(from, &myParticles_[0], start, finish);
      int end = Utilities::binary_search_ge(to, &myParticles_[0], begin, finish);

      int np = end - begin;

      if(np > 0){ 
        ParticleMsg *msg = new (np) ParticleMsg(&myParticles_[begin], np);
        treePieceProxy[i].particles(msg);
      }

      nFlushedParticles += np;

      start = end;
    }
  }
  
  if(nFlushedParticles != myParticles_.size()){
    CkPrintf("[%d] TreePiece::flush held %d flushed %d\n", thisIndex, myParticles_.size(), nFlushedParticles);

    ostringstream oss;
    for(int i = 0; i < splitters.size(); i++){
      oss << "piece " << i << " " << splitters[i].nParticles << " keys [" << splitters[i].from << ", " << splitters[i].to << ")" << endl;
    }
    CkPrintf("[%d] splitters:\n%s\n", thisIndex, oss.str().c_str());


    CkAbort("TreePiece didn't flush all particles\n");
  }

  callback_ = cb;

  if(myNumParticles_ == nExpectedParticles_){
    copyParticles();
  }

}

void TreePiece::particles(ParticleMsg *msg){
  particleMsgs_.push_back(ParticleMsgWrapper(msg));
  myNumParticles_ += msg->nParticles;

  if(myNumParticles_ == nExpectedParticles_){
    copyParticles();
  }
}

void TreePiece::copyParticles(){
  myParticles_.resize(myNumParticles_);

  int target = 0;
  for(int i = 0; i < particleMsgs_.size(); i++){
    ParticleMsg *msg = particleMsgs_[i].msg;
    memcpy(&myParticles_[target], msg->particles, msg->nParticles * sizeof(Particle));
    target += msg->nParticles;

    delete msg;
  }
  CkAssert(target == myNumParticles_);
  myParticles_.quickSort();

  // reset
  myNumParticles_ = 0;
  nExpectedParticles_ = -1;
  particleMsgs_.resize(0);
  myBox_.reset();
  myLeaves_.resize(0);

  contribute(callback_);
}

void TreePiece::printParticlesToFile(){
  ostringstream oss;
  oss << "plot/particles." << thisIndex << ".pos";
  ofstream particleOutFile(oss.str().c_str());
  CkAssert(particleOutFile.is_open());

  for(int i = 0; i < myParticles_.size(); i++){
    const Particle &p = myParticles_[i];
    particleOutFile << (*p.getSphCore()) << endl;
  }

  particleOutFile.close();
}

void TreePiece::build(const CkCallback &cb){
  callback_ = cb;

  // try to build a tree with your particles only
  // if you have been given splitters by the 
  // decomposer. Note that a TreePiece may have 
  // 0 particles even though it has been assigned
  // a splitter range. However, if it hasn't been 
  // assigned a splitter range, it cannot have any
  // particles assigned to it, and therefore there
  // is no need to attempt a local tree build
  if(thisIndex < nPieces_){
    Particle *particles = NULL;
    if(myParticles_.size() > 0){
      particles = &myParticles_[0];
    }

    // arguments to constructor of MyLocalNodeType are:
    // key, depth, particle start index, num particles 
    root_ = new MyLocalNodeType(Key(1), 0, particles, myParticles_.size());
    // global root is shared by all (useful) tree pieces
    root_->setOwnerStart(0);
    root_->setOwnerEnd(nPieces_);

    root_->setParent(NULL);

    // save this for orb 3d load balancing
    myLbCentroid_ = root_->getPayload().moments().com;

    const CkVec<Splitter> &splitters = splitterGroupProxy.ckLocalBranch()->getSplitters();

    // do only local build; don't ask for remote nodes'
    // moments, but ensure they are correctly labeled;
    // submit tree to manager when done

    build(splitters, root_, false);

    //print(root_);
    treeHandle.registration(root_, this);
  }
  else{
    root_ = NULL;
    CkAssert(myParticles_.size() == 0);
  }

  contribute(cb);
}

void TreePiece::print(MyLocalNodeType *root){
  ostringstream oss;
  oss << "tree." << thisIndex << ".dot";
  ofstream out(oss.str().c_str());
  CkAssert(out.is_open());
  out << "digraph tree" << thisIndex << "{" << endl;
  root->dot(out);
  out << "}" << endl;
  out.close();
}

// returns whether node and everything below it is entirely local
bool TreePiece::build(const CkVec<Splitter> &splitters, MyLocalNodeType *node, bool rootLiesOnPath){
  node->tp() = thisIndex;
  // if we haven't seen the root of the tp yet, check whether
  // current node is the root; otherwise, don't perform this check:
  // once the tp root has been seen along a path, it stays that way
  if(!rootLiesOnPath){
    rootLiesOnPath = (node->getKey() == tpRootKey_);
  }

  MultipoleMoments &moments = node->getPayload().moments();

  int nodeNumParticles = node->getNumParticles();
  bool isLight = (nodeNumParticles <= BUCKET_TOLERANCE * Real(parameters.maxppb));
  if(rootLiesOnPath && isLight){
    // we have seen the tp root along the path to this node, and
    // this node is light enough to be made a leaf
    // therefore, local leaf
    if(nodeNumParticles == 0) node->setType(MultiphaseDistree::EmptyLeaf);
    else node->setType(MultiphaseDistree::Leaf);

    // OK to calculate moments of local leaves
    if(node->getNumParticles() > 0){
      moments.fromParticles(node->getParticles(), node->getNumParticles());
    }

    return true;
  }
  else if(isLight && !Utilities::isPrefix(node->getKey(), node->getDepth(),
              tpRootKey_, Utilities::getDepthFromKey(tpRootKey_))){
    // if we haven't seen the tp root yet, but encounter a light node,
    // this node could still be on the path from the global root to the
    // tp root; therefore, we do the Utilities::isPrefix check. if in 
    // addition to the previous conditions, the node is off-path, then
    // it is a remote node.
    CkAssert(nodeNumParticles == 0);
    // we can distinguish between
    // remote nodes and remote leaves at this point,
    // since we have the number of particles 
    // assigned to all the tree pieces, and in 
    // particular, the number of particles assigned
    // to the tree piece with this node as
    // its root
    int ownerStart = node->getOwnerStart();
    int ownerEnd = node->getOwnerEnd();
    int askWhom = -1;

    CkAssert(node->getNumChildren() == 0);
    CkAssert(node->getChildren() == NULL);

    if(ownerStart + 1 == ownerEnd){
      // single owner, i.e. is a tree piece root
      // how many particles were assigned to this tree piece?
      int nParticlesAssigned = splitters[ownerStart].nParticles;
      if(nParticlesAssigned == 0){
        node->setType(MultiphaseDistree::RemoteEmptyLeaf);
      }
      else if(nParticlesAssigned <= BUCKET_TOLERANCE * Real(parameters.maxppb)){
        node->setType(MultiphaseDistree::RemoteLeaf);
      }
      else{
        node->setType(MultiphaseDistree::Remote);
        // we know that this node is not a leaf, so 
        // it must have children in the tree piece that
        // owns it.
        node->setNumChildren(BRANCH_FACTOR);
      }
    }
    else{
      node->setType(MultiphaseDistree::Remote);
      // again, this node is not a leaf, so it must
      // have children in its owner
      node->setNumChildren(BRANCH_FACTOR);
    }


    return false;
  }

  // we reset moments here, and accumulate the moments of 
  // children that are ready, below. 
  moments.begin();

  Key childKey = (node->getKey() << parameters.nNodeChildrenLg); 

  Particle *particles = node->getParticles();

  int start = 0;
  int finish = start + node->getNumParticles();

  int ownerStart = node->getOwnerStart();
  int ownerEnd = node->getOwnerEnd();
  bool singleOwnerTp = (ownerStart + 1 == ownerEnd);

  node->allocateChildren(parameters.nNodeChildren); 

  int nonLocalChildren = 0;
  for(int i = 0; i < parameters.nNodeChildren; i++){
    int firstGeIdx;

    Key siblingSplitterKey = Utilities::toSplitter(childKey + 1);
    if(i < parameters.nNodeChildren-1){
      firstGeIdx = Utilities::binary_search_ge(siblingSplitterKey, particles, start, finish);
    }
    else{
      firstGeIdx = finish; 
    }
    int np = firstGeIdx - start;

    MyLocalNodeType *child;
    child = new (node->getChild(i)) MyLocalNodeType(childKey, node->getDepth() + 1, particles + start, np);
    child->setParent(node);
    
    // set owner tree pieces of child
    child->setOwnerStart(ownerStart);
    if(singleOwnerTp){
      child->setOwnerEnd(ownerEnd);
    }
    else{
      if(i < parameters.nNodeChildren-1){
        // if you have a right sibling, search through array of assigned
        // splitters and find the  index of the splitter that has the first key
        // GEQ the right sibling's, i.e. the first key not underneath this child
        int firstGeTp = Utilities::binary_search_ge(siblingSplitterKey, &splitters[0], ownerStart, ownerEnd); 
        child->setOwnerEnd(firstGeTp);
        ownerStart = firstGeTp;
      }
      else{
        // for last child, there is no right sibling, and its
        // owner end is same as parent's
        child->setOwnerEnd(node->getOwnerEnd());
      }
    }

    bool local = build(splitters, child, rootLiesOnPath);

    if(local){
      // OK to accumulate moments in node, of those children that are local
      moments.accumulate(child->getPayload().moments());
    }
    else{
      nonLocalChildren++;
    }

    start = firstGeIdx;
    childKey++;
  }

  if(nonLocalChildren == 0){
    node->setType(MultiphaseDistree::Internal);
    // OK to calculate moments of internal nodes
    moments.end();
  }
  else{
    node->setType(MultiphaseDistree::Boundary);
  }

  return (nonLocalChildren == 0);
}

void TreePiece::addLeaf(MyLocalNodeType *leaf){
  myLeaves_.push_back(leaf);
}

void TreePiece::gravity(const CkCallback &cb){
  callback_ = cb;

  gravityOutstanding() = myLeaves_.size();
  nLeavesInitiated_ = 0;
#ifdef DEBUG_TRAVERSAL
  checker_.reset();
#endif

  gravityVisitors_.resize(myLeaves_.size());

  checkDoneGravity();

  myProxy()[thisIndex].doGravity();
}

void TreePiece::doGravity(){
  // get global tree's root
  MyLocalNodeType *root = treeHandle.root();
  GravityVisitor::TraversalType *traversal = gravityTraversalHandle.local();

  int yield = parameters.yieldPeriod;
  // traverse the global tree, starting at its root,
  // for every bucket that you own
  for(int i = nLeavesInitiated_; (i < myLeaves_.size()) && (yield > 0); i++){
    MyLocalNodeType *leaf = myLeaves_[i];
    gravityVisitors_[i] = GravityVisitor(leaf, this);
    walk(traversal, &gravityVisitors_[i], root);

    yield--;
    gravityOutstanding()--;
    nLeavesInitiated_++;
    checkDoneGravity();
  }


  if(nLeavesInitiated_ < myLeaves_.size()){
    CkAssert(gravityOutstanding() > 0);
    myProxy()[thisIndex].doGravity();
  }
}

void TreePiece::sph(const CkCallback &cb){
  sphCallback_ = cb;
  sphOutstanding() = myLeaves_.size();

  //CkPrintf("[%d] TreePiece::sph out %d\n", thisIndex, sphOutstanding());
  nSphLeavesInitiated_ = 0;
  densitySphVisitors_.resize(myLeaves_.size());
  sphLeafData_.resize(myLeaves_.size());

  // in case i don't have any particles at all
  checkDoneSph((DensitySphVisitor *) NULL);

  //printParticlesToFile();

  myProxy()[thisIndex].doSph();
}

void TreePiece::doSph(){
  DensitySphVisitor::TraversalType *traversal = sphDensityTraversalHandle.local();

  int yield = parameters.yieldPeriod;
  // traverse the global tree, starting at its root,
  // for every bucket that you own
  for(int i = nSphLeavesInitiated_; (i < myLeaves_.size()) && (yield > 0); i++){
    MyLocalNodeType *leaf = myLeaves_[i];

    //CkPrintf("[%d] TreePiece::doSph leaf %d outstanding %d\n", thisIndex, i, sphOutstanding());
    // placement new operator to reinitialize sph leaf data
    //SphLeafData *sphLeafData = new (&sphLeafData_[i]) SphLeafData(leaf->getNumParticles(), i);
    SphLeafData *sphLeafData = &sphLeafData_[i];
    //*sphLeafData = SphLeafData(leaf->getNumParticles(), i);
    sphLeafData->initialize(leaf->getNumParticles(), i);

    // for now, we treat all particles as SPH particles
    // so that the bounding box of the leaf is also
    // the bounding box to be used for the SPH comptn.
    OrientedBox<Real> &box = leaf->getPayload().moments().box;
    sphLeafData->center() = box.center();
    sphLeafData->size() = 0.5 * box.size().lengthSquared(); 
    // this signifies that the local walk for this 
    // leaf is still outstanding; this is balanced by
    // the hit() call right after the local walk for
    // the leaf is finished, inside doSph()
    sphLeafData->nOutstanding() = 1;

    //DensitySphVisitor *visitor = new (&densitySphVisitors_[i]) DensitySphVisitor(this, leaf, sphLeafData);
    DensitySphVisitor *visitor = &densitySphVisitors_[i];
    *visitor = DensitySphVisitor(this, leaf, sphLeafData);

    visitor->localLeaf(leaf->getKey(), leaf->getParticles(), leaf->getNumParticles());


    // start traversal with visitor (which is 
    // doing work at the behest of 'leaf'. traversal
    // starts at global tree (source) node 'leaf' 
    // itself
    walk(traversal, visitor, leaf);

    yield--;
    nSphLeavesInitiated_++;
    visitor->hit(leaf->getKey(), true);
    checkDoneSph(visitor);
  }

  if(nSphLeavesInitiated_ < myLeaves_.size()){
    CkAssert(sphOutstanding() > 0);
    myProxy()[thisIndex].doSph();
  }
}

void TreePiece::ballSph(const CkCallback &cb){
  sphCallback_ = cb;
  sphOutstanding() = myLeaves_.size();

  //CkPrintf("[%d] TreePiece::ballSph out %d\n", thisIndex, sphOutstanding());
  nSphLeavesInitiated_ = 0;
  ballSphVisitors_.resize(myLeaves_.size());
  CkAssert(sphLeafData_.size() == ballSphVisitors_.size());

  // in case i don't have any particles at all
  checkDoneSph((BallSphVisitor *) NULL);

  myProxy()[thisIndex].doBallSph();
}

void TreePiece::doBallSph(){
  BallSphVisitor::TraversalType *traversal = sphBallTraversalHandle.local();
  MyLocalNodeType *root = treeHandle.root();

  int yield = parameters.yieldPeriod;
  // traverse the global tree, starting at its root,
  // for every bucket that you own
  for(int i = nSphLeavesInitiated_; (i < myLeaves_.size()) && (yield > 0); i++){
    MyLocalNodeType *leaf = myLeaves_[i];
    SphLeafData *sphLeafData = &sphLeafData_[i];
    sphLeafData->nOutstanding() = 1;

    BallSphVisitor *visitor = new (&ballSphVisitors_[i]) BallSphVisitor(this, leaf, sphLeafData);

    // output ball radius for 0-th leaf
    // of 0-th tree piece
    if(thisIndex == 0 && i == 0){
      ostringstream oss;
      oss << sphLeafData->center();
      oss << sphLeafData->size() + sphLeafData->maxDist();
    }


    // start traversal with visitor (which is 
    // doing work at the behest of 'leaf'. 
    walk(traversal, visitor, root);

    yield--;
    nSphLeavesInitiated_++;
    visitor->hit(leaf->getKey(), true);
    checkDoneSph(visitor);
  }

  if(nSphLeavesInitiated_ < myLeaves_.size()){
    CkAssert(sphOutstanding() > 0);
    myProxy()[thisIndex].doBallSph();
  }
}



void TreePiece::integrate(const CkCallback &cb){
  callback_ = cb;
  IntegrateVisitor integrator(iteration_);
  for(int i = 0; i < myLeaves_.size(); i++){
    integrator.localLeaf(myLeaves_[i]);
  }
  contribute(sizeof(IntegrateVisitor), 
             &integrator, 
             MultiphaseDistree::Reduce<IntegrateVisitor>::reducer(), 
             callback_);

  iteration_++;
}

void TreePiece::pup(PUP::er &p){
  MyPieceType::pup(p);
  p | myNumParticles_;
  p | myParticles_;
  if(p.isPacking()){
    CkAssert(particleMsgs_.size() == 0);
  }
  p | myBox_;
  p | iteration_;
  p | callback_;
  if(p.isUnpacking()){
    gravityOutstanding_ = 0;
    nExpectedParticles_ = -1;
    findOrb3dLb();
  }
  CkAssert(nExpectedParticles_ == -1);
}

void TreePiece::resumeGravityNode(const NodeResumption<GravityVisitor> &r){
  r.visitor.data->hit(r.node.data->getKey());
  //CkPrintf("[%d] TreePiece::Gravity leaf %llu RESUME NODE %llu outstanding %d\n", thisIndex, leaf.data->getKey(), node.data->getKey(), outstanding_);
  GravityVisitor::TraversalType *traversal = gravityTraversalHandle.local();
  const GravityTraversalDataInterface::RemoteNodeType *children = r.node.data->getChildren();
  CkAssert(children != NULL);
  for(int i = 0; i < r.node.data->getNumChildren(); i++){
    walk(traversal, r.visitor.data, &children[i]);
  }

  checkDoneGravity();
}

void TreePiece::resumeGravityLeaf(const LeafResumption<GravityVisitor> &r){
  r.visitor.data->hit(r.visitor.data->leaf()->getKey());
  //CkPrintf("[%d] TreePiece::Gravity leaf %llu RESUME LEAF %llu outstanding %d\n", thisIndex, leaf.data->getKey(), srcKey, outstanding_);
  gravityTraversalData_.pp() += Physics::Gravity::forces<GravityTraversalDataInterface::RemoteParticleType>(r.visitor.data->leaf(), r.particles.data, r.nParticles);
#ifdef DEBUG_TRAVERSAL
  checker_.insertComputedLeaf(r.visitor.data->leaf()->getKey(), srcKey);
  r.visitor.data->leaf()->getUserData().pp()++;
#endif

  checkDoneGravity();
}

void TreePiece::resumeDensitySphNode(const NodeResumption<DensitySphVisitor> &r){
  r.visitor.data->hit(r.node.data->getKey());
  //CkPrintf("[%d] TreePiece::SPH leaf %llu RESUME NODE %llu outstanding %d\n", thisIndex, leaf.data->getKey(), node.data->getKey(), outstanding_);
  DensitySphVisitor::TraversalType *traversal = sphDensityTraversalHandle.local();
  const SphDensityTraversalDataInterface::RemoteNodeType *children = r.node.data->getChildren();
  CkAssert(children != NULL);
  for(int i = 0; i < r.node.data->getNumChildren(); i++){
    walk(traversal, r.visitor.data, &children[i]);
  }

  checkDoneSph(r.visitor.data);
}

void TreePiece::resumeDensitySphLeaf(const LeafResumption<DensitySphVisitor> &r){
  r.visitor.data->hit(r.key);
  //CkPrintf("[%d] TreePiece::SPH leaf %llu RESUME LEAF %llu outstanding %d\n", thisIndex, leaf.data->getKey(), srcKey, outstanding_);
  r.visitor.data->remoteLeaf(r.key, r.particles.data, r.nParticles);

  checkDoneSph(r.visitor.data);
}

void TreePiece::resumeBallSphNode(const NodeResumption<BallSphVisitor> &r){
  r.visitor.data->hit(r.node.data->getKey());
  //CkPrintf("[%d] TreePiece::SPH leaf %llu RESUME NODE %llu outstanding %d\n", thisIndex, leaf.data->getKey(), node.data->getKey(), outstanding_);
  BallSphVisitor::TraversalType *traversal = sphBallTraversalHandle.local();
  const SphBallTraversalDataInterface::RemoteNodeType *children = r.node.data->getChildren();
  CkAssert(children != NULL);
  for(int i = 0; i < r.node.data->getNumChildren(); i++){
    walk(traversal, r.visitor.data, &children[i]);
  }

  checkDoneSph(r.visitor.data);
}

void TreePiece::resumeBallSphLeaf(const LeafResumption<BallSphVisitor> &r){
  r.visitor.data->hit(r.key);
  //CkPrintf("[%d] TreePiece::SPH leaf %llu RESUME LEAF %llu outstanding %d\n", thisIndex, leaf.data->getKey(), srcKey, outstanding_);
  r.visitor.data->remoteLeaf(r.key, r.particles.data, r.nParticles);

  checkDoneSph(r.visitor.data);
}

void TreePiece::checkDoneGravity(){
  if(gravityOutstanding() > 0) return;
  doneGravity();
}

void TreePiece::doneGravity(){
  CkAssert(nLeavesInitiated_ == myLeaves_.size());

#ifdef DEBUG_TRAVERSAL
  gravityTraversalData_.good() = checker_.done();
#endif

  contribute(sizeof(TraversalData), 
             &gravityTraversalData_, 
             MultiphaseDistree::Reduce<TraversalData>::reducer(), 
             callback_);
  gravityTraversalData_.reset();
}

void TreePiece::doneSph(TraversalData &data){
  // sphTraversalData_.good_ is set by 
  // SphLeafData for each leaf in its SphLeafData::check() 
  // method
  sphShadowProxy[thisIndex].ckLocal()->contribute(sizeof(TraversalData), 
                                                  &data, 
                                                  MultiphaseDistree::Reduce<TraversalData>::reducer(), 
                                                  sphCallback_);
  data.reset();
}

// miscellaneous book-keeping functions called by
// BarnesHutVisitor2

#ifdef DEBUG_TRAVERSAL
// FIXME - these should handle both MyLocalNodeType and 
// GravityTraversalDataInterface::RemoteDataType types for 'n'
void TreePiece::preOpen(const MyNodeType *leaf, const MyNodeType *n){
  checker_.preOpen(leaf, n);
}

void TreePiece::insertComputedNode(Key leafKey, Key nodeKey){
  checker_.insertComputedNode(leafKey, nodeKey);
}

void TreePiece::insertComputedLeaf(Key leafKey, Key nodeKey){
  checker_.insertComputedLeaf(leafKey, nodeKey);
}
#endif

CmiUInt8 &TreePiece::open(){
  return gravityTraversalData_.open();
}

CmiUInt8 &TreePiece::pn(){
  return gravityTraversalData_.pn();
}

CmiUInt8 &TreePiece::pp(){
  return gravityTraversalData_.pp();
}

int &TreePiece::gravityOutstanding(){
  return gravityOutstanding_;
}

CmiUInt8 &TreePiece::sphOpen(){
  return sphTraversalData_.open();
}

CmiUInt8 &TreePiece::sphPn(){
  return sphTraversalData_.pn();
}

CmiUInt8 &TreePiece::sphPp(){
  return sphTraversalData_.pp();
}

CmiUInt8 &TreePiece::ballSphOpen(){
  return ballSphTraversalData_.open();
}

int &TreePiece::sphOutstanding(){
  return sphOutstanding_;
}

TraversalData &TreePiece::getGravityTraversalData(){
  return gravityTraversalData_; 
}

TraversalData &TreePiece::getDensitySphTraversalData(){
  return sphTraversalData_; 
}

TraversalData &TreePiece::getBallSphTraversalData(){
  return ballSphTraversalData_; 
}

void TreePiece::findOrb3dLb(){
  LBDatabase *lbdb = LBDatabaseObj();
  int numLB = lbdb->getNLoadBalancers();
  BaseLB **lb = lbdb->getLoadBalancers();

  orb3dLbFound_ = false;
  for(int i = 0; i < numLB; i++){
    if(string(lb[i]->lbName()) == "Orb3dLB"){
      orb3dLbGid_ = lb[i]->getGroupID();
      orb3dLbFound_ = true;
      break;
    }
  }
}

void TreePiece::balanceLoad(const CkCallback &cb){
  callback_ = cb;
  if(thisIndex != 0) return;

  if(orb3dLbFound_){
    CProxy_Orb3dLB(orb3dLbGid_)[0].start(CkCallback(CkIndex_TreePiece::callAtSync(), thisProxy), CkCallbackResumeThread());
    CProxy_TreePiece(thisProxy).sendCentroids();
  }
  else{
    CProxy_TreePiece(thisProxy).callAtSync();
  }
}

void TreePiece::sendCentroids(){
  LDObjHandle handle = myRec->getLdHandle();
  CkAssert(orb3dLbFound_);
  TaggedVector3D tv(myLbCentroid_, handle, 0, myParticles_.size(), 0, 0);
  tv.tag = thisIndex;
  CkCallback lbcb(CkIndex_Orb3dLB::receiveCentroids(NULL), 0, CProxy_Orb3dLB(orb3dLbGid_));
  contribute(sizeof(TaggedVector3D), (char *) &tv, CkReduction::concat, lbcb);
}

void TreePiece::callAtSync(){
  AtSync();
}

void TreePiece::ResumeFromSync(){
  //CkPrintf("[%d] TreePiece::ResumeFromSync on node %d\n", thisIndex, CkMyNode());
  contribute(callback_);
}


