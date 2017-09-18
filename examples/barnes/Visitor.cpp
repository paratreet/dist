#include "mdt.h"
#include "Visitor.h"
#include "TemplateInstanceTypedefs.h"
#include "Physics.h"
#include "TreePiece.h"
#include <sstream>
#include "Visitor_def.h"

// VISITOR STATE

#define VISITOR_VERBOSE CkPrintf

extern CProxy_TreePiece treePieceProxy;

/*
BarnesHutVisitor::State::State(){
  reset();
}

void BarnesHutVisitor::State::reset(){
  leaf_ = NULL;
}

// TREE VISITOR

BarnesHutVisitor::BarnesHutVisitor(const BarnesHutVisitor &other){
  state_ = other.state_;
  outstanding_ = other.outstanding_;
  pn_ = other.pn_;
  pp_ = other.pp_;
  open_ = other.open_;
}


BarnesHutVisitor::BarnesHutVisitor(int tp){
#ifdef DEBUG_TRAVERSAL
  checkKey_ = 0;
  noFilter_ = false;
#endif
  tp_ = tp;
  reset();
}

BarnesHutVisitor::BarnesHutVisitor(int tp, Key check){
#ifdef DEBUG_TRAVERSAL
  checkKey_ = check;
  noFilter_ = false;
#endif
  tp_ = tp;
  reset();
}

void BarnesHutVisitor::reset(){
  state_.reset();
  pn_ = 0;
  pp_ = 0;
  open_ = 0;
  outstanding_ = 0;
#ifdef DEBUG_TRAVERSAL
  computed_.clear();
  out_.str("");
#endif
}

#ifdef DEBUG_TRAVERSAL
void BarnesHutVisitor::preOpen(const MyNodeType *leaf, const MyNodeType *node){
  if(noFilter_ || leaf->getKey() == checkKey_){
    out_ << "[" << tp_ 
         << "] BarnesHutVisitor::preOpen leaf " 
         << leaf->getKey() 
         << " (" 
         //<< leaf->getUserData().moments().com
         << ") node " 
         << node->getKey() 
         << " (" 
         << node->getUserData().moments().com
         << ")" << endl;
    CkPrintf("%s", out_.str().c_str());
    out_.str("");
  }
}
#endif

bool BarnesHutVisitor::node(const MyNodeType *n){
  MyNodeType *leaf = state_.leaf();
  CkAssert(n->getType() != MultiphaseDistree::Invalid);

#ifdef DEBUG_TRAVERSAL
  if(noFilter_ || leaf->getKey() == checkKey_) VISITOR_VERBOSE("[%d] BarnesHutVisitor::node leaf %llu node %llu\n", tp_, leaf->getKey(), n->getKey());
#endif

  if(Physics::Gravity::open(leaf, n)){
    open()++;
    return true;
  }
  else{
    int nInter = Physics::Gravity::forces(leaf, n);
    pn() += nInter;
#ifdef DEBUG_TRAVERSAL
    insertComputedNode(leaf->getKey(), n->getKey());
    leaf->getUserData().pn() += nInter;
#endif
    return false;
  }
}

void BarnesHutVisitor::localLeaf(Key srcKey, const Particle *sources, int nSources){
  MyNodeType *leaf = state_.leaf();
#ifdef DEBUG_TRAVERSAL
  if(noFilter_ || leaf->getKey() == checkKey_) VISITOR_VERBOSE("[%d] BarnesHutVisitor::local leaf %llu node %llu\n", tp_, leaf->getKey(), srcKey);
#endif
  int n = Physics::Gravity::forces(leaf, sources, nSources);
  pp() += n;
#ifdef DEBUG_TRAVERSAL
  insertComputedLeaf(leaf->getKey(), srcKey);
  leaf->getUserData().pp() += n; 
#endif
}

void BarnesHutVisitor::remoteLeaf(Key srcKey, const RemoteParticle *sources, int nSources){
  MyNodeType *leaf = state_.leaf();
#ifdef DEBUG_TRAVERSAL
  if(noFilter_ || leaf->getKey() == checkKey_) VISITOR_VERBOSE("[%d] BarnesHutVisitor::remote leaf %llu node %llu\n", tp_, leaf->getKey(), srcKey);
#endif
  int n = Physics::Gravity::forces(leaf, sources, nSources);
  pp() += n;
#ifdef DEBUG_TRAVERSAL
  insertComputedLeaf(leaf->getKey(), srcKey);
  leaf->getUserData().pp() += n;
#endif
}

#ifdef DEBUG_TRAVERSAL
void BarnesHutVisitor::insertComputedNode(Key leaf, Key node){
  if(noFilter_ || leaf == checkKey_){
    out_ << "[" << tp_ << "] BarnesHutVisitor::computed leaf " << leaf << " node " << node << endl;
    CkPrintf("%s", out_.str().c_str());
    out_.str("");
  }
  insertComputed(leaf, node);
}

void BarnesHutVisitor::insertComputedLeaf(Key leaf, Key node){
  if(noFilter_ || leaf == checkKey_){
    out_ << "[" << tp_ << "] BarnesHutVisitor::computed leaf " << leaf << " leaf " << node << endl;
    CkPrintf("%s", out_.str().c_str());
    out_.str("");
  }
  insertComputed(leaf, node);
}

void BarnesHutVisitor::insertComputed(Key leaf, Key node){
#ifndef LOW_OVERHEAD_DEBUG_TRAVERSAL
  std::set<Key> &leafKeys = computed_[leaf];
  std::set<Key>::iterator it;
  do{
    node = (node >> 1);
    it = leafKeys.find(node);
    if(noFilter_ || leaf == checkKey_){
      out_ << "[" << tp_ << "] Set::find " << node <<  " size: " << leafKeys.size() <<  endl;
      CkPrintf("%s", out_.str().c_str());
      out_.str("");
    }
    if(it == leafKeys.end()){
      leafKeys.insert(node);
      break;
    }
    else{
      leafKeys.erase(it);
    }
  }
  while(node > Key(0));
#endif
}
#endif



void BarnesHutVisitor::hit(Key key){
#ifdef DEBUG_TRAVERSAL
  if(noFilter_ || state_.leaf()->getKey() == checkKey_) VISITOR_VERBOSE("[%d] BarnesHutVisitor::hit leaf %llu key %llu outstanding %d\n", tp_, state_.leaf()->getKey(), key, outstanding());
#endif
}

void BarnesHutVisitor::miss(Key key){
#ifdef DEBUG_TRAVERSAL
  if(noFilter_ || state_.leaf()->getKey() == checkKey_) VISITOR_VERBOSE("[%d] BarnesHutVisitor::miss leaf %llu key %llu outstanding %d\n", tp_, state_.leaf()->getKey(), key, outstanding());
#endif
}

#ifdef DEBUG_TRAVERSAL
bool BarnesHutVisitor::checkTraversal(){
  std::map<Key, set<Key> >::iterator it;
  bool good = true;
  for(it = computed_.begin(); it != computed_.end(); ++it){
    std::set<Key> &bkeys = it->second;
    std::set<Key>::iterator jt;
    if(bkeys.size() != 1 || bkeys.find(Key(0)) == bkeys.end()){
      std::ostringstream oss;
      for(jt = bkeys.begin(); jt != bkeys.end(); ++jt){
        oss << *jt << ", ";
      }
      CkPrintf("BarnesHutVisitor::checkTraversal leaf %llu BAD traversal: [%s]\n", 
                 it->first, oss.str().c_str());
      good = false;
    }
  }

  computed_.clear();
  return good;
}
#endif


const BarnesHutVisitor::State &BarnesHutVisitor::state() const {
  return state_;
};

void BarnesHutVisitor::state(const BarnesHutVisitor::State &c){
  state_ = c;
};

void BarnesHutVisitor::pup(PUP::er &p){
  p | state_;
  p | pn_;
  p | pp_;
  p | open_;

}

bool BarnesHutVisitor::done(){
#ifdef DEBUG_TRAVERSAL
  //CkPrintf("[%d] BarnesHutVisitor::done:\n%s", tp_, out_.str().c_str());
#ifndef LOW_OVERHEAD_DEBUG_TRAVERSAL
  return checkTraversal();
#endif
#endif
}

void InterCheckVisitor::localLeaf(const MyNodeType *leaf) const {
  const NodePayload &pl = leaf->getUserData();
#ifdef DEBUG_TRAVERSAL
  CkPrintf("[%d] leaf %llu pn %llu pp %llu\n", CkMyPe(), leaf->getKey(), pl.pn(), pl.pp());
#endif
}
*/

// ********* PARTICLE TRAJECTORY INTEGRATION *********

void IntegrateVisitor::localLeaf(MyLocalNodeType *leaf){
  Particle *particles = leaf->getParticles();
  for(int i = 0; i < leaf->getNumParticles(); i++){
    Particle &p = particles[i];
#ifndef STATIC
    Physics::integrate(p, parameters.dtime, iteration_);
#endif
    box().grow(p.position());
    box().pe += p.mass() * p.potential;
    box().ke += 0.5 * p.mass() * p.velocity().lengthSquared();
    box().mass += p.mass();

    p.acceleration = 0.0;
    p.potential = 0.0;
  }

  box().numParticles += leaf->getNumParticles();
}

void IntegrateVisitor::done() const {
  ostringstream oss;
  oss << box(); 
  //CkPrintf("[%d] IntegrateVisitor: %s\n", CkMyPe(), oss.str().c_str());
}

// ********* GOOD GRAVITY *********
void GravityVisitor::localLeaf(Key sourceKey, const Particle *sources, int nSources){
  tp()->pp() += Physics::Gravity::forces<Particle>(leaf_, sources, nSources);
#ifdef DEBUG_TRAVERSAL
  tp()->insertComputedLeaf(leaf_->getKey(), sourceKey);
  leaf_->getUserData().pp()++;
#endif
}

void GravityVisitor::remoteLeaf(Key sourceKey, const RemoteParticleType *sources, int nSources){
  tp()->pp() += Physics::Gravity::forces<RemoteParticleType>(leaf_, sources, nSources);
#ifdef DEBUG_TRAVERSAL
  tp()->insertComputedLeaf(leaf_->getKey(), sourceKey);
  leaf_->getUserData().pp()++;
#endif
}

void GravityVisitor::miss(Key sourceKey){
  tp()->gravityOutstanding()++;
  getTraversalData().misses()++;
}

void GravityVisitor::hit(Key sourceKey){
  tp()->gravityOutstanding()--;
  getTraversalData().hits()++;
}

void GravityVisitor::deliver(RemoteNodeType *node, int chunk){
  treePieceProxy[tp()->getIndex()].resumeGravityNode(NodeResumption<GravityVisitor>(node, chunk, this));
}

void GravityVisitor::deliver(Key srcKey, RemoteParticleType *particles, int nParticles, int chunk){
  treePieceProxy[tp()->getIndex()].resumeGravityLeaf(LeafResumption<GravityVisitor>(srcKey, particles, nParticles, chunk, this));
}

// ********* SPH STUFF *********
void DensitySphVisitor::localLeaf(Key sourceKey, const Particle *sources, int nSources){
  for(int i = 0; i < nSources; i++){
    Physics::Sph::compare(&sources[i], leaf(), data());
  }
#ifdef DEBUG_TRAVERSAL
  data()->insertComputedLeaf(sourceKey, checkKey_);
#endif
}

void DensitySphVisitor::remoteLeaf(Key sourceKey, const RemoteParticleType *sources, int nSources){
  for(int i = 0; i < nSources; i++){
    Physics::Sph::compare(&sources[i], leaf(), data());
  }
#ifdef DEBUG_TRAVERSAL
  data()->insertComputedRemoteLeaf(sourceKey, checkKey_);
#endif
}

bool DensitySphVisitor::checkDone(){
  if(data()->nOutstanding() == 0){
    CkAssert(leaf()->getNumParticles() == data()->nParticles());
#ifdef DEBUG_TRAVERSAL
    CkAssert(data()->ownerLeafKey() == leaf()->getKey());
    if(data()->ownerLeafKey() == checkKey_){
      std::ostringstream oss;
      for(int i = 0; i < data()->nParticles(); i++){
        const SphParticleData<ParticleCore> &targetSphData = data()->particleData(i);
        if(targetSphData.numNeighbors() > 0){
          CkPrintf("leaf %llu rmax %f\n", data()->ownerLeafKey(), targetSphData.maxDist2());
        }
      }
    }

    tp()->sphGood() &= data()->check(checkKey_);
#endif
    return true;
  }

  return false;
}

void DensitySphVisitor::kernel(){
  Particle *targets = leaf()->getParticles();
  for(int i = 0; i < leaf()->getNumParticles(); i++){
    Particle &t = targets[i];
    const SphParticleData::StorageType &sources = data()->particleData(i).neighbors().storage();
    for(int j = 0; j < sources.size(); j++){
      const SphParticle &jt = sources[j];
      Physics::Sph::density(t, jt.core(), sqrt(jt.dist2()));
    }

    t.pressure() = parameters.soundSpeed2 * (t.density() - parameters.density0);
    // no need to keep around pointers to 
    // neighbors around once we've used them 
    data()->particleData(i).free();
  }
}

void DensitySphVisitor::deliver(RemoteNodeType *node, int chunk){
  treePieceProxy[tp()->getIndex()].resumeDensitySphNode(NodeResumption<DensitySphVisitor>(node, chunk, this));
}

void DensitySphVisitor::deliver(Key srcKey, RemoteParticleType *particles, int nParticles, int chunk){
  treePieceProxy[tp()->getIndex()].resumeDensitySphLeaf(LeafResumption<DensitySphVisitor>(srcKey, particles, nParticles, chunk, this));
}


void BallSphVisitor::localLeaf(Key sourceKey, const Particle *sources, int nSources){
  for(int i = 0; i < nSources; i++){
    Physics::BallSph::compare(&sources[i], leaf(), data());
  }
#ifdef DEBUG_TRAVERSAL
  data()->insertComputedLeaf(sourceKey, checkKey_);
#endif
}

void BallSphVisitor::remoteLeaf(Key sourceKey, const RemoteParticleType *sources, int nSources){
  for(int i = 0; i < nSources; i++){
    Physics::BallSph::compare(&sources[i], leaf(), data());
  }
#ifdef DEBUG_TRAVERSAL
  data()->insertComputedRemoteLeaf(sourceKey, checkKey_);
#endif
}

bool BallSphVisitor::checkDone(){
  if(data()->nOutstanding() == 0){
    CkAssert(leaf()->getNumParticles() == data()->nParticles());
    return true;
  }

  return false;
}

void BallSphVisitor::kernel(){
  Particle *targets = leaf()->getParticles();
  for(int i = 0; i < leaf()->getNumParticles(); i++){
    Particle &t = targets[i];
    const SphParticleData::StorageType &sources = data()->particleData(i).neighbors().storage();
    for(int j = 0; j < sources.size(); j++){
      const SphParticle &jt = sources[j];
      //Physics::BallSph::compute(t, jt.core(), sqrt(jt.dist2()));
    }

    //t.pressure() = parameters.soundSpeed2 * (t.density() - parameters.density0);

    // no need to keep around pointers to 
    // neighbors around once we've used them 
    data()->particleData(i).free();
  }

}

void BallSphVisitor::deliver(RemoteNodeType *node, int chunk){
  treePieceProxy[tp()->getIndex()].resumeBallSphNode(NodeResumption<BallSphVisitor>(node, chunk, this));
}

void BallSphVisitor::deliver(Key srcKey, RemoteParticleType *particles, int nParticles, int chunk){
  treePieceProxy[tp()->getIndex()].resumeBallSphLeaf(LeafResumption<BallSphVisitor>(srcKey, particles, nParticles, chunk, this));
}


void BaseSphVisitor::hit(Key sourceKey, bool fake){
  // see comments above
  tp()->sphOutstanding()--;
  data()->nOutstanding()--;
  //CkPrintf("[%d] BaseSphVisitor::hit leaf %llu out %d leafout %d\n", tp()->getIndex(), leaf()->getKey(), tp()->sphOutstanding(), data()->nOutstanding());
  if(!fake) getTraversalData().hits()++;
}

void BaseSphVisitor::miss(Key sourceKey){
  // this is the counter across all leaves 
  // of the TreePiece
  tp()->sphOutstanding()++;
  // this is the counter for the particular leaf
  // to which this visitor currently belongs
  data()->nOutstanding()++;
  //CkPrintf("[%d] BaseSphVisitor::miss leaf %llu out %d leafout %d\n", tp()->getIndex(), leaf()->getKey(), tp()->sphOutstanding(), data()->nOutstanding());
  getTraversalData().misses()++;
}


TraversalData &GravityVisitor::getTraversalData(){
  return tp()->getGravityTraversalData();
}

TraversalData &DensitySphVisitor::getTraversalData(){
  return tp()->getDensitySphTraversalData();
}

TraversalData &BallSphVisitor::getTraversalData(){
  return tp()->getBallSphTraversalData();
}
