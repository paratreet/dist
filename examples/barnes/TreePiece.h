#ifndef TREE_PIECE_H
#define TREE_PIECE_H

#include <map>

#include "Particle.h"
#include "BoundingBox.h"
#include "ParticleMsg.h"

#include "barnes.decl.h"

#include "mdt.h"
#include "DataInterface.h"
#include "TemplateInstanceTypedefs.h"
#include "Visitor.h"
#include "TraversalData.h"
#include "PointerContainer.h"


class TreePiece : public CBase_TreePiece {
  int myNumParticles_;
  int nExpectedParticles_;
  CkVec<Particle> myParticles_;
  CkVec<ParticleMsgWrapper> particleMsgs_;
  BoundingBox myBox_;

  CkVec<MyLocalNodeType *> myLeaves_;

  // number of pieces actually used; obtained 
  // during flush() 
  int nPieces_;

  // saved callback, in case a phase spans several 
  // entry method invocations
  CkCallback callback_;
  CkCallback sphCallback_;

  MyLocalNodeType *root_;
  // needed while building tree local tree to see
  // which nodes are off path from global root to 
  // local root (and therefore remote) 
  Key tpRootKey_;

  int iteration_;
  // number of outstanding remote data requests 
  int gravityOutstanding_;
  int sphOutstanding_;
  int ballSphOutstanding_;

  // number of leaves for which I have already
  // started the top-down traversal
  int nLeavesInitiated_;
  // same as above, for sph
  int nSphLeavesInitiated_;
  // same as above, for sph
  int nBallSphLeavesInitiated_;

  // visitors and stats associated with 
  // their traversals
  TraversalData gravityTraversalData_;
  CkVec<GravityVisitor> gravityVisitors_;

  TraversalData sphTraversalData_;
  CkVec<DensitySphVisitor> densitySphVisitors_;

  TraversalData ballSphTraversalData_;
  CkVec<BallSphVisitor> ballSphVisitors_;

  // the two sph phases use this data structure
  // to store neighbor information; the search radii
  // set by the density sph visitor are used by the 
  // ball sph visitor.
  CkVec<SphLeafData> sphLeafData_;

#ifdef DEBUG_TRAVERSAL
  // this object contains code to check the
  // correctness of traversals
  BarnesHutVisitor checker_;
#endif

  Vector3D<Real> myLbCentroid_;
  bool orb3dLbFound_;
  CkGroupID orb3dLbGid_;

  public:
  TreePiece();
  TreePiece(CkMigrateMessage *m);

  void load(const CkCallback &cb);
  void assignKeys(BoundingBox &box, const CkCallback &cb);
  void count(CkVec<Key> &keys, const CkCallback &cb);
  void flush(const CkCallback &cb);
  void particles(ParticleMsg *msg);

  void printParticlesToFile();
  void build(const CkCallback &cb);

  void gravity(const CkCallback &cb);
  void doGravity();
  void sph(const CkCallback &cb);
  void doSph();
  void ballSph(const CkCallback &cb);
  void doBallSph();
   
  void integrate(const CkCallback &cb);
  void balanceLoad(const CkCallback &cb);
  void sendCentroids();
  void callAtSync();
  void ResumeFromSync();

  public:
  // called by Manager to tell TreePiece that 
  // it has anothoer leaf
  void addLeaf(MyLocalNodeType *leaf);

  void pup(PUP::er &p);

  // called by BarnesHutVisitor2 in order to do 
  // some book-keeping
#ifdef DEBUG_TRAVERSAL
  void preOpen(const MyLocalNodeType *leaf, const MyLocalNodeType *node); 

  void insertComputedNode(Key leafKey, Key nodeKey);
  void insertComputedLeaf(Key leafKey, Key nodeKey);

  // for sph
  bool &sphGood(){
    return sphTraversalData_.good();
  }
#endif

  CmiUInt8 &open();
  CmiUInt8 &pn();
  CmiUInt8 &pp();

  CmiUInt8 &sphOpen();
  CmiUInt8 &sphPn();
  CmiUInt8 &sphPp();

  CmiUInt8 &ballSphOpen();

  int &gravityOutstanding();
  int &sphOutstanding();

  TraversalData &getGravityTraversalData();
  TraversalData &getDensitySphTraversalData();
  TraversalData &getBallSphTraversalData();

  CProxyElement_ArrayElement &getArrayElement();

  // RESUMPTION ENTRY METHODS
  void resumeGravityNode(const NodeResumption<GravityVisitor> &resumption);
  void resumeGravityLeaf(const LeafResumption<GravityVisitor> &resumption);

  // density sph traversal
  void resumeDensitySphNode(const NodeResumption<DensitySphVisitor> &resumption);
  void resumeDensitySphLeaf(const LeafResumption<DensitySphVisitor> &resumption);

  // ball sph traversal
  void resumeBallSphNode(const NodeResumption<BallSphVisitor> &resumption);
  void resumeBallSphLeaf(const LeafResumption<BallSphVisitor> &resumption);

  template<typename SphVisitorType>
  void checkDoneSph(SphVisitorType *v){
    // visitor is NULL when checkDoneSph is called
    // from sph(). this happens when the TreePiece
    // has no leaves
    if(v != NULL){
      if(v->checkDone()){
        //CkPrintf("[%d] Visitor for leaf %llu done!\n", thisIndex, v->leaf()->getKey());
        v->kernel();
        if(thisIndex == 0 && v->data()->getLeafNum() == 0){
          outputParticleCloud(v);
        }
      }
    }

    if(sphOutstanding() > 0) return;

    CkAssert(nSphLeavesInitiated_ <= myLeaves_.size());
    if(v == NULL){
      TraversalData dummy;
      doneSph(dummy);
    }
    else{
      doneSph(v->getTraversalData());
    }
  }


 
  private:
  void loadTipsy();
  void convertKeysToSplitters(CkVec<Key> &keys);
  void copyParticles();
  bool build(const CkVec<Splitter> &splitters, MyLocalNodeType *root, bool rootLiesOnPath);

  template<typename VisitorType, typename AnyNodeType>
  void walk(typename VisitorType::TraversalType *traversal, VisitorType *visitor, const AnyNodeType *node){
    /*
       if(node->getKey() == Key(1) && leaf->getKey() == Key(143)){
       CkPrintf("here!\n");
       }
     */
    // check whether we should open this node
    traversal->go(node, visitor);
  }


  void checkDoneGravity();
  void doneGravity();

  void doneSphLeaf();
  void doneSph(TraversalData &data);


  template<typename SphVisitorType>
  void outputParticleCloud(SphVisitorType *v){
    return;
    // open particle file
    ostringstream oss;
    oss << "plot/" << v->name() << ".particle." << thisIndex << ".pos";
    ofstream particleFile(oss.str().c_str());
    CkAssert(particleFile.is_open());

    // open neighbor file
    ostringstream noss;
    noss << "plot/" << v->name() << ".neighbors." << thisIndex << ".pos";
    ofstream neighborFile(noss.str().c_str());
    CkAssert(neighborFile.is_open());

    for(int i = 0; i < v->leaf()->getNumParticles(); i++){
      const Particle &p = v->leaf()->getParticles()[i];
      particleFile << *(p.getSphCore()) << endl;
      const SphParticleData::StorageType &neighbors = v->data()->particleData(i).neighbors().storage();
      for(int i = 0; i < neighbors.size(); i++){
        neighborFile << *neighbors[i].core() << endl;
      }
    }

    // close particle file
    particleFile.close();
    neighborFile.close();
  }

  void print(MyLocalNodeType *root);
  CProxy_TreePiece myProxy() const;

  private:
  void findOrb3dLb();
};

#define CK_TEMPLATES_ONLY
#include "barnes.def.h"
#undef CK_TEMPLATES_ONLY

#endif // TREE_PIECE_H
