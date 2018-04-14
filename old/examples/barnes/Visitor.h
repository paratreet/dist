#ifndef VISITOR_H
#define VISITOR_H

#include "defines.h"
#include "TemplateInstanceTypedefs.h"
#include "BoundingBox.h"
#include "SphData.h"
#include "Physics.h"
#include "DataInterface.h"
#include "Resumption.h"

#ifdef DEBUG_TRAVERSAL
#include <map>
#include <set>
#include <sstream>
#endif

class CkArrayID;
class CkArrayIndex;
class CProxyElement_ArrayElement;
class TreePiece;

/*
class BarnesHutVisitor {
  public:

  // the part of the state of this visitor that
  // cannot be shared across several visitors

  class State {
    MyNodeType *leaf_;

    public:
    State();

    void reset();

    MyNodeType *&leaf(){
      return leaf_;
    }

    void pup(PUP::er &p){
    }
  };


  private:
  State state_;

  // the part of the state of this visitor that can
  // be shared across several visitor instances

  int outstanding_;

  CmiUInt8 pn_;
  CmiUInt8 pp_;
  CmiUInt8 open_;

  int tp_;

#ifdef DEBUG_TRAVERSAL
  Key checkKey_;
  bool noFilter_;
  std::ostringstream out_;

  // keep track of frontier of computed
  // nodes; these should collapse to root key when  
  // traversal is done.
  std::map<Key, std::set<Key> > computed_; 
#endif

  public:
  BarnesHutVisitor() {}
  BarnesHutVisitor(const BarnesHutVisitor &other);
  BarnesHutVisitor(int tp);
  BarnesHutVisitor(int tp, Key check);
  void pup(PUP::er &p);

  public:
  // will be called by traversal manager
  void reset();
  bool node(const MyNodeType *node);
  void localLeaf(Key srcKey, const Particle *sources, int nSources);
  void remoteLeaf(Key srcKey, const RemoteParticle *sources, int nSources);

  // called by traversal manager when traversal
  // is suspended pending receipt of remote data
  void miss(Key key);
  int &outstanding(){
    return outstanding_;
  }

  MyNodeType *&leaf(){
    return state_.leaf();
  }
  // called by traversal manager when traversal
  // is resumed, nullifying the effect of a previous
  // miss() event
  void hit(Key key);

  // this method is called when traversal finishes
  bool done();

  const State &state() const;
  void state(const State &s);
  
  BarnesHutVisitor &operator+=(const BarnesHutVisitor &other){
    pn() += other.pn();
    pp() += other.pp();
    open() += other.open();
    return *this;
  }

  CmiUInt8 &open(){
    return open_;
  }

  CmiUInt8 &pn(){
    return pn_;
  }

  CmiUInt8 &pp(){
    return pp_;
  }

  // const versions, used in operator+=
  const CmiUInt8 &open() const {
    return open_;
  }

  const CmiUInt8 &pn() const {
    return pn_;
  }

  const CmiUInt8 &pp() const {
    return pp_;
  }

#ifdef DEBUG_TRAVERSAL
  bool checkTraversal();

  void preOpen(const MyNodeType *leaf, const MyNodeType *node);
  void insertComputedNode(Key bucket, Key node);
  void insertComputedLeaf(Key bucket, Key node);

  void insertComputed(Key bucket, Key node);
#endif

};

class InterCheckVisitor {
  public:
  void localLeaf(const MyNodeType *leaf) const;
  void done() const {}
  InterCheckVisitor &operator+=(const InterCheckVisitor &other){
    return *this;
  }

  void pup(PUP::er &p) {}
};
*/

class IntegrateVisitor {
  int iteration_;
  BoundingBox box_;

  public:
  IntegrateVisitor() : 
    iteration_(-1)
  {}

  IntegrateVisitor(int iteration) : 
    iteration_(iteration)
  {}

  void localLeaf(MyLocalNodeType *leaf);
  void done() const;
  IntegrateVisitor &operator+=(const IntegrateVisitor &other){
    CkAssert(iteration() == other.iteration());
    box().grow(other.box());
    return *this;
  }

  void pup(PUP::er &p){
    p | iteration_;
  }

  public:
  int iteration() const { return iteration_; }
  const BoundingBox &box() const { return box_; }
  BoundingBox &box() { return box_; }

};

class GravityVisitor {
  public:
  typedef MyGravityTraversalType TraversalType;
  typedef MyGravityTraversalHandleType TraversalHandleType;
  typedef GravityTraversalDataInterface::RemoteParticleType RemoteParticleType;
  typedef GravityTraversalDataInterface::RemoteNodeType RemoteNodeType;

  private:
  TreePiece *tp_;
  MyLocalNodeType *leaf_;

  public:
  GravityVisitor() :
    leaf_(NULL),
    tp_(NULL)
  {}

  GravityVisitor(MyLocalNodeType *leaf, TreePiece *tp) : 
    leaf_(leaf),
    tp_(tp)
  {}

  // will be called by traversal manager.
  // has the same behavior for both nodes in 
  // this address space, and remotely fetched
  // ones, so can template on node type
  template<typename AnyNodeType>
  bool node(const AnyNodeType *node);

  // should only be called with particles from this
  // address space, so not templated on particle type
  void localLeaf(Key srcKey, const Particle *sources, int nSources);

  // only called with particles from some remote address
  // space, therefore expects remote particles (typedef at 
  // start of class)
  void remoteLeaf(Key srcKey, const RemoteParticleType *sources, int nSources);

  // called by traversal manager when traversal
  // is suspended pending receipt of remote data
  void miss(Key key);

  void deliver(RemoteNodeType *node, int chunk);
  void deliver(Key srcKey, RemoteParticleType *particles, int nParticles, int chunk);


  TraversalData &getTraversalData();

  MyLocalNodeType *&leaf(){
    return leaf_;
  }

  TreePiece *&tp(){
    return tp_;
  }

  // called by traversal manager when traversal
  // is resumed, nullifying the effect of a previous
  // miss() event
  void hit(Key key);
};

class BaseSphVisitor {
  public:

  TreePiece *tp_;
  MyLocalNodeType *leaf_;
  SphLeafData *data_;
  std::string name_;

  protected:
  BaseSphVisitor() :
    leaf_(NULL),
    tp_(NULL),
    data_(NULL),
    name_("BAD")
  {}

  BaseSphVisitor(TreePiece *tp, MyLocalNodeType *leaf, SphLeafData *data, const std::string &name) : 
    leaf_(leaf),
    tp_(tp),
    data_(data),
    name_(name)
  {}


  public:
  // called by traversal manager when traversal
  // is suspended pending receipt of remote data
  void miss(Key sourceKey);
  // called by traversal manager when traversal
  // is resumed, nullifying the effect of a previous
  // miss() event
  void hit(Key sourceKey, bool fake=false);

  virtual TraversalData &getTraversalData() = 0;

  MyLocalNodeType *&leaf(){
    return leaf_;
  }

  TreePiece *&tp(){
    return tp_;
  }

  SphLeafData *&data(){
    return data_;
  }

#ifdef DEBUG_TRAVERSAL
  Key &checkKey(){
    return checkKey_;
  }
#endif

  std::string &name(){
    return name_;
  }
};

class DensitySphVisitor : public BaseSphVisitor {
  public:
  typedef MySphDensityTraversalType TraversalType;
  typedef MySphDensityTraversalHandleType TraversalHandleType;
  typedef SphDensityTraversalDataInterface::RemoteNodeType RemoteNodeType;
  typedef SphDensityTraversalDataInterface::RemoteParticleType RemoteParticleType;

  private:
#ifdef DEBUG_TRAVERSAL
  Key checkKey_;
#endif

  public:
  DensitySphVisitor() : 
    BaseSphVisitor()
  {
#ifdef DEBUG_TRAVERSAL
    checkKey_ = Key(0);
#endif
  }

  DensitySphVisitor(TreePiece *tp, MyLocalNodeType *leaf, SphLeafData *data) : 
    BaseSphVisitor(tp, leaf, data, "density")
  {
#ifdef DEBUG_TRAVERSAL
    checkKey_ = Key(0);
#endif
  }


  // will be called by traversal manager
  template<typename AnyNodeType>
  bool node(const AnyNodeType *n);

  void localLeaf(Key sourceKey, const Particle *sources, int nSources);
  void remoteLeaf(Key sourceKey, const RemoteParticleType *sources, int nSources);

  // called from TreePiece::checkDoneSph() to check whether
  // we have finished the tree traversal for the leaf to 
  // which this SphVisitor currently belongs.
  bool checkDone();
  void kernel();

  // for resumption of traversals
  void deliver(RemoteNodeType *node, int chunk);
  void deliver(Key srcKey, RemoteParticleType *particles, int nParticles, int chunk);

  TraversalData &getTraversalData();
};

// useful for debugging sph: get rmax as defined by
// sph walk, and check (visually or otherwise) whether 
// particles within rmax ball are the same as the 
// sph particles
class BallSphVisitor : public BaseSphVisitor {
  public:
  typedef MySphBallTraversalType TraversalType;
  typedef MySphBallTraversalHandleType TraversalHandleType;
  typedef SphBallTraversalDataInterface::RemoteNodeType RemoteNodeType;
  typedef SphBallTraversalDataInterface::RemoteParticleType RemoteParticleType;

  public:
  BallSphVisitor() : 
    BaseSphVisitor()
  {}

  BallSphVisitor(TreePiece *tp, MyLocalNodeType *leaf, SphLeafData *data) : 
    BaseSphVisitor(tp, leaf, data, "ball")
  {}

  public:
  // will be called by traversal manager
  template<typename AnyNodeType>
  bool node(const AnyNodeType *n);

  void localLeaf(Key sourceKey, const Particle *sources, int nSources);
  void remoteLeaf(Key sourceKey, const RemoteParticleType *sources, int nSources);

  // called to check whether the traversal 
  // for the owner leaf has finished, so that it
  // is ok to do computation
  bool checkDone();
  void kernel();

  // for resumption of traversals
  void deliver(RemoteNodeType *node, int chunk);
  void deliver(Key srcKey, RemoteParticleType *particles, int nParticles, int chunk);

  TraversalData &getTraversalData();
};


#endif // VISITOR_H
