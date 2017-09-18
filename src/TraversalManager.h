#ifndef MULTIPHASE_DISTREE_TRAVERSAL_MANAGER_H
#define MULTIPHASE_DISTREE_TRAVERSAL_MANAGER_H

#include "mdt.decl.h"
#include "Node.h"
#include "Utility.h"
#include "CkCache.h"
#include "Handle.h"
#include "CachedData.h"
#include "Reduce.h"
#include "TraversalHandle.h"
#include "TraversalStatistics.h"

#define MDT_TRAVERSAL_VERBOSE //CkPrintf

/*
 * To assign priorities to different tasks.
 */
#define NUM_PRIORITY_BITS (sizeof(int)*CHAR_BIT)

#define REQUEST_MOMENTS_PRIORITY (-9)
#define RECV_MOMENTS_PRIORITY (-8)
#define REQUEST_NODE_PRIORITY (-7)
#define RECV_NODE_PRIORITY (-6)
#define REQUEST_PARTICLES_PRIORITY (-5)
#define REMOTE_GRAVITY_PRIORITY (-4)
#define START_TRAVERSAL_PRIORITY (-3)
#define RECV_PARTICLES_PRIORITY (-2)
#define LOCAL_GRAVITY_PRIORITY (-1)

namespace MultiphaseDistree {

template<typename VisitorType, typename StateType>
struct Resumption {
  VisitorType *visitor;
  StateType state; 

  Resumption(VisitorType *v, const StateType &s) : 
    visitor(v),
    state(s)
  {}
};


template<typename TraversalDataInterface>
class TraversalManager : public CBase_TraversalManager<TraversalDataInterface> {
  public:
  typedef typename TraversalDataInterface::TreeBuildDataInterfaceType::LocalNodeType LocalNodeType;
  typedef typename LocalNodeType::KeyType KeyType;
  typedef typename TraversalDataInterface::RemoteNodeType RemoteNodeType;
  typedef typename TraversalDataInterface::TreeBuildDataInterfaceType TreeBuildDataInterface;

  typedef Manager<TreeBuildDataInterface> TreeManagerType;
  typedef CProxy_Manager<TraversalDataInterface> TreeManagerProxyType;
  typedef CProxy_Piece<TreeBuildDataInterface> ClientProxyType;
  typedef CProxy_TraversalManager<TraversalDataInterface> ProxyType;

  typedef CachedLeafContents<TraversalDataInterface> CachedLeafContentsType;
  typedef CachedNodes<TraversalDataInterface> CachedNodesType;
  typedef CkCacheManager<KeyType> CacheType;
  typedef CProxy_CkCacheManager<KeyType> CacheProxyType;
  typedef NodeCacheInterface<TraversalDataInterface> NodeCacheInterfaceType;
  typedef LeafCacheInterface<TraversalDataInterface> LeafCacheInterfaceType;
  typedef typename CkCacheRequestorData<KeyType>::CkCacheCallback CallbackFnType;



  private:

  ClientProxyType clientProxy_;
  CacheProxyType nodeCacheProxy_;
  CacheProxyType leafCacheProxy_;
  CacheType *nodeCache_;
  CacheType *leafCache_;
  TreeManagerType *treeManager_;

  CProxyElement_ArrayElement thisElement_;

  // for CkCache interface
  NodeCacheInterfaceType *nodeCacheInterface_;
  LeafCacheInterfaceType *leafCacheInterface_;

  TraversalStatistics traversalStatistics_;

  protected:
  ClientProxyType &clientProxy() { return clientProxy_; }
  CacheProxyType &nodeCacheProxy() { return nodeCacheProxy_; }
  CacheProxyType &leafCacheProxy() { return leafCacheProxy_; }
  CacheType *&nodeCache() { return nodeCache_; }
  CacheType *&leafCache() { return leafCache_; }
  TreeManagerType *&treeManager() { return treeManager_; }
  NodeCacheInterfaceType *nodeCacheInterface() { return nodeCacheInterface_; }
  LeafCacheInterfaceType *leafCacheInterface() { return leafCacheInterface_; }
  CProxyElement_ArrayElement &thisElement() { return thisElement_; }

  TraversalStatistics &traversalStatistics() { return traversalStatistics_; }

  public:
  TraversalManager() : 
    nodeCache_(NULL),
    leafCache_(NULL),
    thisElement_(this->thisProxy[this->thisIndex])
  {
    this->setMigratable(false);
  }

  TraversalManager(bool writebackNodes, bool writebackLeaves) : 
    nodeCache_(NULL),
    leafCache_(NULL),
    treeManager_(NULL),
    thisElement_(this->thisProxy[this->thisIndex])
  {
    this->setMigratable(false);

    if(writebackNodes){
      CkAbort("TraversalManager::TraversalManager writebackNodes not currently supported\n");
    }

    /*
    if(writebackNodes){
      nodeCacheInterface_ = new WritebackNodeCacheInterface<TraversalDataInterface>;
    }
    else{
      nodeCacheInterface_ = new NodeCacheInterface<TraversalDataInterface>;
    }
    */
    nodeCacheInterface_ = new NodeCacheInterface<TraversalDataInterface>;

    if(writebackLeaves){
      leafCacheInterface_ = new WritebackLeafCacheInterface<TraversalDataInterface>;
    }
    else{
      leafCacheInterface_ = new LeafCacheInterface<TraversalDataInterface>;
    }
  }

  TraversalManager(CkMigrateMessage *) {}

  void ckAboutToMigrate(){
    CkPrintf("[%d] TraversalManager migrating!\n", this->thisIndex);
    CkAbort("TraversalManagers shouldn't migrate\n");
  }

  // for statistics
  void nodeCacheMiss(){
    traversalStatistics().nodeCacheMiss();
  }

  void leafCacheMiss(){
    traversalStatistics().leafCacheMiss();
  }

};

template<typename TraversalDataInterface>
class LeafNodeTopDownTraversalManager : public TraversalManager<TraversalDataInterface> {
  public:
  typedef typename TraversalDataInterface::TreeBuildDataInterfaceType::LocalNodeType LocalNodeType;
  typedef typename LocalNodeType::KeyType KeyType;
  typedef typename TraversalDataInterface::RemoteNodeType RemoteNodeType;
  typedef typename TraversalDataInterface::TreeBuildDataInterfaceType TreeBuildDataInterface;

  typedef Manager<TreeBuildDataInterface> TreeManagerType;
  typedef CProxy_Manager<TreeBuildDataInterface> TreeManagerProxyType;
  typedef CProxy_Piece<TreeBuildDataInterface> ClientProxyType;
  typedef TraversalHandle<LeafNodeTopDownTraversalManager<TraversalDataInterface> > TraversalHandleType;

  typedef CachedLeafContents<TraversalDataInterface> CachedLeafContentsType;
  typedef CachedNodes<TraversalDataInterface> CachedNodesType;
  typedef CkCacheManager<KeyType> CacheType;
  typedef CProxy_CkCacheManager<KeyType> CacheProxyType;
  typedef CProxy_LeafNodeTopDownTraversalManager<TraversalDataInterface> ProxyType;
  typedef NodeCacheInterface<TraversalDataInterface> NodeCacheInterfaceType;
  typedef LeafCacheInterface<TraversalDataInterface> LeafCacheInterfaceType;
  typedef typename CkCacheRequestorData<KeyType>::CkCacheCallback CallbackFnType;

  protected: 

  /*
  template<typename VisitorType>
  class MapReduceCallbacks {
    public:
    CallbackFnType node() const {
      return &LeafNodeTopDownTraversalManager::template NodeCallback<VisitorType>;
    }

    CallbackFnType leaf() const {
      return &LeafNodeTopDownTraversalManager::template LeafCallback<VisitorType>;
    }
  };
  */

  template<typename VisitorType>
  class IteratorCallbacks {
    public:
    CallbackFnType node() const {
      return &LeafNodeTopDownTraversalManager::template IteratorIfaceNodeCallback<VisitorType>;
    }

    CallbackFnType leaf() const {
      return &LeafNodeTopDownTraversalManager::template IteratorIfaceLeafCallback<VisitorType>;
    }
  };


  private:
  CkCallback callback_;


  public:
  LeafNodeTopDownTraversalManager() {}
  LeafNodeTopDownTraversalManager(bool doNodeWriteback, bool doLeafWriteback) : 
    TraversalManager<TraversalDataInterface>(doNodeWriteback, doLeafWriteback)
  {}

  LeafNodeTopDownTraversalManager(CkMigrateMessage *) {}

  /*
  template<typename TreeHandleType, typename VisitorType>
  void map(const TreeHandleType &handle, const VisitorType &visitor, const CkCallback &cb);
  */

  template<typename VisitorType>
  void finished(VisitorType &v);


  // callback after data is received.
  template<typename VisitorType>
  static void NodeCallback(CkArrayID, CkArrayIndex&, KeyType, CkCacheUserData &, void *data, int chunk);

  template<typename VisitorType>
  static void LeafCallback(CkArrayID, CkArrayIndex&, KeyType, CkCacheUserData &, void *data, int chunk);

  // iterator-like interface
  template<typename TreeHandleType>
  void initialize(const TreeHandleType &handle, const CacheProxyType &nodeCacheProxy, const CacheProxyType &leafCacheProxy, const CkCallback &cb);

  void synchronize(const CkCallback &cb);
  void done(const CkCallback &cb);
  void statistics(const CkCallback &cb);


  // entry methods, used internally
  template<typename VisitorType, typename CallbacksType>
  const RemoteNodeType *fetchRemoteNode(const RemoteNodeType *node, VisitorType *context, const CallbacksType &callbacks);

  // a local leaf node
  // can never give us cached contents
  template<typename VisitorType, typename CallbacksType>
  CachedLeafContentsType *fetchRemoteLeaf(const RemoteNodeType *node, VisitorType *context, const CallbacksType &callbacks);


  // when data is received by the iterator interface, we have 
  // different behavior from the map-reduce interface
  template<typename VisitorType>
  static void IteratorIfaceNodeCallback(CkArrayID, CkArrayIndex&, KeyType, CkCacheUserData &, void *data, int chunk);

  template<typename VisitorType>
  static void IteratorIfaceLeafCallback(CkArrayID, CkArrayIndex&, KeyType, CkCacheUserData &, void *data, int chunk);

  // this can be called by user code; we will supply an object of type
  // IteratorCallbacks as the 'callbacks' argument to it
  // the version of 'go' used by 'map' supplies a MapReduceCallbacks arg to 
  // the method.
  template<typename AnyNodeType, typename VisitorType>
  void go(const AnyNodeType *node, VisitorType *v);

  protected:
  template<typename AnyNodeType, typename VisitorType, typename CallbacksType>
  void go(const AnyNodeType *node, VisitorType *v, const CallbacksType &callbacks);

  private:
  void innerSynchronize();

  void innerDone(int chunk);

  public:
  // factory method
  static TraversalHandleType instantiate(const CacheAccessType &nodeCacheAccessType, const CacheAccessType &leafCacheAccessType){
    TraversalHandleType handle;

    bool doNodeCacheWriteback = (nodeCacheAccessType == CacheAccessType::Accumulate);
    CkAssert(!doNodeCacheWriteback);

    bool doLeafCacheWriteback = (leafCacheAccessType == CacheAccessType::Accumulate);

    handle.proxy() = ProxyType::ckNew(doNodeCacheWriteback, doLeafCacheWriteback, CkNumPes());

    handle.nodeCacheAccessType() = nodeCacheAccessType; 
    handle.leafCacheAccessType() = leafCacheAccessType; 

    // cache proxies are set by user later on
    return handle;
  }

  // to get a callback that can be used by my
  // traversal handle to initiate done() calls
  // once completion detection library has been 
  // set up. we use this function to keep from 
  // exposing to the handle the CkIndex_... type.
  static CkCallback getCallDoneCallback(const ProxyType &proxy){
    return CkCallback(CkIndex_LeafNodeTopDownTraversalManager<TraversalDataInterface>::callDone(), proxy);
  }
  
  void callDone();

};

template<typename TraversalDataInterface>
class LeafNodeBottomUpTraversalManager : public LeafNodeTopDownTraversalManager<TraversalDataInterface> {
  public:
  typedef typename TraversalDataInterface::TreeBuildDataInterfaceType::LocalNodeType LocalNodeType;
  typedef typename LocalNodeType::KeyType KeyType;
  typedef typename TraversalDataInterface::RemoteNodeType RemoteNodeType;
  typedef typename TraversalDataInterface::TreeBuildDataInterfaceType TreeBuildDataInterface;

  typedef Manager<TreeBuildDataInterface> TreeManagerType;
  typedef CProxy_Manager<TreeBuildDataInterface> TreeManagerProxyType;
  typedef CProxy_Piece<TreeBuildDataInterface> ClientProxyType;
  typedef CProxy_LeafNodeBottomUpTraversalManager<TraversalDataInterface> ProxyType;
  typedef TraversalHandle<LeafNodeBottomUpTraversalManager<TraversalDataInterface> > TraversalHandleType;

  typedef CachedLeafContents<TraversalDataInterface> CachedLeafContentsType;
  typedef CachedNodes<TraversalDataInterface> CachedNodesType;
  typedef CkCacheManager<KeyType> CacheType;
  typedef CProxy_CkCacheManager<KeyType> CacheProxyType;
  typedef NodeCacheInterface<TraversalDataInterface> NodeCacheInterfaceType;
  typedef LeafCacheInterface<TraversalDataInterface> LeafCacheInterfaceType;
  typedef typename CkCacheRequestorData<KeyType>::CkCacheCallback CallbackFnType;


  public:
  LeafNodeBottomUpTraversalManager() {}
  LeafNodeBottomUpTraversalManager(bool doNodeWriteback, bool doLeafWriteback) : 
    LeafNodeTopDownTraversalManager<TraversalDataInterface>(doNodeWriteback, doLeafWriteback)
  {}

  LeafNodeBottomUpTraversalManager(CkMigrateMessage *) {}

  // this can be called by user code; we will supply an object of type
  // IteratorCallbacks as the 'callbacks' argument to it
  // the version of 'go' used by 'map' supplies a MapReduceCallbacks arg to 
  // the method.
  template<typename AnyNodeType, typename VisitorType>
  void go(const AnyNodeType *node, VisitorType *v);

  private:
  template<typename AnyNodeType, typename VisitorType, typename CallbacksType>
  void goDownOnly(const AnyNodeType *node, VisitorType *v, const CallbacksType &callbacks);

  template<typename AnyNodeType, typename VisitorType, typename CallbacksType>
  void goUpThenDown(const AnyNodeType *node, VisitorType *v, const CallbacksType &callbacks);


  // factory methods
  public:

  static TraversalHandleType instantiate(const CacheAccessType &nodeCacheAccessType, const CacheAccessType &leafCacheAccessType){
    TraversalHandleType handle;

    bool doNodeCacheWriteback = (nodeCacheAccessType == CacheAccessType::Accumulate);
    CkAssert(!doNodeCacheWriteback);

    bool doLeafCacheWriteback = (leafCacheAccessType == CacheAccessType::Accumulate);

    handle.proxy() = ProxyType::ckNew(doNodeCacheWriteback, doLeafCacheWriteback, CkNumPes());

    handle.nodeCacheAccessType() = nodeCacheAccessType; 
    handle.leafCacheAccessType() = leafCacheAccessType; 

    // cache proxies are set by user later on
    return handle;
  }

};


/*
template<typename KeyType,
         typename UserType,
         typename ParticleType,
         typename RemoteParticleType>
class LeafTraversalManager : public TraversalManager<KeyType, UserType, ParticleType, RemoteParticleType> {
  typedef Node<KeyType, UserType, ParticleType> MyNodeType;
  typedef Manager<KeyType, UserType, ParticleType, RemoteParticleType> TreeManagerType;
  typedef CProxy_Manager<KeyType, UserType, ParticleType, RemoteParticleType> TreeManagerProxyType;
  typedef CProxy_Piece<KeyType, UserType, ParticleType, RemoteParticleType> ClientProxyType;
  typedef Handle<KeyType, UserType, ParticleType, RemoteParticleType> TreeHandleType;
  typedef CkCacheManager<KeyType> CacheType;
  typedef CProxy_CkCacheManager<KeyType> CacheProxyType;

  private:
  CkCallback callback_;

  public:
  LeafTraversalManager() {}
  LeafTraversalManager(CkMigrateMessage *) {}

  template<typename VisitorType>
  void map(TreeHandleType &handle, const VisitorType &visitor, const CkCallback &cb);
};
*/


/*
// method definitions
template<typename KeyType,
         typename UserType,
         typename ParticleType,
         typename RemoteParticleType>
template<typename VisitorType>
void
LeafNodeTopDownTraversalManager<KeyType, UserType, ParticleType, RemoteParticleType>::
map(TreeHandleType &handle, const VisitorType &visitor, const CkCallback &cb){
  // save callback
  callback_ = cb;

  // make copy so as to make the visitor object persistent
  // across entry method invocations; we might yield to the other activities
  // (such as responding to remote data requests) in the middle of traversal,
  // so this has to be done
  VisitorType *v = new VisitorType(visitor);

  // disabled, because it is not used, and in order
  // to not get a compile-time error, we need to 
  // pass a traversal handle into this call
  //innrSynchronize(handle);

  CkVec<MyNodeType *> &leaves = this->treeManager()->leaves();
  MyNodeType *root = this->treeManager()->root();

  // start traversal
  v->outstanding() = leaves.size(); 
  MapReduceCallbacks<VisitorType> mapCallbacks;
  if(leaves.size() > 0){
    for(int i = 0; i < leaves.size(); i++){
      v->leaf() = leaves[i];
      go(root, v, mapCallbacks);
      v->outstanding()--;
    }
  }

  this->savedVisitor() = (void *) v;

  if(v->outstanding() == 0){
    CProxyElement_LeafNodeTopDownTraversalManager<KeyType, UserType, ParticleType, RemoteParticleType> elem(this->thisProxy.ckGetArrayID(), CkArrayIndex1D(this->thisIndex));
    elem.finished(*v);
  }

}
*/

template<typename TraversalDataInterface>
void
LeafNodeTopDownTraversalManager<TraversalDataInterface>::
innerSynchronize(){
  CkArrayIndex idx = CkArrayIndex1D(this->thisIndex);
  int localIdx = -1;
  int numChunks = 1;

  // must cacheSync on every iteration

  //CkPrintf("[%d] TopDownTraversalManager::cacheSync\n", this->thisIndex);
  this->nodeCache()->cacheSync(numChunks, idx, localIdx);
  this->leafCache()->cacheSync(numChunks, idx, localIdx);

  CkAssert(localIdx >= 0);
}

template<typename TraversalDataInterface>
template<typename AnyNodeType, typename VisitorType>
void
LeafNodeTopDownTraversalManager<TraversalDataInterface>::
go(const AnyNodeType *n, VisitorType *v){
  go(n, v, IteratorCallbacks<VisitorType>());
}

template<typename TraversalDataInterface>
template<typename AnyNodeType, typename VisitorType>
void
LeafNodeBottomUpTraversalManager<TraversalDataInterface>::
go(const AnyNodeType *n, VisitorType *v){
  if(AnyNodeType::IsCached(n->getType())){
    // user wants to start a traversal at a cached
    // node. in the bottom-up traversal, this can only
    // happen when we are resuming a previously halted 
    // traversal after fetching the required data; since 
    // we never travel upwards after descending into an 
    // aunt's subtree, we must continue descending into
    // this subtree.
    goDownOnly(n, v, typename LeafNodeTopDownTraversalManager<TraversalDataInterface>::template IteratorCallbacks<VisitorType>());
  }
  else{
    goUpThenDown(n, v, typename LeafNodeTopDownTraversalManager<TraversalDataInterface>::template IteratorCallbacks<VisitorType>());
  }
}

template<typename TraversalDataInterface>
template<typename AnyNodeType, typename VisitorType, typename CallbacksType>
void
LeafNodeBottomUpTraversalManager<TraversalDataInterface>::
goDownOnly(const AnyNodeType *n, VisitorType *v, const CallbacksType &callbacks){
  // use the TopDownTraversalManager's recursive descent method
  LeafNodeTopDownTraversalManager<TraversalDataInterface>::go(n, v, callbacks);
}
 
template<typename TraversalDataInterface>
template<typename AnyNodeType, typename VisitorType, typename CallbacksType>
void
LeafNodeBottomUpTraversalManager<TraversalDataInterface>::
goUpThenDown(const AnyNodeType *n, VisitorType *v, const CallbacksType &callbacks){
  // have we reached the root?
  const AnyNodeType *parent = n->getParent();
  if(parent == NULL) return;
  
  // not yet reached the root. recursively descend into aunts' subtrees
  CkAssert(parent->getChildren() != NULL);
  for(int i = 0; i < parent->getNumChildren(); i++){
    const AnyNodeType *child = parent->getChild(i);
    // don't descend into own subtree; that's where we came up from
    if(n == child) continue;

    goDownOnly(child, v, callbacks);
  }

  goUpThenDown(parent, v, callbacks);
}
 

template<typename TraversalDataInterface>
template<typename AnyNodeType, typename VisitorType, typename CallbacksType>
void
LeafNodeTopDownTraversalManager<TraversalDataInterface>::
go(const AnyNodeType *n, VisitorType *v, const CallbacksType &callbacks){
  if(!v->node(n)) return;

  const AnyNodeType *children = NULL;

  NodeType type = n->getType(); 
  // take action based on type of node
  if(AnyNodeType::IsLocalLeaf(n->getType())){
    v->localLeaf(n->getKey(), n->getParticles(), n->getNumParticles());
  }
  else if(AnyNodeType::IsRemoteLeaf(n->getType()) ||
          AnyNodeType::IsCachedLeaf(n->getType())){
    // remote leaf
    const RemoteNodeType *remoteN = (const RemoteNodeType *) n;
    const CachedLeafContentsType *cached = fetchRemoteLeaf(remoteN, v, callbacks);
    if(cached != NULL){
      v->remoteLeaf(remoteN->getKey(), &cached->particles[0], cached->nParticles);
    }
    else{
      v->miss(remoteN->getKey());
      //CkPrintf("[%d] TreePiece::Gravity leaf %llu MISS LEAF %llu outstanding %d\n", thisIndex, leaf->getKey(), n->getKey(), outstanding_);
    }
  }
  else if(AnyNodeType::IsBoundary(n->getType()) || 
          AnyNodeType::IsInternal(n->getType())){
    const LocalNodeType *localN = (const LocalNodeType *) n;
    CkAssert(localN->getNumChildren() > 0);

    const LocalNodeType *children = localN->getChildren();
    CkAssert(children != NULL);
    for(int i = 0; i < localN->getNumChildren(); i++){
      go(&children[i], v, callbacks);
    }
  }
  else{
    const RemoteNodeType *remoteN = (const RemoteNodeType *) n;
    if(AnyNodeType::IsRemote(remoteN->getType())){
      CkAssert(remoteN->getChildren() == NULL);
    }
    else{
      CkAssert(AnyNodeType::IsCached(remoteN->getType()));
    }

    const RemoteNodeType *children = fetchRemoteNode(remoteN, v, callbacks);
    if(children != NULL){
      for(int i = 0; i < n->getNumChildren(); i++){
        go(&children[i], v, callbacks);
      }
    }
    else{
      v->miss(n->getKey());
      //CkPrintf("[%d] TreePiece::Gravity leaf %llu MISS NODE %llu outstanding %d\n", thisIndex, leaf->getKey(), n->getKey(), outstanding_);
    }
  }
}

/*
template<typename KeyType,
template<typename VisitorType>
void
LeafNodeTopDownTraversalManager<KeyType, UserType, ParticleType, RemoteParticleType>::
finished(VisitorType &v){
  int chunk = 0;

  innerDone(chunk);

  v.done();

  this->contribute(sizeof(VisitorType), &v, Reduce<VisitorType>::reducer(), callback_);

  VisitorType *newdVisitor = (VisitorType *) this->savedVisitor();
  CkAssert(newdVisitor != NULL);
  delete newdVisitor;
  this->savedVisitor() = NULL;
}
*/

template<typename TraversalDataInterface>
void
LeafNodeTopDownTraversalManager<TraversalDataInterface>::
innerDone(int chunk){
  //CkPrintf("[%d] TopDownTraversalManager::innerDone finishedChunk\n", this->thisIndex);
  this->nodeCacheInterface()->done(chunk);
  this->leafCacheInterface()->done(chunk);
}



/*
template<typename KeyType,
         typename UserType,
         typename ParticleType,
         typename RemoteParticleType>
template<typename VisitorType>
void
LeafNodeTopDownTraversalManager<KeyType, UserType, ParticleType, RemoteParticleType>::
NodeCallback(CkArrayID requestorArray, CkArrayIndex &requestorIndex, KeyType key, CkCacheUserData &context, void *data, int chunk){ 
  CachedNodesType *c = (CachedNodesType *) data;
  MyNodeType *node = &c->nodes[0];

  Resumption<VisitorType, typename VisitorType::State> *r;
  r = (Resumption<VisitorType, typename VisitorType::State> *) context.d0;

  CkArrayIndex1D idx(requestorIndex.data()[0]);
  CProxyElement_LeafNodeTopDownTraversalManager<KeyType, UserType, ParticleType, RemoteParticleType> elem(requestorArray, idx);
  LeafNodeTopDownTraversalManager *traversal = elem.ckLocal();

  // save current state of visitor, since client is 
  // allowed to reuse the visitor object
  const typename VisitorType::State &saved = r->visitor->state();
  
  // first recover state of suspended traversal
  r->visitor->state(r->state);
  r->visitor->outstanding()--;
  r->visitor->hit(key);
  // then, start with the children
  CkAssert(node->getChildren() != NULL);
  MapReduceCallbacks<VisitorType> mapCallbacks;
  for(int i = 0; i < node->getNumChildren(); i++){
    traversal->go(node->getChild(i), r->visitor, mapCallbacks);
  }

  CkAssert(key == node->getKey());

  if(r->visitor->outstanding() ==0){
    elem.finished(*r->visitor);
  }


  // restore state to what it was before this method 
  // was invoked
  r->visitor->state(saved);

  delete r;
}

template<typename KeyType,
         typename UserType,
         typename ParticleType,
         typename RemoteParticleType>
template <typename VisitorType>
void
LeafNodeTopDownTraversalManager<KeyType, UserType, ParticleType, RemoteParticleType>::
LeafCallback(CkArrayID requestorArray, CkArrayIndex &requestorIndex, KeyType key, CkCacheUserData &context, void *data, int chunk){ 
  Resumption<VisitorType, typename VisitorType::State> *r = (Resumption<VisitorType, typename VisitorType::State> *) context.d0;
  CachedLeafContentsType *c = (CachedLeafContentsType *) data;

  CkArrayIndex1D idx(requestorIndex.data()[0]);
  CProxyElement_LeafNodeTopDownTraversalManager<KeyType, UserType, ParticleType, RemoteParticleType> elem(requestorArray, idx);
  LeafNodeTopDownTraversalManager *traversal = elem.ckLocal();

  // save current state of traversal
  const typename VisitorType::State &saved = r->visitor->state();
  
  // recover state of suspended traversal
  r->visitor->state(r->state);
  r->visitor->outstanding()--;
  r->visitor->hit(key);
  r->visitor->remoteLeaf(key, c->particles, c->nParticles);

  if(r->visitor->outstanding() ==0){
    elem.finished(*r->visitor);
  }

  // restore state of visitor to what it was
  // before this method was invoked
  r->visitor->state(saved);

  delete r;
}
*/

template<typename TraversalDataInterface>
template<typename TreeHandleType>
void
LeafNodeTopDownTraversalManager<TraversalDataInterface>::
initialize(const TreeHandleType &treeHandle, const CacheProxyType &ncp, const CacheProxyType &lcp, const CkCallback &cb){
  // so that the traversal can query the cache for
  // required nodes when it is doing a fetchNext()
  this->nodeCacheProxy() = ncp;
  this->leafCacheProxy() = lcp;

  this->nodeCache() = this->nodeCacheProxy().ckLocalBranch();
  this->leafCache() = this->leafCacheProxy().ckLocalBranch();
  // XXX this is used only in the map method, which
  // we plan to discard 
  this->treeManager() = treeHandle.managerProxy().ckLocalBranch();

  // so that the cache interface can do stats 
  // collection
  this->nodeCacheInterface()->setTraversalManager(this);
  this->leafCacheInterface()->setTraversalManager(this);

  // so that the interface's unpack() method can use
  // the functionality provided by the tree manager
  this->nodeCacheInterface()->setTreeManager(this->treeManager());
  this->leafCacheInterface()->setTreeManager(this->treeManager());

  // give interface objects access to completion detection library
  this->nodeCacheInterface()->setCompletionDetector(treeHandle.completionDetectorProxy());
  this->leafCacheInterface()->setCompletionDetector(treeHandle.completionDetectorProxy());

  // so that the interface can send a message to the
  // appropriate chare array member when it wants data
  this->nodeCacheInterface()->setClientProxy(treeHandle.pieceProxy());
  this->leafCacheInterface()->setClientProxy(treeHandle.pieceProxy());

  // when the interface requests data from a tree piece
  // that owns it, it must tell the tree piece which
  // cache manager it should send the reply to.
  this->nodeCacheInterface()->setCacheProxy(this->nodeCacheProxy());
  this->leafCacheInterface()->setCacheProxy(this->leafCacheProxy());

  this->contribute(cb);
}


template<typename TraversalDataInterface>
void
LeafNodeTopDownTraversalManager<TraversalDataInterface>::
synchronize(const CkCallback &cb){
  innerSynchronize();
  this->traversalStatistics().reset();
  this->contribute(cb);
}

template<typename TraversalDataInterface>
void
LeafNodeTopDownTraversalManager<TraversalDataInterface>::
done(const CkCallback &cb){
  int chunk = 0;
  innerDone(chunk);
  this->contribute(cb);
}

template<typename TraversalDataInterface>
void
LeafNodeTopDownTraversalManager<TraversalDataInterface>::
callDone(){
  done(CkCallback(CkCallback::ignore));
}

template<typename TraversalDataInterface>
void
LeafNodeTopDownTraversalManager<TraversalDataInterface>::
statistics(const CkCallback &cb){
  this->contribute(sizeof(TraversalStatistics), &this->traversalStatistics(), Reduce<TraversalStatistics>::reducer(), cb);
}




template<typename TraversalDataInterface>
template <typename VisitorType, typename CallbacksType>
const typename TraversalDataInterface::RemoteNodeType *
LeafNodeTopDownTraversalManager<TraversalDataInterface>::
fetchRemoteNode(const RemoteNodeType *node, VisitorType *v, const CallbacksType &callbacks){
  // the node we are given can only be a non-home node
  // (i.e. fetched from a remote source); moreover, we know that this 
  // function is only called for non-leaf nodes. 
  const RemoteNodeType *children = node->getChildren();
  // if children found, return them 
  if(children != NULL){
    // since this method is only called for Remote nodes, 
    // and Cached nodes, and Remote nodes do not have children,
    // it follows that this node must be cached
    CkAssert(RemoteNodeType::IsCached(node->getType()));
    return children;
  }
  
  // otherwise, try the cache
  KeyType k = node->getKey();
  int owner = node->home();

  CkCacheUserData user;
  // assume that the user has supplied a persistent context
  VisitorType *context = v;
  //VisitorType *context = new VisitorType(*v); 
  // make a persistent copy of the user-supplied context
  user.d0 = (CmiUInt8) context;

  CkCacheRequestorData<KeyType> req(this->thisElement(), callbacks.node(), user);

  // we assume that owner TreePieces are always one-dimensional array members 
  CkArrayIndex1D ownerIdx(owner);
  // we can use the nodeCacheIterface_ object as the cache interface, since
  // all the basic operations (request, unpack, free, size) remain the same 
  // it is just that this method gives us the flexibility of separating the 
  // requestor (which could be, say, a Compute) from the holder, which is always
  // a one dimensional Piece
  void *data = this->nodeCache()->requestData(k, ownerIdx, 0, this->nodeCacheInterface(), req);
  if(data != NULL){
    //delete context;
    CachedNodesType *cached = (CachedNodesType *) data;
    typename TraversalDataInterface::RemoteNodeType *cachedParent = &cached->nodes[0];
    // sanity check: if children pointer is non-null, node must be known to 
    // have children
    CkAssert(cachedParent->getNumChildren() > 0);
    children = cachedParent->getChildren();
  }

  return children;
}

template<typename TraversalDataInterface>
template <typename VisitorType, typename CallbacksType>
CachedLeafContents<TraversalDataInterface> *
LeafNodeTopDownTraversalManager<TraversalDataInterface>::
fetchRemoteLeaf(const RemoteNodeType *leaf, VisitorType *v, const CallbacksType &callbacks){
  KeyType k = leaf->getKey();

  CkCacheUserData user;
  VisitorType *context = v;
  //VisitorType *context = new VisitorType(*v);
  user.d0 = (CmiUInt8) context; 
  CkCacheRequestorData<KeyType> req(this->thisElement(), callbacks.leaf(), user);

  // should have only one owner
  CkAssert(leaf->getOwnerStart() + 1 == leaf->getOwnerEnd());
  // can't assert this anymore. we could have received this leaf
  // as part of a cached subtree that was fetched from some remote
  // source. if it was a shared node on the PE of its home tree piece
  // then it would have received a home of '-1' there. while unpacking 
  // the received subtree on this PE, we would have set all -1 homes to
  // the index of the tree piece that we believe to be the home of the
  // ancestor for which the request was originally made. but this is not
  // to bad an error, because even if we do ask the tree piece that is the home
  // but not the actual owner of the node, it will always forward the request
  // to its manager, which will have the node in question (since it was 
  // considered a shared node on that PE). We can reason about this as follows:
  // during unpacking, we set the homes (perhaps incorrectly) of only those nodes 
  // that have a -1 home. By construction, each such node be considered shared
  // on the PE that hosts its home, as well as its owner. If we send a request for
  // this node to any of those tree pieces (i.e. the home or the owner), it will 
  // be forwarded to the same manager, and we will get the node regardless.
  //CkAssert(leaf->home() == leaf->getOwnerStart());
  int owner = leaf->home();
  CkArrayIndex1D ownerIdx(owner);

  // assumes that when we want the particles
  // in a leaf we request the leaf itself.
  void *data = this->leafCache()->requestData(k, ownerIdx, 0, this->leafCacheInterface(), req);

  if(data != NULL){
    //delete context;
    CachedLeafContentsType *cached = (CachedLeafContentsType *) data;
    //CkAssert(leaf->getNumParticles() == cached->nParticles);

    return cached;
  }

  return NULL;
}

template<typename TraversalDataInterface>
template <typename VisitorType>
void
LeafNodeTopDownTraversalManager<TraversalDataInterface>::
IteratorIfaceNodeCallback(CkArrayID requestorArray, CkArrayIndex &requestorIndex, KeyType key, CkCacheUserData &context, void *data, int chunk){ 
  VisitorType *v = (VisitorType *) context.d0;
  CachedNodesType *cached = (CachedNodesType *) data;

  // let the context provided by the user translate the
  // saved requestor array and index into a proper array element
  v->deliver(&cached->nodes[0], chunk);
  // don't del the context anymore, since we assume that the
  // user passed in persistent memory that she herself manages.

  // since we allocated this memory (using the copy constructor), we will
  // delete it
  //delete v;
}

template<typename TraversalDataInterface>
template <typename VisitorType>
void
LeafNodeTopDownTraversalManager<TraversalDataInterface>::
IteratorIfaceLeafCallback(CkArrayID requestorArray, CkArrayIndex &requestorIndex, KeyType key, CkCacheUserData &context, void *data, int chunk){ 
  VisitorType *v = (VisitorType *) context.d0;
  CachedLeafContentsType *cached = (CachedLeafContentsType *) data;

  // let the context provided by the user translate the
  // saved requestor array and index into a proper array element
  v->deliver(key, &cached->particles[0], cached->nParticles, chunk);

  // ditto 
  // since we allocated this memory (using the copy constructor), we will
  // delete it
  //delete v;
}

/*
template<typename KeyType,
         typename UserType,
         typename ParticleType,
         typename RemoteParticleType>
template<typename VisitorType>
void
LeafTraversalManager<KeyType, UserType, ParticleType, RemoteParticleType>::
map(TreeHandleType &handle, const VisitorType &visitor, const CkCallback &cb){
  // save callback
  callback_ = cb;
  // so that we don't need to worry about constness
  VisitorType *v = new VisitorType(visitor);

  this->treeManager() = handle.managerProxy().ckLocalBranch();
  CkVec<MyNodeType *> &leaves = this->treeManager()->leaves();

  // all we do in this traversal is iterate over the leaves
  for(int i = 0; i < leaves.size(); i++){
    v->localLeaf(leaves[i]);
  }

  v->done();
  this->contribute(sizeof(VisitorType), v, Reduce<VisitorType>::reducer(), callback_);

  delete v;
}
*/

}; // namespace MultiphaseDistree


#endif // MULTIPHASE_DISTREE_TRAVERSAL_MANAGER_H
