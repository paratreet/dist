#ifndef MULTIPHASE_DISTREE_PIECE_H
#define MULTIPHASE_DISTREE_PIECE_H

#include "mdt.decl.h"
#include "Manager.h"
#include "CkCache.h"
#include <map>

namespace MultiphaseDistree {

template<typename DataInterface>
class Piece : public CBase_Piece<DataInterface> {
  typedef typename DataInterface::LocalNodeType LocalNodeType;
  typedef typename DataInterface::LocalParticleType LocalParticleType;
  typedef typename LocalNodeType::KeyType KeyType;
  typedef Manager<DataInterface> ManagerType;
  typedef CProxy_Manager<DataInterface> ManagerProxyType;
  typedef Handle<DataInterface> HandleType;

  ManagerProxyType treeManagerProxy_;
  ManagerType *treeManager_;
  LocalNodeType *root_;
  // table containing pointers to all of the nodes
  // owned by this piece. the idea is that instead
  // of forwarding every request to the tree manager,
  // we have each piece maintain key-indexed pointers
  // to the nodes that it owns. only for nodes shared
  // by some pieces on this node do we look up the
  // node manager's table. 
  std::map<KeyType, LocalNodeType *> ownedNodes_;

  CProxy_CompletionDetector completionDetectorProxy_;
  bool doneRegistration_;

  public:
  Piece();
  Piece(CkMigrateMessage *msg);

  int getIndex();
  void initialize(const HandleType &handle, const CkCallback &cb);

  void registration(const PointerContainer<LocalNodeType> &n);
  // used during tree construction
  void requestPayload(KeyType key, int replyTo);
  // invoked during traversal
  template<typename TraversalDataInterface>
  void requestNode(KeyType key, const CProxy_CkCacheManager<KeyType> &cacheProxy, int replyTo);
  template<typename TraversalDataInterface>
  void requestLeafContents(KeyType key, const CProxy_CkCacheManager<KeyType> &cacheProxy, int replyTo);

  // invoked when data is being written back
  // NOT SUPPORTED CURRENTLY
  /*
  template<typename TraversalDataInterface>
  void accumulateNodes(CkCacheFillMsg<KeyType> *msg);
  */
  template<typename TraversalDataInterface>
  void accumulateLeafContents(CkCacheFillMsg<KeyType> *msg);

  void deleteLocalTree(const PointerContainer<LocalNodeType> &r);


  void pup(PUP::er &p){
    CBase_Piece<DataInterface>::pup(p);
    p | doneRegistration_;
    CkAssert(!doneRegistration_);
    p | treeManagerProxy_;
    if(p.isUnpacking()){
      treeManager_ = treeManagerProxy_.ckLocalBranch();
      root_ = NULL;
      CkAssert(ownedNodes_.empty());
    }
    p | completionDetectorProxy_;
  }

  public:
  // to be defined in user code, by subclass of Piece
  virtual void addLeaf(LocalNodeType *leaf) {}

  private:
  void searchForLeaves(LocalNodeType *n);
  ManagerType *getTreeManager();


  // debugging
  public:

  void quiescence(){
    int n = doneRegistration_;
    this->contribute(sizeof(int), &n, CkReduction::sum_int, CkCallback(CkIndex_Manager<DataInterface>::quiescence2(NULL), 0, treeManagerProxy_));
  }
};


// METHOD DEFINITIONS
template<typename DataInterface> 
Piece<DataInterface>::
Piece() : 
  root_(NULL),
  doneRegistration_(false)
{}

template<typename DataInterface> 
Piece<DataInterface>::
Piece(CkMigrateMessage *m) :
  root_(NULL),
  doneRegistration_(false)
{}

template<typename DataInterface> 
int Piece<DataInterface>::
getIndex(){
  return this->thisIndex;
}

template<typename DataInterface> 
void Piece<DataInterface>::
initialize(const HandleType &handle, const CkCallback &cb){
  treeManagerProxy_ = handle.managerProxy();
  treeManager_ = treeManagerProxy_.ckLocalBranch();
  completionDetectorProxy_ = handle.completionDetectorProxy();
  this->contribute(cb);
}

template<typename DataInterface> 
void Piece<DataInterface>::
registration(const PointerContainer<LocalNodeType> &n){
  // register all the nodes in the subtree rooted
  // at 'n' as belonging to self.
  // if a request for such a registered node is 
  // received during traversal, the tree piece itself
  // serves the request. otherwise, it forwards it
  // to the manager nodegroup 
  CkAssert(!doneRegistration_);
  root_ = n.data;
  searchForLeaves(root_);
  doneRegistration_ = true;
  treeManagerProxy_[CkMyNode()].doneRegistration(n);
}


template<typename DataInterface>
void Piece<DataInterface>::searchForLeaves(LocalNodeType *node){
  CkAssert(LocalNodeType::IsLocal(node->getType()));
  ownedNodes_[node->getKey()] = node;
  //CkPrintf("[%d] Piece::registration key %llu\n", this->thisIndex, node->getKey());
  if(LocalNodeType::IsLeaf(node->getType())){
    CkAssert(!LocalNodeType::IsRemote(node->getType()));
    addLeaf(node);
  }
  else{
    CkAssert(node->getChildren() != NULL);
    for(int i = 0; i < node->getNumChildren(); i++){
      searchForLeaves(node->getChild(i));
    }
  }
}

template<typename DataInterface> 
void Piece<DataInterface>::
requestPayload(KeyType key, int replyTo){
  // can't do this for nodegroup
  //getTreeManager()->requestPayload(key, replyTo);
  // have to forward the request through a message instead.
  // note that we could have looked to check whether the node
  // is the root of this piece, and if so, serve the request
  // here itself. however, then we have to do more book-keeping:
  // each tree will have to maintain requests for its root node
  // until it has received its root node from the manager. not
  // worth the effort.
  treeManagerProxy_[CkMyNode()].requestPayload(key, replyTo);
}

template<typename DataInterface> 
template<typename TraversalDataInterface> 
void Piece<DataInterface>::
requestNode(KeyType key, const CProxy_CkCacheManager<KeyType> &cacheProxy, int replyTo){
  typename std::map<KeyType, LocalNodeType *>::iterator it;
  it = ownedNodes_.find(key);
  if(it == ownedNodes_.end()){
    // since this is a read-only operation, we can pass on the 
    // request message to the manager through a local method call. 
    // can't do the same for requestPayload above, since that method
    // modifies data structures internal to the tree manager.
    getTreeManager()->template requestNode<TraversalDataInterface>(key, cacheProxy, replyTo);
  }
  else{
    LocalNodeType *node = it->second;
    getTreeManager()->template packAndSendNode<TraversalDataInterface>(node, cacheProxy, replyTo);
  }
}

template<typename DataInterface> 
template<typename TraversalDataInterface> 
void Piece<DataInterface>::
requestLeafContents(KeyType key, const CProxy_CkCacheManager<KeyType> &cacheProxy, int replyTo){
  typename std::map<KeyType, LocalNodeType *>::iterator it;
  it = ownedNodes_.find(key);
  CkAssert(it != ownedNodes_.end());
  LocalNodeType *leaf = it->second; 
  getTreeManager()->template packAndSendLeaf<TraversalDataInterface>(leaf, cacheProxy, replyTo);
  //getTreeManager()->template requestLeafContents<TraversalDataInterface>(key, cacheProxy, replyTo);
}

/*
template<typename DataInterface>
template<typename TraversalDataInterface>
void Piece<DataInterface>::
accumulateNodes(CkCacheFillMsg<KeyType> *msg){
  CacheNodes<TraversalDataInterface> *c = (CachedNodes<TraversalDataInterface> *) msg->data;
  typename TraversalDataInterface::RemoteNodeType *sourceRoot = &c->nodes[i];
  std::map<Key, LocalNodeType *>::iterator it;
  it = ownedNodes_.find(sourceRoot->getKey());
  if(it != ownedNodes_.end()){
    // this node is internal to the tree piece, and can be
    // updated here itself; no need to forward to tree manager
    for(int i = 0; i < c->nNodes; i++){
      const typename TraversalDataInterface::RemoteNodeType *source = &c->nodes[i];
      LocalNodeType *target = ownedNodes_[source->getKey()];
      CkAssert(target != NULL);
      TraversalDataInterface::accumulate(*target, *source);
    }

    getCompletionDetector()->consume();
    delete msg;
  }
  else{
    treeManagerProxy_[CkMyNode()].accumulateNodes<TraversalDataInterface>(msg); 
  }
}
*/

template<typename DataInterface>
template<typename TraversalDataInterface>
void Piece<DataInterface>::
accumulateLeafContents(CkCacheFillMsg<KeyType> *msg){
  //getTreeManager()->template accumulateLeafContents<TraversalDataInterface>(msg); 
  CachedLeafContents<TraversalDataInterface> *c = (CachedLeafContents<TraversalDataInterface> *) msg->data;
  LocalNodeType *leaf = ownedNodes_[msg->key];
  CkAssert(leaf != NULL);
  LocalParticleType *targets = leaf->getParticles(); 
  for(int i = 0; i < c->nParticles; i++){
    const typename TraversalDataInterface::RemoteParticleType *source = &c->particles[i];
    TraversalDataInterface::accumulate(targets[i], *source);
  }

  completionDetectorProxy_.ckLocalBranch()->consume();

  delete msg;
}

template<typename DataInterface>
void Piece<DataInterface>::deleteLocalTree(const PointerContainer<LocalNodeType> &r){
  CkAssert(root_ == r.data);
  //CkPrintf("[%d] Piece::deleteLocalTree\n", getIndex());
  root_->deleteBeneath();
  ownedNodes_.clear();
  doneRegistration_ = false;
  treeManagerProxy_[CkMyNode()].doneDelete();
}

template<typename DataInterface>
Manager<DataInterface> *Piece<DataInterface>::getTreeManager(){
  return treeManagerProxy_.ckLocalBranch(); 
}


/*
template<DataInterface> 
void Piece<DataInterface>::receivePayload(KeyType key, const UserType &Payload){
  getTreeManager()->receivePayload(key, Payload);
}
*/

}; // namespace MultiphaseDistree

#endif // MULTIPHASE_DISTREE_PIECE_H
