#ifndef MULTIPHASE_DISTREE_MANAGER_H
#define MULTIPHASE_DISTREE_MANAGER_H

#include <map>
#include <sstream>
#include <fstream>


#include "mdt.decl.h"
#include "Node.h"
#include "CkCache.h"
#include "Utility.h"
#include "CacheInterface.h"
#include "CachedData.h"

#include <cstdarg>

extern void MDT_MANAGER_VERBOSE_FN(const char *fmt, ...);

namespace MultiphaseDistree {

template<typename TreeBuildDataInterfaceType> 
class Manager : public CBase_Manager<TreeBuildDataInterfaceType> {

  typedef typename TreeBuildDataInterfaceType::LocalNodeType LocalNodeType;
  typedef typename TreeBuildDataInterfaceType::LocalParticleType LocalParticleType;
  typedef typename LocalNodeType::KeyType KeyType;

  typedef Piece<TreeBuildDataInterfaceType> PieceType;
  typedef CProxy_Piece<TreeBuildDataInterfaceType> PieceProxyType;
  typedef CProxy_Manager<TreeBuildDataInterfaceType> ManagerProxyType;

  private:
  struct RegistrationEntry {
    LocalNodeType *root;
    PieceType *tp;

    RegistrationEntry() : 
      root(NULL),
      tp(NULL)
    {}

    RegistrationEntry(LocalNodeType *r, PieceType *t) : 
      root(r),
      tp(t)
    {}
  };

  private:
  CkVec<RegistrationEntry> registrationEntries_;
  CkVec<int> registrationMsgsSentTo_;
  LocalNodeType *mergedRoot_;
  PieceProxyType clientProxy_;

  // these are the requests received for nodes local 
  // to this pe during global tree construction
  std::map<KeyType, CkVec<int> > remoteRequests_;
  // pointers to nodes requested from remote
  // sources during global tree construction
  std::map<KeyType, LocalNodeType *> requestedNodes_;
  // the nodes local to this pe
  std::map<KeyType, LocalNodeType *> localNodes_;

  // saved callback
  CkCallback callback_;

  // depth of chunk to be sent when a node 
  // is requested
  int cacheChunkDepth_;

  // number of tree pieces to which we have sent
  // requests to register/delete the subtrees that they own.
  // we must wait to hear from them before we finish
  // the tree build / start
  // to delete the top-level tree
  int numTpAcks_;
  // whether we have all the moments required for the merged tree.
  bool mergeDone_;

  // lock for registration
  CmiNodeLock registrationLock_;

  public:
  Manager(int cacheChunkDepth);

  void initialize(const Handle<TreeBuildDataInterfaceType> &hdl, const CkCallback &cb);
  void registration(const PointerContainer<LocalNodeType> &root, const PointerContainer<PieceType> &tp);
  void syncToMerge(const CkCallback &cb);
  void syncToDelete(const CkCallback &cb);

  // received from tree piece, as ack that it has
  // deleted the subtree beneath its root
  void doneDelete();

  public:

  // called by a piece when it has finished registering the nodes
  // in the subtree assigned to it
  void doneRegistration(const PointerContainer<LocalNodeType> &n);

  // when the library on another pe asks user code for 
  // some data, the user code redirects the request to
  // the library on its pe via this local method call
  void requestPayload(KeyType k, int replyTo);
  void receivePayload(KeyType k, const typename TreeBuildDataInterfaceType::PayloadType &payload);

  // data interface used during traversal can be 
  // different from the one used during tree building
  template<typename TraversalDataInterface>
  void requestNode(KeyType k, const CProxy_CkCacheManager<KeyType> &cacheProxy, int replyTo);

  /*
  template<typename TraversalDataInterface>
  void requestLeafContents(KeyType k, const CProxy_CkCacheManager<KeyType> &cacheProxy, int replyTo);
  */

  // NOT SUPPORTED CURRENTLY
  /*
  template<typename TraversalDataInterface>
  void accumulateNodes(CkCacheFillMsg<KeyType> *m);
  */

  // not needed since leaves are always private to 
  // tree pieces, so manager should never receive a forwarded
  // request to accumulate to a leaf
  /*
  template<typename TraversalDataInterface>
  void accumulateLeafContents(CkCacheFillMsg<KeyType> *m);
  */



  // called by tree pieces in order to reuse functionality encapsulated
  // within tree manager
  template<typename TraversalDataInterface>
  void packAndSendNode(LocalNodeType *node, const CProxy_CkCacheManager<KeyType> &cacheProxy, int replyTo);

  template<typename TraversalDataInterface>
  void packAndSendLeaf(LocalNodeType *leaf, const CProxy_CkCacheManager<KeyType> &cacheProxy, int replyTo);


  // called by traversal manager when map() is called
  // on it
  LocalNodeType *root();

  private:
  void mergePe(CkVec<RegistrationEntry> &vec, LocalNodeType *dest);
  LocalNodeType *keyToNode(KeyType k);
  void payloadReady(LocalNodeType *node);
  void passUpwards(LocalNodeType *node);
  void finishMerge();
  void trulyFinishMerge();
  void printOutstanding();

  // actions to take when you encounter leaves
  // called from mergePe or noMergeDescendPe
  void handleLocal(LocalNodeType *node, PieceType *ownerTp);
  void searchForLeaves(LocalNodeType *node, PieceType *ownerTp);
  void handleRemote(LocalNodeType *node);
  // actions to take when you encounter a Boundary
  // node that does not have to be merged with
  // other similar nodes. 'ownerTp' is the tree piece 
  // to which the tree rooted at 'node' belongs.
  void noMergeDescend(LocalNodeType *node, PieceType *ownerTp);

  void sendDeleteRequests(LocalNodeType *node);
  void deleteTopLevelTree();

  // 'first' is the new contender for the best deal.
  // the current best deal is 'second'
  bool isBetterDeal(const LocalNodeType *first, const LocalNodeType *second){
    if(LocalNodeType::IsLocal(second->getType())){
      // only a local node with more particles under
      // it will be a better deal than a local node.
      // however, since these nodes are from two 
      // different tree pieces, they cannot both be
      // internal
      CkAssert(!LocalNodeType::IsLocal(first->getType()));
      return false;
    }
    else if(LocalNodeType::IsBoundary(second->getType())){
      // a node can't be a boundary node in one tree piece,
      // and an internal one in another.
      CkAssert(!LocalNodeType::IsLocal(first->getType()));
      if(first->getNumParticles() > second->getNumParticles()){
        return true;
      }

      return false;
    }

    // can't have cached nodes yet
    CkAssert(LocalNodeType::IsRemote(second->getType()));
    // any node (internal, boundary or remote) can displace a 
    // remote node as the new best deal
    return true;
  }

  public:
  // for packing nodes
  template<typename TraversalDataInterface>
  CkCacheFillMsg<KeyType> *pack(const LocalNodeType *root, int depth);
  int countNodesBeneath(const LocalNodeType *root, int depth);

  template<typename TraversalDataInterface>
  void recursivePack(const LocalNodeType *root, typename TraversalDataInterface::RemoteNodeType *copy, typename TraversalDataInterface::RemoteNodeType *&buf, int depth);

  // for packing leaves
  template<typename TraversalDataInterface>
  CkCacheFillMsg<KeyType> *pack(const LocalNodeType *leaf);

  // for unpacking nodes
  template<typename TraversalDataInterface>
  void *unpackNodes(CkCacheFillMsg<KeyType> *msg, int chunk, CkArrayIndex &owner);

  template<typename TraversalDataInterface>
  int recursiveUnpack(typename TraversalDataInterface::RemoteNodeType *node, int home);

  // for unpacking leaves
  template<typename TraversalDataInterface>
  void *unpackLeafContents(CkCacheFillMsg<KeyType> *msg, int chunk, CkArrayIndex &owner);

  int getCacheChunkDepth() const {
    return cacheChunkDepth_;
  }

  void print(LocalNodeType *root);

  // factory methods
  public:
  static ManagerProxyType instantiate(int cacheChunkDepth){
    return ManagerProxyType::ckNew(cacheChunkDepth);
  }

  // debugging
  public:
  void quiescence(){
    //CkPrintf("[%d] Manager::quiescence1 numTpAcks %d mergeDone %d reg %d outstanding %d %d\n", CkMyNode(), numTpAcks_, mergeDone_, registrationEntries_.size(), requestedNodes_.size(), remoteRequests_.size());
    this->contribute(CkCallback(CkIndex_Piece<TreeBuildDataInterfaceType>::quiescence(), clientProxy_));
  }

  void quiescence2(CkReductionMsg *m){
    int *n = (int *) m->getData();
    CkPrintf("[%d] Manager::quiescence2 tp reduction %d\n", CkMyNode(), *n);
    delete m;
    CkExit();
  }
};


template<typename TreeBuildDataInterfaceType>
Manager<TreeBuildDataInterfaceType>::Manager(int cacheChunkDepth) : 
  mergedRoot_(NULL),
  cacheChunkDepth_(cacheChunkDepth)
{}

template<typename TreeBuildDataInterfaceType>
void Manager<TreeBuildDataInterfaceType>::initialize(const Handle<TreeBuildDataInterfaceType> &hdl, const CkCallback &cb){
  clientProxy_ = hdl.pieceProxy();
  registrationLock_ = CmiCreateLock();
  this->contribute(cb);
}

template<typename TreeBuildDataInterfaceType>
void Manager<TreeBuildDataInterfaceType>::registration(const PointerContainer<LocalNodeType> &root, 
                                                       const PointerContainer<PieceType> &tp){
  CmiLock(registrationLock_);
  registrationEntries_.push_back_v(RegistrationEntry(root.data, tp.data));
  CmiUnlock(registrationLock_);
}

template<typename TreeBuildDataInterfaceType> 
void Manager<TreeBuildDataInterfaceType>::syncToMerge(const CkCallback &cb){
  callback_ = cb;
  numTpAcks_ = 0;
  mergeDone_ = false;
  if(registrationEntries_.size() > 0){
    mergedRoot_ = new LocalNodeType;
    mergedRoot_->setParent(NULL);
    mergePe(registrationEntries_, mergedRoot_);

    if(registrationEntries_.size() != numTpAcks_){
      CkPrintf("[%d] Manager::syncToMerge nRegistered %d numTpMsgsSent %d\n", CkMyNode(), registrationEntries_.size(), numTpAcks_);
      std::ostringstream registeredPieces;
      std::ostringstream msgPieces;
      for(int i = 0; i < registrationEntries_.size(); i++){
        registeredPieces << registrationEntries_[i].tp->getIndex() << ", ";
      }
      for(int i = 0; i < registrationMsgsSentTo_.size(); i++){
        msgPieces << registrationMsgsSentTo_[i] << ", ";
      }


      CkPrintf("[%d] Manager::syncToMerge reg: %s \n", CkMyNode(), registeredPieces.str().c_str());
      CkPrintf("[%d] Manager::syncToMerge msg: %s \n", CkMyNode(), msgPieces.str().c_str());

      
      print(mergedRoot_);
    }

    registrationMsgsSentTo_.resize(0);

    // this node was allocated by the PE itself,
    // so it doesn't belong to any particular tree piece
    mergedRoot_->tp() = -1;
    if(LocalNodeType::IsLocal(mergedRoot_->getType())){
      finishMerge();
    }
  }
  else{
    mergedRoot_ = NULL;
    finishMerge();
  }
}



template<typename TreeBuildDataInterfaceType>
void Manager<TreeBuildDataInterfaceType>::mergePe(CkVec<RegistrationEntry> &entries, LocalNodeType *dest){
  CkAssert(entries.size() > 0);
  // there are several nodes to merge; pick that one
  // to serve as the prototype, which has the most particles
  // underneath itself
  LocalNodeType *prototype = entries[0].root;
  CkAssert(prototype != NULL);
  PieceType *ownerTp = entries[0].tp;

  int nChildren = 0;
  int numNodesWithChildren = 0;
  if(prototype->getChildren() != NULL){
    nChildren = prototype->getNumChildren();
    numNodesWithChildren = 1;
  }
  int maxNumParticles = prototype->getNumParticles();

  for(int i = 1; i < entries.size(); i++){
    LocalNodeType *node = entries[i].root;
    CkAssert(node != NULL);
    if(isBetterDeal(node, prototype)){
      prototype = node;
      ownerTp = entries[i].tp;
      maxNumParticles = node->getNumParticles();
    }

    if(node->getChildren() != NULL){
      numNodesWithChildren++;
      nChildren = node->getNumChildren();
    }
  }
  *dest = *prototype;

  // terminate (base case) when there is one node, 
  // or only leaves
  if(entries.size() == 1 || numNodesWithChildren <= 1){
    // if there is only one node to merge, or if there are
    // several nodes to merge but only one node with children,
    // then we can hijack its children
    dest->setChildren(prototype->getChildren());
    if(dest->getChildren() != NULL){
      for(int i = 0; i < dest->getNumChildren(); i++){
        dest->getChild(i)->setParent(dest);
      }
      prototype->setChildren(NULL);
      prototype->setNumChildren(0);
    }

    // we are interested in descending into the tree only in certain locations.
    // for local leaves (Leaf, EmptyLeaf) and Internal nodes, there is no
    // merging to do, and no remote requests to make.  therefore, we don't
    // descend into the tree at these points.  for non-local nodes (Remote,
    // RemoteLeaf, RemoteEmptyLeaf), we just make requests for the nodes, and
    // return. however, we must descend into Boundary nodes even though there
    // may be no merging work to be done, since we have to make requests for the 
    // non-local nodes beneath them.
    noMergeDescend(dest, ownerTp);
    /*
    if(LocalNodeType::IsRemote(dest->getType())){
      // if this node is remote, all other candidates
      // (there may not be any other) should also
      // have no children, otherwise we would have picked
      // the candidate with children in the first place.
      CkAssert(numNodesWithChildren == 0);
      handleRemote(dest);
    }
    else{
      CkAssert(LocalNodeType::IsLocal(dest->getType()) || 
               LocalNodeType::IsBoundary(dest->getType()));
      // we've already usurped the tree beneath 'prototype', so 
      // we can descend into it.
      noMergeDescend(dest, ownerTp);
    }
    */
    
    return;
  }

  CkVec<RegistrationEntry> newEntries;

  dest->allocateChildren(nChildren);
  dest->getPayload().mergeStart();


  bool allChildrenLocal = true;
  bool anyChildLocal = false;

  int nChildrenPending = 0;
  int nChildrenParticles = 0;

  for(int i = 0; i < nChildren; i++){
    for(int j = 0; j < entries.size(); j++){
      LocalNodeType *node = entries[j].root;
      PieceType *tp = entries[j].tp;

      if(node->getChildren() != NULL){
        newEntries.push_back(RegistrationEntry(node->getChild(i), tp));
      }
    }

    LocalNodeType *newDest = dest->getChild(i);
    newDest->setParent(dest);
    mergePe(newEntries, newDest);

    // this node was allocated by the PE itself (in allocateChildren, above)
    // so it doesn't belong to any particular tree piece
    newDest->tp() = -1;

    bool childLocal = LocalNodeType::IsLocal(newDest->getType());

    if(childLocal){
      dest->getPayload().mergeAccumulate(newDest->getPayload());
    }
    else{
      nChildrenPending++;
    }

    allChildrenLocal &= childLocal;
    anyChildLocal |= childLocal;

    nChildrenParticles += newDest->getNumParticles();

    newEntries.resize(0);
  }

  dest->pending() = nChildrenPending;
  dest->setNumParticles(nChildrenParticles);

  MDT_MANAGER_VERBOSE_FN("[%d] Manager<>::mergePe halfready key %llu pending %d\n", CkMyNode(), dest->getKey(), nChildrenPending);

  if(allChildrenLocal){
    CkAssert(nChildrenPending == 0);
    dest->getPayload().mergeFinish();
    // since 'dest' has children and all of them
    // are local, it must itself be internal
    dest->setType(Internal);
    payloadReady(dest);
  }
  else if(anyChildLocal){
    // not all children are local, but at least one is
    CkAssert(nChildrenPending > 0);
    dest->setType(Boundary);
  }
  // otherwise, don't do anything; type of prototype
  // will have been copied (correctly) to dest

  return;
}

template<typename TreeBuildDataInterfaceType>
void Manager<TreeBuildDataInterfaceType>::noMergeDescend(LocalNodeType *node, PieceType *ownerTp){
  if(LocalNodeType::IsLocal(node->getType())){
    handleLocal(node, ownerTp);
    return;
  }
  else if(LocalNodeType::IsRemote(node->getType())){
    handleRemote(node);
    return;
  }

  int nChildrenPending = 0;
  node->getPayload().mergeStart();

  CkAssert(LocalNodeType::IsBoundary(node->getType()));
  CkAssert(node->getChildren() != NULL);
  for(int i = 0; i < node->getNumChildren(); i++){
    LocalNodeType *child = node->getChild(i);

    noMergeDescend(child, ownerTp);

    if(LocalNodeType::IsLocal(child->getType())){
      node->getPayload().mergeAccumulate(child->getPayload());
    }
    else{
      nChildrenPending++;
    }
  }

  CkAssert(nChildrenPending > 0);
  node->pending() = nChildrenPending;
}

template<typename TreeBuildDataInterfaceType>
void Manager<TreeBuildDataInterfaceType>::handleLocal(LocalNodeType *node, PieceType *ownerTp){
  node->pending() = 0;
  
  CkAssert(node->getOwnerStart() == ownerTp->getIndex());
  clientProxy_[ownerTp->getIndex()].registration(PointerContainer<LocalNodeType>(node));
  registrationMsgsSentTo_.push_back(ownerTp->getIndex());


  MDT_MANAGER_VERBOSE_FN("[%d] Manager<>::handleLocal ready local key %llu\n", CkMyNode(), node->getKey());
  payloadReady(node);
  // so that we wait for this tp to ack that it has
  // registered its subtree, before terminating the 
  // tree merge procedure
  numTpAcks_++;
  // this is done by the tree piece itself
  //searchForLeaves(node, ownerTp);
}

template<typename TreeBuildDataInterfaceType>
void Manager<TreeBuildDataInterfaceType>::doneRegistration(const PointerContainer<LocalNodeType> &n){
  LocalNodeType *node = n.data;
  CkAssert(node->pending() == 0);
  CkAssert(numTpAcks_ > 0);
  numTpAcks_--;
  if(numTpAcks_ == 0 && mergeDone_){
    trulyFinishMerge();
  }
}


template<typename TreeBuildDataInterfaceType>
void Manager<TreeBuildDataInterfaceType>::handleRemote(LocalNodeType *node){
  // if this node is remote, request its user data
  // from its owner tree piece
  KeyType key = node->getKey();
  int owner = Utility::pickRandom(node->getOwnerStart(), node->getOwnerEnd());
  // set a single owner that will be used by all traversals on this PE
  node->home() = owner;
  clientProxy_[owner].requestPayload(key, CkMyNode());
  requestedNodes_[key] = node;

  // leaves are always ready since they have no children
  node->pending() = 0;

  MDT_MANAGER_VERBOSE_FN("[%d] Manager<>::handleRemote REQUEST remote key %llu from %d\n", CkMyNode(), node->getKey(), owner);
}

template<typename TreeBuildDataInterfaceType>
typename TreeBuildDataInterfaceType::LocalNodeType *Manager<TreeBuildDataInterfaceType>::keyToNode(KeyType k){
  typename std::map<KeyType, LocalNodeType *>::iterator it;
  it = localNodes_.find(k);
  if(it == localNodes_.end()) return NULL;
  else return it->second;
}

template<typename TreeBuildDataInterfaceType>
void Manager<TreeBuildDataInterfaceType>::requestPayload(KeyType k, int replyTo){
  LocalNodeType *node = keyToNode(k);
  MDT_MANAGER_VERBOSE_FN("[%d] Manager<>::requestPayload key %llu replyTo %d\n", CkMyNode(), k, replyTo);
  if(node != NULL){
    // since nodes are added to this table only once   
    // their moments have been calculated, this node's 
    // moments must be ready: send to requestor.       
    // Moreover, there can be no other pending         
    // requestors whose requests have already been     
    // received, save the current one, because when a  
    // node is added to the node table, its moments are
    // sent to all pending requestors.                   
    MDT_MANAGER_VERBOSE_FN("[%d] Manager<>::requestPayload READY sending key %llu replyTo %d\n", CkMyNode(), k, replyTo);
    this->thisProxy[replyTo].receivePayload(k, node->getPayload());
  }
  else{
    CkVec<int> &requestors = remoteRequests_[k];
    // node not yet ready; enqueue request
    requestors.push_back(replyTo);
  }
}

template<typename TreeBuildDataInterfaceType>
void Manager<TreeBuildDataInterfaceType>::receivePayload(KeyType k, const typename TreeBuildDataInterfaceType::PayloadType &payload){
  MDT_MANAGER_VERBOSE_FN("[%d] Manager<>::receivePayload key %llu\n", CkMyNode(), k);
  typename std::map<KeyType, LocalNodeType *>::iterator it;
  it = requestedNodes_.find(k);
  CkAssert(it != requestedNodes_.end());

  LocalNodeType *node = it->second;
  requestedNodes_.erase(it);

  node->getPayload() = payload;
  passUpwards(node);
}

template<typename TreeBuildDataInterfaceType>
void Manager<TreeBuildDataInterfaceType>::payloadReady(LocalNodeType *node){
  KeyType k = node->getKey();
  localNodes_[k] = node;
  MDT_MANAGER_VERBOSE_FN("[%d] Manager<>::payloadReady READY key %llu\n", CkMyNode(), k);
  // check whether this node has been requested by someone else
  typename std::map<KeyType, CkVec<int> >::iterator it;
  it = remoteRequests_.find(k);
  if(it != remoteRequests_.end()){
    // if so, send it to all who have requested it
    CkVec<int> &requestors = it->second;
    for(int i = 0; i < requestors.size(); i++){
      int replyTo = requestors[i];
      MDT_MANAGER_VERBOSE_FN("[%d] Manager<>::payloadReady READY sending key %llu replyTo %d\n", CkMyNode(), k, replyTo);
      this->thisProxy[replyTo].receivePayload(k, node->getPayload());
    }
    remoteRequests_.erase(it);
  }
}


template<typename TreeBuildDataInterfaceType>
void Manager<TreeBuildDataInterfaceType>::passUpwards(LocalNodeType *node){
  payloadReady(node);
  LocalNodeType *parent = node->getParent();
  if(parent != NULL){
    CkAssert(parent->pending() > 0);
    parent->getPayload().mergeAccumulate(node->getPayload());
    parent->pending()--;
    MDT_MANAGER_VERBOSE_FN("[%d] Manager<>::passUpwards key %llu parent %llu pending %d\n", CkMyNode(), node->getKey(), parent->getKey(), parent->pending());
    if(parent->pending() == 0){
      parent->getPayload().mergeFinish();
      passUpwards(parent);
    }
  }
  else{
    finishMerge();
  }
}

template<typename TreeBuildDataInterfaceType>
void Manager<TreeBuildDataInterfaceType>::finishMerge(){
  // we have all the moments required 
  mergeDone_ = true;
  // have all the local tp's ack'ed?
  if(numTpAcks_ == 0){
    trulyFinishMerge();
  }
}

template<typename TreeBuildDataInterfaceType>
void Manager<TreeBuildDataInterfaceType>::trulyFinishMerge(){
  mergeDone_ = false;
  numTpAcks_ = -1;
  // error checking
  printOutstanding();
  // diagnostic
  //if(mergedRoot_ != NULL) print(mergedRoot_);

  // look at clients' submissions 
  for(int i = 0; i < registrationEntries_.size(); i++){
    PieceType *tp = registrationEntries_[i].tp;
    LocalNodeType *tproot = registrationEntries_[i].root;
    tproot->deleteBeneath();
    delete tproot;
  }
  registrationEntries_.resize(0);


  // done
  this->contribute(callback_);
}

template<typename TreeBuildDataInterfaceType>
void Manager<TreeBuildDataInterfaceType>::printOutstanding(){
  if(requestedNodes_.size() != 0){
    std::ostringstream oss;
    typename std::map<KeyType, LocalNodeType *>::iterator it;
    for(it = requestedNodes_.begin(); it != requestedNodes_.end(); ++it){
      oss << it->first << ", ";
    }
    CkPrintf("[%d] Distree::Manager::finishMerge outstanding requested nodes: %s\n", CkMyNode(), oss.str().c_str());
    CkAbort("Didn't receive all requested nodes\n");
  }

  // error checking
  if(remoteRequests_.size() != 0){
    std::ostringstream oss;
    typename std::map<KeyType, CkVec<int> >::iterator it;
    for(it = remoteRequests_.begin(); it != remoteRequests_.end(); ++it){
      oss << it->first << ", ";
    }
    CkPrintf("[%d] Distree::Manager::finishMerge outstanding moment requests: %s\n", CkMyNode(), oss.str().c_str());
    CkAbort("Didn't fulfill all moment requests\n");
  }
}



// used during traversal

template<typename TreeBuildDataInterfaceType>
template<typename TraversalDataInterface>
void Manager<TreeBuildDataInterfaceType>::requestNode(KeyType key, const CProxy_CkCacheManager<KeyType> &cacheProxy, int replyTo){
  LocalNodeType *node = keyToNode(key);
  packAndSendNode<TraversalDataInterface>(node, cacheProxy, replyTo);
}

template<typename TreeBuildDataInterfaceType>
template<typename TraversalDataInterface>
void Manager<TreeBuildDataInterfaceType>::
packAndSendNode(LocalNodeType *node, const CProxy_CkCacheManager<KeyType> &cacheProxy, int replyTo){
  CkAssert(node != NULL);
  CkAssert(LocalNodeType::IsLocal(node->getType()) || LocalNodeType::IsBoundary(node->getType()));

  int depth = cacheChunkDepth_;
  CkCacheFillMsg<KeyType> *msg = pack<TraversalDataInterface>(node, depth);
  cacheProxy[replyTo].recvData(msg);
}

template<typename TreeBuildDataInterfaceType>
template<typename TraversalDataInterface>
CkCacheFillMsg<typename TreeBuildDataInterfaceType::LocalNodeType::KeyType> *Manager<TreeBuildDataInterfaceType>::
pack(const LocalNodeType *root, int depth){
  // this count doesn't include 'root' itself
  int nNodes = countNodesBeneath(root, depth);
  int msgSize = sizeof(CachedNodes<TraversalDataInterface>) + 
                nNodes * sizeof(typename TraversalDataInterface::RemoteNodeType);
  int prioBits = 8 * sizeof(int);

  CkCacheFillMsg<KeyType> *msg = new (msgSize, prioBits) CkCacheFillMsg<KeyType>(root->getKey());
  CachedNodes<TraversalDataInterface> *c = (CachedNodes<TraversalDataInterface> *) msg->data;
  // no need to set the message pointer right now; it'll be lost upon 
  // message transmission anyway
  typename TraversalDataInterface::RemoteNodeType *copy = &c->nodes[0];
  typename TraversalDataInterface::RemoteNodeType *buf = copy + 1;
  *copy = *root;
  recursivePack<TraversalDataInterface>(root, copy, buf, depth);

  long offset = 1;
  copy->setChildren((typename TraversalDataInterface::RemoteNodeType *) offset);
  c->nNodes = nNodes;

  *((int*) CkPriorityPtr(msg)) = TraversalDataInterface::NodeReceivePriority;
  CkSetQueueing(msg, CK_QUEUEING_IFIFO);

  return msg;
}

template<typename TreeBuildDataInterfaceType>
int Manager<TreeBuildDataInterfaceType>::
countNodesBeneath(const LocalNodeType *root, int depth){
  if(root->getChildren() == NULL || depth == 0) return 0;
  else{
    int num = 0;
    CkAssert(root->getChildren() != NULL);
    for(int i = 0; i < root->getNumChildren(); i++){
      num += countNodesBeneath(root->getChild(i), depth-1);
      num++;
    }
    return num;
  }
}

template<typename TreeBuildDataInterfaceType>
template<typename TraversalDataInterface>
void 
Manager<TreeBuildDataInterfaceType>::
recursivePack(const LocalNodeType *root, typename TraversalDataInterface::RemoteNodeType *copy, typename TraversalDataInterface::RemoteNodeType *&buf, int depth){
  if(depth == 0 || LocalNodeType::IsRemote(copy->getType())){
    // this will be read as offset '0' while unpacking
    copy->setChildren(NULL);
    return;
  }

  copy->setChildren(buf);
  // since children of a node are stored together, we 
  // copy en masse
  //memcpy(buf, root->getChildren(), root->getNumChildren() * sizeof(MyNodeType));
  const LocalNodeType *rootChildren = root->getChildren();
  if(rootChildren != NULL){
    for(int i = 0; i < root->getNumChildren(); i++){
      // data interface class should have defined:
      // this is done so that vptr is set appropriately
      typename TraversalDataInterface::RemoteNodeType *child = new (&buf[i]) typename TraversalDataInterface::RemoteNodeType; 
      TraversalDataInterface::copy(*child, rootChildren[i]);
    }
  }

  // FIXME - fuse these two if's?
  if(root->getChildren() != NULL){
    buf += root->getNumChildren(); 
    depth--;
    for(int i = 0; i < root->getNumChildren(); i++){
      recursivePack<TraversalDataInterface>(root->getChild(i), copy->getChild(i), buf, depth); 
    }
  }

  // now, change children pointers to offsets
  long offset = (copy->getChildren() - copy);
  copy->setChildren((typename TraversalDataInterface::RemoteNodeType *)offset);
}

// this method shouldn't be called, because leaves should always be
// internal to treepieces, and requests for them shouldn't get past the
// Piece::requestLeafCntents method
/*
template<typename TreeBuildDataInterfaceType>
template<typename TraversalDataInterface>
void Manager<TreeBuildDataInterfaceType>::
requestLeafContents(KeyType key, const CProxy_CkCacheManager<KeyType> &cacheProxy, int replyTo){
  LocalNodeType *leaf = keyToNode(key);
  packAndSendLeaf(leaf, cacheProxy, replyTo);
}
*/

template<typename TreeBuildDataInterfaceType>
template<typename TraversalDataInterface>
void Manager<TreeBuildDataInterfaceType>::
packAndSendLeaf(LocalNodeType *leaf, const CProxy_CkCacheManager<KeyType> &cacheProxy, int replyTo){
  CkAssert(leaf != NULL);
  CkAssert(LocalNodeType::IsLeaf(leaf->getType()));
  CkCacheFillMsg<KeyType> *msg = pack<TraversalDataInterface>(leaf);
  cacheProxy[replyTo].recvData(msg);
}


template<typename TreeBuildDataInterfaceType>
template<typename TraversalDataInterface>
CkCacheFillMsg<typename TreeBuildDataInterfaceType::LocalNodeType::KeyType> *
Manager<TreeBuildDataInterfaceType>::
pack(const LocalNodeType *leaf){
  // first particle is copied to 'CachedLeafContents' struct
  int msgSize = sizeof(CachedLeafContents<TraversalDataInterface>) + 
                (leaf->getNumParticles() - 1) * sizeof(typename TraversalDataInterface::RemoteParticleType);
  int prioBits = 8 * sizeof(int);

  CkCacheFillMsg<KeyType> *msg = NULL;
  msg = new (msgSize, prioBits) CkCacheFillMsg<KeyType>(leaf->getKey());
  CachedLeafContents<TraversalDataInterface> *c = (CachedLeafContents<TraversalDataInterface> *) msg->data;
  // can't set the message pointer right now; 
  for(int i = 0; i < leaf->getNumParticles(); i++){
    typename TraversalDataInterface::RemoteParticleType *p = new (&c->particles[i]) typename TraversalDataInterface::RemoteParticleType;
    TraversalDataInterface::copy(*p, leaf->getParticles()[i]);
  }
  c->nParticles = leaf->getNumParticles();

  *((int*) CkPriorityPtr(msg)) = TraversalDataInterface::LeafReceivePriority;
  CkSetQueueing(msg, CK_QUEUEING_IFIFO);

  return msg;
}

template<typename TreeBuildDataInterfaceType>
template<typename TraversalDataInterface>
void *Manager<TreeBuildDataInterfaceType>::
unpackNodes(CkCacheFillMsg<KeyType> *msg, int chunk, CkArrayIndex &owner){
  // extract struct to save message for entry
  CachedNodes<TraversalDataInterface> *c = (CachedNodes<TraversalDataInterface> *) msg->data;
  // save pointer to message, to be deleted when free() is called
  c->msg = msg;
  // home of this node is given by owner
  int home = owner.data()[0];

  typename TraversalDataInterface::RemoteNodeType *root = &c->nodes[0];
  // recunpack
  int nUnpacked = recursiveUnpack<TraversalDataInterface>(root, home);
  CkAssert(c->nNodes == nUnpacked);

  return c;
}

template<typename TreeBuildDataInterfaceType>
template<typename TraversalDataInterface>
int Manager<TreeBuildDataInterfaceType>::
recursiveUnpack(typename TraversalDataInterface::RemoteNodeType *node, int home){
  // we can't have cached nodes here, since they
  // aren't connected to the local trees, and
  // therefore couldn't have been packed into this
  // response.
  CkAssert(!TraversalDataInterface::RemoteNodeType::IsCached(node->getType()));

  // newly created nodes have -1 home. during tree building,
  // we set homes of remote nodes to non-negative values. so, only
  // remote nodes should have non-negative homes. 
  // note that after consolidation (merging), pieces on an SMP
  // node are connected by a top-level tree structure. So,
  // when we make a request for a node, given that we will fetch a
  // depth-limited subtree underneath that node, we could get portions
  // of several tree pieces in that subtree. However, it is incorrect
  // to treat all nodes received here in this manner as coming from the 
  // same home. so, we need to check whether a node has a >= 0 home. if 
  // so, it was a remote node on the SMP processor from which we fetched it, 
  // and we preserve its home. 
  // otherwise, it was a local node to that SMP processor. in this case, we must
  // distinguish between top-level nodes (whose home can be any of the 
  // 
  if(node->home() < 0){
    if(node->getOwnerStart() + 1 == node->getOwnerEnd()){
      // only one owner, so it must be internal to that
      // owner tree piece. therefore, home has to be set
      // to that owner tree piece
      node->home() = node->getOwnerStart();
    }
    else{
      // several owners, including the tree piece that we just fetched
      // this chunk from; set home to that tree piece.
      node->home() = home;
    }
  }
  else{
    CkAssert(TraversalDataInterface::RemoteNodeType::IsRemote(node->getType()));
  }

  node->setType(TraversalDataInterface::RemoteNodeType::GetCached(node->getType()));

  int underChildren = 0;
  int nChildren = 0;
  // compiler complains about loss of precision if this is cast to an int
  CmiUInt8 childOffset = (CmiUInt8) (node->getChildren());
  if(childOffset != 0){
    // in recursivePack(), offsets were calculated relative to each node 
    node->setChildren(node + childOffset);
    for(int i = 0; i < node->getNumChildren(); i++){
      underChildren += recursiveUnpack<TraversalDataInterface>(node->getChild(i), home);
      nChildren++;
    }
  }
  else{
    node->setChildren(NULL);
  }
  return nChildren + underChildren;
}

template<typename TreeBuildDataInterfaceType>
template<typename TraversalDataInterface>
void *Manager<TreeBuildDataInterfaceType>::
unpackLeafContents(CkCacheFillMsg<KeyType> *msg, int chunk, CkArrayIndex &owner){
  CachedLeafContents<TraversalDataInterface> *c = (CachedLeafContents<TraversalDataInterface> *) msg->data;
  c->msg = msg;
  return c;
}

/*
template<typename TreeBuildDataInterfaceType>
template<typename TraversalDataInterface>
void Manager<TreeBuildDataInterfaceType>::accumulateNodes(CkCacheFillMsg<KeyType> *msg){
  CachedNodes<TraversalDataInterface> *c = (CachedNodes<TraversalDataInterface> *) msg->data;
  for(int i = 1; i < c->nNodes; i++){
    const typename TraversalDataInterface::RemoteNodeType *source = &c->nodes[i];
    LocalNodeType *target = keyToNode(source->getKey()); 
    CkAssert(target != NULL);
    TraversalDataInterface::accumulate(*target, *source);
  }

  CkAbort("Currently unsafe to call accumulateNodes with nodegroup TreeManager\n");
  cmpletionDetectorProxy.ckLocalBranch()->consume();

  delete msg;
}
*/

/*
template<typename TreeBuildDataInterfaceType>
template<typename TraversalDataInterface>
void Manager<TreeBuildDataInterfaceType>::accumulateLeafContents(CkCacheFillMsg<KeyType> *msg){
  CachedLeafContents<TraversalDataInterface> *c = (CachedLeafContents<TraversalDataInterface> *) msg->data;
  LocalNodeType *leaf = keyToNode(msg->key);
  CkAssert(leaf != NULL);
  LocalParticleType *targets = leaf->getParticles(); 
  for(int i = 0; i < c->nParticles; i++){
    typename TraversalDataInterface::RemoteParticleType *source = &c->particles[i];
    TraversalDataInterface::accumulate(targets[i], *source);
  }

  cmpletionDetectorProxy.ckLocalBranch()->consume();

  delete msg;
}
*/


template<typename TreeBuildDataInterfaceType>
typename TreeBuildDataInterfaceType::LocalNodeType *Manager<TreeBuildDataInterfaceType>::
root(){
  return mergedRoot_;
}

template<typename TreeBuildDataInterfaceType>
void Manager<TreeBuildDataInterfaceType>::
print(LocalNodeType *root){
  std::ostringstream oss;
  oss << "tree.merged." << CkMyNode() << ".dot";
  std::ofstream out(oss.str().c_str());
  CkAssert(out.is_open());
  out << "digraph mergedtree" << CkMyNode() << "{" << std::endl;
  root->dot(out);
  out << "}" << std::endl;
  out.close();
}

template<typename TreeBuildDataInterfaceType>
void Manager<TreeBuildDataInterfaceType>::
syncToDelete(const CkCallback &cb){
  //CkPrintf("[%d] Manager::syncToDelete\n", CkMyNode());
  callback_ = cb;
  // how many deletion acks must i await 
  // before deleting the top level of the
  // tree?
  numTpAcks_ = 0;
  if(mergedRoot_ != NULL){
    sendDeleteRequests(mergedRoot_);
  }
  //CkPrintf("[%d] Manager::syncToDelete numTpAcks: %d\n", CkMyNode(), numTpAcks_);
  if(numTpAcks_ == 0){
    deleteTopLevelTree();
  }
}

template<typename TreeBuildDataInterfaceType>
void Manager<TreeBuildDataInterfaceType>::
sendDeleteRequests(LocalNodeType *node){
  // if a subtree is owned by a single piece that is 
  // local to this node, send a message to it, asking it
  // to delete the subtree beneath that node.
  if((node->getOwnerEnd() - node->getOwnerStart()) == 1 && 
      LocalNodeType::IsLocal(node->getType())){
    clientProxy_[node->getOwnerStart()].deleteLocalTree(PointerContainer<LocalNodeType>(node));
    //CkPrintf("[%d] Manager::sendDeleteRequests tp: %d numTpAcks_: %d\n", CkMyNode(), node->getOwnerStart(), numTpAcks_);
    numTpAcks_++;
  }
  else if(node->getChildren() == NULL){
    // since the node is not local and it 
    // has no children, it must be a remote node
    CkAssert(LocalNodeType::IsRemote(node->getType()));
  }
  else{
    for(int i = 0; i < node->getNumChildren(); i++){
      sendDeleteRequests(node->getChild(i));
    }
  }
}

template<typename TreeBuildDataInterfaceType>
void Manager<TreeBuildDataInterfaceType>::
doneDelete(){
  CkAssert(numTpAcks_ > 0);
  numTpAcks_--;
  //CkPrintf("[%d] Manager::doneDelete numTpAcks_: %d\n", CkMyNode(), numTpAcks_);

  if(numTpAcks_ == 0){
    deleteTopLevelTree();
  }
}

template<typename TreeBuildDataInterfaceType>
void Manager<TreeBuildDataInterfaceType>::
deleteTopLevelTree(){
  if(mergedRoot_ != NULL){
    mergedRoot_->deleteBeneath();
    delete mergedRoot_;
    mergedRoot_ = NULL;

    localNodes_.clear();
  }

  CkAssert(remoteRequests_.size() == 0);
  CkAssert(requestedNodes_.size() == 0);

  this->contribute(callback_);
}
 

}; // namespace MultiphaseDistree

#endif // MULTIPHASE_DISTREE_MANAGER_H
