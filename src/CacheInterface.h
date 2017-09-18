#ifndef MULTIPHASE_DISTREE_CACHE_INTERFACE_H
#define MULTIPHASE_DISTREE_CACHE_INTERFACE_H

#include "mdt.decl.h"
#include "Node.h"
#include "Utility.h"
#include "CkCache.h"
#include "Handle.h"
#include "CachedData.h"

#include "completion.h"

namespace MultiphaseDistree {

template<typename TraversalDataInterface>
class NodeCacheInterface : public CkCacheEntryType<typename TraversalDataInterface::TreeBuildDataInterfaceType::KeyType> {
  protected:
  typedef typename TraversalDataInterface::TreeBuildDataInterfaceType TreeBuildDataInterface;
  typedef typename TreeBuildDataInterface::KeyType KeyType;
  typedef CProxy_CkCacheManager<KeyType> CacheProxyType;
  typedef CProxy_Piece<TreeBuildDataInterface> ClientProxyType;
  typedef CachedNodes<TraversalDataInterface> CachedNodesType;
  typedef typename TraversalDataInterface::RemoteNodeType RemoteNodeType;

  typedef Manager<TreeBuildDataInterface> TreeManagerType;
  typedef CProxy_Manager<TreeBuildDataInterface> TreeManagerProxyType;
  typedef TraversalManager<TraversalDataInterface> TraversalManagerType;

  private:
  CacheProxyType cacheProxy_; 
  ClientProxyType clientProxy_;
  TreeManagerType *treeManager_;
  CProxy_CompletionDetector completionDetectorProxy_;
  TraversalManagerType *traversal_;


  public:
  NodeCacheInterface() : 
    treeManager_(NULL),
    traversal_(NULL)
  {}

  void setTraversalManager(TraversalManagerType *traversal){
    traversal_ = traversal;
  }

  void setClientProxy(const ClientProxyType &clientProxy){
    clientProxy_ = clientProxy;
  }

  ClientProxyType &clientProxy(){
    return clientProxy_;
  }

  void setCacheProxy(const CacheProxyType &cacheProxy){
    cacheProxy_ = cacheProxy;
  }

  CacheProxyType &getCacheProxy(){
    return cacheProxy_;
  }

  void setTreeManager(TreeManagerType *treeManager){
    treeManager_ = treeManager;
  }

  void setCompletionDetector(const CProxy_CompletionDetector &completionDetectorProxy){
    completionDetectorProxy_ = completionDetectorProxy;
  }

  CProxy_CompletionDetector &getCompletionDetectorProxy(){
    return completionDetectorProxy_;
  }

  void *request(CkArrayIndex &idx, KeyType key){
    int home = *idx.data();
    CkAssert(home >= 0);
    CkEntryOptions opts;
    opts.setPriority(TraversalDataInterface::NodeRequestPriority);
    clientProxy_[home].template requestNode<TraversalDataInterface>(key, cacheProxy_, CkMyPe(), &opts);

    traversal_->nodeCacheMiss();

    return NULL;
  }

  // Return data from fufilled cache request.
  void *unpack(CkCacheFillMsg<KeyType> *msg, int chunk, CkArrayIndex &owner){
    return treeManager_->template unpackNodes<TraversalDataInterface>(msg, chunk, owner);
  }

  virtual void writeback(CkArrayIndex &idx, KeyType key, void *data){
  }

  /// @brief free cached data.
  virtual void free(void *data){
    CachedNodesType *c = (CachedNodesType *) data;
    delete c->msg;
  }

  /// @brief return size of cached data.
  int size(void *data){
    CachedNodesType *c = (CachedNodesType *) data;
    return sizeof(CachedNodesType) + c->nNodes * sizeof(RemoteNodeType);
  }

  virtual void done(int chunk){
    cacheProxy_[CkMyPe()].ckLocalBranch()->finishedChunk(chunk, 0);
  }
};

/*
template<typename TraversalDataInterface>
class WritebackNodeCacheInterface : public NodeCacheInterface<TraversalDataInterface> {
  typedef NodeCacheInterface<TraversalDataInterface> ParentType;
  public:
  void writeback(CkArrayIndex &idx, typename ParentType::KeyType key, void *data){
    int home = *idx.data();
    CkAssert(home >= 0);
    // we can get the msg in which we received the data,
    // from the data itself, and send that message back to the
    // home piece
    typename ParentType::CachedNodesType *c = (typename ParentType::CachedNodesType *) data;
    this->clientProxy()[home].template accumulateNodes<TraversalDataInterface>(c->msg);
    this->getCompletionDetectorProxy().ckLocalBranch()->produce();
  }


  // don't do anything: message will have been
  // reused in the sending back of accumulation data
  void free(void *data){
  }

  void done(int chunk){
    this->getCacheProxy()[CkMyPe()].ckLocalBranch()->writebackChunk(chunk);
    this->getCompletionDetectorProxy().ckLocalBranch()->done();
    this->getCacheProxy()[CkMyPe()].ckLocalBranch()->finishedChunk(chunk, 0);
  }
};
*/

template<typename TraversalDataInterface>
class LeafCacheInterface : public CkCacheEntryType<typename TraversalDataInterface::TreeBuildDataInterfaceType::KeyType> {
  protected:
  typedef typename TraversalDataInterface::TreeBuildDataInterfaceType TreeBuildDataInterface;
  typedef typename TreeBuildDataInterface::KeyType KeyType;
  typedef CProxy_Piece<TreeBuildDataInterface> ClientProxyType;
  typedef CachedLeafContents<TraversalDataInterface> CachedLeafContentsType;
  typedef CProxy_CkCacheManager<KeyType> CacheProxyType;
  typedef typename TraversalDataInterface::RemoteParticleType RemoteParticleType;

  typedef Manager<TreeBuildDataInterface> TreeManagerType;
  typedef CProxy_Manager<TreeBuildDataInterface> TreeManagerProxyType;
  typedef TraversalManager<TraversalDataInterface> TraversalManagerType;

  private:
  CacheProxyType cacheProxy_;
  ClientProxyType clientProxy_;
  TreeManagerType *treeManager_;
  CProxy_CompletionDetector completionDetectorProxy_;
  TraversalManagerType *traversal_;

  public:
  LeafCacheInterface() :
    treeManager_(NULL),
    traversal_(NULL)
  {}

  void setTraversalManager(TraversalManagerType *traversal){
    traversal_ = traversal;
  }

  void setClientProxy(const ClientProxyType &clientProxy){
    clientProxy_ = clientProxy;
  }

  ClientProxyType &clientProxy(){
    return clientProxy_;
  }

  void setCacheProxy(const CacheProxyType &cacheProxy){
    cacheProxy_ = cacheProxy;
  }

  CacheProxyType &getCacheProxy(){
    return cacheProxy_;
  }

  void setTreeManager(TreeManagerType *treeManager){
    treeManager_ = treeManager;
  }

  void setCompletionDetector(const CProxy_CompletionDetector &completionDetectorProxy){
    completionDetectorProxy_ = completionDetectorProxy;
  }

  CProxy_CompletionDetector &getCompletionDetectorProxy(){
    return completionDetectorProxy_;
  }

  void *request(CkArrayIndex &idx, KeyType key){
    int home = *idx.data();
    CkAssert(home >= 0);
    CkEntryOptions opts;
    opts.setPriority(TraversalDataInterface::LeafRequestPriority);
    clientProxy_[home].template requestLeafContents<TraversalDataInterface>(key, cacheProxy_, CkMyPe(), &opts);

    traversal_->leafCacheMiss();

    return NULL;
  }

  // Return data from fufilled cache request.
  void *unpack(CkCacheFillMsg<KeyType> *msg, int chunk, CkArrayIndex &owner){
    return treeManager_->template unpackLeafContents<TraversalDataInterface>(msg, chunk, owner);
  }

  virtual void writeback(CkArrayIndex &idx, KeyType key, void *data){
  }

  // free cached data.
  virtual void free(void *data){
    // here, the data pointer points to the first particle
    // go back to the start of the msg pointer
    CachedLeafContentsType *c = (CachedLeafContentsType *) data;
    delete c->msg;
  }

  // return size of cached data.
  int size(void *data){
    CachedLeafContentsType *c = (CachedLeafContentsType *) data;
    return sizeof(CachedLeafContentsType) + (c->nParticles * sizeof(RemoteParticleType));
  }

  virtual void done(int chunk){
    cacheProxy_[CkMyPe()].ckLocalBranch()->finishedChunk(chunk, 0);
  }
};

template<typename TraversalDataInterface>
class WritebackLeafCacheInterface : public LeafCacheInterface<TraversalDataInterface> {
  typedef LeafCacheInterface<TraversalDataInterface> ParentType;

  public:
  // don't do anything: will have reused the message
  // associated with this entry in sending back the
  // data to accumulate on the home piece
  void free(void *data){
  }

  void writeback(CkArrayIndex &idx, typename ParentType::KeyType key, void *data){
    int home = *idx.data();
    CkAssert(home >= 0);
    // we can get the msg in which we received the data,
    // from the data itself, and send that message back to the
    // home piece
    typename ParentType::CachedLeafContentsType *c = (typename ParentType::CachedLeafContentsType *) data;
    this->clientProxy()[home].template accumulateLeafContents<TraversalDataInterface>(c->msg);
    this->getCompletionDetectorProxy().ckLocalBranch()->produce();
  }

  void done(int chunk){
    this->getCacheProxy()[CkMyPe()].ckLocalBranch()->writebackChunk(chunk);
    this->getCompletionDetectorProxy().ckLocalBranch()->done();
    this->getCacheProxy()[CkMyPe()].ckLocalBranch()->finishedChunk(chunk, 0);
  }
};

}; // namespace MultiphaseDistree


#endif // MULTIPHASE_DISTREE_CACHE_INTERFACE_H
