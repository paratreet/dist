#ifndef MULTIPHASE_DISTREE_TRAVERSAL_HANDLE_H
#define MULTIPHASE_DISTREE_TRAVERSAL_HANDLE_H

#include "mdt.decl.h"
#include "Node.h"
#include "Utility.h"
#include "CkCache.h"
#include "CacheAccessType.h"
#include "completion.h"
#include "TraversalStatistics.h"

namespace MultiphaseDistree {

template<typename TraversalType>
class TraversalHandle {
  typename TraversalType::ProxyType proxy_;
  typename TraversalType::CacheProxyType nodeCacheProxy_;
  typename TraversalType::CacheProxyType leafCacheProxy_;
  CProxy_CompletionDetector completionDetectorProxy_;

  CacheAccessType nodeCacheAccessType_;
  CacheAccessType leafCacheAccessType_;

  public:

  TraversalHandle() :
    nodeCacheAccessType_(CacheAccessType::Invalid),
    leafCacheAccessType_(CacheAccessType::Invalid)
  {}

  // const versions
  const typename TraversalType::CacheProxyType &nodeCacheProxy() const { return nodeCacheProxy_; }
  const typename TraversalType::CacheProxyType &leafCacheProxy() const { return leafCacheProxy_; }
  const typename TraversalType::ProxyType &proxy() const { return proxy_; }

  // non-const versions
  typename TraversalType::CacheProxyType &nodeCacheProxy() { return nodeCacheProxy_; }
  typename TraversalType::CacheProxyType &leafCacheProxy() { return leafCacheProxy_; }
  typename TraversalType::ProxyType &proxy() { return proxy_; }


  void pup(PUP::er &p){
    p | nodeCacheAccessType_;
    p | leafCacheAccessType_;
    p | proxy_;
    p | nodeCacheProxy_;
    p | leafCacheProxy_;
  }

  TraversalType *local() {
    return proxy()[CkMyPe()].ckLocal();
  }

  template<typename TreeHandleType>
  void initialize(const TreeHandleType &treeHandle, const CkCallback &cb) {
    proxy().initialize(treeHandle, nodeCacheProxy(), leafCacheProxy(), cb);
    completionDetectorProxy_ = treeHandle.completionDetectorProxy();
  }

  void synchronize(){
    proxy().synchronize(CkCallbackResumeThread());
  }

  void done(){
    bool doNodeCacheWriteback = (nodeCacheAccessType_ == CacheAccessType::Accumulate); 
    bool doLeafCacheWriteback = (leafCacheAccessType_ == CacheAccessType::Accumulate); 
    //CkPrintf("TraversalHandle::done node: %d leaf %d\n", doNodeCacheWriteback, doLeafCacheWriteback);

    if(doNodeCacheWriteback || doLeafCacheWriteback){
      CkCallback initCb(TraversalType::getCallDoneCallback(proxy()));
      CkCallback prodCb(CkCallback::ignore);
      completionDetectorProxy_.start_detection(CkNumPes(), 
                                               initCb, 
                                               prodCb, 
                                               CkCallbackResumeThread(), 0);
    }
    else{
      proxy().done(CkCallbackResumeThread());
    }
  }


  TraversalStatistics statistics(){
    CkReductionMsg *msg = NULL;
    proxy().statistics(CkCallbackResumeThread((void *&) msg));
    CkAssert(msg->getSize() == sizeof(TraversalStatistics));

    TraversalStatistics *stats = static_cast<TraversalStatistics *>(msg->getData());
    return *stats;
    delete msg;
  }

  // called by traversalmanager::instantiate
  CacheAccessType &nodeCacheAccessType(){
    return nodeCacheAccessType_;
  }

  CacheAccessType &leafCacheAccessType(){
    return leafCacheAccessType_;
  }

  // called by cacheclienttraversals when a 
  // traversal is added as a client for some cache
  const CacheAccessType &nodeCacheAccessType() const {
    return nodeCacheAccessType_;
  }

  const CacheAccessType &leafCacheAccessType() const {
    return leafCacheAccessType_;
  }

};

}; // namespace MultiphaseDistree

#endif // MULTIPHASE_DISTREE_TRAVERSAL_HANDLE_H
