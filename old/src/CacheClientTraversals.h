#ifndef MULTIPHASE_DISTREE_CACHE_CLIENT_TRAVERSALS_H
#define MULTIPHASE_DISTREE_CACHE_CLIENT_TRAVERSALS_H

#include "mdt.decl.h"
#include "Node.h"
#include "Utility.h"
#include "CkCache.h"
#include "TraversalHandle.h"

namespace MultiphaseDistree {

class CacheClientTraversals {
  CkVec<CkGroupID> readonly_;
  CkVec<CkGroupID> accumulate_;

  public:

  template<typename TraversalType>
  void addReadonly(const TraversalHandle<TraversalType> &traversalHandle){
    readonly_.push_back(traversalHandle.proxy().ckLocMgr()->getGroupID());
  }

  template<typename TraversalType>
  void addAccumulate(const TraversalHandle<TraversalType> &traversalHandle){
    accumulate_.push_back(traversalHandle.proxy().ckLocMgr()->getGroupID());
  }

  void pup(PUP::er &p){
    p | readonly_;
    p | accumulate_;
  }

  const CkVec<CkGroupID> &readonly() const {
    return readonly_;
  }

  const CkVec<CkGroupID> &accumulate() const {
    return accumulate_;
  }
};

class NodeCacheClientTraversals : public CacheClientTraversals {
  public:
  template<typename TraversalType>
  void add(const TraversalHandle<TraversalType> &traversalHandle){
    if(traversalHandle.nodeCacheAccessType() == CacheAccessType::Readonly){
      addReadonly(traversalHandle);
    }
    else if(traversalHandle.nodeCacheAccessType() == CacheAccessType::Accumulate){
      addAccumulate(traversalHandle);
    }
    else{
      CkAbort("Bad TraversalHandle CacheAccessType\n");
    }
  }

  void pup(PUP::er &p){
    CacheClientTraversals::pup(p);
  }
};

class LeafCacheClientTraversals : public CacheClientTraversals {
  public:
  template<typename TraversalType>
  void add(const TraversalHandle<TraversalType> &traversalHandle){
    if(traversalHandle.leafCacheAccessType() == CacheAccessType::Readonly){
      addReadonly(traversalHandle);
    }
    else if(traversalHandle.leafCacheAccessType() == CacheAccessType::Accumulate){
      addAccumulate(traversalHandle);
    }
    else{
      CkAbort("Bad TraversalHandle CacheAccessType\n");
    }
  }

  void pup(PUP::er &p){
    CacheClientTraversals::pup(p);
  }
};



}; // namespace MultiphaseDistree

#endif // MULTIPHASE_DISTREE_CACHE_CLIENT_TRAVERSALS_H
