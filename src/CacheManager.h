#ifndef MULTIPHASE_DISTREE_CACHE_MANAGER_H
#define MULTIPHASE_DISTREE_CACHE_MANAGER_H

#include "mdt.decl.h"
#include "CkCache.h"
#include "CacheClientTraversals.h"

namespace MultiphaseDistree {

template<typename TraversalDataInterface>
class CacheManager {
  public:
  typedef typename TraversalDataInterface::TreeBuildDataInterfaceType TreeBuildDataInterface;
  typedef CProxy_CkCacheManager<typename TreeBuildDataInterface::LocalNodeType::KeyType> CacheProxyType;

  class Handle {
    typename CacheManager::CacheProxyType proxy_;

    public:
    CacheProxyType &proxy() { return proxy_; }
    void pup(PUP::er &p){
      p | proxy_;
    }
  };

  static Handle instantiate(unsigned int size, const CacheClientTraversals &clients){
    Handle handle;
    handle.proxy() = CacheProxyType::ckNew(size, 
                                           clients.readonly().size(), 
                                           clients.readonly().getVec(),
                                           clients.accumulate().size(), 
                                           clients.accumulate().getVec()
                                           );
    return handle;
  }

};

};

#endif // MULTIPHASE_DISTREE_CACHE_MANAGER_H
