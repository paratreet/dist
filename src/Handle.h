#ifndef MULTIPHASE_DISTREE_HANDLE_H
#define MULTIPHASE_DISTREE_HANDLE_H

#include "mdt.decl.h"
#include "CkCache.h"

#include "completion.h"

namespace MultiphaseDistree {

struct ManagerQuiescence {
  int numTpAcks;
  int mergeDone;
  int registrations;
};

/**
 A Handle holds proxy objects for a (Tree)Piece, a Manager, 
 and a completion detector.
*/
template<typename DataInterface> 
struct Handle {
  typedef typename DataInterface::LocalNodeType LocalNodeType;
  typedef CProxy_Piece<DataInterface> PieceProxyType;
  typedef CProxy_Manager<DataInterface> ManagerProxyType;
  typedef Manager<DataInterface> ManagerType;

  typedef CProxy_CkCacheManager<typename LocalNodeType::KeyType> CacheProxyType;
  typedef CProxy_CompletionDetector CompletionDetectorProxyType;


  PieceProxyType pieceProxy_;
  ManagerProxyType managerProxy_;
  CompletionDetectorProxyType completionDetectorProxy_;

  // const versions
  const PieceProxyType &pieceProxy() const {return pieceProxy_;}
  const ManagerProxyType &managerProxy() const {return managerProxy_;}
  const CompletionDetectorProxyType &completionDetectorProxy() const { return completionDetectorProxy_; }

  // non-const versions
  PieceProxyType &pieceProxy() {return pieceProxy_;}
  ManagerProxyType &managerProxy() {return managerProxy_;}
  CompletionDetectorProxyType &completionDetectorProxy() { return completionDetectorProxy_; }

  void initialize();
  void registration(LocalNodeType *root, Piece<DataInterface> *piece);
  void syncToMerge(const CkCallback &cb);
  void syncToDelete(const CkCallback &cb);

  ManagerType *manager() { return managerProxy().ckLocalBranch(); }
  LocalNodeType *root() { return manager()->root(); }

  void pup(PUP::er &p){
    p | pieceProxy_;
    p | managerProxy_;
    p | completionDetectorProxy_;
  }

  // factory method
  public:
  static Handle instantiate(const PieceProxyType &proxy, int cacheChunkDepth){
    Handle handle;
    handle.pieceProxy() = proxy;
    handle.managerProxy() = ManagerType::instantiate(cacheChunkDepth);
    handle.completionDetectorProxy() = CompletionDetectorProxyType::ckNew();
    return handle;
  }
};

// this method should be called from a [threaded] entry method

template<typename DataInterface>
void Handle<DataInterface>::
initialize(){
  managerProxy().initialize(*this, CkCallbackResumeThread());
  pieceProxy().initialize(*this, CkCallbackResumeThread());
  // FIXME - might have to initialize cache/register pieces 
  // with it when new version of cache is rolled out
}

template<typename DataInterface>
void Handle<DataInterface>::
registration(LocalNodeType *node, Piece<DataInterface> *piece){
  managerProxy().ckLocalBranch()->registration(node, piece);
}

template<typename DataInterface>
void Handle<DataInterface>::syncToMerge(const CkCallback &cb){
  CkStartQD(CkCallback(CkIndex_Manager<DataInterface>::quiescence(), managerProxy()));
  managerProxy().syncToMerge(cb);
}

template<typename DataInterface> 
void Handle<DataInterface>::syncToDelete(const CkCallback &cb){
  managerProxy().syncToDelete(cb);
}

}; // end namespace MultiphaseDistree

#endif // MULTIPHASE_DISTREE_HANDLE_H
