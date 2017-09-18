#ifndef MULTIPHASE_DISTREE_REDUCE_H
#define MULTIPHASE_DISTREE_REDUCE_H

#include "mdt.decl.h"
#include "Node.h"
#include "Utility.h"
#include "CkCache.h"
#include "Handle.h"

namespace MultiphaseDistree {

template<typename VisitorType>
class Reduce {
  static CkReduction::reducerType VisitorReductionType;

  public:
  // to be called by user in initnode function
  static void registration(){
    VisitorReductionType = CkReduction::addReducer(VisitorReducer);
  }

  static CkReductionMsg *VisitorReducer(int nMsgs, CkReductionMsg **msgs) {
    VisitorType* v = static_cast<VisitorType *>(msgs[0]->getData());
    VisitorType* msgv;
    for(int i = 1; i < nMsgs; i++) {
      msgv = static_cast<VisitorType *>(msgs[i]->getData());
      *v += *msgv;
    }

    return CkReductionMsg::buildNew(sizeof(VisitorType), v);
  }

  static CkReduction::reducerType reducer(){
    return VisitorReductionType;
  }
};

template<typename VisitorType>
CkReduction::reducerType Reduce<VisitorType>::VisitorReductionType;

}; // namespace MultiphaseDistree


#endif // MULTIPHASE_DISTREE_CACHE_REDUCE_H
