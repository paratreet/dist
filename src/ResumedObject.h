#ifndef MULTIPHASE_DISTREE_RESUMED_OBJECT_H
#define MULTIPHASE_DISTREE_RESUMED_OBJECT_H

#include "mdt.decl.h"
#include "Node.h"
#include "Utility.h"
#include "CkCache.h"

namespace MultiphaseDistree {

#define GetProxy1D(UserPieceType, objId) CProxyElement_#UserPieceType(objId.arrayId, objId.index.data()[0])
#define GetProxy2D(UserPieceType, objId) CProxyElement_#UserPieceType(objId.arrayId, objId.index.data()[0], objId.index.data()[1])
#define GetProxy3D(UserPieceType, objId) CProxyElement_#UserPieceType(objId.arrayId, objId.index.data()[0], objId.index.data()[1], objId.index.data()[2])

struct ResumedObject {
  CkArrayID arrayId;
  CkArrayIndex index;

  ResumedObject(const CkArrayID &aid, const CkArrayIndex &idx) : 
    arrayId(aid),
    index(idx)
  {}
};

}; // namespace MultiphaseDistree


#endif // MULTIPHASE_DISTREE_RESUMED_OBJECT_H
