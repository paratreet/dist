#ifndef MULTIPHASE_DISTREE_CACHED_DATA_H
#define MULTIPHASE_DISTREE_CACHED_DATA_H

#include "mdt.decl.h"
#include "Node.h"
#include "Utility.h"
#include "CkCache.h"

namespace MultiphaseDistree {

template<typename TraversalDataInterface>
struct CachedLeafContents {
  CkCacheFillMsg<typename TraversalDataInterface::TreeBuildDataInterfaceType::LocalNodeType::KeyType> *msg;
  int nParticles;
  typename TraversalDataInterface::RemoteParticleType particles[1];
};

// Hack courtesy Filippo Gioachin (ChaNGa); since 
// last field must be aligned to sizeof(MyNodeType)
// boundary, if the CachedNodes struct adjoins
// a run of MyNodeType elements, then we can 
// overrun the 'nodes' array and start to address 
// elements in adjoining MyNodeType array.
// However, this works only if the compiler lays out
// the 'nodes' field as the last in the run of bytes
// allocated for CachedNodes.

template<typename TraversalDataInterface>
struct CachedNodes {
  CkCacheFillMsg<typename TraversalDataInterface::TreeBuildDataInterfaceType::LocalNodeType::KeyType> *msg;
  int nNodes;
  typename TraversalDataInterface::RemoteNodeType nodes[1];
};

}; // namespace MultiphaseDistree


#endif // MULTIPHASE_DISTREE_CACHED_DATA_H
