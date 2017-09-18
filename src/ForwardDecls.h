#ifndef MULTIPHASE_DISTREE_FORWARD_DECLS_H
#define MULTIPHASE_DISTREE_FORWARD_DECLS_H

#include "PointerContainer.h"

namespace MultiphaseDistree {

// these are needed because charmxi doesn't like
// the keyword 'template' preceding "dependenent scope"
// e.g. 'typename DataInterface::LocalNodeType::KeyType'
#define DATA_INTERFACE_KEY_TYPE_MACRO typename DataInterface::KeyType
#define DATA_INTERFACE_PAYLOAD_TYPE_MACRO typename DataInterface::PayloadType
#define DATA_INTERFACE_NODE_TYPE_MACRO typename DataInterface::LocalNodeType
#define DATA_INTERFACE_PIECE_TYPE_MACRO Piece<DataInterface>

#define TRAVERSAL_DATA_INTERFACE_TREE_BUILD_DATA_INTERFACE_TYPE_MACRO typename TraversalDataInterface::TreeBuildDataInterfaceType
#define TRAVERSAL_DATA_INTERFACE_KEY_TYPE_MACRO typename TraversalDataInterface::TreeBuildDataInterfaceType::KeyType

template <typename DataInterface> class Handle;

};

#endif // MULTIPHASE_DISTREE_FORWARD_DECLS_H
