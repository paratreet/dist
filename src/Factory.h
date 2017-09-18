#ifndef MULTIPHASE_DISTREE_FACTORY_H
#define MULTIPHASE_DISTREE_FACTORY_H

#include "mdt.decl.h"
#include "Node.h"
#include "Utility.h"
#include "Handle.h"

namespace MultiphaseDistree {

template<typename DataInterface>
class Factory {
  typedef Node<DataInterface> MyNodeType;
  typedef Piece<DataInterface> PieceType;
  typedef Manager<DataInterface> ManagerType;
  //typedef Traversal<DataInterface> TraversalType;

  typedef CProxy_Piece<DataInterface> PieceProxyType;
  typedef CProxy_Manager<DataInterface> ManagerProxyType;
  //typedef CProxy_Traversal<DataInterface> TraversalProxyType;

  typedef Handle<DataInterface> HandleType;

  public:
};


}; // namespace MultiphaseDistree

#endif // MULTIPHASE_DISTREE_FACTORY_H
