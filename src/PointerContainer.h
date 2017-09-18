#ifndef MULTIPHASE_DISTREE_POINTER_CONTAINER_H
#define MULTIPHASE_DISTREE_POINTER_CONTAINER_H

#include "charm++.h"

namespace MultiphaseDistree {

template<typename PointedType>
struct PointerContainer {
  PointedType *data;

  PointerContainer() : 
    data(0)
  {}

  PointerContainer(PointedType *d) : 
    data(d)
  {}

  void pup(PUP::er &p){
    p((char *) this, sizeof(PointerContainer<PointedType>));
  }
};

}; // namespace MultiphaseDistree

#endif // MULTIPHASE_DISTREE_POINTER_CONTAINER_H
