#ifndef RESUMPTION_H
#define RESUMPTION_H

#include "PointerContainer.h"
using namespace MultiphaseDistree;

template<typename VisitorType>
struct NodeResumption {
  PointerContainer<typename VisitorType::RemoteNodeType> node;
  int chunk;
  PointerContainer<VisitorType> visitor;

  NodeResumption(typename VisitorType::RemoteNodeType *n, int c, VisitorType *v) : 
    node(n),
    chunk(c),
    visitor(v)
  {}

  void pup(PUP::er &p){
    p | node;
    p | chunk;
    p | visitor;
  }
};

template<typename VisitorType>
struct LeafResumption {
  Key key;
  PointerContainer<typename VisitorType::RemoteParticleType> particles;
  int nParticles;
  int chunk;
  PointerContainer<VisitorType> visitor;

  LeafResumption(Key k, typename VisitorType::RemoteParticleType *p, int n, int c, VisitorType *v) : 
    key(k),
    particles(p),
    nParticles(n),
    chunk(c),
    visitor(v)
  {}

  void pup(PUP::er &p){
    p | key;
    p | particles;
    p | nParticles;
    p | chunk;
    p | visitor;
  }
};

#endif // RESUMPTION_H
