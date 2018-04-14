#ifndef NODE_PAYLOAD_H
#define NODE_PAYLOAD_H

#include <iostream>

#include "MultipoleMoments.h"

struct NodePayload {
  MultipoleMoments moments_;
#ifdef DEBUG_TRAVERSAL
  CmiUInt8 open_;
  CmiUInt8 pp_;
  CmiUInt8 pn_;
#endif

  NodePayload() : 
    moments_()
  {
#ifdef DEBUG_TRAVERSAL
    open() = 0;
    pp() = 0;
    pn() = 0;
#endif
  }

  MultipoleMoments &moments(){
    return moments_;
  }

  const MultipoleMoments &moments() const {
    return moments_;
  }

  void mergeStart(){
    moments_.begin();
  }

  void mergeAccumulate(const NodePayload &child){
    moments_.accumulate(child.moments_);
  }

  void mergeFinish(){
    moments_.end();
  }

#ifdef DEBUG_TRAVERSAL
  CmiUInt8 &open(){
    return open_;
  }

  CmiUInt8 &pn(){
    return pn_;
  }

  CmiUInt8 &pp(){
    return pp_;
  }

  const CmiUInt8 &open() const {
    return open_;
  }

  const CmiUInt8 &pn() const {
    return pn_;
  }

  const CmiUInt8 &pp() const {
    return pp_;
  }
#endif

  void pup(PUP::er &p);
};

std::ostream &operator<<(std::ostream &out, const NodePayload &pl);


struct BallSphPayload {
  // no need for node mass in sph traversal
  Real rsq;
  Vector3D<Real> com;
  OrientedBox<Real> box;

  BallSphPayload() {}
  BallSphPayload(const NodePayload &other){
    *this = other;
  }

  BallSphPayload &operator=(const NodePayload &other){
    rsq = other.moments().rsq;
    com = other.moments().com;
    box = other.moments().box;
    return *this;
  }

  void pup(PUP::er &p);
};

#endif // NODE_PAYLOAD_H

