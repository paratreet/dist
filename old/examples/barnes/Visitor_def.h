#include "TreePiece.h"

template<typename AnyNodeType>
bool GravityVisitor::node(const AnyNodeType *n){
#ifdef DEBUG_TRAVERSAL
  tp()->preOpen(leaf_, n);
#endif
  bool doOpen = Physics::Gravity::open(leaf_, n);
  if(AnyNodeType::IsEmpty(n->getType()) || !doOpen){
    tp()->pn() += Physics::Gravity::forces(leaf_, n);
#ifdef DEBUG_TRAVERSAL
    tp()->insertComputedNode(leaf_->getKey(), n->getKey());
    leaf_->getUserData().pn()++;
#endif
    return false;
  }

  tp()->open()++;
#ifdef DEBUG_TRAVERSAL
  leaf_->getUserData().open()++;
#endif
  return true;
}

template<typename AnyNodeType>
bool DensitySphVisitor::node(const AnyNodeType *n){
  int status = -1;
  bool doOpen = Physics::Sph::open(leaf(), data(), n, status);
  if(AnyNodeType::IsEmpty(n->getType()) || !doOpen){
    return false;
  }

  tp()->sphOpen()++;
  return true;
}

template<typename AnyNodeType>
bool BallSphVisitor::node(const AnyNodeType *n){
  bool doOpen = Physics::BallSph::open(data(), n);
  //CkPrintf("[%d] BallSphVisitor::node leaf %llu source %llu open: %d\n", tp()->getIndex(), leaf()->getKey(), n->getKey(), doOpen);

  if(AnyNodeType::IsEmpty(n->getType()) || !doOpen){
    return false;
  }

  tp()->ballSphOpen()++;
  return true;
}



