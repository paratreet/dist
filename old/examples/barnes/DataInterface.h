#ifndef DATA_INTERFACE_H
#define DATA_INTERFACE_H

#include "defines.h"
#include "NodePayload.h"
#include "mdt.h"
#include "pup.h"

// TREE DATA INTERFACES
// These are the interfaces from which we derive
// different types of nodes. Basically, these are
// passed as template parameters to 
// MultiphaseDistree::Node<> to get different types
// of nodes (i.e. Nodes with different payloads)

class TreeDataInterface {
  public:
  typedef Key KeyType;
  typedef NodePayload PayloadType;
  typedef Particle LocalParticleType;
};

// for now we have the same interface for 
// the gravity traversal tree as we do for the 
// built tree
typedef TreeDataInterface GravityRemoteDataInterface;

class RemoteBallSphDataInterface {
  public:
  typedef Key KeyType;
  typedef BallSphPayload PayloadType;
  typedef Particle LocalParticleType;
};

/*
class TreeNode : public MultiphaseDistree::Node<TreeDataInterface> {
  public:
  const Real &rsq() const {
    return getPayload().moments().rsq;
  }

  const Vector3D<Real> &com() const {
    return getPayload().moments().com;
  }

  const OrientedBox<Real> &box() const {
    return getPayload().moments().box;
  }

  const Real &mass() const {
    return getPayload().moments().mass;
  }
};

class SphBallRemoteNode : public MultiphaseDistree::Node<RemoteBallSphDataInterface> {
  public:
  SphBallRemoteNode(const TreeNode &other) : 
    MultiphaseDistree::Node<RemoteBallSphDataInterface>(other)
  {}

  const Real &rsq() const {
    return getPayload().rsq;
  }

  const Vector3D<Real> &com() const {
    return getPayload().com;
  }

  const OrientedBox<Real> &box() const {
    return getPayload().box;
  }

};
*/

typedef MultiphaseDistree::Node<TreeDataInterface> TreeNode;
typedef MultiphaseDistree::Node<RemoteBallSphDataInterface> SphBallRemoteNode;

// TreeBuild and Traversal DataInterfaces. These are
// the ones that are passed as template parameters to
// Piece, Manager (TreeBuildDataInterface) 
// and Traversals (TraversalDataInterface). The
// TreeBuildDataInterface provides the local node
// type used by the program, and the various traversal
// data inferfaces give the remote node types and
// remote node particle types.

class TreeBuildDataInterface {
  public:
  typedef TreeNode LocalNodeType; 
  typedef TreeDataInterface::KeyType KeyType;
  typedef TreeDataInterface::PayloadType PayloadType;
  typedef TreeDataInterface::LocalParticleType LocalParticleType;


  static const Real &rsq(const LocalNodeType *n){
    return n->getPayload().moments().rsq;
  }

  static const Vector3D<Real> &com(const LocalNodeType *n){
    return n->getPayload().moments().com;
  }

  static const OrientedBox<Real> &box(const LocalNodeType *n){
    return n->getPayload().moments().box;
  }

  static const Real &mass(const LocalNodeType *n){
    return n->getPayload().moments().mass;
  }

};

class GravityTraversalDataInterface {
  public: 
  typedef TreeBuildDataInterface TreeBuildDataInterfaceType;
  typedef TreeBuildDataInterface::LocalNodeType LocalNodeType;
  typedef LocalNodeType RemoteNodeType;
  typedef GravityRemoteParticle RemoteParticleType;

  static int NodeRequestPriority;
  static int LeafRequestPriority;
  static int NodeReceivePriority;
  static int LeafReceivePriority;

  public:
  static void copy(RemoteNodeType &to, const TreeBuildDataInterface::LocalNodeType &from){
    to = from;
  }

  static void copy(RemoteParticleType &to, const TreeBuildDataInterface::LocalParticleType &from){
    to = from;
  }

  static void accumulate(TreeBuildDataInterface::LocalNodeType &to, const RemoteNodeType &from){
  }

  static void accumulate(TreeBuildDataInterface::LocalParticleType &to, const RemoteParticleType &from){
  }

  // for gravity traversal, home and away nodes are
  // the same. therefore, we can have the same 
  // field-fetcher functions for both types of 
  // traversal
  template<typename AnyNodeType>
  static const Real &rsq(const AnyNodeType *n){
    return TreeBuildDataInterfaceType::rsq(n); 
  }

  template<typename AnyNodeType>
  static const Vector3D<Real> &com(const AnyNodeType *n){
    return TreeBuildDataInterfaceType::com(n); 
  }

  template<typename AnyNodeType>
  static const OrientedBox<Real> &box(const AnyNodeType *n){
    return TreeBuildDataInterfaceType::box(n); 
  }

  template<typename AnyNodeType>
  static const Real &mass(const AnyNodeType *n){
    return TreeBuildDataInterfaceType::mass(n); 
  }
};

// density is calculated together with gravity, and
// since sph density traversal and gravity traversal
// share caches, their interfaces must have the same
// underlying data types.
typedef GravityTraversalDataInterface SphDensityTraversalDataInterface;

class SphBallTraversalDataInterface {
  public:
  typedef TreeBuildDataInterface TreeBuildDataInterfaceType;
  typedef TreeBuildDataInterface::LocalNodeType LocalNodeType;
  typedef SphBallRemoteNode RemoteNodeType;
  typedef BallSphRemoteParticle RemoteParticleType;

  static int NodeRequestPriority;
  static int LeafRequestPriority;
  static int NodeReceivePriority;
  static int LeafReceivePriority;

  public:
  static void copy(RemoteNodeType &to, const TreeBuildDataInterface::LocalNodeType &from){
    to = from;
  }

  static void copy(RemoteParticleType &to, const TreeBuildDataInterface::LocalParticleType &from){
    to = from;
  }

  static void accumulate(TreeBuildDataInterface::LocalNodeType &to, const RemoteNodeType &from){
  }

  static void accumulate(TreeBuildDataInterface::LocalParticleType &to, const RemoteParticleType &from){
    // FIXME - do some accumulation here
  }

  // for ball sph traversal, the away node type is actually
  // different from the home node type, and consequently, we
  // have different property fetchers for the two types:
  static const Real &rsq(const LocalNodeType *n){
    return TreeBuildDataInterfaceType::rsq(n);
  }

  static const Vector3D<Real> &com(const LocalNodeType *n){
    return TreeBuildDataInterfaceType::com(n);
  }

  static const OrientedBox<Real> &box(const LocalNodeType *n){
    return TreeBuildDataInterfaceType::box(n);
  }

  static const Real &rsq(const RemoteNodeType *n){
    return n->getPayload().rsq;
  }

  static const Vector3D<Real> &com(const RemoteNodeType *n){
    return n->getPayload().com;
  }

  static const OrientedBox<Real> &box(const RemoteNodeType *n){
    return n->getPayload().box;
  }
};

// for sph, there are two walks: density, which happens
// concurrently with gravity, and a ball walk, which
// happens afterwards. 
// after the density walk, we have density and pressure
// values at every point. there is no writeback for this
// phase, i.e. only a readonly cache is used.

// next is the ball sph walk, which uses the ball 
// defined for each particle by the previous traversal,
// and also the pressure and density values (for
// pressure gradient) and density and velocity values
// (for the viscosity term). There is only one 
// traversal performed, the ball sph traversal. At the
// end of this traversal, we invoke two kernels,
// one for pressure gradient calculation, and the other
// for viscosity calculation.

#endif // DATA_INTERFACE_H
