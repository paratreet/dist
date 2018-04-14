#ifndef MULTIPHASE_DISTREE_TREE_NODE_H
#define MULTIPHASE_DISTREE_TREE_NODE_H
#include <iostream>
#include <string>

namespace MultiphaseDistree {

enum NodeType {
  Invalid = 0,

  // All particles under the node are on this PE
  Internal,
  // This is a leaf
  Leaf,
  // This is a leaf with no particles in it
  EmptyLeaf,

  /* 
   * Some of the particles underneath this node are
   * owned by TreePieces that are not on this PE
  */
  Boundary,

  /*
   * None of the particles underneath this node are owned 
   * by TreePieces on this PE, although they may have been
   * fetched from a remote source.
   */
  Remote,
  // Same as above, but this is a leaf
  RemoteLeaf,
  // Same as above, but this is a leaf with no particles in it
  RemoteEmptyLeaf,

  // cached versions of above
  CachedInternal,
  CachedLeaf,
  CachedEmptyLeaf,
  CachedBoundary,
  // note that when a tree piece asks for the subtree
  // underneath a node, it can get nodes that are remote
  // to the tree piece that was queried. 
  CachedRemote,
  CachedRemoteLeaf,
  CachedRemoteEmptyLeaf
};

template<typename TreeDataInterface>
class Node {
  public:
  typedef typename TreeDataInterface::KeyType KeyType;
  typedef typename TreeDataInterface::PayloadType PayloadType;
  typedef typename TreeDataInterface::LocalParticleType ParticleType;

  private:
  KeyType key_;
  NodeType type_;
  int depth_;
  int ownerStart_;
  int ownerEnd_;
  int home_;

  ParticleType *particles_;
  int numParticles_;
 
  int numChildren_;
  Node *children_;
  Node *parent_;

  PayloadType payload_;
  int pending_;

  // useful while debugging
  // tracks the original tree piece that this node
  // belonged to, when it was submitted for merging 
  // -1 means that it was created by the merging PE
  // during the merge procedure. 0 and above gives the
  // original owner tree piece
  int tp_;

  // constructors
  public:
  Node(KeyType key, int depth, ParticleType *particles, int numParticles);
  Node();

  template<typename OtherTreeDataInterface>
  Node(const Node<OtherTreeDataInterface> &other);

  // get-set methods
  public:
  template<typename OtherTreeDataInterface>
  Node &operator=(const Node<OtherTreeDataInterface> &other);
  template<typename OtherTreeDataInterface>
  Node &operator=(Node<OtherTreeDataInterface> &other);

  template<typename OtherTreeDataInterface>
  void copyEverythingButParticles(const Node<OtherTreeDataInterface> &other);

  KeyType getKey() const;
  NodeType getType() const;
  int getDepth() const;
  int getOwnerStart() const;
  int getOwnerEnd() const;
  int &home();
  const int &home() const;
  ParticleType *getParticles();
  const ParticleType *getParticles() const;
  int getNumParticles() const;
  int getNumChildren() const;
  Node *getChildren();
  Node *getChild(int i);
  const Node *getChild(int i) const;
  const Node *getChildren() const;
  Node *getParent();
  const Node *getParent() const;
  PayloadType &getPayload();
  const PayloadType &getPayload() const;
  const int &pending() const;
  int &pending();
  const int &tp() const;
  int &tp();

  void setKey(KeyType k);
  void setType(NodeType t);
  void setDepth(int depth);
  void setOwnerStart(int start);
  void setOwnerEnd(int end);
  void setParticles(ParticleType *particles);
  void setNumParticles(int np);
  void setChildren(Node *children);
  void setNumChildren(int nc);
  void setParent(Node *parent);
  void setPayload(const PayloadType &data);

  void allocateChildren(int nChildren);
  void dot(std::ostream &oss) const;

  static NodeType makeRemote(NodeType type);

  // operations on tree nodes
  public:
  void deleteBeneath();

  private:
  // Initialize child fields based on those of parent.
  void initChild(int i, int *splitters, KeyType childKey, int childDepth);
  void reset();

  public:
  // to flatten to transport over network
  void serialize(Node *placeInBuf, Node *&emptyBuf, int subtreeDepth);
  // to revive pointers from flattened representation
  void deserialize(Node *start, int nn);

  // utilities
  public:
  static std::string TypeString(NodeType type){
    switch(type){
      case Invalid:                     return "Invalid";
      case Internal:                    return "Internal";
      case Leaf:                        return "Leaf";
      case EmptyLeaf:                   return "EmptyLeaf";
      case Boundary:                    return "Boundary";
      case Remote:                      return "Remote";
      case RemoteLeaf:                  return "RemoteLeaf";
      case RemoteEmptyLeaf:             return "RemoteEmptyLeaf";
      case CachedInternal:              return "CachedInternal";
      case CachedLeaf:                  return "CachedLeaf";
      case CachedEmptyLeaf:             return "CachedEmptyLeaf";
      case CachedBoundary:              return "CachedBoundary";
      case CachedRemote:                return "CachedRemote";
      case CachedRemoteLeaf:            return "CachedRemoteLeaf";
      case CachedRemoteEmptyLeaf:       return "CachedRemoteEmptyLeaf";
      default:                          return "BadNode";
    }
  }

  static std::string TypeDotColor(NodeType type){
    switch(type){
      case Invalid:                     return "firebrick1";
      case Internal:                    return "darkolivegreen1";
      case Leaf:                        return "darkolivegreen3";
      case EmptyLeaf:                   return "darksalmon";
      case Boundary:                    return "darkkhaki";
      case Remote:                      return "deepskyblue1";
      case RemoteLeaf:                  return "dodgerblue4";
      case RemoteEmptyLeaf:             return "deeppink";
      default:                          return "black";
    }
  }

  static bool IsInvalid(NodeType type){
    return type == Invalid;
  }

  static bool IsLeaf(NodeType type){
    return (type == Leaf || 
            type == EmptyLeaf || 
            type == RemoteLeaf || 
            type == RemoteEmptyLeaf || 
            type == CachedLeaf || 
            type == CachedEmptyLeaf ||
            type == CachedRemoteLeaf ||
            type == CachedRemoteEmptyLeaf);
  }

  static bool IsLocal(NodeType type){
    return (type == Leaf || type == EmptyLeaf || type == Internal);
  }

  static bool IsLocalLeaf(NodeType type){
    return (type == Leaf || type == EmptyLeaf);
  }

  static bool IsInternal(NodeType type){
    return (type == Internal);
  }

  static bool IsBoundary(NodeType type){
    return (type == Boundary);
  }

  static bool IsRemote(NodeType type){
    return (type == Remote || type == RemoteLeaf || type == RemoteEmptyLeaf);
  }

  static bool IsRemoteLeaf(NodeType type){
    return (type == RemoteLeaf || type == RemoteEmptyLeaf);
  }

  static bool IsCachedLeaf(NodeType type){
    return (type == CachedLeaf || type == CachedEmptyLeaf || type == CachedRemoteLeaf || type == CachedRemoteEmptyLeaf);
  }

  static bool IsCached(NodeType type){
    return (!IsInvalid(type) &&
            !IsLocal(type) &&
            !IsBoundary(type) &&
            !IsRemote(type));
  }

  static bool IsEmpty(NodeType type){
    return (type == EmptyLeaf || 
            type == RemoteEmptyLeaf || 
            type == CachedEmptyLeaf || 
            type == CachedRemoteEmptyLeaf);
  }

  static NodeType GetCached(NodeType type){
    switch(type){
      case Internal:            return CachedInternal;
      case Leaf:                return CachedLeaf;
      case EmptyLeaf:           return CachedEmptyLeaf;
      case Boundary:            return CachedBoundary;
      case Remote:              return CachedRemote;
      case RemoteLeaf:          return CachedRemoteLeaf;
      case RemoteEmptyLeaf:     return CachedRemoteEmptyLeaf;
      default:                  return Invalid;
    }
  }

};


template<typename TreeDataInterface>
Node<TreeDataInterface>::Node(KeyType key, int depth, ParticleType *particles, int numParticles) : 
  key_(key),
  depth_(depth),
  particles_(particles),
  numParticles_(numParticles),
  type_(Invalid),
  ownerStart_(-1),
  ownerEnd_(-1),
  home_(-1),
  numChildren_(0),
  children_(NULL),
  parent_(NULL),
  pending_(-1),
  tp_(-1)
{
}

template<typename TreeDataInterface>
Node<TreeDataInterface>::Node() : 
  key_(KeyType(0)),
  depth_(-1),
  particles_(NULL),
  numParticles_(-1),
  type_(Invalid),
  ownerStart_(-1),
  ownerEnd_(-1),
  home_(-1),
  numChildren_(0),
  children_(NULL),
  parent_(NULL),
  pending_(-1),
  tp_(-1)
{
}

// copy everything except for:
// 1. children
// 2. parent
template<typename TreeDataInterface>
template<typename OtherTreeDataInterface>
Node<TreeDataInterface>::Node(const Node<OtherTreeDataInterface> &other) : 
  children_(NULL),
  parent_(NULL)
{
  *this = other;
}

// we can't promise to do nothing to particles
// in const version, so don't copy them
template<typename TreeDataInterface>
template<typename OtherTreeDataInterface>
Node<TreeDataInterface> &Node<TreeDataInterface>::operator=(const Node<OtherTreeDataInterface> &other){
  copyEverythingButParticles(other);
  return *this;
}

// in non-const version, since we make no 
// promises about doing nothing to particles, 
// we copy everything, including particles
template<typename TreeDataInterface>
template<typename OtherTreeDataInterface>
Node<TreeDataInterface> &Node<TreeDataInterface>::operator=(Node<OtherTreeDataInterface> &other){
  copyEverythingButParticles(other);
  particles_ = other.getParticles();
  return *this;
}

template<typename TreeDataInterface>
template<typename OtherTreeDataInterface>
void Node<TreeDataInterface>::copyEverythingButParticles(const Node<OtherTreeDataInterface> &other){
  key_ = other.getKey();
  type_ = other.getType();
  depth_ = other.getDepth();
  ownerStart_ = other.getOwnerStart();
  ownerEnd_ = other.getOwnerEnd();

  // have to copy this, in case it is not -1
  // note that the home is -1 for local nodes, and
  // non-negative for remote ones (since 'home' is 
  // initialized to -1, and we have only reset it to
  // a non-negative value for remote nodes fetched
  // from their owners during tree building).
  // moreover, while unpacking, we must check whether
  // the home is >= 0, and if so retain its value
  home_ = other.home();

  //particles_ = other.getParticles();
  numParticles_ = other.getNumParticles();

  numChildren_ = other.getNumChildren();
  // user must define an appropriate
  // operator= for this to work
  payload_ = other.getPayload();

  pending_ = other.pending();
  tp_ = other.tp();
}




template<typename TreeDataInterface>
int Node<TreeDataInterface>::getNumParticles() const {
  return numParticles_;
}

template<typename TreeDataInterface>
typename TreeDataInterface::KeyType Node<TreeDataInterface>::getKey() const {
  return key_;
}

template<typename TreeDataInterface>
typename TreeDataInterface::LocalParticleType *Node<TreeDataInterface>::getParticles(){
  return particles_;
}

template<typename TreeDataInterface>
const typename TreeDataInterface::LocalParticleType *Node<TreeDataInterface>::getParticles() const {
  return particles_;
}



/*
template<typename TreeDataInterface>
int Node<TreeDataInterface>::getParticleStartIndex() const {
  return particleStartIndex_;
}
*/

template<typename TreeDataInterface>
Node<TreeDataInterface> *Node<TreeDataInterface>::getChild(int i){
  CkAssert(i < numChildren_);
  return children_ + i;
}

template<typename TreeDataInterface>
Node<TreeDataInterface> *Node<TreeDataInterface>::getParent(){
  return parent_;
}

template<typename TreeDataInterface>
const Node<TreeDataInterface> *Node<TreeDataInterface>::getParent() const {
  return parent_;
}

template<typename TreeDataInterface>
Node<TreeDataInterface> *Node<TreeDataInterface>::getChildren(){
  return children_;
}

template<typename TreeDataInterface>
int Node<TreeDataInterface>::getNumChildren() const {
  return numChildren_;
}

template<typename TreeDataInterface>
const Node<TreeDataInterface> *Node<TreeDataInterface>::getChild(int i) const {
  CkAssert(i < numChildren_);
  return children_ + i;
}

template<typename TreeDataInterface>
const Node<TreeDataInterface> *Node<TreeDataInterface>::getChildren() const {
  return children_;
}

template<typename TreeDataInterface>
NodeType Node<TreeDataInterface>::getType() const {
  return type_;
}

template<typename TreeDataInterface>
int Node<TreeDataInterface>::getDepth() const {
  return depth_;
}

template<typename TreeDataInterface>
int Node<TreeDataInterface>::getOwnerStart() const {
  return ownerStart_;
}

template<typename TreeDataInterface>
int Node<TreeDataInterface>::getOwnerEnd() const {
  return ownerEnd_;
}

template<typename TreeDataInterface>
int &Node<TreeDataInterface>::home(){
  return home_;
}

template<typename TreeDataInterface>
const int &Node<TreeDataInterface>::home() const {
  return home_;
}

template<typename TreeDataInterface>
typename TreeDataInterface::PayloadType &Node<TreeDataInterface>::getPayload(){
  return payload_;
}

template<typename TreeDataInterface>
const typename TreeDataInterface::PayloadType &Node<TreeDataInterface>::getPayload() const {
  return payload_;
}

template<typename TreeDataInterface>
const int &Node<TreeDataInterface>::pending() const {
  return pending_;
}

template<typename TreeDataInterface>
int &Node<TreeDataInterface>::pending(){
  return pending_;
}

template<typename TreeDataInterface>
const int &Node<TreeDataInterface>::tp() const {
  return tp_;
}

template<typename TreeDataInterface>
int &Node<TreeDataInterface>::tp(){
  return tp_;
}

template<typename TreeDataInterface>
void Node<TreeDataInterface>::setKey(KeyType key){
  key_ = key;
}

template<typename TreeDataInterface>
void Node<TreeDataInterface>::setParent(Node *parent){
  parent_ = parent;
}

template<typename TreeDataInterface>
void Node<TreeDataInterface>::setChildren(Node *children){
  children_ = children;
}

template<typename TreeDataInterface>
void Node<TreeDataInterface>::setNumChildren(int numChildren){
  numChildren_ = numChildren;
}

template<typename TreeDataInterface>
void Node<TreeDataInterface>::setType(NodeType type){
  type_ = type;
}

template<typename TreeDataInterface>
void Node<TreeDataInterface>::setOwnerStart(int ownerStart){
  ownerStart_ = ownerStart;
}

template<typename TreeDataInterface>
void Node<TreeDataInterface>::setOwnerEnd(int ownerEnd){
  ownerEnd_ = ownerEnd;
}

template<typename TreeDataInterface>
void Node<TreeDataInterface>::setNumParticles(int np){
  numParticles_ = np;
}

template<typename TreeDataInterface>
void Node<TreeDataInterface>::allocateChildren(int nChildren){
  children_ = new Node[nChildren];
  numChildren_ = nChildren;
}

template<typename TreeDataInterface>
void Node<TreeDataInterface>::deleteBeneath(){
  if(getChildren() == NULL){
    return;
  }

  CkAssert(getNumChildren() > 0);
  for(int i = 0; i < getNumChildren(); i++){
    getChild(i)->deleteBeneath();
  }

  delete[] getChildren();
  setChildren(NULL);
}



template<typename TreeDataInterface>
void Node<TreeDataInterface>::dot(std::ostream &out) const {
  out << key_ << " [";

  out << "label=\"";
  out << key_ << ", ";
  out << numParticles_ << ", ";
  out << "[" << ownerStart_ << ", " << ownerEnd_ << ")";
  out << "\\n" << payload_; 
  out << "\\n" << tp_; 
  out << "\",";

  out << "color=\"" << TypeDotColor(type_) << "\", ";
  out << "style=\"filled\"";

  out << "];" << std::endl;

  if(IsLocal(getType())) return;

  if(getChildren() == NULL) return;

  for(int i = 0; i < numChildren_; i++){
    const Node *child = getChild(i);
    out << key_ << " -> " << child->getKey() << ";" << std::endl; 
    child->dot(out);
  }
}

// gives user access to the leaves assigned to a 
// particular tree piece that submitted a tree for
// merging
template<typename TreeDataInterface>
class LeafIterator {
  typedef Node<TreeDataInterface> MyNodeType;

  CkVec<MyNodeType *> *allLeaves_;
  int firstLeafIndex_;
  int nLeaves_;

  int currentLeafOffset_;

  public:
  LeafIterator() : 
    allLeaves_(NULL),
    firstLeafIndex_(-1),
    nLeaves_(0),
    currentLeafOffset_(0)
  {}

  LeafIterator(CkVec<MyNodeType> *allLeaves, int firstLeafIndex) : 
    allLeaves_(allLeaves),
    firstLeafIndex_(firstLeafIndex),
    nLeaves_(0),
    currentLeafOffset_(0)
  {}

  void reset(){
    currentLeafOffset_ = 0;
  }

  MyNodeType *next(){
    if(currentLeafOffset_ == nLeaves()) return NULL;

    return (*allLeaves_)[firstLeafIndex_ + currentLeafOffset_++];
  }

  // should be used only by Manager to set number 
  // of leaves for tree piece
  int &nLeaves(){
    return nLeaves_;
  }
};

}; // namespace MultiphaseDistree

#endif // MULTIPHASE_DISTREE_TREE_NODE_H
