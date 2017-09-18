// User must define the data interface classes of the 
// following pattern. Framework requires one for the tree
// structure itself (used in tree building), and one per
// traversal

class DataInterface {
  public:
  typedef X1 LocalNodeType;
  typedef X2 RemoteNodeType;
  typedef X3 LocalParticleType;
  typedef X4 RemoteParticleType;
  typedef X5 KeyType;

  static RemoteNodeType &extract(const LocalNodeType &other);
  static RemoteParticleType &extract(const LocalParticleType &other);
};

