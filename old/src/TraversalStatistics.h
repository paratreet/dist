#ifndef MULTIPHASE_DISTREE_TRAVERSAL_STATISTICS_H
#define MULTIPHASE_DISTREE_TRAVERSAL_STATISTICS_H

namespace MultiphaseDistree {

class TraversalStatistics {
  typedef CmiUInt8 CounterType;

  CounterType nNodeCacheMisses_;
  CounterType nLeafCacheMisses_;

  public:
  TraversalStatistics(){
    reset();
  }

  void pup(PUP::er &p){
    p | nNodeCacheMisses_;
    p | nLeafCacheMisses_;
  }

  void nodeCacheMiss(){
    nNodeCacheMisses_++;
  }

  void leafCacheMiss(){
    nLeafCacheMisses_++;
  }

  const CounterType &nNodeCacheMisses() const {
    return nNodeCacheMisses_;
  }

  const CounterType &nLeafCacheMisses() const {
    return nLeafCacheMisses_;
  }


  TraversalStatistics &operator+=(const TraversalStatistics &other) {
    nNodeCacheMisses_ += other.nNodeCacheMisses_;
    nLeafCacheMisses_ += other.nLeafCacheMisses_;
    return *this;
  }

  void reset(){
    nNodeCacheMisses_ = 0;
    nLeafCacheMisses_ = 0;
  }

};


std::ostream &operator<<(std::ostream &out, const TraversalStatistics &stats);

}; // namespace MultiphaseDistree

#endif // MULTIPHASE_DISTREE_TRAVERSAL_STATISTICS_H
