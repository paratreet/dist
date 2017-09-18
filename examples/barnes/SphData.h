#ifndef SPH_DATA_H
#define SPH_DATA_H

#include <set>
#include "defines.h"
#include "Vector3D.h"
#include "MyPriorityQueue.h"
#include "TraversalData.h"
#include "DataInterface.h"

using namespace std;


class SphParticle {
  Real dist2_;
  const SphCore *core_;

  public:
  SphParticle() : 
    core_(NULL)
  {}

  
  SphParticle(Real dist2, const SphCore *core) : 
    dist2_(dist2),
    core_(core)
  {
  }

  Real &dist2() { return dist2_; }

  const Real &dist2() const { return dist2_; }
  const SphCore *core() const { return core_; }

  bool operator<(const SphParticle &other) const {
    return dist2() < other.dist2();
  }
};

class SphParticleData {
  public:
  typedef MyPriorityQueue<SphParticle> QueueType;
  typedef QueueType::StorageType StorageType;

  private:
  QueueType neighbors_;
  // this is only required for the ball-sph walk
  Real rCutoff2_;

  public:
  void addNeighbor(SphParticle &neighbor);
  QueueType &neighbors() { return neighbors_; }
  const QueueType &neighbors() const { return neighbors_; }

  Real &maxDist2() {
    return neighbors().top().dist2();
  }

  int numNeighbors(){
    return neighbors().size();
  }

  void addWithReplacement(const SphParticle &p){
    neighbors().pop();
    neighbors().push(p);
    rCutoff2() = maxDist2();
  }

  void add(const SphParticle &p){
    neighbors().push(p);
    rCutoff2() = maxDist2();
  }

  // used in ball sph walk: doesn't depend on
  // current list of neighbors, but is determined
  // by preceding sph walk. all particles within
  // the cutoff sphere are included as neighbors
  Real &rCutoff2(){
    return rCutoff2_;
  }

  void free(){
    neighbors_.free();
  }
};

class SphLeafData {
  int leafNum_;
  Real smoothSize_;
  Vector3D<Real> smoothCenter_;
  Real maxDist_;

  int nParticles_;
  CkVec<SphParticleData> particleData_;
  int nOutstanding_;

  public:

  SphLeafData() : 
    nParticles_(0),
    maxDist_(REAL_MAX),
    nOutstanding_(0),
    leafNum_(-1)
  {}


  void initialize(int nParticles, int leafNum){
    nParticles_ = nParticles;
    maxDist_ = REAL_MAX;
    nOutstanding_ = 0;
    leafNum_ = leafNum;
    particleData_.resize(nParticles);
    free();
  }

  Real &size() { return smoothSize_; }
  Vector3D<Real> &center() { return smoothCenter_; }
  Real &maxDist() { return maxDist_; }
  int &nParticles() { return nParticles_; }
  const int &nParticles() const { return nParticles_; }
  SphParticleData &particleData(int i) { return particleData_[i]; }
  const SphParticleData &particleData(int i) const { return particleData_[i]; }
  int &nOutstanding() { return nOutstanding_; }

  int getLeafNum() const {
    return leafNum_;
  }

  void free(){
    for(int i = 0; i < particleData_.size(); i++){
      particleData_[i].free();
    }
  }

};

#endif // SPH_DATA_H
