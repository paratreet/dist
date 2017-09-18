#ifndef ORB3DLB_H
#define ORB3DLB_H

#include "CentralLB.h"
#include "Orb3dLB.decl.h"
#include "TaggedVector3D.h"
#include "OrientedBox.h"

#define XDIM 0
#define YDIM 1
#define ZDIM 2
#define NDIMS 3

#define LEFT_PARTITION 0
#define RIGHT_PARTITION 1
#define INVALID_PARTITION -1

void CreateOrb3dLB();
BaseLB * AllocateOrb3dLB();

struct OrbObject {
  int partition;
  int lbindex;
  Vector3D<float> centroid;
  int numParticles;
  OrbObject() : partition(-1), lbindex(-1), numParticles(0) {}
  OrbObject(int tag, int np) : partition(-1), lbindex(tag), numParticles(np) {}
};

struct Event {
  int owner;
  float load;
  float position;

  Event(float pos, float ld, int o) : 
    position(pos),
    load(ld),
    owner(o)
  {
  }

  Event() : 
    owner(-1),
    load(0.0),
    position(0.0)
  {
  }

  bool operator<=(const Event &e){
    return position <= e.position;
  }

  bool operator>=(const Event &e){
    return position >= e.position;
  }
};



class Orb3dLB : public CentralLB {
private:
  CkVec<int> *mapping;
  CkVec<int> *from;

  CkVec<OrbObject> tps;
  CkVec<float> procload;
  CkVec<OrientedBox<double> > procbox;

  // things are stored in here before work
  // is ever called.
  TaggedVector3D *tpCentroids;
  CkReductionMsg *tpmsg;
  int nrecvd;
  bool haveTPCentroids;

  int nextProc;

  int numRefinementsDone_;

  CkCallback callback_;

  bool QueryBalanceNow(int step);
  void printData(BaseLB::LDStats &stats, int phase, int *revObjMap);
  void orbPartition(CkVec<Event> *events, OrientedBox<double> &box, int procs);
  int partitionRatioLoad(CkVec<Event> &events, float ratio);

public:
  Orb3dLB(const CkLBOptions &);
  Orb3dLB(CkMigrateMessage *m):CentralLB(m) { lbname = "Orb3dLB"; }
  
  void start(const CkCallback &cbDoneReceiveCentroids, const CkCallback &cb);
  void work(BaseLB::LDStats* stats);
  void receiveCentroids(CkReductionMsg *msg);
};

#endif /* ORB3DLB_H */
