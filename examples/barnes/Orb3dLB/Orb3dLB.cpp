#include <charm++.h>
#include "Orb3dLB.h"
#include "Refiner.h"
#include "Vector3D.h"

#define ORB3DLB_DEBUG 

using namespace std;

CreateLBFunc_Def(Orb3dLB, "3d ORB mapping of tree piece space onto linear set of processors");

Orb3dLB::Orb3dLB(const CkLBOptions &opt): CentralLB(opt)
{
  lbname = "Orb3dLB";
  if (CkMyPe() == 0){
    CkPrintf("[%d] Orb3dLB created\n",CkMyPe());
  }
  haveTPCentroids = false;
}

void Orb3dLB::start(const CkCallback &cbDoneReceiveCentroids, const CkCallback &cb){
  callback_ = cbDoneReceiveCentroids;
  cb.send();
}

void Orb3dLB::receiveCentroids(CkReductionMsg *msg){
  if(haveTPCentroids){
    delete tpmsg;
  }
  tpCentroids = (TaggedVector3D *)msg->getData();
  nrecvd = msg->getSize()/sizeof(TaggedVector3D);
  tpmsg = msg;
  haveTPCentroids = true;
  CkPrintf("Orb3dLB: receiveCentroids %d elements, msg length: %d\n", nrecvd, msg->getLength()); 
  callback_.send();
}

bool Orb3dLB::QueryBalanceNow(int step){
  return true;
}

void Orb3dLB::work(BaseLB::LDStats* stats)
{
  int numobjs = stats->n_objs;
  int nmig = stats->n_migrateobjs;

  stats->makeCommHash();
  CkAssert(nrecvd == nmig);

  CkVec<Event> tpEvents[NDIMS];
  for(int i = 0; i < NDIMS; i++){
    tpEvents[i].reserve(nmig);
  }
  tps.resize(nmig);

  OrientedBox<double> box;

  int numProcessed = 0;

  for(int i = 0; i < nmig; i++){
    TaggedVector3D *data = tpCentroids+i;
    LDObjHandle &handle = data->handle;
    int tag = stats->getHash(handle.id, handle.omhandle.id);

    float load = stats->objData[tag].wallTime;
    //float load = (float) data->numActiveParticles;

    tps[i] = OrbObject(tag, data->myNumParticles);
    tps[i].centroid = data->vec;
    
    tpEvents[XDIM].push_back(Event(data->vec.x, load, i));
    tpEvents[YDIM].push_back(Event(data->vec.y, load, i));
    tpEvents[ZDIM].push_back(Event(data->vec.z, load, i));

      /*
    CkPrintf("[Orb3dLB_notopo] tree piece %d particles %d load %f centroid %f %f %f\n", 
                                      data->tag,
                                      data->myNumParticles,
                                      load,
                                      data->vec.x,
                                      data->vec.y,
                                      data->vec.z
                                      );
    */
    /*
    if(step() == 1){
      CkPrintf("[tpinfo] %f %f %f %f %d %d\n",
          tp[tag].centroid.x,
          tp[tag].centroid.y,
          tp[tag].centroid.z,
          tp[tag].load,
          tpCentroids[i].numActiveParticles,
          tp[tag].lbindex
          );
    }
    */

    numProcessed++;
  }

  CkAssert(numProcessed == nmig);
  CkAssert(tpEvents[XDIM].length() == nmig);
  CkAssert(tpEvents[YDIM].length() == nmig);
  CkAssert(tpEvents[ZDIM].length() == nmig);

  mapping = &stats->to_proc;
  from = &stats->from_proc;
  int dim = 0;

  CkPrintf("[Orb3dLB] sorting\n");
  for(int i = 0; i < NDIMS; i++){
    tpEvents[i].quickSort();
  }

  box.lesser_corner.x = tpEvents[XDIM][0].position;
  box.lesser_corner.y = tpEvents[YDIM][0].position;
  box.lesser_corner.z = tpEvents[ZDIM][0].position;

  box.greater_corner.x = tpEvents[XDIM][nmig-1].position;
  box.greater_corner.y = tpEvents[YDIM][nmig-1].position;
  box.greater_corner.z = tpEvents[ZDIM][nmig-1].position;

  nextProc = 0;

  procload.resize(stats->count);
  procbox.resize(stats->count);
  for(int i = 0; i < stats->count; i++){
    procload[i] = 0.0;
  }

  orbPartition(tpEvents, box, stats->count);

  //int *from_procs = Refiner::AllocProcs(stats->count, stats);
  int migr = 0;
  for(int i = 0; i < numobjs; i++){
    if(stats->to_proc[i] != stats->from_proc[i]) migr++;
    //int pe = stats->to_proc[i];
    //from_procs[i] = pe;
  }

  /*
  int *to_procs = Refiner::AllocProcs(stats->count, stats);

  Refiner refiner(1.010);

  refiner.Refine(stats->count,stats,from_procs,to_procs);

  int numRefineMigrated = 0;
  for(int i = 0; i < numobjs; i++){
    if(to_procs[i] != from_procs[i]) numRefineMigrated++;
    stats->to_proc[i] = to_procs[i];
  }
  */


  double minWall = 0.0;
  double maxWall = 0.0;
  double avgWall = 0.0;

  double minObj = 0.0;
  double maxObj = 0.0;
  double avgObj = 0.0;

  CkPrintf("***************************\n");
  CkPrintf("Before LB step %d\n", step());
  CkPrintf("***************************\n");
  CkPrintf("i pe wall idle bg_wall objload\n");
  for(int i = 0; i < stats->count; i++){
    double wallTime = stats->procs[i].total_walltime;
    double idleTime = stats->procs[i].idletime;
    double bgTime = stats->procs[i].bg_walltime;
    double objTime = wallTime-(idleTime+bgTime);
    /*
    CkPrintf("[pestats] %d %d %f %f %f %f\n", 
        i,
        stats->procs[i].pe, 
        wallTime,
        idleTime,
        bgTime,
        objTime);
        */

    avgWall += wallTime; 
    avgObj += objTime; 

    if(i==0 || minWall > wallTime) minWall = wallTime;
    if(i==0 || maxWall < wallTime) maxWall = wallTime;

    if(i==0 || minObj > objTime) minObj = objTime;
    if(i==0 || maxObj < objTime) maxObj = objTime;

  }

  avgWall /= stats->count;
  avgObj /= stats->count;

  float minload, maxload, avgload;
  minload = maxload = procload[0];
  avgload = 0.0;
  for(int i = 0; i < stats->count; i++){
    /*
    CkPrintf("pe %d load %f box %f %f %f %f %f %f\n", i, procload[i], 
                                procbox[i].lesser_corner.x,
                                procbox[i].lesser_corner.y,
                                procbox[i].lesser_corner.z,
                                procbox[i].greater_corner.x,
                                procbox[i].greater_corner.y,
                                procbox[i].greater_corner.z
                                );
    */
    avgload += procload[i];
    if(minload > procload[i]) minload = procload[i];
    if(maxload < procload[i]) maxload = procload[i];
  }

  avgload /= stats->count;

 
  CkPrintf("Orb3dLB stats: minWall %f maxWall %f avgWall %f maxWall/avgWall %f\n", minWall, maxWall, avgWall, maxWall/avgWall);
  //CkPrintf("Orb3dLB stats: minObj %f maxObj %f avgObj %f maxObj/avgObj %f\n", minObj, maxObj, avgObj, maxObj/avgObj);
  CkPrintf("Orb3dLB stats: orb migrated %d objects\n", migr);
  CkPrintf("Orb3dLB stats: EXPECTED: \n");
  CkPrintf("Orb3dLB stats: min %f max %f avg %f max/avg %f\n", minload, maxload, avgload, maxload/avgload);
  //CkPrintf("Orb3dLB stats: orb migrated %d refine migrated %d objects\n", migr, numRefineMigrated);

  /*
  // Free the refine buffers
  Refiner::FreeProcs(from_procs);
  Refiner::FreeProcs(to_procs);
  */
}

#define ZERO_THRESHOLD 0.001

void Orb3dLB::orbPartition(CkVec<Event> *events, OrientedBox<double> &box, int nprocs){

  ORB3DLB_DEBUG("partition events %d %d %d nprocs %d\n", 
            events[XDIM].length(),
            events[YDIM].length(),
            events[ZDIM].length(),
            nprocs
            );
  int numEvents = events[XDIM].length();
  CkAssert(numEvents == events[YDIM].length());
  CkAssert(numEvents == events[ZDIM].length());

  if(nprocs == 1){
    ORB3DLB_DEBUG("base: assign %d tps to proc %d\n", numEvents, nextProc);
    // direct assignment of tree pieces to processors
    //if(numEvents > 0) CkAssert(nprocs != 0);
    float totalLoad = 0.0;
    for(int i = 0; i < events[XDIM].length(); i++){
      Event &ev = events[XDIM][i];
      OrbObject &orb = tps[ev.owner];
      if(orb.numParticles > 0){
        (*mapping)[orb.lbindex] = nextProc;
        totalLoad += ev.load;
        procbox[nextProc].grow(orb.centroid);
      }
      else{
        int fromPE = (*from)[orb.lbindex];
        procload[fromPE] += ev.load;
      }
      // in order to control the number of objects moved
    }
    procload[nextProc] += totalLoad;

    if(numEvents > 0) nextProc++;
    return;
  }

  // find longest dimension

  int longestDim = XDIM;
  float longestDimLength = box.greater_corner[longestDim] - box.lesser_corner[longestDim];
  for(int i = YDIM; i <= ZDIM; i++){
    float thisDimLength = box.greater_corner[i]-box.lesser_corner[i];
    if(thisDimLength > longestDimLength){
      longestDimLength = thisDimLength;
      longestDim = i;
    }
  }

  ORB3DLB_DEBUG("dimensions %f %f %f longest %d\n", 
            box.greater_corner[XDIM]-box.lesser_corner[XDIM],
            box.greater_corner[YDIM]-box.lesser_corner[YDIM],
            box.greater_corner[ZDIM]-box.lesser_corner[ZDIM],
            longestDim
          );

  int nlprocs = nprocs/2;
  int nrprocs = nprocs-nlprocs;

  float ratio = (1.0*nlprocs)/(1.0*nrprocs);

  ORB3DLB_DEBUG("nlprocs %d nrprocs %d ratio %f\n", nlprocs, nrprocs, ratio);

  int splitIndex = partitionRatioLoad(events[longestDim],ratio);
  int nleft = splitIndex;
  int nright = numEvents-nleft;

  OrientedBox<double> leftBox;
  OrientedBox<double> rightBox;

  leftBox = rightBox = box;
  float splitPosition = events[longestDim][splitIndex].position;
  leftBox.greater_corner[longestDim] = splitPosition;
  rightBox.lesser_corner[longestDim] = splitPosition;

  // classify events
  for(int i = 0; i < splitIndex; i++){
    Event &ev = events[longestDim][i];
    CkAssert(ev.owner >= 0);
    CkAssert(tps[ev.owner].partition == INVALID_PARTITION);
    tps[ev.owner].partition = LEFT_PARTITION;
  }
  for(int i = splitIndex; i < numEvents; i++){
    Event &ev = events[longestDim][i];
    CkAssert(ev.owner >= 0);
    CkAssert(tps[ev.owner].partition == INVALID_PARTITION);
    tps[ev.owner].partition = RIGHT_PARTITION;
  }

  CkVec<Event> leftEvents[NDIMS];
  CkVec<Event> rightEvents[NDIMS];

  for(int i = 0; i < NDIMS; i++){
    if(i == longestDim){ 
      leftEvents[i].resize(nleft);
      rightEvents[i].resize(nright);
    }
    else{
      leftEvents[i].reserve(nleft);
      rightEvents[i].reserve(nright);
    }
  }

  // copy events of split dimension
  memcpy(leftEvents[longestDim].getVec(),events[longestDim].getVec(),sizeof(Event)*nleft);
  memcpy(rightEvents[longestDim].getVec(),events[longestDim].getVec()+splitIndex,sizeof(Event)*nright);
  
  // copy events of other dimensions
  for(int i = XDIM; i <= ZDIM; i++){
    if(i == longestDim) continue;
    for(int j = 0; j < numEvents; j++){
      Event &ev = events[i][j];
      CkAssert(ev.owner >= 0);
      OrbObject &orb = tps[ev.owner];
      CkAssert(orb.partition != INVALID_PARTITION);
      if(orb.partition == LEFT_PARTITION) leftEvents[i].push_back(ev);
      else if(orb.partition == RIGHT_PARTITION) rightEvents[i].push_back(ev);
    }
  }

  // cleanup
  // next, reset the ownership information in the
  // OrbObjects, so that the next invocation may use
  // the same locations for its book-keeping
  CkVec<Event> &eraseVec = events[longestDim];
  for(int i = 0; i < numEvents; i++){
    Event &ev = eraseVec[i];
    CkAssert(ev.owner >= 0);
    OrbObject &orb = tps[ev.owner];
    CkAssert(orb.partition != INVALID_PARTITION);
    orb.partition = INVALID_PARTITION;
  }

  // free events from parent node,
  // since they are not needed anymore
  // (we have partition all events into the
  // left and right event subsets)
  for(int i = 0; i < NDIMS; i++){
    events[i].free();
  }

  orbPartition(leftEvents,leftBox,nlprocs);
  orbPartition(rightEvents,rightBox,nrprocs);
}

#if 0
void Orb3dLB_notopo::directMap(int tpstart, int tpend, int nodestart, int nodeend){
  //CkPrintf("[Orb3dLB_notopo] mapping %d objects to Node (%d,%d,%d)\n", ntp, nodes[0].x, nodes[0].y, nodes[0].z);

  std::priority_queue<TPObject> pq_obj;
  std::priority_queue<Processor> pq_proc;

  float load = 0.0;
  CkAssert(nodestart==(nodeend-1));
  for(int i = tpstart; i < tpend; i++){
    //CkPrintf("obj %d thisindex %d %d %f %f %f %f to node %d %d %d\n", tp[i].lbindex, tp[i].index, tp[i].nparticles, tp[i].load, tp[i].centroid.x, tp[i].centroid.y, tp[i].centroid.z, nodes[0].x, nodes[0].y, nodes[0].z);
    load += tps[i].load;
    pq_obj.push(tps[i]);
  }
  //CkPrintf("node %d %d %d total load %f\n", nodes[0].x, nodes[0].y, nodes[0].z, load);
  
  for(int i = 0; i < procsPerNode; i++){
    Processor p;
    p.load = 0.0;
    p.t = i;
    pq_proc.push(p);
  }

  int currentZeroProc = 0;
  while(!pq_obj.empty()){
    TPObject tp = pq_obj.top();
    pq_obj.pop();

    // spread around zero-load objects
    // disabled to reduce the number of migrations, and
    // check whether this might solve the BG/P crash
    if(tp.load < ZERO_THRESHOLD){
      /*
      (*mapping)[tp.lbindex] = nodes[0].procRanks[currentZeroProc];
      currentZeroProc = currentZeroProc+1;
      if(currentZeroProc == procsPerNode){
        currentZeroProc = 0;
      }
      */
    }
    else{
      // if object has some non-zero load, assign it to a proc greedily
      Processor p = pq_proc.top();
      pq_proc.pop();

      //CkPrintf("proc %d load %f gets obj %d load %f\n", p.t, p.load, tp.lbindex, tp.load);

      p.load += tp.load;
      (*mapping)[tp.lbindex] = nodes[nodestart].procRanks[p.t];
#ifdef PRINT_BOUNDING_BOXES
      nodes[nodestart].box.grow(tp.centroid);
#endif

      pq_proc.push(p);
    }
  }
}
#endif

#define LOAD_EQUAL_TOLERANCE 1.02

int Orb3dLB::partitionRatioLoad(CkVec<Event> &events, float ratio){
  float totalLoad = 0.0;
  for(int i = 0; i < events.length(); i++){
    totalLoad += events[i].load;
  }
  //CkPrintf("************************************************************\n");
  //CkPrintf("partitionEvenLoad start %d end %d total %f\n", tpstart, tpend, totalLoad);
  float lload = 0.0;
  float rload = totalLoad;
  float prevDiff = lload-ratio*rload;
  if(prevDiff < 0.0){
    prevDiff = -prevDiff;
  }

  int consider;
  for(consider = 0; consider < events.length();){
    float newll = lload + events[consider].load;
    float newrl = rload - events[consider].load;

    float newdiff = newll-ratio*newrl;
    if(newdiff < 0.0){
      newdiff = -newdiff;
    }

    ORB3DLB_DEBUG("consider load %f newdiff %f prevdiff %f\n", events[consider].load, newdiff, prevDiff);

    if(newdiff > prevDiff){
      break;
    }
    else{
      consider++;
      lload = newll;
      rload = newrl;
      prevDiff = newdiff;
    }
  }

  ORB3DLB_DEBUG("partitionEvenLoad mid %d lload %f rload %f ratio %f\n", consider, lload, rload, lload/rload);
  return consider;
}

#include "Orb3dLB.def.h"

/*@}*/
