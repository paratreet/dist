#include <sstream>

#include "Main.h"
#include "Params.h"
#include "ParamStorage.h"

#include "TreePiece.h"
#include "SphShadow.h"
#include "SplitterGroup.h"
#include "LbOnOffGroup.h"
#include "Decomposer.h"
#include "ParticleMsg.h"

#include "mdt.h"
using namespace MultiphaseDistree;

// readonly
ParamStorage parameters;
CProxy_TreePiece treePieceProxy;
CProxy_SphShadow sphShadowProxy;
CProxy_Decomposer decomposerProxy;
CProxy_SplitterGroup splitterGroupProxy;
CProxy_LbOnOffGroup lbOnOffGroupProxy;

MyTreeHandleType treeHandle;
GravityVisitor::TraversalHandleType gravityTraversalHandle;
DensitySphVisitor::TraversalHandleType sphDensityTraversalHandle;
BallSphVisitor::TraversalHandleType sphBallTraversalHandle;

Main::Main(CkArgMsg *m){
  ParamCollection args;
  // compulsory
  std::string fname;
  args.add("inFile", &fname, true);
  // max num tree pieces
  args.add("nPieces", &parameters.nPieces, 1, true);
  // max particles per chare (tree piece)
  args.add("maxppc", &parameters.maxppc, 1000, true);
  // max particles per leaf 
  args.add("maxppb", &parameters.maxppb, 10, true);
  // whether to print tree after merging
  args.add("printTree", &parameters.doPrintTree, true, false);
  // opening angle
  args.add("theta", &parameters.theta, Real(0.5), true);
  // integration timestep
  args.add("dtime", &parameters.dtime, Real(0.0001), true);
  // num iterations
  args.add("nIterations", &parameters.nIterations, 10, true);
  // num leaves to traverse before yielding
  args.add("yieldPeriod", &parameters.yieldPeriod, 10, true);
  args.add("cacheChunkDepth", &parameters.cacheChunkDepth, 5, true);
  args.add("lbPeriod", &parameters.lbPeriod, 10, true);

  // softening
  args.add("epssq", &parameters.epssq, Real(0.05), true);

  // sph
  args.add("enforceSphHminLimit", &parameters.enforceSphHminLimit, true, false);
  args.add("maxSphNeighbors", &parameters.maxSphNeighbors, 32, true);
  args.add("ball2OverSoft2", &parameters.ball2OverSoft2, Real(1.0), true);

  // FIXME - sound speed
  args.add("soundSpeed2", &parameters.soundSpeed2, Real(0), true);
  // FIXME - density0
  args.add("density0", &parameters.density0, Real(0), true);

  args.add("noGravity", &parameters.noGravity, false, false);
  args.add("noSphKnn", &parameters.noSphKnn, false, false);
  args.add("noSphBall", &parameters.noSphBall, false, false);

  if(!args.process(m->argc, m->argv)){
    CkPrintf("[main] compulsory arguments missing; pgm not started\n");
    CkExit();
  }

  CkAssert(fname.size() < MAX_FILENAME_CHARS);
  parameters.nFilenameChars = fname.size();
  fname.copy(parameters.inFile, fname.size());
  parameters.inFile[fname.size()] = '\0';

  CkPrintf("[main] Main::Main\n");

  // legacy from SPLASH barnes
  Real tol = parameters.theta;
  parameters.tolsq = tol * tol;

  // branching factor
  parameters.nNodeChildren = BRANCH_FACTOR;
  parameters.nNodeChildrenLg = LOG_BRANCH_FACTOR;

  parameters.invRootTwoPi = 1.0/sqrt(0.5 * PI);

  CkArrayOptions opts(parameters.nPieces);
  treePieceProxy = CProxy_TreePiece::ckNew(opts);
  opts.bindTo(treePieceProxy);

  sphShadowProxy = CProxy_SphShadow::ckNew(opts);

  decomposerProxy = CProxy_Decomposer::ckNew();
  splitterGroupProxy = CProxy_SplitterGroup::ckNew();
  lbOnOffGroupProxy = CProxy_LbOnOffGroup::ckNew();

  treeHandle = MyTreeHandleType::instantiate(treePieceProxy, parameters.cacheChunkDepth);

  // CREATE TRAVERSAL MANAGERS
  // gravity
  CacheAccessType ro = CacheAccessType::GetReadonly();
  if(!parameters.noGravity){
    gravityTraversalHandle = GravityVisitor::TraversalType::instantiate(ro, ro);
  }

  // sph density
  if(!parameters.noSphKnn){
    sphDensityTraversalHandle = DensitySphVisitor::TraversalType::instantiate(ro, ro);
  }

  // sph ball
  CacheAccessType acc = CacheAccessType::GetAccumulate();
  if(!parameters.noSphBall){
    sphBallTraversalHandle = BallSphVisitor::TraversalType::instantiate(ro, acc);
  }


  // CREATE CACHES
  // GRAVITY + DENSITY-SPH
  // nodes
  NodeCacheClientsType nodeCacheClients;
  if(!parameters.noGravity) nodeCacheClients.add(gravityTraversalHandle);
  if(!parameters.noSphKnn) nodeCacheClients.add(sphDensityTraversalHandle);
  MyGravityCacheManagerType::Handle nodeCacheHandle = MyGravityCacheManagerType::instantiate(0, nodeCacheClients);

  // particles
  LeafCacheClientsType particleCacheClients;
  if(!parameters.noGravity) particleCacheClients.add(gravityTraversalHandle);
  if(!parameters.noSphKnn) particleCacheClients.add(sphDensityTraversalHandle);
  MyGravityCacheManagerType::Handle particleCacheHandle = MyGravityCacheManagerType::instantiate(0, particleCacheClients);

  // BALL-SPH
  // for ball sph phase, we have only one (writeback) traversal, one node cache,
  // one particle cache, and one client (the lone traversal) for each cache.
  // nodes
  NodeCacheClientsType sphBallNodeCacheClients;
  if(!parameters.noSphBall) sphBallNodeCacheClients.add(sphBallTraversalHandle);
  MySphBallCacheManagerType::Handle sphBallNodeCacheHandle = MySphBallCacheManagerType::instantiate(0, sphBallNodeCacheClients);

  // particles
  LeafCacheClientsType sphBallParticleCacheClients;
  if(!parameters.noSphBall) sphBallParticleCacheClients.add(sphBallTraversalHandle);
  MySphBallCacheManagerType::Handle sphBallParticleCacheHandle = MySphBallCacheManagerType::instantiate(0, sphBallParticleCacheClients);


  // SETUP CACHES FOR TRAVERSALS
  // gravity
  if(!parameters.noGravity){
    gravityTraversalHandle.nodeCacheProxy() = nodeCacheHandle.proxy();
    gravityTraversalHandle.leafCacheProxy() = particleCacheHandle.proxy();
  }

  // density-sph
  if(!parameters.noSphKnn){
    sphDensityTraversalHandle.nodeCacheProxy() = nodeCacheHandle.proxy();
    sphDensityTraversalHandle.leafCacheProxy() = particleCacheHandle.proxy();
  }

  // ball-sph
  if(!parameters.noSphBall){
    sphBallTraversalHandle.nodeCacheProxy() = sphBallNodeCacheHandle.proxy();
    sphBallTraversalHandle.leafCacheProxy() = sphBallParticleCacheHandle.proxy();
  }
  
  /*
  CkPrintf("[main] created node cache: gid %d\n", nodeCacheHandle.proxy().ckGetGroupID().idx);
  CkPrintf("[main] created particle cache: gid %d\n", particleCacheHandle.proxy().ckGetGroupID().idx);
  CkPrintf("[main] created sph-ball node cache: gid %d\n", sphBallNodeCacheHandle.proxy().ckGetGroupID().idx);
  CkPrintf("[main] created sph-ball particle cache: gid %d\n", sphBallParticleCacheHandle.proxy().ckGetGroupID().idx);
  */

  thisProxy.commence();
  delete m;
}

void Main::commence(){
  double startTime;

  startTime = CkWallTimer();
  // initialize tree 
  treeHandle.initialize();
  CkPrintf("[main] framework initialized in %g sec\n", CkWallTimer() - startTime);

  lbOnOffGroupProxy.off(CkCallbackResumeThread());

  CkReductionMsg *result = NULL;
  startTime = CkWallTimer();
  // load particles from file
  treePieceProxy.load(CkCallbackResumeThread((void *&) result));
  CkPrintf("[main] particles loaded in %g sec\n", CkWallTimer() - startTime);

  BoundingBox universe = *((BoundingBox *) result->getData());
  delete result;

  traversalsDone_ = 0;

  // initialize traversals; this is required so that the traversals
  // can get access to tree manager
  if(!parameters.noGravity) gravityTraversalHandle.initialize(treeHandle, CkCallbackResumeThread());
  if(!parameters.noSphKnn) sphDensityTraversalHandle.initialize(treeHandle, CkCallbackResumeThread());
  if(!parameters.noSphBall) sphBallTraversalHandle.initialize(treeHandle, CkCallbackResumeThread());

  numConcurrentTraversals_ = 0;
  // only gravity and sph knn are done concurrently; ball sph, if done at all,
  // happens after both have finished
  if(!parameters.noGravity) numConcurrentTraversals_++;
  if(!parameters.noSphKnn) numConcurrentTraversals_++;

  for(int iteration = 0; iteration < parameters.nIterations; iteration++){
    Real pad = 0.0001;
    universe.expand(pad);
    cout << endl << "[main] Iteration " << iteration << " bb: " << universe << endl;

    // DO LOAD BALANCING?
    if(iteration % parameters.lbPeriod == 0 && iteration > 0){
      cout << endl << "[main] Balancing load iteration " << iteration << endl;
      treePieceProxy.balanceLoad(CkCallbackResumeThread());
    }

    startTime = CkWallTimer();
    // assign a 64-bit position-dependent key to each particle
    treePieceProxy.assignKeys(universe, CkCallbackResumeThread());
    CkPrintf("[main] keys assigned/sorted in %g sec\n", CkWallTimer() - startTime);

    startTime = CkWallTimer();
    // spatially decompose the particles onto tree pieces
    decomposerProxy.decompose(universe, CkCallbackResumeThread());
    CkPrintf("[main] decomposition done in %g sec\n", CkWallTimer() - startTime);

    startTime = CkWallTimer();
    // build local trees from the particles assigned to pieces;
    // when done, submit local tree to framework for merging
    treePieceProxy.build(CkCallbackResumeThread());
    CkPrintf("[main] tree build done in %g sec\n", CkWallTimer() - startTime);

    startTime = CkWallTimer();
    // merge trees to create shared copy (one per PE/ or one per node)
    // depending on level of support requested
    treeHandle.syncToMerge(CkCallbackResumeThread());
    CkPrintf("[main] tree merge done in %g sec\n", CkWallTimer() - startTime);

    startTime = CkWallTimer();
    // DO GRAVITY AND DENSITY-SPH TOGETHER
    // first synchronize traversals
    if(!parameters.noGravity) gravityTraversalHandle.synchronize();
    if(!parameters.noSphKnn) sphDensityTraversalHandle.synchronize();

    // TURN ON LB INSTRUMENTATION
    lbOnOffGroupProxy.on(CkCallbackResumeThread());

    // initiate gravity
    if(!parameters.noGravity) treePieceProxy.gravity(CkCallback(CkIndex_Main::gravityDone(NULL), thisProxy));
    // initiate density sph
    if(!parameters.noSphKnn) treePieceProxy.sph(CkCallback(CkIndex_Main::sphDone(NULL), thisProxy));

    if(numConcurrentTraversals_ > 0) setTraversalCallback(CkCallbackResumeThread()); 

    CkPrintf("[main] Gravity + Density-SPH done in %g sec\n", CkWallTimer() - startTime);

    // done with gravity traversal
    //CkPrintf("[main] Gravity::done\n");
    if(!parameters.noGravity) gravityTraversalHandle.done();
    // done with gravity traversal
    //CkPrintf("[main] DensitySph::done\n");
    if(!parameters.noSphKnn) sphDensityTraversalHandle.done();

    // NOW INITIATE BALL-SPH TRAVERSAL
    if(!parameters.noSphBall){
      // first synchronize traversal
      startTime = CkWallTimer();
      sphBallTraversalHandle.synchronize();

      treePieceProxy.ballSph(CkCallbackResumeThread((void *&) result));

      CkPrintf("[main] Ball-SPH done in %g sec\n", CkWallTimer() - startTime);
      TraversalData *d = (TraversalData *) result->getData();
      std::ostringstream oss;
      oss << *d;
      CkPrintf("[main] Ball SPH traversal stats: %s\n", oss.str().c_str());
      delete result;


      // done with ball-sph traversal
      sphBallTraversalHandle.done();
    }

    // TURN OFF LB INSTRUMENTATION
    lbOnOffGroupProxy.off(CkCallbackResumeThread());

    // INTEGRATE PARTICLE TRAJECTORIES
    startTime = CkWallTimer();
    treePieceProxy.integrate(CkCallbackResumeThread((void *&)result));
    CkPrintf("[main] integration done in %g sec\n", CkWallTimer() - startTime);
    universe = ((IntegrateVisitor *) result->getData())->box();
    delete result;

    if(!parameters.noGravity) cout << "gravity: " << gravityTraversalHandle.statistics() << endl;
    if(!parameters.noSphKnn) cout << "sph-knn: " << sphDensityTraversalHandle.statistics() << endl;
    if(!parameters.noSphBall) cout << "sph-ball: " << sphBallTraversalHandle.statistics() << endl;

    treeHandle.syncToDelete(CkCallbackResumeThread());
  }

  CkPrintf("[main] exit.\n");
  CkExit();
}

void Main::gravityDone(CkReductionMsg *result){
  TraversalData *d = (TraversalData *) result->getData();
  std::ostringstream oss;
  oss << *d;
  CkPrintf("[main] gravity traversal stats: %s\n", oss.str().c_str());
  CkAssert(d->good());
  delete result;

  checkTraversalsDone();
}

void Main::sphDone(CkReductionMsg *result){
  TraversalData *d = (TraversalData *) result->getData();
  std::ostringstream oss;
  oss << *d;
  CkPrintf("[main] SPH traversal stats: %s\n", oss.str().c_str());
  CkAssert(d->good());
  delete result;

  checkTraversalsDone();
}

void Main::setTraversalCallback(const CkCallback &cb){
  traversalCallback_ = cb;
}

void Main::checkTraversalsDone(){
  if(++traversalsDone_ == numConcurrentTraversals_){
    traversalsDone_ = 0;
    traversalCallback_.send();
  }
}

namespace MultiphaseDistree {

void userInitializeReducers(){
  //Reduce<BarnesHutVisitor>::registration();
  //Reduce<InterCheckVisitor>::registration();
  Reduce<IntegrateVisitor>::registration();
  Reduce<TraversalData>::registration();
  Reduce<BoundingBox>::registration();
}

}

int GravityTraversalDataInterface::NodeRequestPriority = -130000000;
int GravityTraversalDataInterface::LeafRequestPriority = -110000000;
int GravityTraversalDataInterface::NodeReceivePriority = -120000000;
int GravityTraversalDataInterface::LeafReceivePriority = -100000000;

int SphBallTraversalDataInterface::NodeRequestPriority = -130000000;
int SphBallTraversalDataInterface::LeafRequestPriority = -110000000;
int SphBallTraversalDataInterface::NodeReceivePriority = -120000000;
int SphBallTraversalDataInterface::LeafReceivePriority = -100000000;


Real Physics::openingGeometryFactor = 2.0/sqrt(3);

#include "TemplateInstanceTypedefs.h"
#include "barnes.def.h"
