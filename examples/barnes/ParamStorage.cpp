#include "ParamStorage.h"

void ParamStorage::pup(PUP::er &p){
  p | nPieces;
  p | nFilenameChars;
  PUParray(p, inFile, nFilenameChars);

  p | maxppc;
  p | maxppb;
  p | nNodeChildren;
  p | nNodeChildrenLg;
  p | doPrintTree;
  p | nIterations;

  p | yieldPeriod;
  p | cacheChunkDepth;
  p | lbPeriod;

  p | theta;
  p | dtime;
  p | epssq;
  p | tolsq;

  p | ball2OverSoft2;
  p | maxSphNeighbors;
  p | enforceSphHminLimit;

  p | invRootTwoPi;
  p | soundSpeed2;
  p | density0;

  p | noGravity;
  p | noSphKnn;
  p | noSphBall;

}
