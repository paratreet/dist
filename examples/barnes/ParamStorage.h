#ifndef PARAM_STORAGE_H
#define PARAM_STORAGE_H

#include "pup_stl.h"
#include <string>

using namespace std;

#include "defines.h"

#define MAX_FILENAME_CHARS 256

struct ParamStorage {
  int nPieces;

  int nFilenameChars;
  char inFile[MAX_FILENAME_CHARS]; 
  int maxppc;
  int maxppb;
  int nNodeChildren;
  int nNodeChildrenLg;
  bool doPrintTree;
  int nIterations;

  int yieldPeriod;
  int cacheChunkDepth;
  int lbPeriod;

  bool noGravity;
  bool noSphKnn;
  bool noSphBall;

  Real theta;
  Real dtime;
  Real epssq;
  Real tolsq;

  // sph
  Real ball2OverSoft2;
  int maxSphNeighbors;
  bool enforceSphHminLimit;

  Real invRootTwoPi;
  Real soundSpeed2;
  Real density0;

  void pup(PUP::er &p);
};

#endif // PARAM_STORAGE_H
