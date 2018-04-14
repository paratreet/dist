#ifndef DECOMPOSER_H
#define DECOMPOSER_H

#include "barnes.decl.h"
#include "Splitter.h"

class BoundingBox;

template<typename T>
class DoubleBufferedVec;

class Decomposer : public CBase_Decomposer {
  CkVec<Splitter> splitters_;

  public:
  Decomposer();

  void decompose(BoundingBox &box, const CkCallback &cb);

  private:
  void processCounts(int *counts, int nCounts, DoubleBufferedVec<Key> &keys);
  void printKeys(DoubleBufferedVec<Key> &keys);
};


#endif // DECOMPOSER_H
