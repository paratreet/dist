#ifndef SPLITTER_GROUP_H
#define SPLITTER_GROUP_H

#include "Splitter.h"
#include "barnes.decl.h"

class SplitterGroup : public CBase_SplitterGroup {
  CkVec<Splitter> splitters_;

  public:
  SplitterGroup();

  void splitters(CkVec<Splitter> &splitters, const CkCallback &cb);
  const CkVec<Splitter> &getSplitters();
};

#endif // SPLITTER_GROUP_H
