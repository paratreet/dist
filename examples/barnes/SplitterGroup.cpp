#include "SplitterGroup.h"

SplitterGroup::SplitterGroup(){
}

void SplitterGroup::splitters(CkVec<Splitter> &splitters, const CkCallback &cb){
  splitters_ = splitters;
  contribute(cb);
}

const CkVec<Splitter> &SplitterGroup::getSplitters(){
  return splitters_;
}
