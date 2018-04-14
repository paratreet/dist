#include "mdt.h"

namespace MultiphaseDistree {

extern void userInitializeReducers();

void systemInitializeReducers(){
  Reduce<TraversalStatistics>::registration();
}

/* initnode */ 
void InitializeReducers(){
  systemInitializeReducers();
  userInitializeReducers();
}

// statics
int CacheAccessType::Invalid = 0;
int CacheAccessType::Readonly = 1;
int CacheAccessType::Accumulate = 2;

std::ostream &operator<<(std::ostream &out, const TraversalStatistics &stats){
  out << "[";
  out << "nodeCacheMisses: " << stats.nNodeCacheMisses() << ", ";
  out << "leafCacheMisses: " << stats.nLeafCacheMisses();
  out << "]";
  return out;
}

}; // namespace MultiphaseDistree

void MDT_MANAGER_VERBOSE_FN(const char *fmt, ...)
{
#ifdef MDT_MANAGER_VERBOSE_FN
#endif
}



#include "mdt.def.h"
