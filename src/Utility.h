#ifndef MULTIPHASE_DISTREE_UTILITY_H
#define MULTIPHASE_DISTREE_UTILITY_H

#include <cstdlib>

namespace MultiphaseDistree {

class Utility {
  public:
  static int pickRandom(int begin, int end){
    int nObjects = end - begin;
    return ((rand() % nObjects) + begin);
  }
};

}; // namespace MultiphaseDistree

#endif // MULTIPHASE_DISTREE_UTILITY_H
