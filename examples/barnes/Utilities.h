#ifndef UTILITIES_H
#define UTILITIES_H

#include "defines.h"

class Utilities {

  public:

  static int mssb64_pos(Key x) {
    int n;

    if (x == Key(0)) return -1;
    n = 0;
    if (x > Key(0x00000000FFFFFFFF)) {n += 32; x >>= 32;}
    if (x > Key(0x000000000000FFFF)) {n += 16; x >>= 16;}
    if (x > Key(0x00000000000000FF)) {n += 8; x >>= 8;}
    if (x > Key(0x000000000000000F)) {n += 4; x >>= 4;}
    if (x > Key(0x0000000000000003)) {n += 2; x >>= 2;}
    if (x > Key(0x0000000000000001)) {n += 1;}
    return n;
  }

  static Key mssb64(Key x){
    x |= (x >> 1);
    x |= (x >> 2);
    x |= (x >> 4);
    x |= (x >> 8);
    x |= (x >> 16);
    x |= (x >> 32);
    return(x & ~(x >> 1));
  }

  template<typename KEY_TYPE, typename OBJ_TYPE>
  static int binary_search_ge(const KEY_TYPE &check, const OBJ_TYPE *particles, int start, int end){
    int lo = start;
    int hi = end;
    int mid;
    while(lo < hi){
      mid = lo+((hi-lo)>>1);
      if(particles[mid] >= check){
        hi = mid;
      }
      else{
        lo = mid+1;
      }
    }
    return lo;
  }

  static Key getParticleLevelKey(Key k, int depth){
    return (k<<(TREE_KEY_BITS-(LOG_BRANCH_FACTOR*depth+1)));
  }

  static Key getLastParticleLevelKey(Key k, int depth){
    int nshift = TREE_KEY_BITS-(LOG_BRANCH_FACTOR*depth+1);
    k <<= nshift;
    Key mask = (Key(1) << nshift);
    mask--;
    k |= mask;
    return k;
  }

  // Find depth of node from position of prepended '1' bit in node
  static int getDepthFromKey(Key k){
    int nUsedBits = mssb64_pos(k);
    return (nUsedBits/LOG_BRANCH_FACTOR);
  }

  static Key toSplitter(Key k){
    int depth = Utilities::getDepthFromKey(k);
    k = Utilities::getParticleLevelKey(k, depth);
    return k;
  }

  static int completeTreeSize(int levels){
    return numLeaves(levels+1)-1;
  }

  static int numLeaves(int levels){
    return (1<<(LOG_BRANCH_FACTOR*levels));
  }

  // is k1 a prefix of k2?
  static bool isPrefix(Key k1, int d1, Key k2, int d2){
    if(d1 > d2) return false;
    k2 = (k2 >> (d2 - d1));
    return k1 == k2;
  }

};

#endif // UTILITIES_H
