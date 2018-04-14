#ifndef DOUBLE_BUFFERED_VEC_H
#define DOUBLE_BUFFERED_VEC_H

#include "cklists.h"

template<typename T>
class DoubleBufferedVec {
  CkVec<T> *one;
  CkVec<T> *two;

  public:
  DoubleBufferedVec(){
    one = new CkVec<T>;
    two = new CkVec<T>;
  }

  void add(T &t){
    one->push_back(t);
  }

  void add(T t){
    one->push_back(t);
  }

  void sync(){
    CkVec<T> *tmp = one;
    one = two;
    two = tmp;
    one->resize(0);
  }

  T &get(int i){
    return (*two)[i];
  }

  CkVec<T> &get(){
    return *two;
  }

  int size(){
    return two->size();
  }

  ~DoubleBufferedVec(){
    delete one;
    delete two;
  }
};

#endif // DOUBLE_BUFFERED_VEC_H
