#include <algorithm>
#include <vector>

template<typename T>
class Heap {
  CkVec<T> storage_;

  public:
  typedef T DataType;
  typedef CkVec<T> StorageType;

  Heap()
  {}

  void push(const T &t){
    try{
      storage_.push_back(t);
    }
    catch(std::bad_alloc &e){
      CkPrintf("Heap::push memory allocation error %d\n", CmiMemoryUsage());
      CkExit();
    }
    up(storage_.size()-1);
  }

  void pop(){
    int i = storage_.size()-1;
    swap(0,i);
    //storage_.pop_back();
    //storage_.length()--;
    storage_.resize(storage_.size()-1);
    down(0);
  }

  T &top(){
    CkAssert(!empty());
    return storage_[0];
  }

  bool empty() const {
    return storage_.size() == 0;
  }

  const CkVec<T> &storage() const {
    return storage_;
  }

  int size() const {
    return storage_.size();
  }

  void free() {
    storage_.free();
  }

  private:
  void up(int i){
    if(i == 0) return;

    int pi = (i-1)/2;
    CkAssert(pi >= 0);

    if(storage_[pi] < storage_[i]){
      swap(pi, i);
    }

    up(pi);
  }

  void down(int j){
    int li = -1;
    T compare = storage_[j];

    for(int i = 0; i < 2; i++){
      int ci = 2 * j + i + 1;
      if(ci < storage_.size() && compare < storage_[ci]){ 
        li = ci;
        compare = storage_[ci];
      }
    }

    CkAssert(li != 0);
    if(li > 0){
      swap(j, li);
      down(li);
    }
  }

  void swap(int i, int j){
    T tmp = storage_[i];
    storage_[i] = storage_[j];
    storage_[j] = tmp;
  }
};

template<typename T>
std::ostream &operator<<(std::ostream &out, const Heap<T> &h){
  out << "[";
  for(int i = 0; i < h.storage().size(); i++){
    out << h.storage()[i] << ", ";
  }
  out << "]";
  return out;
}



template<typename T>
class MyPriorityQueue {
  std::vector<T> storage_;

  public:
  typedef std::vector<T> StorageType;

  void free(){
    std::vector<T>().swap(storage_);
  }

  void push(const T &t){
    storage_.push_back(t);
    std::push_heap(storage_.begin(), storage_.end());
  }

  T &top(){
    return storage_[0];
  }

  void pop(){
    CkAssert(storage_.size() > 0);
    std::pop_heap(storage_.begin(), storage_.end());
    storage_.pop_back();
  }

  int size() const {
    return storage_.size();
  }

  const std::vector<T> &storage() const {
    return storage_;
  }
};
