#ifndef MULTIPHASE_DISTREE_CACHE_ACCESS_TYPE_H
#define MULTIPHASE_DISTREE_CACHE_ACCESS_TYPE_H

namespace MultiphaseDistree {

class CacheAccessType {
  int code_;

  public:
  static int Invalid;
  static int Readonly;
  static int Accumulate;

  CacheAccessType(int code=Invalid) : 
    code_(code)
  {}

  void pup(PUP::er &p){
    p | code_;
  }

  bool operator==(const int &code) const {
    return code_ == code; 
  }

  static CacheAccessType GetInvalid(){
    return CacheAccessType(Invalid);
  }

  static CacheAccessType GetReadonly(){
    return CacheAccessType(Readonly);
  }

  static CacheAccessType GetAccumulate(){
    return CacheAccessType(Accumulate);
  }

};


}; // namespace MultiphaseDistree

#endif // MULTIPHASE_DISTREE_CACHE_ACCESS_TYPE_H
