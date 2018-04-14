#ifndef TRAVERSAL_DATA_H
#define TRAVERSAL_DATA_H

#include "pup.h"
#include <iostream>

class TraversalData {
  CmiUInt8 open_;
  CmiUInt8 pn_;
  CmiUInt8 pp_;
  CmiUInt8 misses_;
  CmiUInt8 hits_;

  bool good_;

  public:
  TraversalData(){
    reset();
  }

  public:
  void reset(){
    open_ = 0;
    pn_ = 0;
    pp_ = 0;
    misses_ = 0;
    hits_ = 0;
    good_ = true;
  }

  TraversalData &operator+=(const TraversalData &other){
    pn() += other.pn();
    pp() += other.pp();
    open() += other.open();
    misses() += other.misses();
    hits() += other.hits();
    good() &= other.good();
    return *this;
  }

  CmiUInt8 &open(){
    return open_;
  }

  CmiUInt8 &pn(){
    return pn_;
  }

  CmiUInt8 &pp(){
    return pp_;
  }

  CmiUInt8 &misses(){
    return misses_;
  }

  CmiUInt8 &hits(){
    return hits_;
  }

  bool &good(){
    return good_;
  }

  // const versions, used in operator+=
  const CmiUInt8 &open() const {
    return open_;
  }

  const CmiUInt8 &pn() const {
    return pn_;
  }

  const CmiUInt8 &pp() const {
    return pp_;
  }

  const CmiUInt8 &misses() const {
    return misses_;
  }

  const CmiUInt8 &hits() const {
    return hits_;
  }

  const bool &good() const {
    return good_;
  }

  void pup(PUP::er &p){
    p | open_;
    p | pn_;
    p | pp_;
    p | misses_;
    p | hits_;
    p | good_;
  }
};

std::ostream &operator<<(std::ostream &out, const TraversalData &data);

#endif // TRAVERSAL_DATA_H
