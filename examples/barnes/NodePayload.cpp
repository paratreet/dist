#include "NodePayload.h"

std::ostream &operator<<(std::ostream &out, const NodePayload &pl){
  out << pl.moments().com << "\\n" << pl.moments().rsq;
  return out;
}

void NodePayload::pup(PUP::er &p){
  p | moments_;
}

void BallSphPayload::pup(PUP::er &p){
  p | rsq;
  p | com;
  p | box;
}
