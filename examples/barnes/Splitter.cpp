#include "Splitter.h"

Splitter::Splitter() : 
  from(~Key(0)),
  to(Key(0)),
  nParticles(-1)
{}

Splitter::Splitter(Key f, Key t, Key r, int n) :
  from(f),
  to(t),
  tpRootKey(r),
  nParticles(n)
{
}

void Splitter::pup(PUP::er &p){
  p | from;
  p | to;
  p | tpRootKey;
  p | nParticles;
}

bool Splitter::operator<=(const Splitter &other) const {
  return from <= other.from;
}

bool Splitter::operator>=(const Splitter &other) const {
  return from >= other.from;
}

bool Splitter::operator>=(const Key &k) const {
  return from >= k;
}
