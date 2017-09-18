#ifndef SPLITTER_H
#define SPLITTER_H

#include "defines.h"

struct Splitter {
  Key from;
  Key to;
  Key tpRootKey;
  int nParticles;

  Splitter();
  Splitter(Key fron, Key to, Key tpRootKey, int n);

  void pup(PUP::er &p);

  bool operator<=(const Splitter &other) const;
  bool operator>=(const Splitter &other) const;
  bool operator>=(const Key &k) const;
};

#endif // SPLITTER_H
