#ifndef PARTICLE_MSG_H
#define PARTICLE_MSG_H

#include "Particle.h"
#include "barnes.decl.h"

struct ParticleMsg : public CMessage_ParticleMsg {
  Particle *particles;
  int nParticles;

  ParticleMsg();
  ParticleMsg(Particle *p, int n);
};

struct ParticleMsgWrapper {
  ParticleMsg *msg;
  Key sortKey;

  ParticleMsgWrapper();
  ParticleMsgWrapper(ParticleMsg *m);

  bool operator<=(const ParticleMsgWrapper &other) const;
  bool operator>=(const ParticleMsgWrapper &other) const;
};

#endif // PARTICLE_MSG_H
