#include "ParticleMsg.h"
ParticleMsg::ParticleMsg() : 
  particles(NULL),
  nParticles(0)
{}

ParticleMsg::ParticleMsg(Particle *p, int n) : 
  nParticles(n)
{
  memcpy(particles, p, n * sizeof(Particle));
}

ParticleMsgWrapper::ParticleMsgWrapper() : 
  msg(NULL),
  sortKey(~Key(0))
{}

ParticleMsgWrapper::ParticleMsgWrapper(ParticleMsg *m) : 
  msg(m),
  sortKey(m->particles[0].key)
{}

bool ParticleMsgWrapper::operator<=(const ParticleMsgWrapper &other) const {
  return sortKey <= other.sortKey;
}

bool ParticleMsgWrapper::operator>=(const ParticleMsgWrapper &other) const {
  return sortKey >= other.sortKey;
}
