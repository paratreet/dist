#include "Particle.h"

void ParticleCore::pup(PUP::er &p){
  p|position_;
  p|mass_;
}

void BallSphParticleCore::pup(PUP::er &p){
  p | position_;
  p | mass_;
  p | velocity_;
  p | density_;
  p | pressure_;
}

void GravityRemoteParticle::pup(PUP::er &p){
  p|core_;
}

void BallSphRemoteParticle::pup(PUP::er &p){
  p | ballSphCore;
}

Particle::Particle() : 
  key(Key(0))
{ 
  reset();
}

void Particle::reset(){
  potential = 0.0;
  acceleration = 0.0;
  density() = 0.0;
  pressure() = 0.0;
}


bool Particle::operator<=(const Particle &other) const { 
  return key <= other.key; 
}

bool Particle::operator>=(const Particle &other) const {
  return key >= other.key; 
}

bool Particle::operator>=(const Key &k) const {
  return key >= k; 
}

void Particle::pup(PUP::er &p){
  p|key;
  p|potential;
  p|order;
  p|ballSphCore;
}

std::ostream &operator<<(std::ostream &out, const SphCore &p){
  out << p.position();
  return out;
}


Real SphCore::dummyDensity_;
Real SphCore::dummyPressure_;
Vector3D<Real> SphCore::dummyVelocity_;

