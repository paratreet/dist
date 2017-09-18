#ifndef MULTIPOLE_MOMENTS_H
#define MULTIPOLE_MOMENTS_H

#include "defines.h"
#include "Particle.h"
#include "OrientedBox.h"

struct MultipoleMoments {
  Real rsq;
  Real mass;
  Vector3D<Real> com;
  OrientedBox<Real> box;

  int nParticles;

  MultipoleMoments(){
    reset();
  }

  void reset(){
    mass = 0.0;
    com.x = com.y = com.z = 0.0;
    nParticles = 0;
    box.reset();
  }

  void begin(){
    reset();
  }

  void accumulate(const MultipoleMoments &other){
    com += other.mass * other.com;
    mass += other.mass;
    nParticles += other.nParticles;
    box.grow(other.box);
  }

  void end(){
    if(mass > 0.0){
      com = com / mass;
    }

    // Radius: distance between center of mass and the corner that is
    // farthest from it.
    Vector3D<Real> delta1 = com - box.lesser_corner;
    Vector3D<Real> delta2 = box.greater_corner - com;
    delta1.x = (delta1.x > delta2.x ? delta1.x : delta2.x);
    delta1.y = (delta1.y > delta2.y ? delta1.y : delta2.y);
    delta1.z = (delta1.z > delta2.z ? delta1.z : delta2.z);
    rsq = delta1.lengthSquared();
  }

  void fromParticles(const Particle *particles, int numParticles){
    for(int i = 0; i < numParticles; i++){
      const Particle &p = particles[i];
      mass += p.mass();
      com += p.mass() * p.position();
      box.grow(p.position());
    }

    if(mass > 0.0){
      com = com / mass;
    }

    nParticles = numParticles;

    Real d;
    // Radius: distance between center of mass and particle farthest from it
    rsq = 0.0;
    for(const Particle *p = particles; p != particles + numParticles; ++p) {
      d = (com - p->position()).lengthSquared();         
      if(d > rsq) rsq = d;  
    }                                                 
  }

  void pup(PUP::er &p){
    p | mass;
    p | com;
    p | nParticles;
    p | box;
    p | rsq;
  }
};

#endif // MULTIPOLE_MOMENTS_H
