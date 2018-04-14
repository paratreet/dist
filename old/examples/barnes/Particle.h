#ifndef PARTICLE_H
#define PARTICLE_H

#include "defines.h"
#include "Vector3D.h"

#ifdef __CHARMC__
#include "pup.h"
#endif

class SphCore {
  static Real dummyDensity_;
  static Real dummyPressure_;
  static Vector3D<Real> dummyVelocity_;

  public:
  virtual Real &mass() = 0;
  virtual Vector3D<Real> &position() = 0;

  virtual Real &density() {
    CkAbort("SphCore::density called!");
    return dummyDensity_;
  }

  virtual Real &pressure() {
    CkAbort("SphCore::pressure called!");
    return dummyPressure_;
  }

  virtual Vector3D<Real> &velocity(){
    CkAbort("SphCore::velocity called!");
    return dummyVelocity_;
  }

  virtual const Real &mass() const = 0;
  virtual const Vector3D<Real> &position() const  = 0;

  virtual const Real &density() const {
    CkAbort("SphCore::density(const) called!");
    return dummyDensity_;
  }

  virtual const Real &pressure() const {
    CkAbort("SphCore::pressure(const) called!");
    return dummyPressure_;
  }

  virtual const Vector3D<Real> &velocity() const {
    CkAbort("SphCore::velocity(const) called!");
    return dummyVelocity_;
  }
};

class ParticleCore : public SphCore {
  Vector3D<Real> position_;
  Real mass_;

  public:
  Vector3D<Real> &position(){
    return position_;
  }

  Real &mass(){
    return mass_;
  }

  const Vector3D<Real> &position() const {
    return position_;
  }

  const Real &mass() const {
    return mass_;
  }

  void pup(PUP::er &p);
};

class BallSphParticleCore : public SphCore {
  Real mass_;
  Vector3D<Real> position_;
  Vector3D<Real> velocity_;
  Real density_;
  Real pressure_;

  public:
  Real &mass(){
    return mass_;
  }

  Vector3D<Real> &position(){
    return position_;
  }

  Real &density() {
    return density_;
  }

  Real &pressure() {
    return pressure_;
  }

  Vector3D<Real> &velocity(){
    return velocity_;
  }

  const Real &mass() const {
    return mass_;
  }

  const Vector3D<Real> &position() const {
    return position_;
  }

  const Real &density() const {
    return density_;
  }

  const Real &pressure() const {
    return pressure_;
  }

  const Vector3D<Real> &velocity() const {
    return velocity_;
  }

  void pup(PUP::er &p);
};

struct Particle {
  Key key;
  Vector3D<Real> acceleration;
  Real potential;
  int order;

  // for SPH
  BallSphParticleCore ballSphCore;

#ifdef CHECK_INTER
  Real interMass;
#endif

  Particle();

  bool operator<=(const Particle &other) const;
  bool operator>=(const Particle &other) const;
  bool operator>=(const Key &key) const;

  const Vector3D<Real> &position() const {
    return ballSphCore.position();
  }

  Vector3D<Real> &position() {
    return ballSphCore.position();
  }

  Vector3D<Real> &velocity() {
    return ballSphCore.velocity();
  }

  const Real &mass() const {
    return ballSphCore.mass();
  }

  Real &mass() {
    return ballSphCore.mass();
  }

  Real &pressure() {
    return ballSphCore.pressure();
  }

  Real &density() {
    return ballSphCore.density();
  }

  SphCore *getSphCore() {
    return &ballSphCore;
  }

  const SphCore *getSphCore() const {
    return &ballSphCore;
  }

#ifdef __CHARMC__
  void pup(PUP::er &p);
#endif

  void reset();
};

struct GravityRemoteParticle {
  ParticleCore core_;

#ifdef __CHARMC__
  void pup(PUP::er &p);
#endif

  GravityRemoteParticle() {}

  GravityRemoteParticle(const Particle &other){
    *this = other;
  }

  GravityRemoteParticle &operator=(const Particle &other){
    const SphCore *otherSphCore = other.getSphCore();
    core_.mass() = otherSphCore->mass(); 
    core_.position() = otherSphCore->position(); 
    return *this;
  }

  const Vector3D<Real> &position() const{
    return core_.position();
  }

  const Real &mass() const {
    return core_.mass();
  }

  SphCore *getSphCore() {
    return &core_;
  }

  const SphCore *getSphCore() const {
    return &core_;
  }

};

std::ostream &operator<<(std::ostream &out, const SphCore &p);

struct BallSphRemoteParticle {
  BallSphParticleCore ballSphCore;

#ifdef __CHARMC__
  void pup(PUP::er &p);
#endif

  BallSphRemoteParticle(){}

  BallSphRemoteParticle(const Particle &other){
    *this = other;
  }

  BallSphRemoteParticle &operator=(const Particle &other){
    const SphCore *otherSphCore = other.getSphCore();
    ballSphCore.mass() = otherSphCore->mass(); 
    ballSphCore.position() = otherSphCore->position(); 
    ballSphCore.velocity() = otherSphCore->velocity(); 
    ballSphCore.density() = otherSphCore->density(); 
    ballSphCore.pressure() = otherSphCore->pressure(); 
    return *this;
  }

  const Vector3D<Real> &position() const{
    return ballSphCore.position();
  }

  Vector3D<Real> &velocity() {
    return ballSphCore.velocity();
  }

  const Real &mass() const {
    return ballSphCore.mass();
  }

  Real &density() {
    return ballSphCore.density();
  }

  Real &pressure() {
    return ballSphCore.pressure();
  }

  SphCore *getSphCore() {
    return &ballSphCore;
  }

  const SphCore *getSphCore() const {
    return &ballSphCore;
  }


};

#endif // PARTICLE_H
