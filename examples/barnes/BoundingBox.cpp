#include "BoundingBox.h"

BoundingBox::BoundingBox(){
  reset();
}

void BoundingBox::reset(){
  numParticles = 0;
  box.reset();
  pe = 0.0;
  ke = 0.0;
  mass = 0.0;
}

void BoundingBox::grow(const Vector3D<Real> &v){
  box.grow(v);
}

  /*
   * This method is called when performing a reduction over 
   * BoundingBox's. It subsumes the bounding box of the 'other'
   * and accumulates its energy in its own. If a PE has no
   * particles, its contributions are not counted.
   */
void BoundingBox::grow(const BoundingBox &other){
  if(other.numParticles == 0) return;
  if(numParticles == 0){
    *this = other;
  }
  else{
    box.grow(other.box);
    numParticles += other.numParticles;
    pe += other.pe;
    ke += other.ke;
    mass += other.mass;
  }
}

void BoundingBox::expand(Real pad){
  box.greater_corner = box.greater_corner*pad+box.greater_corner;
  box.lesser_corner = box.lesser_corner-pad*box.lesser_corner;
}

void BoundingBox::pup(PUP::er &p){
  p | box;
  p | numParticles;
  p | pe;
  p | ke;
  p | mass;
}


ostream &operator<<(ostream &os, const BoundingBox &bb){
  os << "<"
     << bb.numParticles << ", "
     << bb.mass << ", "
     << bb.pe + bb.ke << ", "
     << bb.box
     << ">";
      
  return os;
}
