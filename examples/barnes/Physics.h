#ifndef PHYSICS_H
#define PHYSICS_H

#include "defines.h"
#include "MultipoleMoments.h"
#include "ParamStorage.h"
#include "Vector3D.h"
#include "Sphere.h"

#define SPH_TRACKED_LEAF 121
#define SPH_TRACKED_PARTICLE 66

extern ParamStorage parameters;

class Physics {
  static Real openingGeometryFactor;
  public:

  // adapted from ChaNGa, courtesy trq@astro.washington.edu
  static 
  bool intersect(const OrientedBox<Real> &box, const Sphere<Real> &sphere){
    Real dsq = 0.0;
    Real rsq = sphere.radius * sphere.radius;
    Real delta;

    if((delta = box.lesser_corner.x - sphere.origin.x) > 0)
      dsq += delta * delta;
    else if((delta = sphere.origin.x - box.greater_corner.x) > 0)
      dsq += delta * delta;
    if(rsq < dsq)
      return false;
    if((delta = box.lesser_corner.y - sphere.origin.y) > 0)
      dsq += delta * delta;
    else if((delta = sphere.origin.y - box.greater_corner.y) > 0)
      dsq += delta * delta;
    if(rsq < dsq)
      return false;
    if((delta = box.lesser_corner.z - sphere.origin.z) > 0)
      dsq += delta * delta;
    else if((delta = sphere.origin.z - box.greater_corner.z) > 0)
      dsq += delta * delta;
    return (dsq <= rsq);
  }



  class Gravity {
    public:
    // has to work for both home- and away-nodes
    template<typename AnyNodeType>
    static 
    bool open(const MyLocalNodeType *leaf, const AnyNodeType *node){
      Real nodeMomentsRadius = sqrt(GravityTraversalDataInterface::rsq(node));

      Real radius = openingGeometryFactor * nodeMomentsRadius / parameters.theta;
      if(radius < nodeMomentsRadius) radius = nodeMomentsRadius;

      Sphere<Real> nodeSphere(GravityTraversalDataInterface::com(node), radius);

      bool intersection = intersect(TreeBuildDataInterface::box(leaf), nodeSphere);

      return intersection; 
      //Vector3D<Real> dr = node->getUserData().moments().com - leaf->getUserData().moments().com;
      //return (parameters.tolsq * dr.lengthSquared() < node->getUserData().moments().rsq);
    }

    template<typename AnyNodeType>
    static int forces(MyLocalNodeType *leaf, const AnyNodeType *node){
      Particle *targets = leaf->getParticles();
      int nTargets = leaf->getNumParticles();
      for(Particle *t = targets; t != targets + nTargets; t++){
        grav(t, 
             GravityTraversalDataInterface::mass(node), 
             GravityTraversalDataInterface::com(node));
      }     
      return nTargets;
    }

    template<typename SourceType>
    static 
    int forces(MyLocalNodeType *leaf, const SourceType *sources, int nSources){
      Particle *targets = leaf->getParticles();
      int nTargets = leaf->getNumParticles();
      int nInteractions = 0;
      for(Particle *t = targets; t != targets + nTargets; t++){
        for(const SourceType *s = sources; s != sources + nSources; s++){
          if((void *) t == (void *) s) continue;
          grav(t, s->mass(), s->position());
          nInteractions++;
        }
      }

      return nInteractions;
    }

    static 
    void grav(Particle *p, Real mass, const Vector3D<Real> &position){
      Vector3D<Real> dr;
      Real drsq;
      Real drabs;
      Real phii;
      Real mor3;

      dr = position - p->position();
      drsq = dr.lengthSquared();
      drsq += parameters.epssq;
      drabs = sqrt((double) drsq);
      phii = mass/drabs;
      p->potential -= phii;
      mor3 = phii/drsq;
      p->acceleration += mor3*dr;
#ifdef CHECK_INTER
      p->interMass += mass;
#endif
    }

    // uses half-step method
    static void integrate(Particle &p, Real dt, int step){
      if(step == 0){
        // bootstrap with Euler integration
        p.velocity() += p.acceleration * 0.5 * dt;
      }
      else{
        p.velocity() += p.acceleration * dt; 
      }
      p.position() += p.velocity() * dt;
    }
  };

  class Sph {
    public:

    // no need to template this method on the type of core used:
    // should be invoked only from density sph traversal
    template<typename AnyNodeType>
    static bool open(MyLocalNodeType *leaf, SphLeafData *leafSphData, const AnyNodeType *source, int &status){
      Real rLeaf = leafSphData->size() + leafSphData->maxDist();
      // is it the case that the source's bounding box
      // doesn't come close to any of the particles in
      // this leaf?
      Sphere<Real> leafSphere(leafSphData->center(), rLeaf);


      if(leaf->getKey() == Key(SPH_TRACKED_LEAF)){
        int x = 0;
      }

      bool intersection = intersect(SphDensityTraversalDataInterface::box(source), leafSphere);
      /*
      if(leaf->getKey() == Key(SPH_TRACKED_LEAF)){
        CkPrintf("Sph::open test leaf %llu source %llu INTERSECT %d\n",
                 leaf->getKey(), source->getKey(),
                 intersection);
      }
      */

      if(!intersection){
        status = 0;
        return false;
      }

      // there might be a particle in this node
      // that is closer to some target particle in 
      // the leaf than the farthest neighbor of that 
      // target particle
      Particle *particles = leaf->getParticles();
      for(int i = 0; i < leaf->getNumParticles(); i++){
        SphParticleData &targetSphData = leafSphData->particleData(i);
        // if there is a target particle in this leaf that doesn't
        // have enough neighbors, we must greedily open the source node.
        // if there isn't such a target in the leaf, we open only if at
        // least one particle's smoothing sphere intersects with the source
        // node's bounding box
        if(targetSphData.numNeighbors() < parameters.maxSphNeighbors ||
           intersect(SphDensityTraversalDataInterface::box(source),
                     Sphere<Real>(particles[i].position(), sqrt(targetSphData.maxDist2())))){
          status = 1;
          return true;
        }
      }

      status = 2;
      return false;
    }


    // adapted from ChaNGa, courtesy trq@astro.washington.edu
    // check particle source 'p' against the neighbor queues of 
    // all the target particles in 'leaf'. the code i have here is 
    // a little different from ChaNGa's, since i don't do a 
    // boot-strapping of leaves with local particles in tree-order


    // since this is only called from the sph density traversal, we shouldn't
    // template on sph core type, and expect a ParticleCore template argument
    // for SphLeafData
    template<typename ParticleType>
    static void compare(const ParticleType *p, MyLocalNodeType *leaf, SphLeafData *leafSphData){
      // The logic behind this code is the following: initially, each leaf
      // has a maxDist of infty and is greedy. When applying the opening
      // greedy leaves have infinite maxDist; so, they are always opened. 
      // Although this will cause the first several source nodes
      // encountered, to be opened until their respective leaves, once we are
      // at the leaves, we fill the neighbor queues of the target leaf's
      // particles, and check whether all of the targets have (at least) 'k'
      // (=32) nbrs.  if so, the leaf is not greedy anymore, and its
      // maxDist is updated to max_{t \in leaf}{r(t,p) : p \in neighbors(t)}.
      // subsequently, the traversal will be more conservative in opening
      // source nodes.
 
      Vector3D<Real> drLeaf = leafSphData->center() - p->position();
      Real drLeaf2 = drLeaf.lengthSquared();
      Real rExpansion = leafSphData->size() + leafSphData->maxDist();
      Real rExpansion2 = rExpansion * rExpansion;

      /*
      if(leaf->getKey() == Key(SPH_TRACKED_LEAF)){
        CkPrintf("Sph::compare leaf %llu source %d dr2 %g rEx2 %g\n", 
                  leaf->getKey(), p->order, drLeaf2, rExpansion2);
      }
      */

      if(rExpansion2 < drLeaf2){
        // source particle is outside all smoothing radii
        return;
      }

      // this variable stores the maximum non-infinite farthest-neighbor
      // distance for all the targets in this leaf. we use it to update
      // the maxDist of the leaf only if no target particles in the leaf
      // are found to be greedy
      Real leafMaxDist2 = 0.0;
      const Particle *targets = leaf->getParticles();
      // are there still some particles that haven't received their fill of
      // neighbors, and are therefore 'greedy'? note that if even a single
      // target particle in a leaf is greedy, so is the leaf itself.
      bool noTargetsGreedy = true;

      // in original changa code, hmin limit was enforced on a per-particle basis
      // here we assume the force smoothing length to be uniform across all particles
      //Real hMin2 = parameters.ball2OverSoft2 * targets[i].soft * targets[i].soft;
      Real hMin2 = parameters.ball2OverSoft2 * parameters.epssq;

      for(int i = 0; i < leaf->getNumParticles(); i++) {
        SphParticleData &targetSphData = leafSphData->particleData(i); 
        if((void *) p == (void *) &targets[i]){
          noTargetsGreedy &= (targetSphData.numNeighbors() >= parameters.maxSphNeighbors); 
          continue;
        }

        Vector3D<Real> dr = targets[i].position() - p->position();
        Real dr2 = dr.lengthSquared();


        // do i have enough source particles for this target?
        if(targetSphData.numNeighbors() >= parameters.maxSphNeighbors){
          // i can afford to be picky here, or take on neighbors
          // only if forced to by the hmin limit.

          if(parameters.enforceSphHminLimit && dr2 < hMin2){
            // i'm being forced by hmin limit to 
            // accept this particle as a neighbor
            /*
            if(targets[i].order == SPH_TRACKED_PARTICLE){
              CkPrintf("Sph::compare ENFORCE-HMIN target %d dr2 %g hMin2 %g neighbors %d\n", 
                        targets[i].order, dr2, hMin2, targetSphData.numNeighbors());
            }
            */

            targetSphData.add(SphParticle(dr2, p->getSphCore()));
          }
          else if(dr2 < targetSphData.maxDist2()){
            /*
            if(targets[i].order == SPH_TRACKED_PARTICLE){
              CkPrintf("Sph::compare CLOSER THAN FARTHEST target %d dr2 %g farthest %g neighbors %d\n", 
                        targets[i].order, dr2, targetSphData.maxDist2(), targetSphData.numNeighbors());
            }
            */

            // not being forced to accept particle as neighbor, but it 
            // is closer than the neighbor that is currently farthest
            // from the target, so accept it. however, since we already have
            // our fill of particles, first pop the farthest source neighbor
            targetSphData.addWithReplacement(SphParticle(dr2, p->getSphCore()));
          }
        }
        else{
          /*
          if(targets[i].order == SPH_TRACKED_PARTICLE){
            CkPrintf("Sph::compare leaf %llu GREEDY target %d neighbors %d\n", 
                      leaf->getKey(), targets[i].order, targetSphData.numNeighbors());
          }
          */

          // this is a greedy target, it'll accept as a neighbor 
          // whichever source particle it gets
          targetSphData.add(SphParticle(dr2, p->getSphCore()));
          // check whether this target is still greedy
          noTargetsGreedy &= (targetSphData.numNeighbors() >= parameters.maxSphNeighbors);
        }

        if(leafMaxDist2 < targetSphData.maxDist2()){
          leafMaxDist2 = targetSphData.maxDist2();
        }
      }

      // this is needed to ensure that every leaf has an opening
      // radius of at least hMin2
      if(parameters.enforceSphHminLimit && leafMaxDist2 < hMin2){
        leafMaxDist2 = hMin2;
      }

      if(noTargetsGreedy){
        leafSphData->maxDist() = sqrt(leafMaxDist2);
      }
    }

    static void density(Particle &target, const SphCore *sourceCore, Real distance){
      target.density() += sourceCore->mass() * weight(distance);
    }

    static Real weight(Real distance){
      return parameters.invRootTwoPi * exp(-0.5 * distance * distance);
    }

  };

  // BallSph inherits method intersect from Sph
  class BallSph : public Sph {
    public:
    // no need to template this method on type of core data: should 
    // be called only for the BallSph walk
    template<typename AnyNodeType>
    static bool open(SphLeafData *data, const AnyNodeType *source){
      Real rExpansion = data->size() + data->maxDist(); 
      const OrientedBox<Real> &box = SphBallTraversalDataInterface::box(source);
      Sphere<Real> sphere(data->center(), rExpansion);
      //ostringstream oss;
      bool did = intersect(box, sphere);
      //oss << "BallSph::open intersect: " << did << " box: " << box << " sphere: "  << sphere;
      //CkPrintf("%s\n", oss.str().c_str());
      return did;
    }

    template<typename ParticleType>
    static void compare(const ParticleType *source, MyLocalNodeType *leaf, SphLeafData *data){
      Particle *targets = leaf->getParticles();
      int nTargets = leaf->getNumParticles();
      for(int i = 0; i < nTargets; i++){
        if(((const void *) source) == ((const void *) &targets[i])) continue;

        Vector3D<Real> dr = source->position() - targets[i].position(); 
        Real dist2 = dr.lengthSquared();
        SphParticleData &targetSphData = data->particleData(i);
        if(dist2 <= targetSphData.rCutoff2()){
          // the source particle is within the cutoff
          // sphere of the target particle, i.e. it is
          // a neighbor of the target
          ostringstream oss;
          oss << *source->getSphCore();
          targetSphData.add(SphParticle(dist2, source->getSphCore()));
        }
      }
    }
  };
};

#endif // PHYSICS_H
