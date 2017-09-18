#ifndef DATA_INTERFACE_H
#define DATA_INTERFACE_H

#include "defines.h"
#include "mdt.h"
#include "pup.h"



// TREE DATA INTERFACES
// These are the interfaces from which we derive
// different types of nodes. Basically, these are
// passed as template parameters to 
// MultiphaseDistree::Node<> to get different types
// of nodes (i.e. Nodes with different payloads)

class MyParticleType;
class MyPayloadType;

class TreeDataInterface {
  public:
  typedef Key KeyType;
  typedef MultiphaseDistree::Node<TreeDataInterface> LocalNodeType;
  typedef MyParticleType LocalParticleType;
  typedef MyPayloadType PayloadType;
};

class MyParticleType {
public:
    float x,y,z;
    void pup(PUP::er &p) {
        p|x; p|y; p|z;
    }
};

class MyPayloadType {
public:
    float k;
    
    void pup(PUP::er &p) {
        p|k;
    }  
  
  friend std::ostream &operator<<(std::ostream &os,const MyPayloadType &p) {
    os<<p.k<<" ";
  }
  
  // Reduction interface:
  void mergeStart(){
    k=0;
  }

  void mergeAccumulate(const MyPayloadType &child){
    k+=child.k;
  }

  void mergeFinish(){
    
  }
};

typedef MultiphaseDistree::Node<TreeDataInterface> MyNodeType;


#endif // DATA_INTERFACE_H
