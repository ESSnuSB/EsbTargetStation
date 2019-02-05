#ifndef SBPhysicsList_h
#define SBPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class SBPhysicsList: public G4VUserPhysicsList
{
public:
  SBPhysicsList();
  virtual ~SBPhysicsList();

  // Construct particle and physics
  void ConstructParticle();
  void ConstructProcess();
 
  void SetCuts();
   
private:
  // these methods Construct particles 
  void ConstructBosons();
  void ConstructLeptons();
  void ConstructMesons();
  void ConstructBaryons();

  // these methods Construct physics processes and register them
  void ConstructDecay();
  void ConstructEM();
};

#endif



