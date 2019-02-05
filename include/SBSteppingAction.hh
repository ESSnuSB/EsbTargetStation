#ifndef SBSteppingAction_h
#define SBSteppingAction_h 1
#include "globals.hh"
#include "G4UserSteppingAction.hh"
#include "G4ThreeVector.hh"
#include <fstream>

class SBDetectorConstruction;
class SBEventAction;
class SBRunAction;
class SBPrimaryGeneratorAction;

class SBSteppingAction : public G4UserSteppingAction
{

  public:
    SBSteppingAction(SBDetectorConstruction*, SBEventAction*, SBPrimaryGeneratorAction*, SBRunAction*);
    virtual ~SBSteppingAction();

    void UserSteppingAction(const G4Step*);
    
  private:
    SBDetectorConstruction* detector;
    SBEventAction*          evtAct;  
    SBPrimaryGeneratorAction* primaryGen;
    SBRunAction* runAct;
    std::ofstream pOutFileExitTunnel;
    std::ofstream pOutFileExitTarget;
    G4int FLUKA_G4_PID(G4int);
};

#endif
