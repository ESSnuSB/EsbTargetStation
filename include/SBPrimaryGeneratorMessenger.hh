#ifndef SBPrimaryGeneratorMessenger_h
#define SBPrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class SBPrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;

class SBPrimaryGeneratorMessenger: public G4UImessenger
{
public:
  SBPrimaryGeneratorMessenger(SBPrimaryGeneratorAction*);
  virtual ~SBPrimaryGeneratorMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);
    
private:
  SBPrimaryGeneratorAction* SBAction;
  G4UIdirectory* gunDir; 
  G4UIcmdWithAnInteger* ExtGenCmd;
  G4UIcmdWithAString* RndmCmd;
  G4UIdirectory* gunDir1; 
  G4UIcmdWithAString* FileInCmd;
  G4UIcmdWithADoubleAndUnit* EKinCmd;
  G4UIcmdWithADoubleAndUnit* x0_command;
  G4UIcmdWithADoubleAndUnit* y0_command;
  G4UIcmdWithAnInteger* JOBIDCmd;
};

#endif

