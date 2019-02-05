#ifndef SBRunActionMessenger_h
#define SBRunActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class SBRunAction;
class G4UIdirectory;

class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;

class SBRunActionMessenger: public G4UImessenger
{
public:
  SBRunActionMessenger(SBRunAction*);
  virtual ~SBRunActionMessenger();

  void SetNewValue(G4UIcommand*, G4String);

private:
  SBRunAction*     runAction;
  G4UIdirectory* SBrunDir;
  G4UIcmdWithAnInteger* NbinCmd;    
  G4UIcmdWithADoubleAndUnit* EnuMINCmd;
  G4UIcmdWithADoubleAndUnit* EnuMAXCmd;
};
#endif
