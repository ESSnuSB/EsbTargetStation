#ifndef SBEventActionMessenger_h
#define SBEventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class SBEventAction;
class G4UIdirectory;

class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;

class SBEventActionMessenger: public G4UImessenger
{
public:
  SBEventActionMessenger(SBEventAction*);
  virtual ~SBEventActionMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);
    
private:
  SBEventAction*     eventAction;
  G4UIdirectory*        eventDir;   
  G4UIcmdWithAnInteger* PrintCmd;    
  G4UIcmdWithADouble* EKinProtCmd;    
  G4UIcmdWithADouble* PowerCmd;    
  G4UIcmdWithAnInteger* SBVerbosityCmd;
};

#endif
