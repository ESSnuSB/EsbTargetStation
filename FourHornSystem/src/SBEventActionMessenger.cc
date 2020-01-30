#include "SBEventActionMessenger.hh"

#include "SBEventAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"

#include "globals.hh"

SBEventActionMessenger::SBEventActionMessenger(SBEventAction* EvAct)
:eventAction(EvAct)
{
  eventDir = new G4UIdirectory("/SB/event/");
  eventDir->SetGuidance("event control");
   
  PrintCmd = new G4UIcmdWithAnInteger("/SB/event/printModulo",this);
  PrintCmd->SetGuidance("Print events modulo n");
  PrintCmd->SetParameterName("EventNb",false);
  PrintCmd->SetRange("EventNb>0");

  EKinProtCmd = new G4UIcmdWithADouble("/SB/event/SetEKinProt",this);
  EKinProtCmd->SetGuidance("Set incoming protons kinetic energy (for normalization to power only)");
  EKinProtCmd->SetParameterName("EKinProt",false);
  EKinProtCmd->SetRange("EKinProt>0");

  PowerCmd = new G4UIcmdWithADouble("/SB/event/SetPower",this);
  PowerCmd->SetGuidance("Set proton power in MW (for fluxes normalization)");
  PowerCmd->SetParameterName("Power",false);
  PowerCmd->SetRange("Power>0");

  SBVerbosityCmd = new G4UIcmdWithAnInteger("/SB/event/verbosity",this);
  SBVerbosityCmd->SetGuidance("verbosity level");
  SBVerbosityCmd->SetParameterName("SBVerbosity",false);
  SBVerbosityCmd->SetRange("SBVerbosity>=0");

}

SBEventActionMessenger::~SBEventActionMessenger()
{
  delete PrintCmd;
  delete SBVerbosityCmd;
  delete EKinProtCmd;
  delete PowerCmd;
  delete eventDir;
}

void SBEventActionMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if(command == PrintCmd)
    {eventAction->SetPrintModulo(PrintCmd->GetNewIntValue(newValue));}

  if(command == EKinProtCmd)
    {eventAction->SetEKinProt(EKinProtCmd->GetNewDoubleValue(newValue));}

  if(command == PowerCmd)
    {eventAction->SetPower(PowerCmd->GetNewDoubleValue(newValue));}

  if(command == SBVerbosityCmd)
    {eventAction->SetSBVerbosity(PrintCmd->GetNewIntValue(newValue));}
}
