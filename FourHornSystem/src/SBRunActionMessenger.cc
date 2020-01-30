#include "SBRunActionMessenger.hh"

#include "SBRunAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"

#include "globals.hh"

SBRunActionMessenger::SBRunActionMessenger(SBRunAction* RunAct) :runAction(RunAct)
{

  SBrunDir = new G4UIdirectory("/SB/run/");
  SBrunDir->SetGuidance("RunAction control");

  NbinCmd = new G4UIcmdWithAnInteger("/SB/run/SetNbin",this);
  NbinCmd->SetGuidance("Set number of bins in neutrino energy");
  NbinCmd->SetParameterName("Nbin",false);
  NbinCmd->SetRange("Nbin>0");
  NbinCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  EnuMINCmd = new G4UIcmdWithADoubleAndUnit("/SB/run/SetEnuMIN",this);
  EnuMINCmd->SetGuidance("Set MIN neutrino energy in flux histogram");
  EnuMINCmd->SetParameterName("EnuMIN",false);
  EnuMINCmd->SetRange("EnuMIN>=0");
  EnuMINCmd->SetUnitCategory("Energy");
  EnuMINCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  EnuMAXCmd = new G4UIcmdWithADoubleAndUnit("/SB/run/SetEnuMAX",this);
  EnuMAXCmd->SetGuidance("Set MAX neutrino energy in flux histogram");
  EnuMAXCmd->SetParameterName("EnuMAX",false);
  EnuMAXCmd->SetRange("EnuMAX>0");
  EnuMAXCmd->SetUnitCategory("Energy");
  EnuMAXCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

SBRunActionMessenger::~SBRunActionMessenger()
{
  delete SBrunDir;
  delete NbinCmd;
  delete EnuMAXCmd;
  delete EnuMINCmd;
}

void SBRunActionMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if(command == NbinCmd)
    {runAction->SetNbin(NbinCmd->GetNewIntValue(newValue));}
  if(command == EnuMAXCmd)
    {runAction->SetEnuMAX(EnuMAXCmd->GetNewDoubleValue(newValue));}
  if(command == EnuMINCmd)
    {runAction->SetEnuMIN(EnuMINCmd->GetNewDoubleValue(newValue));}
}
