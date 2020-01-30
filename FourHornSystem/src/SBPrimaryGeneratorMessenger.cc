#include "SBPrimaryGeneratorMessenger.hh"
#include "SBPrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

SBPrimaryGeneratorMessenger::SBPrimaryGeneratorMessenger(SBPrimaryGeneratorAction* SBGun)
:SBAction(SBGun)
{
  gunDir = new G4UIdirectory("/SB/gun/");
  gunDir->SetGuidance("PrimaryGenerator control");
   
  RndmCmd = new G4UIcmdWithAString("/SB/gun/rndm",this);
  RndmCmd->SetGuidance("Shoot randomly the incident particle.");
  RndmCmd->SetGuidance("  Choice : on(default), off");
  RndmCmd->SetParameterName("choice",true);
  RndmCmd->SetDefaultValue("on");
  RndmCmd->SetCandidates("on off");
  RndmCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ExtGenCmd = new G4UIcmdWithAnInteger("/SB/det/ExtGen",this);
  ExtGenCmd->SetGuidance("Use external generator if par is 1");
  ExtGenCmd->SetParameterName("ExtGen",true);
  ExtGenCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  gunDir1 = new G4UIdirectory("/SB/file/");
  gunDir1->SetGuidance("PrimaryGenerator input file");
   
  FileInCmd = new G4UIcmdWithAString("/SB/file/myinput",this);
  FileInCmd->SetGuidance("Input file in FLUKA style.");
  FileInCmd->SetGuidance("  Choice : FLUKAcard001_input.dat(default)");
  //FileInCmd->SetParameterName("choice",true);
  FileInCmd->SetDefaultValue("FLUKAcard001_input.dat");
  //FileInCmd->SetCandidates("pippo.dat FLUKAcard001_input.dat");
  FileInCmd->AvailableForStates(G4State_PreInit,G4State_Idle,G4State_Idle,G4State_GeomClosed,G4State_EventProc);

  x0_command = new G4UIcmdWithADoubleAndUnit("/SB/gun/Setx0",this);
  x0_command->SetGuidance("Set gun position x");
  x0_command->SetParameterName("x0",true);
  x0_command->SetUnitCategory("Length");
  x0_command->AvailableForStates(G4State_PreInit,G4State_Idle);

  y0_command = new G4UIcmdWithADoubleAndUnit("/SB/gun/Sety0",this);
  y0_command->SetGuidance("Set gun position y");
  y0_command->SetParameterName("y0",true);
  y0_command->SetUnitCategory("Length");
  y0_command->AvailableForStates(G4State_PreInit,G4State_Idle);

  EKinCmd = new G4UIcmdWithADoubleAndUnit("/SB/gun/SetEKin",this);
  EKinCmd->SetGuidance("Set kinetic energy and unit");
  EKinCmd->SetParameterName("gunEKin",true);
  EKinCmd->SetRange("gunEKin>=0.");
  EKinCmd->SetUnitCategory("Energy");
  EKinCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  JOBIDCmd = new G4UIcmdWithAnInteger("/SB/det/JOBID",this);
  JOBIDCmd->SetGuidance("Job ID");
  JOBIDCmd->SetParameterName("JOBID",true);
  JOBIDCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

SBPrimaryGeneratorMessenger::~SBPrimaryGeneratorMessenger()
{
  delete ExtGenCmd;
  delete RndmCmd;
  delete gunDir;
  delete FileInCmd;
  delete gunDir1;
  delete EKinCmd;
  delete JOBIDCmd;
  delete x0_command;
  delete y0_command;
}

void SBPrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if( command == ExtGenCmd )
    { 
      G4cout << "ciccio GENERATOR "<<ExtGenCmd->GetNewIntValue(newValue)<<G4endl;
      SBAction->SetExtGen(ExtGenCmd->GetNewIntValue(newValue));}

  if( command == RndmCmd )
   { SBAction->SetRndmFlag(newValue);}
  if( command == FileInCmd )
   {
     G4cout << "pippo SetInputFileName("<<newValue<<")"<<G4endl;
     //SBAction->SetInputFileName(FileInCmd->GetNewStringValue(newValue));
     SBAction->SetInputFileName(newValue);
   }
  if( command == EKinCmd )
   {
     SBAction->SetGunEKin(EKinCmd->GetNewDoubleValue(newValue));
   }
  if( command == x0_command )
   {
     SBAction->Setx0(x0_command->GetNewDoubleValue(newValue));
   }
  if( command == y0_command )
   {
     SBAction->Sety0(y0_command->GetNewDoubleValue(newValue));
   }
  if(command == JOBIDCmd)
    {SBAction->SetJOBID(JOBIDCmd->GetNewIntValue(newValue));}
}
