#include "SBDetectorConstructionMessenger.hh"

// *********************************************************************************************************

SBDetectorConstructionMessenger::SBDetectorConstructionMessenger (SBDetectorConstruction* aSBDetectorConstruction)
{

	fSBDetectorConstruction = aSBDetectorConstruction ; 

	fHornDistanceToXAxis_command    = new G4UIcmdWithADoubleAndUnit ("/SB/geometry/horn/distanceToXAxis",this) ;
	fHornDistanceToRAxis_command    = new G4UIcmdWithADoubleAndUnit ("/SB/geometry/horn/distanceToRAxis",this) ;
        fHornDistanceToYAxis_command    = new G4UIcmdWithADoubleAndUnit ("/SB/geometry/horn/distanceToYAxis",this) ;
        fHorncurrent1_command           = new G4UIcmdWithADoubleAndUnit ("/SB/geometry/horn/current1",this) ;
        fHorncurrent2_command           = new G4UIcmdWithADoubleAndUnit ("/SB/geometry/horn/current2",this) ;
	fUpdate_command                 = new G4UIcmdWithoutParameter   ("/SB/geometry/update"              ,this) ;

}


// *********************************************************************************************************


SBDetectorConstructionMessenger::~SBDetectorConstructionMessenger ()
{
	delete fHornDistanceToXAxis_command ; fHornDistanceToXAxis_command = NULL ; 
	delete fHornDistanceToRAxis_command ; fHornDistanceToRAxis_command = NULL ;
        delete fHornDistanceToYAxis_command ; fHornDistanceToYAxis_command = NULL ; 
        delete fHorncurrent1_command        ; fHorncurrent1_command        = NULL ;
        delete fHorncurrent2_command        ; fHorncurrent2_command        = NULL ;
	delete fUpdate_command              ; fUpdate_command              = NULL ; 
}


// *********************************************************************************************************


void SBDetectorConstructionMessenger::SetNewValue (G4UIcommand* command, G4String newValue)
{ 	
	if ( command == fHornDistanceToXAxis_command ) 
    { fSBDetectorConstruction->SetHornDistanceToXAxis(fHornDistanceToXAxis_command->GetNewDoubleValue(newValue)) ; } 
	if ( command == fHornDistanceToYAxis_command ) 
    { fSBDetectorConstruction->SetHornDistanceToYAxis(fHornDistanceToYAxis_command->GetNewDoubleValue(newValue)) ; } 
	if ( command == fHornDistanceToRAxis_command ) 
    { fSBDetectorConstruction->SetHornDistanceToRAxis(fHornDistanceToRAxis_command->GetNewDoubleValue(newValue)) ; }
        if ( command == fHorncurrent1_command ) 
    { fSBDetectorConstruction->SetHorncurrent1(fHorncurrent1_command->GetNewDoubleValue(newValue)) ; }
        if ( command == fHorncurrent2_command ) 
    { fSBDetectorConstruction->SetHorncurrent2(fHorncurrent2_command->GetNewDoubleValue(newValue)) ; }
        if ( command == fUpdate_command                   ) { fSBDetectorConstruction->UpdateGeometry() ; } 

	return ; 	
}
