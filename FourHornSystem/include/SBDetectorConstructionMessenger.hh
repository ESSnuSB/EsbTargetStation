#ifndef SBDetectorConstructionMessenger_hh
#define SBDetectorConstructionMessenger_hh

#include "G4UImessenger.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "globals.hh"

#include "SBDetectorConstruction.hh"

class SBDetectorConstruction;

class SBDetectorConstructionMessenger : public G4UImessenger
{
                
	private :  

		SBDetectorConstruction* fSBDetectorConstruction ;  

		G4UIcmdWithADoubleAndUnit* fHornDistanceToXAxis_command ;  
		G4UIcmdWithADoubleAndUnit* fHornDistanceToYAxis_command ; 
                G4UIcmdWithADoubleAndUnit* fHornDistanceToRAxis_command ;
                G4UIcmdWithADoubleAndUnit* fHorncurrent1_command        ;
                G4UIcmdWithADoubleAndUnit* fHorncurrent2_command        ;
		G4UIcmdWithoutParameter*   fUpdate_command              ; 
	
	
	public :  

		SBDetectorConstructionMessenger (SBDetectorConstruction* aSBDetectorConstruction) ;  
		~SBDetectorConstructionMessenger () ;  

		void SetNewValue (G4UIcommand* command, G4String newValue) ;  
	        
} ;  


#endif
