#ifndef SBPrimaryGeneratorAction_h
#define SBPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
//#include "SBRunAction.hh"

#include "globals.hh"
//#include <CLHEP/Vector/ThreeVector.h>
#include <fstream>

class G4ParticleGun;
class G4Event;
class SBDetectorConstruction;
//class SBRunAction;
class SBPrimaryGeneratorMessenger;
class G4VPrimaryGenerator;
//class G4ThreeVector;
class SBPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  //SBPrimaryGeneratorAction(SBRunAction*,SBDetectorConstruction*);
  SBPrimaryGeneratorAction(SBDetectorConstruction*);
  virtual ~SBPrimaryGeneratorAction();

  void GeneratePrimaries(G4Event*);
  void SetRndmFlag(G4String val) { rndmFlag = val;};
  void SetExtGen(G4int val){ExtGen = val;};
  void SetJOBID(G4int val){JOBID = val;};

  void SetInputFileName(G4String val) {filenamein = val;};
  void SetGunEKin(G4double val){gunEKin = val;G4cout << " gunEKin = " << gunEKin << G4endl;};
  void Setx0(G4double val){x0 = val;G4cout << " x0 = " << x0 << G4endl;};
  void Sety0(G4double val){y0 = val;G4cout << " y0 = " << y0 << G4endl;};
  G4String GetInputFileName() {return filenamein;};

  std::ifstream fin;
  G4bool myEOF;

  G4double XPrimary;
  G4double YPrimary;
  G4double ZPrimary;
  G4double PXPrimary;
  G4double PYPrimary;
  G4double PZPrimary;

  G4double ZPrimary_G4frame;

  G4String particleName;
  G4double x0;
  G4double y0;
  G4double gunEKin;

  G4double GetPrimaryX(){return XPrimary;};
  G4double GetPrimaryY(){return YPrimary;};
  G4double GetPrimaryZ(){return ZPrimary;};
  G4double GetPrimaryPX(){return PXPrimary;};
  G4double GetPrimaryPY(){return PYPrimary;};
  G4double GetPrimaryPZ(){return PZPrimary;};

  G4double GetFLUKAPOTs(){return FLUKApots;};
  //G4int GetPOT_id(){return POT_id;};
  G4int GetJOBID(){return JOBID;};

  // record exiting target kaons for duplication
 
  G4int K_REPL;

  G4int n_kplus_exit;
  G4int n_kminus_exit;
  G4int n_k0s_exit;
  G4int n_k0l_exit;

  G4double x_kplus_exit[10];
  G4double y_kplus_exit[10];
  G4double z_kplus_exit[10];
  G4double px_kplus_exit[10];
  G4double py_kplus_exit[10];
  G4double pz_kplus_exit[10];

  G4double x_kminus_exit[10];
  G4double y_kminus_exit[10];
  G4double z_kminus_exit[10];
  G4double px_kminus_exit[10];
  G4double py_kminus_exit[10];
  G4double pz_kminus_exit[10];

  G4double x_k0s_exit[10];
  G4double y_k0s_exit[10];
  G4double z_k0s_exit[10];
  G4double px_k0s_exit[10];
  G4double py_k0s_exit[10];
  G4double pz_k0s_exit[10];

  G4double x_k0l_exit[10];
  G4double y_k0l_exit[10];
  G4double z_k0l_exit[10];
  G4double px_k0l_exit[10];
  G4double py_k0l_exit[10];
  G4double pz_k0l_exit[10];


private:
  G4ParticleGun*                particleGun;	  //pointer a to G4  class
  //SBRunAction* runAct;
  SBDetectorConstruction*    SBDetector;    //pointer to the geometry
    
  SBPrimaryGeneratorMessenger* gunMessenger;   //m  G4double GetFLUKAPOTs(){return FLUKApots;};essenger of this class
  G4String rndmFlag;	  //flag for a rndm impact point
  G4int ExtGen;
  G4int JOBID;
  G4String filenamein;

  G4bool externalGEN;
  G4VPrimaryGenerator* HEPEvt;

  G4int procEvts;
  G4int procTracks;
  //G4int POT_id;

  G4double OLDidev;
  G4double FLUKApots;
  G4bool InputFileOK;

};

#endif
