#ifndef SBRunAction_h
#define SBRunAction_h 1

#include "G4UserRunAction.hh"
#include "SBPrimaryGeneratorAction.hh"
//#include "HistoManager.hh"
//#include "DumpDataROOT.hh"
#include "globals.hh"

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

class G4Run;
class SBRunActionMessenger;
//class DumpDataROOT;
class SBRunAction : public G4UserRunAction
{
public:
  SBRunAction(SBPrimaryGeneratorAction*);
  virtual ~SBRunAction();

  void BeginOfRunAction(const G4Run*);
  void   EndOfRunAction(const G4Run*);
    
  void CheckEOF();
  void fillPerEvent(G4double, G4double, G4double, G4double);
  void fillPerEvent1(G4double, G4double, G4double, G4double,
		     G4double, G4double, G4double, G4double);
  void fillPerEvent1_W(G4double, G4double, G4double, G4double,
		       G4double, G4double, G4double, G4double,G4int,G4int);
  void fillPerEvent2_W(G4double, G4double, G4double, G4double,
		       G4double, G4double, G4double, G4double,
		       G4double, G4double, G4double, G4double,
		       G4double, G4double, G4double, G4double);
  void fillDIF(G4double,
	       G4double,
	       G4double,
	       G4double,
	       G4double,
	       G4double,
	       G4double,
	       G4double,
	       G4double,
	       G4int,
	       G4int);
  void DIFStat(G4double, G4double, G4double, G4double,G4double, G4double);
  void fillEXIT(G4double, G4double, G4double, G4double, G4double, G4double, G4int);
  void fillTARG(G4double/*G4double, G4double, G4double, G4double, G4double, G4double, G4int*/);

  void fillMU(G4double,G4double,G4double,G4int,G4int);

  void FluxesStats();
  void SumFluxesFlavours();
  
  //void PrepareGLoBESFlux(G4String,G4double);
  void DoFinalOps(G4int);

  void SetNbin(G4int val){ENbins = val;};
  void SetEnuMIN(G4double val){EnuMIN = val;};
  void SetEnuMAX(G4double val){EnuMAX = val;};

  //DumpDataROOT* m_hh;

  TFile *fROOT;
  
  G4double rMAX;
  G4double xMAX;
  G4double yMAX;

  G4double zMIN;
  G4double zMAX;

  G4double NDIF_piplus_TOT;
  G4double NDIF_piminus_TOT;
  G4double NDIF_Kplus_TOT;
  G4double NDIF_Kminus_TOT;
  G4double NDIF_KzeroS_TOT;
  G4double NDIF_KzeroL_TOT;

  G4int ENbins;
  G4double EnuMIN;
  G4double EnuMAX;

  int NENB;//9
  //const double ENULIMS[NENB+1];//change if NENB initialized differently ...
  double ENULIMS[10];//change if NENB initialized differently ...
 
  //TH1F *h_numu;// all
  //TH1F *h_anumu;
  //TH1F *h_nue;
  //TH1F *h_anue;

  //--------------------------------//

  // quantities at the exit of the horn
  

  TH1F *h_pt_exit_piplus;
  TH1F *h_p_exit_piplus;
  TH2F *h_xy_exit_piplus;
  TH1F *h_z_exit_piplus;
  TH1F *h_r_exit_piplus;
  TH2F *h_pVSr_exit_piplus;       //impulsion p versus distance r
  TH2F *h_pVStheta_exit_piplus;   //impulsion p versus angle theta
  TH2F *h_thetaVSr_exit_piplus;   //impulsion theta versus distance r
  TH1F *h_theta_exit_piplus;   
  TH2F *h_pVStheta_exit_piplus_W; //impulsion p versus angle theta

  TH1F *h_pt_exit_piminus;
  TH1F *h_p_exit_piminus;
  TH2F *h_xy_exit_piminus;
  TH1F *h_z_exit_piminus;
  TH1F *h_r_exit_piminus;
  TH2F *h_pVSr_exit_piminus;
  TH2F *h_pVStheta_exit_piminus;
  TH2F *h_thetaVSr_exit_piminus;
  TH1F *h_theta_exit_piminus;

  TH1F *h_pt_exit_kplus;
  TH1F *h_p_exit_kplus;
  TH2F *h_xy_exit_kplus;
  TH1F *h_z_exit_kplus;
  TH1F *h_r_exit_kplus;
  TH2F *h_pVSr_exit_kplus;
  TH2F *h_pVStheta_exit_kplus;
  TH2F *h_thetaVSr_exit_kplus;
  TH1F *h_theta_exit_kplus;

  TH1F *h_pt_exit_kminus;
  TH1F *h_p_exit_kminus;
  TH2F *h_xy_exit_kminus;
  TH1F *h_z_exit_kminus;
  TH1F *h_r_exit_kminus;
  TH2F *h_pVSr_exit_kminus;
  TH2F *h_pVStheta_exit_kminus;
  TH2F *h_thetaVSr_exit_kminus;
  TH1F *h_theta_exit_kminus;

  TH1F *h_pt_exit_muplus;
  TH1F *h_p_exit_muplus;
  TH2F *h_xy_exit_muplus;
  TH1F *h_z_exit_muplus;
  TH1F *h_r_exit_muplus;
  TH2F *h_pVSr_exit_muplus;    
  TH2F *h_pVStheta_exit_muplus;
  TH2F *h_thetaVSr_exit_muplus;
  TH1F *h_theta_exit_muplus;

  TH1F *h_pt_exit_muminus;
  TH1F *h_p_exit_muminus;
  TH2F *h_xy_exit_muminus;
  TH1F *h_z_exit_muminus;
  TH1F *h_r_exit_muminus;
  TH2F *h_pVSr_exit_muminus;
  TH2F *h_pVStheta_exit_muminus;
  TH2F *h_thetaVSr_exit_muminus;
  TH1F *h_theta_exit_muminus;

  TH1F *h_pt_exit_k0L;
  TH1F *h_p_exit_k0L;
  TH2F *h_xy_exit_k0L;
  TH1F *h_z_exit_k0L;
  TH1F *h_r_exit_k0L;
  TH2F *h_pVSr_exit_k0L;
  TH2F *h_pVStheta_exit_k0L;
  TH2F *h_thetaVSr_exit_k0L;
  TH1F *h_theta_exit_k0L;

  TH1F *h_pt_targ_piplus;
  TH1F *h_p_targ_piplus;
  TH1F *h_z_targ_piplus;
  TH1F *h_zabs_targ_piplus;
  TH1F *h_theta_targ_piplus;
  TH2F *h_pVStheta_targ_piplus;
  TH2F *h_pVStheta_targ_piplus_Enubins[9];//NENB
  TH1F *h_p_targ_piplus_Enubins[9];
  TH1F *h_theta_targ_piplus_Enubins[9];

  TH1F *h_pt_targ_piminus;
  TH1F *h_p_targ_piminus;
  TH1F *h_z_targ_piminus;
  TH1F *h_zabs_targ_piminus;
  TH1F *h_theta_targ_piminus;
  TH2F *h_pVStheta_targ_piminus;
  TH2F *h_pVStheta_targ_piminus_Enubins[9];//NENB
  TH1F *h_p_targ_piminus_Enubins[9];
  TH1F *h_theta_targ_piminus_Enubins[9];

  TH1F *h_MUDEC_Enumustar;
  TH1F *h_MUDEC_Enuestar;
  TH1F *h_MUDEC_Eanumustar;
  TH1F *h_MUDEC_Eanuestar;
  //TH1F *h_MUDEC_Eeplusstar;
  //TH1F *h_MUDEC_Eeminusstar;
  //TH1F *h_MUDEC_costhstarmu;
  //TH1F *h_MUDEC_costhstareplus;

  TH1F *h_MUDEC_probaTunnel;
  TH1F *h_MUDEC_pathTunnel;

  TH2F *h_MUDEC_probaTunnel_pmu;
  TH2F *h_MUDEC_probaTunnel_thmu;
  TH2F *h_MUDEC_probaTunnel_rmu;

  TH2F *h_MUDEC_muxy;
  TH2F *h_MUDEC_muzr;

  TH1F *h_MUDEC_cthstarMu;
  TH1F *h_MUDEC_cthstarNu;

  TH1F *h_MUDEC_polaL_muplus;
  TH1F *h_MUDEC_polaL_muminus;
  TH1F *h_MUDEC_thpimu_lab;

  TH2F *h_MUDEC_pmupar;

  TH2F *h_geom_g4;

  TH1F *h_decay_flags;
  TH1F *h_decay_NUM;
  TH1F *h_decay_BR;
  TH1F *h_decay_BR_PDG;

  TH2F *h_prim_match_pi;
  TH2F *h_prim_match_ka;
  TH2F *h_prim_match_k0;

private:

  SBRunActionMessenger* runActMessenger;//messenger of this class
  G4double sumETarg, sum2ETarg;
  G4double sumEGap, sum2EGap;
    
  G4double sumLTarg, sum2LTarg;
  G4double sumLGap, sum2LGap;    

  SBPrimaryGeneratorAction* primaryGen;

  //HistoManager* histoManager;

};

#endif

