#include "SBRunAction.hh"
#include "SBRunActionMessenger.hh"
#include "SBAnalysisManager.hh"
#include "TH1F.h"
#include "TFile.h"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SBRunAction::SBRunAction(SBPrimaryGeneratorAction* primGen):primaryGen(primGen)
{
  //create a messenger for this class
  runActMessenger = new SBRunActionMessenger(this);
  //============================================================  
  EnuMAX=1.5;// SPL standard
  ENbins=75;//SPL standard  (75)

  //EnuMAX=1.5;// SPL prova
  //ENbins=20;//SPL prova

  //EnuMAX=10.;// PS2
  //ENbins=50;//PS2
  //EnuMAX=50.;// CNGS
  //ENbins=50;// CNGS

  //EnuMAX=50.;// NOVA
  //ENbins=1000;// NOVA
  //=============================================================
  rMAX=250;

  xMAX=2.5;
  yMAX=2.5;

  zMIN=-30.; // Tunnel make this depend on tunnel parameters!!
  zMAX=30.;  // Tunnel

  NENB=9;
  for(int j=0;j<NENB;j++){
    ENULIMS[j]=0.8*((double)j/(double)NENB);
  }
  ENULIMS[NENB]=0.8;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SBRunAction::~SBRunAction()
{
  G4cout<<"CALL SBRunAction::~SBRunAction "<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SBRunAction::BeginOfRunAction(const G4Run* aRun)
{ 
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  
  G4cout<<"CALL SBRunAction::BeginOfRunAction "<<G4endl;
  G4cout<<"CALL SBRunAction::BeginOfRunAction Nbins = "<<ENbins<<G4endl;
  G4cout<<"CALL SBRunAction::BeginOfRunAction EnuMIN = "<<EnuMIN/CLHEP::GeV<<" (GeV)"<<G4endl;
  G4cout<<"CALL SBRunAction::BeginOfRunAction EnuMAX = "<<EnuMAX/CLHEP::GeV<<" (GeV)"<<G4endl;

  // Inform the runManager to save random number seed
  // --
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
    
  // Initialize cumulative quantities
  // -- 
  sumETarg = sum2ETarg = sumEGap = sum2EGap = 0.;
  sumLTarg = sum2LTarg = sumLGap = sum2LGap = 0.; 

  NDIF_piplus_TOT 	= 0;
  NDIF_piminus_TOT 	= 0;
  NDIF_Kplus_TOT 	= 0;
  NDIF_Kminus_TOT 	= 0;
  NDIF_KzeroL_TOT 	= 0;
  NDIF_KzeroS_TOT 	= 0;
  
  if(primaryGen->GetJOBID()>=0){
    fROOT = new TFile(Form("output_%06d.root",primaryGen->GetJOBID()),"recreate");
  }else{
    fROOT = new TFile(Form("output_-%06d.root",-primaryGen->GetJOBID()),"recreate");    
  }

 
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
  G4double pmax = 1.5;

  h_pt_exit_piplus 		= new TH1F("h_pt_exit_piplus","h_pt_exit_piplus",100,0,1.);
  h_p_exit_piplus 		= new TH1F("h_p_exit_piplus","h_p_exit_piplus",100,0,pmax);
  //h_p_exit_piplus_W 		= new TH1F("h_p_exit_piplus_W","h_p_exit_piplus_W",100,0,pmax);
  h_xy_exit_piplus 		= new TH2F("h_xy_exit_piplus","h_xy_exit_piplus",100,-100,100,100,-100,100);
  h_z_exit_piplus 		= new TH1F("h_z_exit_piplus","h_z_exit_piplus",500,zMIN,zMAX);
  h_r_exit_piplus 		= new TH1F("h_r_exit_piplus","h_r_exit_piplus",208,0,rMAX);
  h_pVSr_exit_piplus 		= new TH2F("h_pVSr_exit_piplus","h_pVSr_exit_piplus",208,0,rMAX,100,0,1);
  h_thetaVSr_exit_piplus 	= new TH2F("h_thetaVSr_exit_piplus","h_thetaVSr_exit_piplus",208,0,rMAX,100,0,1);
  h_theta_exit_piplus 		= new TH1F("h_theta_exit_piplus","h_theta_exit_piplus",100,0,1.);
  //h_theta_exit_piplus_W 	= new TH1F("h_theta_exit_piplus_W","h_theta_exit_piplus_W",100,0,1.);
  h_pVStheta_exit_piplus 	= new TH2F("h_pVStheta_exit_piplus","h_pVStheta_exit_piplus",100,0.,1.,100,0,pmax);
  //h_pVStheta_exit_piplus_W 	= new TH2F("h_pVStheta_exit_piplus_W","h_pVStheta_exit_piplus_W",100,0.,1.,100,0,pmax);
  
  h_pt_exit_piminus 		= new TH1F("h_pt_exit_piminus","h_pt_exit_piminus",100,0,1.);
  h_p_exit_piminus		= new TH1F("h_p_exit_piminus","h_p_exit_piminus",100,0,pmax);
  //h_p_exit_piminus_W 		= new TH1F("h_p_exit_piminus_W","h_p_exit_piminus_W",100,0,pmax);
  h_xy_exit_piminus 		= new TH2F("h_xy_exit_piminus","h_xy_exit_piminus",100,-100,100,100,-100,100);
  h_z_exit_piminus 		= new TH1F("h_z_exit_piminus","h_z_exit_piminus",500,zMIN,zMAX);
  h_r_exit_piminus 		= new TH1F("h_r_exit_piminus","h_r_exit_piminus",208,0,rMAX);
  h_pVSr_exit_piminus 		= new TH2F("h_pVSr_exit_piminus","h_pVSr_exit_piminus",208,0,rMAX,100,0,1);
  h_thetaVSr_exit_piminus 	= new TH2F("h_thetaVSr_exit_piminus","h_thetaVSr_exit_piminus",208,0,rMAX,100,0,1);
  h_theta_exit_piminus 		= new TH1F("h_theta_exit_piminus","h_theta_exit_piminus",100,0,1.);
  //h_theta_exit_piminus_W 	= new TH1F("h_theta_exit_piminus_W","h_theta_exit_piminus_W",100,0,1.);
  h_pVStheta_exit_piminus 	= new TH2F("h_pVStheta_exit_piminus","h_pVStheta_exit_piminus",100,0.,1.,100,0,pmax);
  //h_pVStheta_exit_piminus_W 	= new TH2F("h_pVStheta_exit_piminus_W","h_pVStheta_exit_piminus_W",100,0.,1.,100,0,pmax);
  
  h_pt_exit_kplus 		= new TH1F("h_pt_exit_kplus","h_pt_exit_kplus",100,0,1.);
  h_p_exit_kplus 		= new TH1F("h_p_exit_kplus","h_p_exit_kplus",100,0,pmax);
  //h_p_exit_kplus_W 		= new TH1F("h_p_exit_kplus_W","h_p_exit_kplus_W",100,0,pmax);
  h_xy_exit_kplus 		= new TH2F("h_xy_exit_kplus","h_xy_exit_kplus",100,-100,100,100,-100,100);
  h_z_exit_kplus 		= new TH1F("h_z_exit_kplus","h_z_exit_kplus",500,zMIN,zMAX);
  h_r_exit_kplus		= new TH1F("h_r_exit_kplus","h_r_exit_kplus",208,0,rMAX);
  h_pVSr_exit_kplus 		= new TH2F("h_pVSr_exit_kplus","h_pVSr_exit_kplus",208,0,rMAX,100,0,1);
  h_thetaVSr_exit_kplus 	= new TH2F("h_thetaVSr_exit_kplus","h_thetaVSr_exit_kplus",208,0,rMAX,100,0,1);
  h_theta_exit_kplus 		= new TH1F("h_theta_exit_kplus","h_theta_exit_kplus",100,0,1.);
  //h_theta_exit_kplus_W 	= new TH1F("h_theta_exit_kplus_W","h_theta_exit_kplus_W",100,0,1.);
  h_pVStheta_exit_kplus 	= new TH2F("h_pVStheta_exit_kplus","h_pVStheta_exit_kplus",100,0.,1.,100,0,pmax);
  
  h_pt_exit_kminus 		= new TH1F("h_pt_exit_kminus","h_pt_exit_kminus",100,0,1.);
  h_p_exit_kminus 		= new TH1F("h_p_exit_kminus","h_p_exit_kminus",100,0,pmax);
  //h_p_exit_kminus_W 		= new TH1F("h_p_exit_kminus_W","h_p_exit_kminus_W",100,0,pmax);
  h_xy_exit_kminus 		= new TH2F("h_xy_exit_kminus","h_xy_exit_kminus",100,-100,100,100,-100,100);
  h_z_exit_kminus 		= new TH1F("h_z_exit_kminus","h_z_exit_kminus",500,zMIN,zMAX);
  h_r_exit_kminus 		= new TH1F("h_r_exit_kminus","h_r_exit_kminus",208,0,rMAX);
  h_pVSr_exit_kminus 		= new TH2F("h_pVSr_exit_kminus","h_pVSr_exit_kminus",208,0,rMAX,100,0,1);
  h_thetaVSr_exit_kminus 	= new TH2F("h_thetaVSr_exit_kminus","h_thetaVSr_exit_kminus",208,0,rMAX,100,0,1);
  h_theta_exit_kminus 		= new TH1F("h_theta_exit_kminus","h_theta_exit_kminus",100,0,1.);
  //h_theta_exit_kminus_W 	= new TH1F("h_theta_exit_kminus_W","h_theta_exit_kminus_W",100,0,1.);
  h_pVStheta_exit_kminus 	= new TH2F("h_pVStheta_exit_kminus","h_pVStheta_exit_kminus",100,0.,1.,100,0,pmax);
  
  h_pt_exit_muplus 		= new TH1F("h_pt_exit_muplus","h_pt_exit_muplus",100,0,1.);
  h_p_exit_muplus 		= new TH1F("h_p_exit_muplus","h_p_exit_muplus",100,0,pmax);
  //h_p_exit_muplus_W 		= new TH1F("h_p_exit_muplus_W","h_p_exit_muplus_W",100,0,pmax);
  h_xy_exit_muplus 		= new TH2F("h_xy_exit_muplus","h_xy_exit_muplus",100,-100,100,100,-100,100);
  h_z_exit_muplus 		= new TH1F("h_z_exit_muplus","h_z_exit_muplus",500,zMIN,zMAX);
  h_r_exit_muplus 		= new TH1F("h_r_exit_muplus","h_r_exit_muplus",208,0,rMAX);
  h_pVSr_exit_muplus 		= new TH2F("h_pVSr_exit_muplus","h_pVSr_exit_muplus",208,0,rMAX,100,0,1);
  h_thetaVSr_exit_muplus 	= new TH2F("h_thetaVSr_exit_muplus","h_thetaVSr_exit_muplus",208,0,rMAX,100,0,1);
  h_theta_exit_muplus 		= new TH1F("h_theta_exit_muplus","h_theta_exit_muplus",100,0,1.);
  //h_theta_exit_muplus_W 	= new TH1F("h_theta_exit_muplus_W","h_theta_exit_muplus_W",100,0,1.);
  h_pVStheta_exit_muplus 	= new TH2F("h_pVStheta_exit_muplus","h_pVStheta_exit_muplus",100,0.,1.,100,0,pmax);
  
  h_pt_exit_muminus 		= new TH1F("h_pt_exit_muminus","h_pt_exit_muminus",100,0,1.);
  h_p_exit_muminus 		= new TH1F("h_p_exit_muminus","h_p_exit_muminus",100,0,pmax);
  //h_p_exit_muminus_W 		= new TH1F("h_p_exit_muminus_W","h_p_exit_muminus_W",100,0,pmax);
  h_xy_exit_muminus 		= new TH2F("h_xy_exit_muminus","h_xy_exit_muminus",100,-100,100,100,-100,100);
  h_z_exit_muminus 		= new TH1F("h_z_exit_muminus","h_z_exit_muminus",500,zMIN,zMAX);
  h_r_exit_muminus 		= new TH1F("h_r_exit_muminus","h_r_exit_muminus",208,0,rMAX);
  h_pVSr_exit_muminus 		= new TH2F("h_pVSr_exit_muminus","h_pVSr_exit_muminus",208,0,rMAX,100,0,1);
  h_thetaVSr_exit_muminus 	= new TH2F("h_thetaVSr_exit_muminus","h_thetaVSr_exit_muminus",208,0,rMAX,100,0,1);
  h_theta_exit_muminus 		= new TH1F("h_theta_exit_muminus","h_theta_exit_muminus",100,0,1.);
  //h_theta_exit_muminus_W 	= new TH1F("h_theta_exit_muminus_W","h_theta_exit_muminus_W",100,0,1.);
  h_pVStheta_exit_muminus 	= new TH2F("h_pVStheta_exit_muminus","h_pVStheta_exit_muminus",100,0.,1.,100,0,pmax);
  
  h_pt_exit_k0L 		= new TH1F("h_pt_exit_k0L","h_pt_exit_k0L",100,0,1.);
  h_p_exit_k0L 			= new TH1F("h_p_exit_k0L","h_p_exit_k0L",100,0,pmax);
  //h_p_exit_k0L_W 		= new TH1F("h_p_exit_k0L_W","h_p_exit_k0L_W",100,0,pmax);
  h_xy_exit_k0L 		= new TH2F("h_xy_exit_k0L","h_xy_exit_k0L",100,-100,100,100,-100,100);
  h_z_exit_k0L 			= new TH1F("h_z_exit_k0L","h_z_exit_k0L",500,zMIN,zMAX);
  h_r_exit_k0L 			= new TH1F("h_r_exit_k0L","h_r_exit_k0L",208,0,rMAX);
  h_pVSr_exit_k0L 		= new TH2F("h_pVSr_exit_k0L","h_pVSr_exit_k0L",208,0,rMAX,100,0,1);
  h_thetaVSr_exit_k0L 		= new TH2F("h_thetaVSr_exit_k0L","h_thetaVSr_exit_k0L",208,0,rMAX,100,0,1);
  h_theta_exit_k0L		= new TH1F("h_theta_exit_k0L","h_theta_exit_k0L",100,0,1.);
  //h_theta_exit_k0L_W 		= new TH1F("h_theta_exit_k0L_W","h_theta_exit_k0L_W",100,0,1.);
  h_pVStheta_exit_k0L 		= new TH2F("h_pVStheta_exit_k0L","h_pVStheta_exit_k0L",100,0.,1.,100,0,pmax);

  // At target level
  // --
  h_pt_targ_piplus 		= new TH1F("h_pt_targ_piplus","h_pt_targ_piplus",100,0,1.);
  h_p_targ_piplus 		= new TH1F("h_p_targ_piplus","h_p_targ_piplus",100,0,pmax);
  h_z_targ_piplus 		= new TH1F("h_z_targ_piplus","h_z_targ_piplus",100,0,90);
  h_zabs_targ_piplus 		= new TH1F("h_zabs_targ_piplus","h_zabs_targ_piplus",500,zMIN,zMAX);
  h_theta_targ_piplus 		= new TH1F("h_theta_targ_piplus","h_theta_targ_piplus",100,0,1.);
  h_pVStheta_targ_piplus 	= new TH2F("h_pVStheta_targ_piplus","h_pVStheta_targ_piplus",100,0.,1.,100,0,pmax);

  for(int i=0;i<NENB;i++){
    h_pVStheta_targ_piplus_Enubins[i] 	= new TH2F(Form("h_pVStheta_targ_piplus_Enubins%d",i),
						Form("h_pVStheta_targ_piplus_Enubins%d",i),100,0.,1.,100,0,3.);
    h_p_targ_piplus_Enubins[i] 		= new TH1F(Form("h_p_targ_piplus_Enubins%d",i),
					  	Form("h_p_targ_piplus_Enubins%d",i),100,0,3.);
    h_theta_targ_piplus_Enubins[i] 	= new TH1F(Form("h_theta_targ_piplus_Enubins%d",i),
					       	Form("h_theta_targ_piplus_Enubins%d",i),100,0,1.);
  }

  h_pt_targ_piminus 		= new TH1F("h_pt_targ_piminus","h_pt_targ_piminus",100,0,1.);
  h_p_targ_piminus 		= new TH1F("h_p_targ_piminus","h_p_targ_piminus",100,0,pmax);
  h_z_targ_piminus		= new TH1F("h_z_targ_piminus","h_z_targ_piminus",100,0,90);
  h_zabs_targ_piminus 		= new TH1F("h_zabs_targ_piminus","h_zabs_targ_piminus",500,zMIN,zMAX);
  h_theta_targ_piminus 		= new TH1F("h_theta_targ_piminus","h_theta_targ_piminus",100,0,1.);
  h_pVStheta_targ_piminus 	= new TH2F("h_pVStheta_targ_piminus","h_pVStheta_targ_piminus",100,0.,1.,100,0,pmax);
  for(int i=0;i<NENB;i++){
    h_pVStheta_targ_piminus_Enubins[i] 	= new TH2F(Form("h_pVStheta_targ_piminus_Enubins%d",i),
						Form("h_pVStheta_targ_piminus_Enubins%d",i),100,0.,1.,100,0,3.);
    h_p_targ_piminus_Enubins[i] 	= new TH1F(Form("h_p_targ_piminus_Enubins%d",i),
					   	Form("h_p_targ_piminus_Enubins%d",i),100,0,3.);
    h_theta_targ_piminus_Enubins[i] 	= new TH1F(Form("h_theta_targ_piminus_Enubins%d",i),
					       	Form("h_theta_targ_piminus_Enubins%d",i),100,0,1.);
  }

  // Muon decay
  // --
  h_MUDEC_Enumustar 		= new TH1F("h_MUDEC_Enumustar","h_MUDEC_Enumustar",100,0.,55.);//MeV
  h_MUDEC_Enuestar 		= new TH1F("h_MUDEC_Enuestar","h_MUDEC_Enuestar",100,0.,55.);//MeV
  h_MUDEC_Eanumustar 		= new TH1F("h_MUDEC_Eanumustar","h_MUDEC_Eanumustar",100,0.,55.);//MeV
  h_MUDEC_Eanuestar 		= new TH1F("h_MUDEC_Eanuestar","h_MUDEC_Eanuestar",100,0.,55.);//MeV

  h_MUDEC_cthstarMu 		= new TH1F("h_MUDEC_cthstarMu","h_MUDEC_cthstarMu",100,-1.,1.);
  h_MUDEC_cthstarNu		= new TH1F("h_MUDEC_cthstarNu","h_MUDEC_cthstarNu",100,-1.,1.);

  h_MUDEC_probaTunnel 		= new TH1F("h_MUDEC_probaTunnel","h_MUDEC_probaTunnel",100,0,0.05);
  h_MUDEC_pathTunnel 		= new TH1F("h_MUDEC_pathTunnel","h_MUDEC_pathTunnel",100,0,50);

  h_MUDEC_probaTunnel_pmu 	= new TH2F("h_MUDEC_probaTunnel_pmu","h_MUDEC_probaTunnel_pmu",100,0,3,100,0,0.05);
  h_MUDEC_probaTunnel_thmu 	= new TH2F("h_MUDEC_probaTunnel_thmu","h_MUDEC_probaTunnel_thmu",100,0,2,100,0,0.05);
  h_MUDEC_probaTunnel_rmu 	= new TH2F("h_MUDEC_probaTunnel_rmu","h_MUDEC_probaTunnel_rmu",100,0,2.5,100,0,0.05);

  h_MUDEC_muxy 			= new TH2F("h_MUDEC_muxy","h_MUDEC_muxy",100,-3.,3.,100,-3.,3.);//m
  h_MUDEC_muzr 			= new TH2F("h_MUDEC_muzr","h_MUDEC_muzr",100,zMIN,zMAX,100,0.,3.);//m

  h_MUDEC_polaL_muplus 		= new TH1F("h_MUDEC_polaL_muplus","h_MUDEC_polaL_muplus",100,-1.,1.);
  h_MUDEC_polaL_muminus 	= new TH1F("h_MUDEC_polaL_muminus","h_MUDEC_polaL_muminus",100,-1.,1.);
  h_MUDEC_thpimu_lab 		= new TH1F("h_MUDEC_thpimu_lab","h_MUDEC_thpimu_lab",100,0,1.);

  h_MUDEC_pmupar 		= new TH2F("h_MUDEC_pmupar","h_MUDEC_pmupar",100,0.,3.,100,0.,3.);//GeV

  // to check the geometry and compare to g3 through hits map
  h_geom_g4 			= new TH2F("h_geom_g4","h_geom_g4",2000,-22,-18,2000,0,2.1);
  //h_geom_g4 			= new TH2F("h_geom_g4","h_geom_g4",2000,-23,-18,2000,0,2.5);
  //h_geom_g4 			= new TH2F("h_geom_g4","h_geom_g4",2000,-23,23,2000,0,2.5);

  //h_decay_flags = new TH1F("h_decay_flags","h_decay_flags",23,-1.5,20.5);
  h_decay_flags 		= new TH1F("h_decay_flags","h_decay_flags",200,0,32);
  h_decay_NUM 			= new TH1F("h_decay_NUM","h_decay_NUM",200,0,32);
  h_decay_BR 			= new TH1F("h_decay_BR","h_decay_BR",200,0,32);
  h_decay_BR_PDG 		= new TH1F("h_decay_BR_PDG","h_decay_BR_PDG",200,0,32);

  h_decay_BR_PDG->Fill(30,99.9877);  
  h_decay_BR_PDG->Fill(31,99.9877);  
  h_decay_BR_PDG->Fill(1,63.51);  
  h_decay_BR_PDG->Fill(2,63.51);  
  h_decay_BR_PDG->Fill(3,21.17);  
  h_decay_BR_PDG->Fill(4,21.17);  
  h_decay_BR_PDG->Fill(5,5.59);  
  h_decay_BR_PDG->Fill(6,5.59);  
  h_decay_BR_PDG->Fill(7,4.82);  
  h_decay_BR_PDG->Fill(8,4.82);  
  h_decay_BR_PDG->Fill(9,3.18);  
  h_decay_BR_PDG->Fill(10,3.18);  
  h_decay_BR_PDG->Fill(11,1.73);  
  h_decay_BR_PDG->Fill(12,1.73);  
  h_decay_BR_PDG->Fill(13,19.35);  
  h_decay_BR_PDG->Fill(14,19.35);  
  h_decay_BR_PDG->Fill(15,13.5);  
  h_decay_BR_PDG->Fill(16,13.5);  
  h_decay_BR_PDG->Fill(17,21.5);  
  h_decay_BR_PDG->Fill(18,12.38);  
  h_decay_BR_PDG->Fill(19,68.61);  
  h_decay_BR_PDG->Fill(20,31.39);  

  //h_prim_match_pi 		= new TH2F("h_prim_match_pi","h_prim_match_pi",200,0.,3.,200,0.,3.);
  //h_prim_match_ka 		= new TH2F("h_prim_match_ka","h_prim_match_ka",200,0.,3.,200,0.,3.);
  //h_prim_match_k0 		= new TH2F("h_prim_match_k0","h_prim_match_k0",200,0.,3.,200,0.,3.);

  h_prim_match_pi 		= new TH2F("h_prim_match_pi","h_prim_match_pi",200,-0.05,0.05,200,0.98,1.02);
  h_prim_match_ka 		= new TH2F("h_prim_match_ka","h_prim_match_ka",200,-0.05,0.05,200,0.98,1.02);
  h_prim_match_k0 		= new TH2F("h_prim_match_k0","h_prim_match_k0",200,-0.05,0.05,200,0.98,1.02);
  
  SBAnalysisManager* analysis = SBAnalysisManager::getInstance();
  analysis->BeginOfRun();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SBRunAction::CheckEOF(){
  if(primaryGen->myEOF){
    G4cout << "********************************************************" << G4endl;
    G4cout << "SBRunAction::CheckEOF() : EOF detected. Run will finish." << G4endl;
    G4cout << "********************************************************" << G4endl;
    
    G4RunManager::GetRunManager()->AbortRun();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SBRunAction::fillPerEvent(G4double ETarg,
			       G4double EGap,
			       G4double LTarg,
			       G4double LGap)
{

  G4cout<<"CALL SBRunAction::fillPerEvent "<<G4endl;

  // Accumulate statistic
  // --
  sumETarg += ETarg;  sum2ETarg += ETarg*ETarg;
  sumEGap += EGap;  sum2EGap += EGap*EGap;
  
  sumLTarg += LTarg;  sum2LTarg += LTarg*LTarg;
  sumLGap += LGap;  sum2LGap += LGap*LGap;  

  SBAnalysisManager* analysis = SBAnalysisManager::getInstance();
  analysis->FillHisto(ETarg,LTarg);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SBRunAction::DIFStat(G4double ndifpiplus,
			  G4double ndifpiminus,
			  G4double ndifkplus,
			  G4double ndifkminus,
			  G4double ndifk0S,
			  G4double ndifk0L)
{

  G4cout<<"CALL SBRunAction::DIFStat "<<G4endl;

  // Accumulate statistic
  // --
  NDIF_piplus_TOT += ndifpiplus;
  NDIF_piminus_TOT += ndifpiminus;
  NDIF_Kplus_TOT += ndifkplus;
  NDIF_Kminus_TOT += ndifkminus;
  NDIF_KzeroS_TOT += ndifk0S;
  NDIF_KzeroL_TOT += ndifk0L;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SBRunAction::fillPerEvent1(G4double p_muplus_PI,
			        G4double p_muminus_PI,
			        G4double p_muplus_K,
			        G4double p_muminus_K,
			        G4double E_numu_PI,
			        G4double E_anumu_PI,
			        G4double E_numu_K,
			        G4double E_anumu_K)
{
  G4cout<<"CALL SBRunAction::fillPerEvent1 "<<G4endl;
 
  SBAnalysisManager* analysis = SBAnalysisManager::getInstance();
  analysis->FillHisto2(p_muplus_PI,p_muminus_PI,p_muplus_K,p_muminus_K,E_numu_PI,E_anumu_PI,E_numu_K,E_anumu_K);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SBRunAction::fillDIF(G4double x,
			  G4double y,
			  G4double z,
			  G4double px,
			  G4double py,
			  G4double pz,
			  G4double pxnu,
			  G4double pynu,
			  G4double pznu,
			  G4int moth,
			  G4int q)
{
  G4cout<<"CALL SBRunAction::fillDIF "<<G4endl;
 
  SBAnalysisManager* analysis = SBAnalysisManager::getInstance();
  analysis->FillGeometryDecay(x,y,z,px,py,pz,pxnu,pynu,pznu,moth,q);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SBRunAction::fillEXIT(G4double x,
			   G4double y,
			   G4double z,
			   G4double px,
			   G4double py,
			   G4double pz,
			   G4int moth)
{
  G4cout<<"CALL SBRunAction::fillEXIT "<<G4endl;
  G4double r = sqrt(x*x+y*y);
  G4double p = sqrt(px*px+py*py+pz*pz);
  G4double pt = sqrt(px*px+py*py);
  G4double th = 0;
  if(p)th=acos(pz/p);

  if((p!=0.)&&(x!=0.)&&(y!=0.)&&(z!=0.)){

    if(r<3/CLHEP::cm)G4cout <<"SBRunAction: EXIT POINT x "<< 
	     x/CLHEP::cm <<" y "<< y/CLHEP::cm <<" r "<< r/CLHEP::cm<<" (cm) "
		  <<" p "<<p/CLHEP::GeV<<" (GeV) z "<<z/CLHEP::m<< " (m) "
		  <<G4endl;

/////////////////////////////////////////////////

    if(moth==1){//piplus
      h_pt_exit_piplus->Fill(pt/CLHEP::GeV);
      h_p_exit_piplus->Fill(p/CLHEP::GeV);
      h_xy_exit_piplus->Fill(x/CLHEP::cm,y/CLHEP::cm);
      h_z_exit_piplus->Fill(z/CLHEP::m);
      //G4cout <<" crash " <<z/m <<G4endl;
      h_r_exit_piplus->Fill(r/CLHEP::cm);
      h_pVSr_exit_piplus->Fill(p/CLHEP::GeV,r/CLHEP::cm);
      h_pVStheta_exit_piplus->Fill(th,p/CLHEP::GeV);
      h_thetaVSr_exit_piplus->Fill(th,r/CLHEP::cm);
      h_theta_exit_piplus->Fill(th);
    } else if (moth==2){//piminus
      h_pt_exit_piminus->Fill(pt/CLHEP::GeV);
      h_p_exit_piminus->Fill(p/CLHEP::GeV);
      h_xy_exit_piminus->Fill(x/CLHEP::cm,y/CLHEP::cm);
      h_z_exit_piminus->Fill(z/CLHEP::m);
      h_r_exit_piminus->Fill(r/CLHEP::cm);
      h_pVSr_exit_piminus->Fill(p/CLHEP::GeV,r/CLHEP::cm);
      h_pVStheta_exit_piminus->Fill(th,p/CLHEP::GeV);
      h_thetaVSr_exit_piminus->Fill(th,r/CLHEP::cm);
      h_theta_exit_piminus->Fill(th);
    }else if(moth==3){//muplus
      h_pt_exit_muplus->Fill(pt/CLHEP::GeV);
      h_p_exit_muplus->Fill(p/CLHEP::GeV);
      h_xy_exit_muplus->Fill(x/CLHEP::cm,y/CLHEP::cm);
      h_z_exit_muplus->Fill(z/CLHEP::m);
      h_r_exit_muplus->Fill(r/CLHEP::cm);
      h_pVSr_exit_muplus->Fill(p/CLHEP::GeV,r/CLHEP::cm);
      h_pVStheta_exit_muplus->Fill(th,p/CLHEP::GeV);
      h_thetaVSr_exit_muplus->Fill(th,r/CLHEP::cm);
      h_theta_exit_muplus->Fill(th);    
    } else if (moth==4){//muminus
      h_pt_exit_muminus->Fill(pt/CLHEP::GeV);
      h_p_exit_muminus->Fill(p/CLHEP::GeV);
      h_xy_exit_muminus->Fill(x/CLHEP::cm,y/CLHEP::cm);
      h_z_exit_muminus->Fill(z/CLHEP::m);
      h_r_exit_muminus->Fill(r/CLHEP::cm);
      h_pVSr_exit_muminus->Fill(p/CLHEP::GeV,r/CLHEP::cm);
      h_pVStheta_exit_muminus->Fill(th,p/CLHEP::GeV);
      h_thetaVSr_exit_muminus->Fill(th,r/CLHEP::cm);
      h_theta_exit_muminus->Fill(th);    
    } else if(moth==5){//kplus
      h_pt_exit_kplus->Fill(pt/CLHEP::GeV);
      h_p_exit_kplus->Fill(p/CLHEP::GeV);
      h_xy_exit_kplus->Fill(x/CLHEP::cm,y/CLHEP::cm);
      h_z_exit_kplus->Fill(z/CLHEP::m);
      h_r_exit_kplus->Fill(r/CLHEP::cm);
      h_pVSr_exit_kplus->Fill(p/CLHEP::GeV,r/CLHEP::cm);
      h_pVStheta_exit_kplus->Fill(th,p/CLHEP::GeV);
      h_thetaVSr_exit_kplus->Fill(th,r/CLHEP::cm);
      h_theta_exit_kplus->Fill(th);    
    } else if (moth==6){//kminus
      h_pt_exit_kminus->Fill(pt/CLHEP::GeV);
      h_p_exit_kminus->Fill(p/CLHEP::GeV);
      h_xy_exit_kminus->Fill(x/CLHEP::cm,y/CLHEP::cm);
      h_z_exit_kminus->Fill(z/CLHEP::m);
      h_r_exit_kminus->Fill(r/CLHEP::cm);
      h_pVSr_exit_kminus->Fill(p/CLHEP::GeV,r/CLHEP::cm);
      h_pVStheta_exit_kminus->Fill(th,p/CLHEP::GeV);
      h_thetaVSr_exit_kminus->Fill(th,r/CLHEP::cm);
      h_theta_exit_kminus->Fill(th);    
    } else  if(moth==7){//k0L
      h_pt_exit_k0L->Fill(pt/CLHEP::GeV);
      h_p_exit_k0L->Fill(p/CLHEP::GeV);
      h_xy_exit_k0L->Fill(x/CLHEP::cm,y/CLHEP::cm);
      h_z_exit_k0L->Fill(z/CLHEP::m);
      h_r_exit_k0L->Fill(r/CLHEP::cm);
      h_pVSr_exit_k0L->Fill(p/CLHEP::GeV,r/CLHEP::cm);
      h_pVStheta_exit_k0L->Fill(th,p/CLHEP::GeV);
      h_thetaVSr_exit_k0L->Fill(th,r/CLHEP::cm);
      h_theta_exit_k0L->Fill(th);    
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SBRunAction::fillTARG(/*G4double x,
			   G4double y,
			   G4double z,
			   G4double px,
			   G4double py,
			   G4double pz,
			   G4int moth*/
			   G4double enu=0)
{
  G4cout<<"CALL SBRunAction::fillTARG "<<G4endl;
  G4double x=primaryGen->XPrimary;
  G4double y=primaryGen->YPrimary;
  G4double z=primaryGen->ZPrimary;
  G4double zabs=primaryGen->ZPrimary_G4frame;
  G4double px=primaryGen->PXPrimary;
  G4double py=primaryGen->PYPrimary;
  G4double pz=primaryGen->PZPrimary;

  G4double p = sqrt((px*px)+(py*py)+(pz*pz));
  G4double pt = sqrt((px*px)+(py*py));
  G4double th = 0;
  if(p)th=acos(pz/p);
  
  if((p!=0.)&&(x!=0.)&&(y!=0.)&&(z!=0.)){
    if(primaryGen->particleName=="pi+"){
      h_pt_targ_piplus->Fill(pt/CLHEP::GeV);
      h_p_targ_piplus->Fill(p/CLHEP::GeV);
      h_z_targ_piplus->Fill(z/CLHEP::cm);
      h_zabs_targ_piplus->Fill(zabs/CLHEP::m);
      h_pVStheta_targ_piplus->Fill(th,p/CLHEP::GeV);
      h_theta_targ_piplus->Fill(th);

      for(int j=0;j<NENB;j++){
	if((enu/CLHEP::GeV>ENULIMS[j])&&(enu/CLHEP::GeV<ENULIMS[j+1])){
	  h_pVStheta_targ_piplus_Enubins[j]->Fill(th,p/CLHEP::GeV);
	  h_theta_targ_piplus_Enubins[j]->Fill(th);
	  h_p_targ_piplus_Enubins[j]->Fill(p/CLHEP::GeV);
	}
      }

    } else if(primaryGen->particleName=="pi-"){
      h_pt_targ_piminus->Fill(pt/CLHEP::GeV);
      h_p_targ_piminus->Fill(p/CLHEP::GeV);
      h_z_targ_piminus->Fill(z/CLHEP::cm);
      h_zabs_targ_piminus->Fill(zabs/CLHEP::m);
      h_pVStheta_targ_piminus->Fill(th,p/CLHEP::GeV);
      h_theta_targ_piminus->Fill(th);

      for(int j=0;j<NENB;j++){
	if((enu/CLHEP::GeV>ENULIMS[j])&&(enu/CLHEP::GeV<ENULIMS[j+1])){
	  h_pVStheta_targ_piminus_Enubins[j]->Fill(th,p/CLHEP::GeV);
	  h_theta_targ_piminus_Enubins[j]->Fill(th);
	  h_p_targ_piminus_Enubins[j]->Fill(p/CLHEP::GeV);
	}
      }

    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SBRunAction::EndOfRunAction(const G4Run* aRun)
{
  G4cout<<"CALL SBRunAction::EndOfRunAction "<<G4endl;
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents == 0) return;
  
  DoFinalOps(NbOfEvents);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SBRunAction::SumFluxesFlavours()
{
  SBAnalysisManager* analysis = SBAnalysisManager::getInstance();
  analysis->FillSumFluxesFlavours();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SBRunAction::DoFinalOps(G4int NbOfEvents){
  
  G4cout<<"SBRunAction::DoFinalOps"<<G4endl;

  // Compute statistics: mean and rms
  // -- 
  sumETarg /= NbOfEvents; sum2ETarg /= NbOfEvents;
  G4double rmsETarg = sum2ETarg - sumETarg*sumETarg;
  if (rmsETarg >0.) rmsETarg = std::sqrt(rmsETarg); else rmsETarg = 0.;
  
  sumEGap /= NbOfEvents; sum2EGap /= NbOfEvents;
  G4double rmsEGap = sum2EGap - sumEGap*sumEGap;
  if (rmsEGap >0.) rmsEGap = std::sqrt(rmsEGap); else rmsEGap = 0.;
  
  sumLTarg /= NbOfEvents; sum2LTarg /= NbOfEvents;
  G4double rmsLTarg = sum2LTarg - sumLTarg*sumLTarg;
  if (rmsLTarg >0.) rmsLTarg = std::sqrt(rmsLTarg); else rmsLTarg = 0.;
  
  sumLGap /= NbOfEvents; sum2LGap /= NbOfEvents;
  G4double rmsLGap = sum2LGap - sumLGap*sumLGap;
  if (rmsLGap >0.) rmsLGap = std::sqrt(rmsLGap); else rmsLGap = 0.;
  
  /**************************************/
  // normalize to number of processed pots !
  G4double ProcEvtsScale = 1.e6/primaryGen->GetFLUKAPOTs();
  G4double K_REP = primaryGen->K_REPL;

  G4cout << "SBRunAction: K replication factor used in fluxes scaling "<< K_REP << G4endl;
  /////////////
  SBAnalysisManager* analysis = SBAnalysisManager::getInstance();
  analysis->FillFinalOps(ProcEvtsScale,K_REP);
  ////////////
    
  G4cout << "Processed FLUKA events is " << primaryGen->GetFLUKAPOTs()<<"."<< G4endl;
  G4cout << "Normalization factor " << ProcEvtsScale <<" applied to fluxes."<< G4endl;  

  G4cout << "decay in flight summary: "<< 
    " pi+ "<<NDIF_piplus_TOT<<
    " pi- "<<NDIF_piminus_TOT<<
    " K+ "<<NDIF_Kplus_TOT<<
    " K- "<<NDIF_Kminus_TOT<<
    " K0L "<<NDIF_KzeroL_TOT<<
    " K0S "<<NDIF_KzeroS_TOT<<G4endl;

  G4cout << "    per/p.o.t.: "<< 
    " pi+ "<<NDIF_piplus_TOT/primaryGen->GetFLUKAPOTs()<<
    " pi- "<<NDIF_piminus_TOT/primaryGen->GetFLUKAPOTs()<<
    " K+ "<<NDIF_Kplus_TOT/primaryGen->GetFLUKAPOTs()<<
    " K- "<<NDIF_Kminus_TOT/primaryGen->GetFLUKAPOTs()<<G4endl;

  /********/

  h_pt_exit_piplus->Scale(ProcEvtsScale);
  h_p_exit_piplus->Scale(ProcEvtsScale);
  h_xy_exit_piplus->Scale(ProcEvtsScale);
  h_z_exit_piplus->Scale(ProcEvtsScale);
  h_r_exit_piplus->Scale(ProcEvtsScale);
  h_pVSr_exit_piplus->Scale(ProcEvtsScale);
  h_pVStheta_exit_piplus->Scale(ProcEvtsScale);
  h_thetaVSr_exit_piplus->Scale(ProcEvtsScale);
  h_theta_exit_piplus->Scale(ProcEvtsScale);
  
  h_pt_exit_piminus->Scale(ProcEvtsScale);
  h_p_exit_piminus->Scale(ProcEvtsScale);
  h_xy_exit_piminus->Scale(ProcEvtsScale);
  h_z_exit_piminus->Scale(ProcEvtsScale);
  h_r_exit_piminus->Scale(ProcEvtsScale);
  h_pVSr_exit_piminus->Scale(ProcEvtsScale);
  h_pVStheta_exit_piminus->Scale(ProcEvtsScale);
  h_thetaVSr_exit_piminus->Scale(ProcEvtsScale);
  h_theta_exit_piminus->Scale(ProcEvtsScale);
  
  h_pt_exit_kplus->Scale(ProcEvtsScale);
  h_p_exit_kplus->Scale(ProcEvtsScale);
  h_xy_exit_kplus->Scale(ProcEvtsScale);
  h_z_exit_kplus->Scale(ProcEvtsScale);
  h_r_exit_kplus->Scale(ProcEvtsScale);
  h_pVSr_exit_kplus->Scale(ProcEvtsScale);
  h_pVStheta_exit_kplus->Scale(ProcEvtsScale);
  h_thetaVSr_exit_kplus->Scale(ProcEvtsScale);
  h_theta_exit_kplus->Scale(ProcEvtsScale);

  h_pt_exit_kminus->Scale(ProcEvtsScale);
  h_p_exit_kminus->Scale(ProcEvtsScale);
  h_xy_exit_kminus->Scale(ProcEvtsScale);
  h_z_exit_kminus->Scale(ProcEvtsScale);
  h_r_exit_kminus->Scale(ProcEvtsScale);
  h_pVSr_exit_kminus->Scale(ProcEvtsScale);
  h_pVStheta_exit_kminus->Scale(ProcEvtsScale);
  h_thetaVSr_exit_kminus->Scale(ProcEvtsScale);
  h_theta_exit_kminus->Scale(ProcEvtsScale);
  
  h_pt_exit_muplus->Scale(ProcEvtsScale);
  h_p_exit_muplus->Scale(ProcEvtsScale);
  h_xy_exit_muplus->Scale(ProcEvtsScale);
  h_z_exit_muplus->Scale(ProcEvtsScale);
  h_r_exit_muplus->Scale(ProcEvtsScale);
  h_pVSr_exit_muplus->Scale(ProcEvtsScale);
  h_pVStheta_exit_muplus->Scale(ProcEvtsScale);
  h_thetaVSr_exit_muplus->Scale(ProcEvtsScale);
  h_theta_exit_muplus->Scale(ProcEvtsScale);
  
  h_pt_exit_muminus->Scale(ProcEvtsScale);
  h_p_exit_muminus->Scale(ProcEvtsScale);
  h_xy_exit_muminus->Scale(ProcEvtsScale);
  h_z_exit_muminus->Scale(ProcEvtsScale);
  h_r_exit_muminus->Scale(ProcEvtsScale);
  h_pVSr_exit_muminus->Scale(ProcEvtsScale);
  h_pVStheta_exit_muminus->Scale(ProcEvtsScale);
  h_thetaVSr_exit_muminus->Scale(ProcEvtsScale);
  h_theta_exit_muminus->Scale(ProcEvtsScale);
  
  h_pt_exit_k0L->Scale(ProcEvtsScale);
  h_p_exit_k0L->Scale(ProcEvtsScale);
  h_xy_exit_k0L->Scale(ProcEvtsScale);
  h_z_exit_k0L->Scale(ProcEvtsScale);
  h_r_exit_k0L->Scale(ProcEvtsScale);
  h_pVSr_exit_k0L->Scale(ProcEvtsScale);
  h_pVStheta_exit_k0L->Scale(ProcEvtsScale);
  h_thetaVSr_exit_k0L->Scale(ProcEvtsScale);
  h_theta_exit_k0L->Scale(ProcEvtsScale);

  // target level

  h_pt_targ_piplus->Scale(ProcEvtsScale);
  h_p_targ_piplus->Scale(ProcEvtsScale);
  h_z_targ_piplus->Scale(ProcEvtsScale);
  h_zabs_targ_piplus->Scale(ProcEvtsScale);
  h_pVStheta_targ_piplus->Scale(ProcEvtsScale);
  h_theta_targ_piplus->Scale(ProcEvtsScale);
  
  h_pt_targ_piminus->Scale(ProcEvtsScale);
  h_p_targ_piminus->Scale(ProcEvtsScale);
  h_z_targ_piminus->Scale(ProcEvtsScale);
  h_zabs_targ_piminus->Scale(ProcEvtsScale);
  h_pVStheta_targ_piminus->Scale(ProcEvtsScale);
  h_theta_targ_piminus->Scale(ProcEvtsScale);

  // DIF                                       //////////////////////////////////////////////////////////////////////////////////////////
/*
  h_piminus_DIF_r->Scale(ProcEvtsScale);
  h_piminus_DIF_z->Scale(ProcEvtsScale);
  h_piminus_DIF_xy->Scale(ProcEvtsScale);
  h_piminus_DIF_zr->Scale(ProcEvtsScale);
  
  h_Kminus_DIF_r->Scale(ProcEvtsScale);
  h_Kminus_DIF_z->Scale(ProcEvtsScale);
  h_Kminus_DIF_xy->Scale(ProcEvtsScale);
  h_Kminus_DIF_zr->Scale(ProcEvtsScale);
  
  h_piplus_DIF_r->Scale(ProcEvtsScale);
  h_piplus_DIF_z->Scale(ProcEvtsScale);
  h_piplus_DIF_xy->Scale(ProcEvtsScale);
  h_piplus_DIF_zr->Scale(ProcEvtsScale);
  
  h_Kplus_DIF_r->Scale(ProcEvtsScale);
  h_Kplus_DIF_z->Scale(ProcEvtsScale);
  h_Kplus_DIF_xy->Scale(ProcEvtsScale);
  h_Kplus_DIF_zr->Scale(ProcEvtsScale);
*/  
  /******************************************************************************************************************************/


  SumFluxesFlavours();

  //h_numu->Scale(ProcEvtsScale);
  //h_anumu->Scale(ProcEvtsScale);
  //h_nue->Scale(ProcEvtsScale);
  //h_anue->Scale(ProcEvtsScale);

  // ------------------- SPL
  G4double BL = 130.;
  G4String fnam="";
  if(primaryGen->GetJOBID()>=0){
    fnam = Form("nufl_%06d_GLOBES_%3.0fKm.txt",primaryGen->GetJOBID(),BL);
  }else{
    fnam = Form("nufl_-%06d_GLOBES_%3.0fKm.txt",-primaryGen->GetJOBID(),BL);
  }
  analysis->PrepareGLoBESFlux(fnam,BL);
  
  //-------------------- SLANIC
  BL = 1544.;
  if(primaryGen->GetJOBID()>=0){
    fnam = Form("nufl_%06d_GLOBES_%4.0fKm.txt",primaryGen->GetJOBID(),BL);
  }else{
    fnam = Form("nufl_-%06d_GLOBES_%4.0fKm.txt",-primaryGen->GetJOBID(),BL);
  }
  analysis->PrepareGLoBESFlux(fnam,BL);

  //-------------------- GENERIC 100Km
  BL = 100.;
  if(primaryGen->GetJOBID()>=0){
    fnam = Form("nufl_%06d_GLOBES_%3.0fKm.txt",primaryGen->GetJOBID(),BL);
  }else{
    fnam = Form("nufl_-%06d_GLOBES_%3.0fKm.txt",-primaryGen->GetJOBID(),BL);
  }
  analysis->PrepareGLoBESFlux(fnam,BL);
 //---------------------


  //FluxesStats();
  
  bool WriteHistogramsToFile = true;
  
  if (WriteHistogramsToFile) {
    
    // fROOT->cd();
   //
   analysis->EndOfRun();
   

    /////////////////////////////////////

    h_pt_exit_piplus->Write();
    h_p_exit_piplus->Write();
    h_xy_exit_piplus->Write();
    h_z_exit_piplus->Write();
    h_r_exit_piplus->Write();
    h_pVSr_exit_piplus->Write();
    h_pVStheta_exit_piplus->Write();
    h_thetaVSr_exit_piplus->Write();
    h_theta_exit_piplus->Write();

    h_pt_exit_piminus->Write();
    h_p_exit_piminus->Write();
    h_xy_exit_piminus->Write();
    h_z_exit_piminus->Write();
    h_r_exit_piminus->Write();
    h_pVSr_exit_piminus->Write();
    h_pVStheta_exit_piminus->Write();
    h_thetaVSr_exit_piminus->Write();
    h_theta_exit_piminus->Write();

    h_pt_exit_kplus->Write();
    h_p_exit_kplus->Write();
    h_xy_exit_kplus->Write();
    h_z_exit_kplus->Write();
    h_r_exit_kplus->Write();
    h_pVSr_exit_kplus->Write();
    h_pVStheta_exit_kplus->Write();
    h_thetaVSr_exit_kplus->Write();
    h_theta_exit_kplus->Write();

    h_pt_exit_kminus->Write();
    h_p_exit_kminus->Write();
    h_xy_exit_kminus->Write();
    h_z_exit_kminus->Write();
    h_r_exit_kminus->Write();
    h_pVSr_exit_kminus->Write();
    h_pVStheta_exit_kminus->Write();
    h_thetaVSr_exit_kminus->Write();
    h_theta_exit_kminus->Write();
 
    h_pt_exit_muplus->Write();
    h_p_exit_muplus->Write();
    h_xy_exit_muplus->Write();
    h_z_exit_muplus->Write();
    h_r_exit_muplus->Write();
    h_pVSr_exit_muplus->Write();
    h_pVStheta_exit_muplus->Write();
    h_thetaVSr_exit_muplus->Write();
    h_theta_exit_muplus->Write();

    h_pt_exit_muminus->Write();
    h_p_exit_muminus->Write();
    h_xy_exit_muminus->Write();
    h_z_exit_muminus->Write();
    h_r_exit_muminus->Write();
    h_pVSr_exit_muminus->Write();
    h_pVStheta_exit_muminus->Write();
    h_thetaVSr_exit_muminus->Write();
    h_theta_exit_muminus->Write();

    h_pt_exit_k0L->Write();
    h_p_exit_k0L->Write();
    h_xy_exit_k0L->Write();
    h_z_exit_k0L->Write();
    h_r_exit_k0L->Write();
    h_pVSr_exit_k0L->Write();
    h_pVStheta_exit_k0L->Write();
    h_thetaVSr_exit_k0L->Write();
    h_theta_exit_k0L->Write();

    h_pt_targ_piplus->Write();
    h_p_targ_piplus->Write();
    h_z_targ_piplus->Write();
    h_zabs_targ_piplus->Write();
    h_theta_targ_piplus->Write();
    h_pVStheta_targ_piplus->Write();
    
    for(int i=0;i<NENB;i++){
      h_pVStheta_targ_piplus_Enubins[i]->Write();
      h_p_targ_piplus_Enubins[i]->Write();
      h_theta_targ_piplus_Enubins[i]->Write();
    }

    h_pt_targ_piminus->Write();
    h_p_targ_piminus->Write();
    h_z_targ_piminus->Write();
    h_zabs_targ_piminus->Write();
    h_theta_targ_piminus->Write();
    h_pVStheta_targ_piminus->Write();

    for(int i=0;i<NENB;i++){
      h_pVStheta_targ_piminus_Enubins[i]->Write();
      h_p_targ_piminus_Enubins[i]->Write();
      h_theta_targ_piminus_Enubins[i]->Write();
    }

    // muon decay
    h_MUDEC_probaTunnel->Write();
    h_MUDEC_pathTunnel->Write();

    h_MUDEC_probaTunnel_pmu->Write();
    h_MUDEC_probaTunnel_thmu->Write();
    h_MUDEC_probaTunnel_rmu->Write();

    h_MUDEC_Enumustar->Write();
    h_MUDEC_Enuestar->Write();
    h_MUDEC_Eanumustar->Write();
    h_MUDEC_Eanuestar->Write();

    h_MUDEC_muxy->Write();
    h_MUDEC_muzr->Write();
    
    h_MUDEC_cthstarMu->Write();
    h_MUDEC_cthstarNu->Write();

    h_MUDEC_polaL_muplus->Write();
    h_MUDEC_polaL_muminus->Write();
    h_MUDEC_thpimu_lab->Write();

    h_MUDEC_pmupar->Write();
    h_geom_g4->Write();

    h_decay_flags->Write();
    h_decay_NUM->Write();
    h_decay_BR->Write();
    h_decay_BR_PDG->Write();

    h_prim_match_pi->Write();
    h_prim_match_ka->Write();
    h_prim_match_k0->Write();

    fROOT->Close();    
  }

  //print
  //
  G4cout
     << "\n--------------------End of Run------------------------------\n"
     //<< "\n mean Energy in Target : " << G4BestUnit(sumETarg,"Energy")
     //<< " +- "                          << G4BestUnit(rmsETarg,"Energy")  
     << "\n mean Energy in Target : " << sumETarg
     << " +- "                          << rmsETarg  
     //<< "\n mean Energy in Gap      : " << G4BestUnit(sumEGap,"Energy")
     //<< " +- "                          << G4BestUnit(rmsEGap,"Energy")
     << G4endl;
     
  G4cout
     //<< "\n mean trackLength in Target : " << G4BestUnit(sumLTarg,"Length")
     //<< " +- "                               << G4BestUnit(rmsLTarg,"Length")  
     << "\n mean trackLength in Target : " << sumLTarg
     << " +- "                               << rmsLTarg  
     //<< "\n mean trackLength in Gap      : " << G4BestUnit(sumLGap,"Length")
     //<< " +- "                               << G4BestUnit(rmsLGap,"Length")
     << "\n------------------------------------------------------------\n"
     << G4endl;
}
