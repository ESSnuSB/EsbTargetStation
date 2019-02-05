#include "SBAnalysisManager.hh"
#include "SBDetectorConstruction.hh"
#include "G4UnitsTable.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"

#include <TH1.h>
#include "TF1.h"
#include <TH2.h>

#include <iostream>
#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SBAnalysisManager* SBAnalysisManager::fManager = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SBAnalysisManager* SBAnalysisManager::getInstance()
{
 if(!fManager) { fManager = new SBAnalysisManager() ;  }
 return fManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SBAnalysisManager::SBAnalysisManager()
{
//============================================================  
  EnuMAX=1.5;// SPL standard
  ENbins=75;//SPL standard
  EnuMIN=0;

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

   zMIN=-30.;//tunnel make this depend on tunnel parameters!!
   zMAX=30.;//tunnel
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SBAnalysisManager::~SBAnalysisManager()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SBAnalysisManager::InitDataStructure()
{
  h_ETarg = new TH1F("h_ETarg","h_ETarg",100,0.,2000.) ;
  h_LTarg = new TH1F("h_LTarg","h_LTarg",100,0.,50.)   ;
  
  //muons
  h_p_muplus_PI = new TH1F("h_p_muplus_PI","h_p_muplus_PI",ENbins,EnuMIN/CLHEP::GeV,0.8);
  h_p_muminus_PI = new TH1F("h_p_muminus_PI","h_p_muminus_PI",ENbins,EnuMIN/CLHEP::GeV,0.8);
  h_p_muplus_K = new TH1F("h_p_muplus_K","h_p_muplus_K",ENbins,EnuMIN/CLHEP::GeV,0.8);
  h_p_muminus_K = new TH1F("h_p_muminus_K","h_p_muminus_K",ENbins,EnuMIN/CLHEP::GeV,0.8);

  // neutrino fluxes

  // unweighted
  pions_numu_du_pi_NW = new TH1F("pions_numu_du_pi_NW","pions_numu_du_pi_NW",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);
  pions_anumu_du_pi_NW = new TH1F("pions_anumu_du_pi_NW","pions_anumu_du_pi_NW",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);
  kaons_numu_du_ka_NW = new TH1F("kaons_numu_du_ka_NW","kaons_numu_du_ka_NW",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);
  kaons_anumu_du_ka_NW = new TH1F("kaons_anumu_du_ka_NW","kaons_anumu_du_ka_NW",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);
 // weighted
  kaons_numu_du_ka = new TH1F("kaons_numu_du_ka","kaons_numu_du_ka",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);
  kaons_anumu_du_ka = new TH1F("kaons_anumu_du_ka","kaons_anumu_du_ka",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);

  pions_numu_du_pi_IDEAL = new TH1F("pions_numu_du_pi_IDEAL","pions_numu_du_pi_IDEAL",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);

  pions_numu_du_pi = new TH1F("pions_numu_du_pi","pions_numu_du_pi",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);
  pions_anumu_du_pi = new TH1F("pions_anumu_du_pi","pions_anumu_du_pi",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);

  kaons_numu_du_pi = new TH1F("kaons_numu_du_pi","kaons_numu_du_pi",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);
  kaons_anumu_du_pi = new TH1F("kaons_anumu_du_pi","kaons_anumu_du_pi",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);

  kzeros_numu_du_pi = new TH1F("kzeros_numu_du_pi","kzeros_numu_du_pi",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);
  kzeros_anumu_du_pi = new TH1F("kzeros_anumu_du_pi","kzeros_anumu_du_pi",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);
  
  kaons_numu_du_3 = new TH1F("kaons_numu_du_3","kaons_numu_du_3",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);
  kaons_anumu_du_3 = new TH1F("kaons_anumu_du_3","kaons_anumu_du_3",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);
  kaons_nue_du_3 = new TH1F("kaons_nue_du_3","kaons_nue_du_3",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);
  kaons_anue_du_3 = new TH1F("kaons_anue_du_3","kaons_anue_du_3",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);

  kzeros_numu_du_3 = new TH1F("kzeros_numu_du_3","kzeros_numu_du_3",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);
  kzeros_anumu_du_3 = new TH1F("kzeros_anumu_du_3","kzeros_anumu_du_3",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);
  kzeros_nue_du_3 = new TH1F("kzeros_nue_du_3","kzeros_nue_du_3",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);
  kzeros_anue_du_3 = new TH1F("kzeros_anue_du_3","kzeros_anue_du_3",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);

  pions_numu_du_mu = new TH1F("pions_numu_du_mu","pions_numu_du_mu",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);
  pions_anumu_du_mu = new TH1F("pions_anumu_du_mu","pions_anumu_du_mu",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);
  pions_nue_du_mu = new TH1F("pions_nue_du_mu","pions_nue_du_mu",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);
  pions_anue_du_mu = new TH1F("pions_anue_du_mu","pions_anue_du_mu",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);

  kaons_nue_du_mu = new TH1F("kaons_nue_du_mu","kaons_nue_du_mu",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);
  kaons_anumu_du_mu = new TH1F("kaons_anumu_du_mu","kaons_anumu_du_mu",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);
  kaons_anue_du_mu = new TH1F("kaons_anue_du_mu","kaons_anue_du_mu",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);
  kaons_numu_du_mu = new TH1F("kaons_numu_du_mu","kaons_numu_du_mu",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);

  kzeros_nue_du_mu = new TH1F("kzeros_nue_du_mu","kzeros_nue_du_mu",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);
  kzeros_anumu_du_mu = new TH1F("kzeros_anumu_du_mu","kzeros_anumu_du_mu",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);
  kzeros_anue_du_mu = new TH1F("kzeros_anue_du_mu","kzeros_anue_du_mu",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);
  kzeros_numu_du_mu = new TH1F("kzeros_numu_du_mu","kzeros_numu_du_mu",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);
  // sum
  h_numu = new TH1F("h_numu","h_numu",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);
  h_anumu = new TH1F("h_anumu","h_anumu",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);
  h_nue = new TH1F("h_nue","h_nue",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);
  h_anue = new TH1F("h_anue","h_anue",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);

  // geometrical distribution of decay in flight of pions and kaons

  h_piminus_DIF_r = new TH1F("h_piminus_DIF_r","h_piminus_DIF_r",208,0.,rMAX);
  h_piminus_DIF_z = new TH1F("h_piminus_DIF_z","h_piminus_DIF_z",500,zMIN,zMAX);
  h_piminus_DIF_xy = new TH2F("h_piminus_DIF_xy","h_piminus_DIF_xy",100,-xMAX,xMAX,100,-yMAX,yMAX);
  h_piminus_DIF_xy->Draw("CONT0Z");
  h_piminus_DIF_zr = new TH2F("h_piminus_DIF_zr","h_piminus_DIF_zr",100,zMIN,zMAX,208,0.,rMAX);
  h_piminus_DIF_thp = new TH2F("h_piminus_DIF_thp","h_piminus_DIF_thp",100,0,1,100,0,1.5);
  h_piminus_DIF_thp_EnuW = new TH2F("h_piminus_DIF_thp_EnuW","h_piminus_DIF_thp_EnuW",100,0,1,100,0,1.5);
  h_piminus_DIF_thpinu = new TH2F("h_piminus_DIF_thpinu","h_piminus_DIF_thpinu",100,0,3.141592,100,0,3.141592);
  h_piminus_DIF_Ethnu = new TH2F("h_piminus_DIF_Ethnu","h_piminus_DIF_Ethnu",100,0,3.141592,100,0,1);
  h_piminus_DIF_EnuForw = new TH1F("h_piminus_DIF_EnuForw","h_piminus_DIF_EnuForw",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);

  h_Kminus_DIF_r =  new TH1F("h_Kminus_DIF_r","h_Kminus_DIF_r",208,0.,rMAX);
  h_Kminus_DIF_z =  new TH1F("h_Kminus_DIF_z","h_Kminus_DIF_z",500,zMIN,zMAX);
  h_Kminus_DIF_xy =  new TH2F("h_Kminus_DIF_xy","h_Kminus_DIF_xy",100,-xMAX,xMAX,100,-yMAX,yMAX);
  h_Kminus_DIF_zr =  new TH2F("h_Kminus_DIF_zr","h_Kminus_DIF_zr",100,zMIN,zMAX,208,0.,rMAX);
  h_Kminus_DIF_thp =  new TH2F("h_Kminus_DIF_thp","h_Kminus_DIF_thp",100,0,1,100,0,1.5);
  h_Kminus_DIF_thp_EnuW =  new TH2F("h_Kminus_DIF_thp_EnuW","h_Kminus_DIF_thp_EnuW",100,0,1,100,0,1.5);
  h_Kminus_DIF_thpinu =  new TH2F("h_Kminus_DIF_thpinu","h_Kminus_DIF_thpinu",100,0,3.141592,100,0,3.141592);
  h_Kminus_DIF_Ethnu =  new TH2F("h_Kminus_DIF_Ethnu","h_Kminus_DIF_Ethnu",100,0,3.141592,100,0,1);
  h_Kminus_DIF_EnuForw = new TH1F("h_Kminus_DIF_EnuForw","h_Kminus_DIF_EnuForw",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);

  h_piplus_DIF_r =  new TH1F("h_piplus_DIF_r","h_piplus_DIF_r",208,0.,rMAX);
  h_piplus_DIF_z =  new TH1F("h_piplus_DIF_z","h_piplus_DIF_z",500,zMIN,zMAX);
  h_piplus_DIF_xy =  new TH2F("h_piplus_DIF_xy","h_piplus_DIF_xy",100,-xMAX,xMAX,100,-yMAX,yMAX);
  h_piplus_DIF_zr =  new TH2F("h_piplus_DIF_zr","h_piplus_DIF_zr",100,zMIN,zMAX,208,0.,rMAX);
  h_piplus_DIF_thp =  new TH2F("h_piplus_DIF_thp","h_piplus_DIF_thp",100,0,1,100,0,1.5);
  h_piplus_DIF_thp_EnuW =  new TH2F("h_piplus_DIF_thp_EnuW","h_piplus_DIF_thp_EnuW",100,0,1,100,0,1.5);
  h_piplus_DIF_thpinu =  new TH2F("h_piplus_DIF_thpinu","h_piplus_DIF_thpinu",100,0,3.141592,100,0,3.141592);
  h_piplus_DIF_Ethnu =  new TH2F("h_piplus_DIF_Ethnu","h_piplus_DIF_Ethnu",100,0,3.141592,100,0,1);
  h_piplus_DIF_EnuForw = new TH1F("h_piplus_DIF_EnuForw","h_piplus_DIF_EnuForw",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);

  h_Kplus_DIF_r =  new TH1F("h_Kplus_DIF_r","h_Kplus_DIF_r",208,0.,rMAX);
  h_Kplus_DIF_z =  new TH1F("h_Kplus_DIF_z","h_Kplus_DIF_z",500,zMIN,zMAX);
  h_Kplus_DIF_xy =  new TH2F("h_Kplus_DIF_xy","h_Kplus_DIF_xy",100,-xMAX,xMAX,100,-yMAX,yMAX);
  h_Kplus_DIF_zr =  new TH2F("h_Kplus_DIF_zr","h_Kplus_DIF_zr",100,zMIN,zMAX,208,0.,rMAX);
  h_Kplus_DIF_thp =  new TH2F("h_Kplus_DIF_thp","h_Kplus_DIF_thp",100,0,1,100,0,1.5);
  h_Kplus_DIF_thp_EnuW =  new TH2F("h_Kplus_DIF_thp_EnuW","h_Kplus_DIF_thp_EnuW",100,0,1,100,0,1.5);
  h_Kplus_DIF_thpinu =  new TH2F("h_Kplus_DIF_thpinu","h_Kplus_DIF_thpinu",100,0,3.141592,100,0,3.141592);
  h_Kplus_DIF_Ethnu =  new TH2F("h_Kplus_DIF_Ethnu","h_Kplus_DIF_Ethnu",100,0,3.141592,100,0,1);
  h_Kplus_DIF_EnuForw = new TH1F("h_Kplus_DIF_EnuForw","h_Kplus_DIF_EnuForw",ENbins,EnuMIN/CLHEP::GeV,EnuMAX);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SBAnalysisManager::FillHisto(Double_t EnerTarg, Double_t LengthTarg)
{
  h_ETarg->Fill(EnerTarg);
  h_LTarg->Fill(LengthTarg);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SBAnalysisManager::FillHisto2(Double_t p_muplus_PI,
			           Double_t p_muminus_PI,
			           Double_t p_muplus_K,
			           Double_t p_muminus_K,
			           Double_t E_numu_PI,
			           Double_t E_anumu_PI,
			           Double_t E_numu_K,
			           Double_t E_anumu_K)
{
  h_p_muplus_PI->Fill(p_muplus_PI/CLHEP::GeV);
  h_p_muminus_PI->Fill(p_muminus_PI/CLHEP::GeV);
  h_p_muplus_K->Fill(p_muplus_K/CLHEP::GeV);
  h_p_muminus_K->Fill(p_muminus_K/CLHEP::GeV);

  pions_numu_du_pi_NW->Fill(E_numu_PI/CLHEP::GeV);
  pions_anumu_du_pi_NW->Fill(E_anumu_PI/CLHEP::GeV);
  kaons_numu_du_ka_NW->Fill(E_numu_K/CLHEP::GeV);
  kaons_anumu_du_ka_NW->Fill(E_anumu_K/CLHEP::GeV);

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SBAnalysisManager::FillGeometryDecay(Double_t x,
			                  Double_t y,
			                  Double_t z,
			                  Double_t px,
			                  Double_t py,
			                  Double_t pz,
			                  Double_t pxnu,
			                  Double_t pynu,
			                  Double_t pznu,
			                  Int_t moth,
			                  Int_t q)
{
  Double_t r=sqrt(x*x+y*y);
  Double_t p=sqrt((px*px)+(py*py)+(pz*pz));
  Double_t enu=sqrt((pxnu*pxnu)+(pynu*pynu)+(pznu*pznu));
  Double_t th=0,thnu=0;
  //Double_t THFOR = 0.01;

  if(p)th=acos(pz/p);

  Double_t OffAxisAngle=3.141592*2.5/180.;
  Double_t OffAxisPhi=0;

  Double_t eOME[3]={0.,0.,0.};
  eOME[0]=sin(OffAxisAngle)*cos(OffAxisPhi);
  eOME[1]=sin(OffAxisAngle)*sin(OffAxisPhi);
  eOME[2]=cos(OffAxisAngle);

  if(enu)thnu=acos(((eOME[0]*pxnu)+(eOME[1]*pynu)+(eOME[2]*pznu))/enu);   //angle theta
  
  if(q==1){//plus

    if(moth==1){//pi
      if(r)h_piplus_DIF_r->Fill(r/CLHEP::cm);
      if(z)h_piplus_DIF_z->Fill(z/CLHEP::m);
      if(x&&y)h_piplus_DIF_xy->Fill(x/CLHEP::m,y/CLHEP::m);
      if(z&&r)h_piplus_DIF_zr->Fill(z/CLHEP::m,r/CLHEP::cm);                
      if(p&&th)h_piplus_DIF_thp->Fill(th,p/CLHEP::GeV);           
      if(p&&th)h_piplus_DIF_thp_EnuW->Fill(th,p/CLHEP::GeV,enu/CLHEP::GeV);
      if(th&&thnu)h_piplus_DIF_thpinu->Fill(th,thnu);
      if(enu&&thnu)h_piplus_DIF_Ethnu->Fill(thnu,enu/CLHEP::GeV);
      if(thnu)h_piplus_DIF_EnuForw->Fill(enu/CLHEP::GeV);   //if(thnu<THFOR)
    } else if (moth==2){//K
      if(r)h_Kplus_DIF_r->Fill(r/CLHEP::cm);
      if(z)h_Kplus_DIF_z->Fill(z/CLHEP::m);
      if(x&&y)h_Kplus_DIF_xy->Fill(x/CLHEP::m,y/CLHEP::m);
      if(z&&r)h_Kplus_DIF_zr->Fill(z/CLHEP::m,r/CLHEP::cm);
      if(p&&th)h_Kplus_DIF_thp->Fill(th,p/CLHEP::GeV);
      if(p&&th)h_Kplus_DIF_thp_EnuW->Fill(th,p/CLHEP::GeV,enu/CLHEP::GeV);
      if(th&&thnu)h_Kplus_DIF_thpinu->Fill(th,thnu);
      if(enu&&thnu)h_Kplus_DIF_Ethnu->Fill(thnu,enu/CLHEP::GeV);
      if(thnu)h_Kplus_DIF_EnuForw->Fill(enu/CLHEP::GeV);   //if(thnu<THFOR)
    }
    
  } else if (q==-1){
    
    if(moth==1){//pi
      if(r)h_piminus_DIF_r->Fill(r/CLHEP::cm);
      if(z)h_piminus_DIF_z->Fill(z/CLHEP::m);
      if(x&&y)h_piminus_DIF_xy->Fill(x/CLHEP::m,y/CLHEP::m);
      if(z&&r)h_piminus_DIF_zr->Fill(z/CLHEP::m,r/CLHEP::cm);
      if(p&&th)h_piminus_DIF_thp->Fill(th,p/CLHEP::GeV);
      if(p&&th)h_piminus_DIF_thp_EnuW->Fill(th,p/CLHEP::GeV,enu/CLHEP::GeV);
      if(th&&thnu)h_piminus_DIF_thpinu->Fill(th,thnu);
      if(enu&&thnu)h_piminus_DIF_Ethnu->Fill(thnu,enu/CLHEP::GeV);
      if(thnu)h_piminus_DIF_EnuForw->Fill(enu/CLHEP::GeV);  //if(thnu<THFOR)
    } else if (moth==2){//K
      if(r)h_Kminus_DIF_r->Fill(r/CLHEP::cm);
      if(z)h_Kminus_DIF_z->Fill(z/CLHEP::m);
      if(x&&y)h_Kminus_DIF_xy->Fill(x/CLHEP::m,y/CLHEP::m);
      if(z&&r)h_Kminus_DIF_zr->Fill(z/CLHEP::m,r/CLHEP::cm);
      if(p&&th)h_Kminus_DIF_thp->Fill(th,p/CLHEP::GeV);
      if(p&&th)h_Kminus_DIF_thp_EnuW->Fill(th,p/CLHEP::GeV,enu/CLHEP::GeV);
      if(th&&thnu)h_Kminus_DIF_thpinu->Fill(th,thnu);
      if(enu&&thnu)h_Kminus_DIF_Ethnu->Fill(thnu,enu/CLHEP::GeV);
      if(thnu)h_Kminus_DIF_EnuForw->Fill(enu/CLHEP::GeV);  //if(thnu<THFOR)
    }
    
  }

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SBAnalysisManager::FillProba2(Int_t q,
                                   G4String partyp,
                                   Int_t parfl,
                                   Double_t enuLAB,
                                   Double_t weight,
                                   Double_t weightIDEAL)
{

    if(partyp=="kaon"){
      if(parfl==0){
	if(q==(+1))kaons_numu_du_ka->Fill(enuLAB/CLHEP::GeV,weight);
	if(q==(-1))kaons_anumu_du_ka->Fill(enuLAB/CLHEP::GeV,weight);
      }
    }

    if(partyp=="pion"){
      if(parfl==0){
	if(q==(+1))pions_numu_du_pi->Fill(enuLAB/CLHEP::GeV,weight);
	if(q==(+1))pions_numu_du_pi_IDEAL->Fill(enuLAB/CLHEP::GeV,weightIDEAL);
	if(q==(-1))pions_anumu_du_pi->Fill(enuLAB/CLHEP::GeV,weight);	
      }else if(parfl==1){
	if(q==(+1))kaons_numu_du_pi->Fill(enuLAB/CLHEP::GeV,weight);
	if(q==(-1))kaons_anumu_du_pi->Fill(enuLAB/CLHEP::GeV,weight);
      }else if(parfl==2){
	if(q==(+1))kzeros_numu_du_pi->Fill(enuLAB/CLHEP::GeV,weight);
	if(q==(-1))kzeros_anumu_du_pi->Fill(enuLAB/CLHEP::GeV,weight);
      }    
    }

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SBAnalysisManager::FillProba3K(Int_t decflag,
                                    Double_t e_nu,
                                    Double_t weight)
{

    if(decflag==7)kaons_nue_du_3->Fill(e_nu/CLHEP::GeV,weight);
    if(decflag==8)kaons_anue_du_3->Fill(e_nu/CLHEP::GeV,weight);
    if(decflag==9)kaons_numu_du_3->Fill(e_nu/CLHEP::GeV,weight);
    if(decflag==10)kaons_anumu_du_3->Fill(e_nu/CLHEP::GeV,weight);
    
    if(decflag==13)kzeros_nue_du_3->Fill(e_nu/CLHEP::GeV,weight);
    if(decflag==14)kzeros_anue_du_3->Fill(e_nu/CLHEP::GeV,weight);
    if(decflag==15)kzeros_numu_du_3->Fill(e_nu/CLHEP::GeV,weight);
    if(decflag==16)kzeros_anumu_du_3->Fill(e_nu/CLHEP::GeV,weight);


}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void SBAnalysisManager::FillprobaMu(Int_t ipart,
                                    Int_t parflag,
                                    Double_t e_nu,
                                    Double_t weight,
                                    Double_t weightnue)
{
      if(ipart==1){//mu+
	if(parflag==1){
	  pions_anumu_du_mu->Fill(e_nu/CLHEP::GeV,weight);
	  pions_nue_du_mu->Fill(e_nu/CLHEP::GeV,weightnue);
	} else if(parflag==2) {
	  kaons_anumu_du_mu->Fill(e_nu/CLHEP::GeV,weight);
	  kaons_nue_du_mu->Fill(e_nu/CLHEP::GeV,weightnue);
	} else if(parflag==3) {
	  kzeros_anumu_du_mu->Fill(e_nu/CLHEP::GeV,weight);
	  kzeros_nue_du_mu->Fill(e_nu/CLHEP::GeV,weightnue);
	}
	
      } else if (ipart==-1){//mu-
	if(parflag==1){
	  pions_numu_du_mu->Fill(e_nu/CLHEP::GeV,weight);
	  pions_anue_du_mu->Fill(e_nu/CLHEP::GeV,weightnue);
	} else if(parflag==2) {
	  kaons_numu_du_mu->Fill(e_nu/CLHEP::GeV,weight);
	  kaons_anue_du_mu->Fill(e_nu/CLHEP::GeV,weightnue);
	} else if(parflag==3) {
	  kzeros_numu_du_mu->Fill(e_nu/CLHEP::GeV,weight);
	  kzeros_anue_du_mu->Fill(e_nu/CLHEP::GeV,weightnue);
	}
      }

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SBAnalysisManager::FillSumFluxesFlavours()
{

  h_numu->Add(pions_numu_du_pi);
  h_numu->Add(pions_numu_du_mu);

  h_numu->Add(kaons_numu_du_ka);
  h_numu->Add(kaons_numu_du_pi);
  h_numu->Add(kaons_numu_du_mu);
  h_numu->Add(kaons_numu_du_3);

  h_numu->Add(kzeros_numu_du_pi);
  h_numu->Add(kzeros_numu_du_mu);
  h_numu->Add(kzeros_numu_du_3);

  //----//

  h_anumu->Add(pions_anumu_du_pi);
  h_anumu->Add(pions_anumu_du_mu);

  h_anumu->Add(kaons_anumu_du_ka);
  h_anumu->Add(kaons_anumu_du_pi);
  h_anumu->Add(kaons_anumu_du_mu);
  h_anumu->Add(kaons_anumu_du_3);

  h_anumu->Add(kzeros_anumu_du_pi);
  h_anumu->Add(kzeros_anumu_du_mu);
  h_anumu->Add(kzeros_anumu_du_3);

  //----//

  h_nue->Add(pions_nue_du_mu);
  h_nue->Add(kaons_nue_du_mu);
  h_nue->Add(kaons_nue_du_3);
  h_nue->Add(kzeros_nue_du_mu);
  h_nue->Add(kzeros_nue_du_3);

  //----//

  h_anue->Add(pions_anue_du_mu);
  h_anue->Add(kaons_anue_du_mu);
  h_anue->Add(kaons_anue_du_3);
  h_anue->Add(kzeros_anue_du_mu);
  h_anue->Add(kzeros_anue_du_3);

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SBAnalysisManager::PrepareGLoBESFlux(G4String fname="nufluxes_GLOBESformat.txt",G4double L=130.)
{
  double scal = pow(100./L,2);
  //double bwscale = 0.02/h_nue->GetBinWidth(1);// norm factor in GLOBES valid for 20MeV bins !        /////////////
  //scal*=bwscale;
  
  std::ofstream fout(fname);
  double nutau=0,antinutau=0;
  for(int i=1;i<h_numu->GetNbinsX();i++){
    fout 
      << h_numu->GetBinCenter(i)+0.5*h_nue->GetBinWidth(i)<< " " 
      << scal*h_nue->GetBinContent(i) << " " 
      << scal*h_numu->GetBinContent(i) << " " 
      << scal*nutau << " "
      << scal*h_anue->GetBinContent(i) << " " 
      << scal*h_anumu->GetBinContent(i) << " " 
      << scal*antinutau
      << G4endl;
    G4cout 
      << h_numu->GetBinCenter(i)+0.5*h_nue->GetBinWidth(i)<< " " 
      << scal*h_nue->GetBinContent(i) << " " 
      << scal*h_numu->GetBinContent(i) << " " 
      << scal*nutau << " "
      << scal*h_anue->GetBinContent(i) << " " 
      << scal*h_anumu->GetBinContent(i) << " " 
      << scal*antinutau
      << G4endl;
  }

  int nleft = 502-ENbins;
  for(int i=0;i<nleft;i++){
    fout 
      << EnuMAX/CLHEP::GeV+h_nue->GetBinWidth(i)*i << " " 
      << 0. << " " 
      << 0. << " " 
      << 0. << " "
      << 0. << " " 
      << 0. << " " 
      << 0.
      << G4endl;
  }
  fout.close();
  //system(Form("cat completamento.dat >> %s",fname));
  G4cout << fname << " created. Scaled by " << scal <<" (100/130)^2. "<< L <<" Km "<< G4endl; 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SBAnalysisManager::FillFinalOps(Double_t ProcEvtsScale,Double_t K_REP)
{
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
  
  //neutrino.........weight
  kaons_numu_du_ka->Scale(ProcEvtsScale/K_REP);
  kaons_anumu_du_ka->Scale(ProcEvtsScale/K_REP);

  pions_numu_du_pi_IDEAL->Scale(ProcEvtsScale);

  pions_numu_du_pi->Scale(ProcEvtsScale);
  pions_anumu_du_pi->Scale(ProcEvtsScale);

  kaons_numu_du_pi->Scale(ProcEvtsScale/K_REP);
  kaons_anumu_du_pi->Scale(ProcEvtsScale/K_REP);

  kzeros_numu_du_pi->Scale(ProcEvtsScale/K_REP);
  kzeros_anumu_du_pi->Scale(ProcEvtsScale/K_REP);

  kaons_numu_du_3->Scale(ProcEvtsScale/K_REP);
  kaons_anumu_du_3->Scale(ProcEvtsScale/K_REP);
  kaons_nue_du_3->Scale(ProcEvtsScale/K_REP);
  kaons_anue_du_3->Scale(ProcEvtsScale/K_REP);

  kzeros_numu_du_3->Scale(ProcEvtsScale/K_REP);
  kzeros_anumu_du_3->Scale(ProcEvtsScale/K_REP);
  kzeros_nue_du_3->Scale(ProcEvtsScale/K_REP);
  kzeros_anue_du_3->Scale(ProcEvtsScale/K_REP);

  pions_numu_du_mu->Scale(ProcEvtsScale);
  pions_anumu_du_mu->Scale(ProcEvtsScale);
  pions_nue_du_mu->Scale(ProcEvtsScale);
  pions_anue_du_mu->Scale(ProcEvtsScale);

  kaons_nue_du_mu->Scale(ProcEvtsScale/K_REP);
  kaons_anumu_du_mu->Scale(ProcEvtsScale/K_REP);
  kaons_anue_du_mu->Scale(ProcEvtsScale/K_REP);
  kaons_numu_du_mu->Scale(ProcEvtsScale/K_REP);

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SBAnalysisManager::BeginOfRun()
{
 SBAnalysisManager();
 InitDataStructure();
 
 /* 
 ROOTFileName = "SB.root";
 //ROOTDirectory = "/scratch9/fkoll/";
 //ROOTFileName = ROOTDirectory+ROOTFileName;
 G4cout << "Opening the output file : " << ROOTFileName << G4endl ;
 
 rootFile = new TFile(ROOTFileName,"RECREATE");
 if (!rootFile) { rootFile = new TFile("SB.root") ; } ;
*/

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SBAnalysisManager::EndOfRun()
{

  h_ETarg->Write();
  h_LTarg->Write();
  //
  h_p_muplus_PI->Write();
  h_p_muminus_PI->Write();
  h_p_muplus_K->Write();
  h_p_muminus_K->Write();
  //neutrino........unweight
  pions_numu_du_pi_NW->Write();
  pions_anumu_du_pi_NW->Write();
  kaons_numu_du_ka_NW->Write();
  kaons_anumu_du_ka_NW->Write();
  //neutrino........weight
  kaons_numu_du_ka->Write();
  kaons_anumu_du_ka->Write();

  pions_numu_du_pi_IDEAL->Write();

  pions_numu_du_pi->Write();
  pions_anumu_du_pi->Write();

  kaons_numu_du_pi->Write();
  kaons_anumu_du_pi->Write();

  kzeros_numu_du_pi->Write();
  kzeros_anumu_du_pi->Write();

  kaons_numu_du_3->Write();
  kaons_anumu_du_3->Write();
  kaons_nue_du_3->Write();
  kaons_anue_du_3->Write();

  kzeros_numu_du_3->Write();
  kzeros_anumu_du_3->Write();
  kzeros_nue_du_3->Write();
  kzeros_anue_du_3->Write(); 

  pions_numu_du_mu->Write();
  pions_anumu_du_mu->Write();
  pions_nue_du_mu->Write();
  pions_anue_du_mu->Write();

  kaons_nue_du_mu->Write();
  kaons_anumu_du_mu->Write();
  kaons_anue_du_mu->Write();
  kaons_numu_du_mu->Write();

  kzeros_nue_du_mu->Write();
  kzeros_anumu_du_mu->Write();
  kzeros_anue_du_mu->Write();
  kzeros_numu_du_mu->Write();
  
  h_numu->Write();
  h_anumu->Write();
  h_nue->Write();
  h_anue->Write();

  //
  h_piminus_DIF_r->Write();
  h_piminus_DIF_z->Write();
  h_piminus_DIF_xy->Write();
  h_piminus_DIF_zr->Write();
  h_piminus_DIF_thp->Write();
  h_piminus_DIF_thp_EnuW->Write();
  h_piminus_DIF_thpinu->Write();
  h_piminus_DIF_Ethnu->Write();
  h_piminus_DIF_EnuForw->Write();

  h_Kminus_DIF_r->Write();
  h_Kminus_DIF_z->Write();
  h_Kminus_DIF_xy->Write();
  h_Kminus_DIF_zr->Write();
  h_Kminus_DIF_thp->Write();
  h_Kminus_DIF_thp_EnuW->Write();
  h_Kminus_DIF_thpinu->Write();
  h_Kminus_DIF_Ethnu->Write();
  h_Kminus_DIF_EnuForw->Write();

  h_piplus_DIF_r->Write();
  h_piplus_DIF_z->Write();
  h_piplus_DIF_xy->Write();
  h_piplus_DIF_zr->Write();
  h_piplus_DIF_thp->Write();
  h_piplus_DIF_thp_EnuW->Write();
  h_piplus_DIF_thpinu->Write();
  h_piplus_DIF_Ethnu->Write();
  h_piplus_DIF_EnuForw->Write();

  h_Kplus_DIF_r->Write();
  h_Kplus_DIF_z->Write();
  h_Kplus_DIF_xy->Write();
  h_Kplus_DIF_zr->Write();
  h_Kplus_DIF_thp->Write();
  h_Kplus_DIF_thp_EnuW->Write();
  h_Kplus_DIF_thpinu->Write();
  h_Kplus_DIF_Ethnu->Write();
  h_Kplus_DIF_EnuForw->Write();
  //







/*
  h_ETarg->Delete();
  h_LTarg->Delete();
  
  rootFile->Close();
  rootFile->Delete();

  G4cout << "Closing root file" << G4endl;
*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SBAnalysisManager::BeginOfEvent()
{
  //G4cout << " Event ID = " << evtID << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SBAnalysisManager::EndOfEvent()
{
}
