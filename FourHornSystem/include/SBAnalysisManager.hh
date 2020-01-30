#ifndef SBAnalysisManager_h
#define SBAnalysisManager_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "globals.hh"
#include <vector>
#include <TObject.h>
#include <TObjString.h>
#include <TCollection.h>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TH2F.h"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include <TH1.h>
#include "TF1.h"
#include <TH2.h>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class SBAnalysisManager
{

  private:

	 SBAnalysisManager();

  static SBAnalysisManager* fManager;

  	 TFile *rootFile; // Root File

	 G4String ROOTDirectory;
         G4String ROOTFileName;

	public:

	G4int evtID;
  	static SBAnalysisManager* getInstance();

	~SBAnalysisManager();

  	void BeginOfRun();
  	void EndOfRun();

  	void BeginOfEvent();
  	void EndOfEvent();

	void InitDataStructure();
	void ClearDataStructure();

        void SetRootFileName(G4String aRootFileName)       { ROOTFileName  = aRootFileName  ; return ; } ;
        void SetRootDirectoryName(G4String aRootDirectory) { ROOTDirectory = aRootDirectory ; return ; } ;
        
        Int_t ENbins;
        Double_t EnuMIN;
        Double_t EnuMAX;

        Double_t rMAX;
        Double_t xMAX;
        Double_t yMAX;

        Double_t zMIN;
        Double_t zMAX;

        // Define Histogram
        TH1F *h_ETarg;
        TH1F *h_LTarg;
        //energy,impulsion (from PI K)
        TH1F *h_p_muplus_PI;
        TH1F *h_p_muminus_PI;
        TH1F *h_p_muplus_K;
        TH1F *h_p_muminus_K;
        
        //Neutrino...................unweight
        TH1F *pions_numu_du_pi_NW;
        TH1F *pions_anumu_du_pi_NW;
        TH1F *kaons_numu_du_ka_NW;
        TH1F *kaons_anumu_du_ka_NW;
       
        //Neutrino...................weight
        TH1F *pions_numu_du_pi_IDEAL;// from direct pi+ decay, IDEAL FOCUSING
        TH1F *pions_numu_du_pi;// from direct pi+ decay
        TH1F *pions_anumu_du_pi;// from direct pi- decay

        TH1F *pions_numu_du_mu;// from decay of mu coming from pi decay
        TH1F *pions_anumu_du_mu;// from decay of mu coming from pi decay
        TH1F *pions_nue_du_mu;// from decay of mu coming from pi decay
        TH1F *pions_anue_du_mu;// from decay of mu coming from pi decay
        //sum
        TH1F *h_numu;
        TH1F *h_anumu;
        TH1F *h_nue;
        TH1F *h_anue;

        // charged kaons chain

        TH1F *kaons_numu_du_ka;// from direct K+ decay
        TH1F *kaons_anumu_du_ka;// from direct K- decay
        
        TH1F *kaons_numu_du_pi;// from pi+/- coming from the K+-
        TH1F *kaons_anumu_du_pi;// from pi+/- coming from the K+-

        TH1F *kaons_nue_du_mu;// from decay of mu coming from K+/- decay
        TH1F *kaons_anumu_du_mu;// from decay of mu coming from K+/- decay
        TH1F *kaons_anue_du_mu;// from decay of mu coming from K+/- decay
        TH1F *kaons_numu_du_mu;// from decay of mu coming from K+/- decay

        // kaons chain

        TH1F *kzeros_numu_du_pi;// from pions coming from the kzeros
        TH1F *kzeros_anumu_du_pi;// from pions coming from the kzeros

        TH1F *kaons_numu_du_3;// from kplus 3 body decays
        TH1F *kaons_anumu_du_3;// from kplus 3 body decays
        TH1F *kaons_nue_du_3;// from kplus 3 body decays
        TH1F *kaons_anue_du_3;// from kplus 3 body decays

        TH1F *kzeros_numu_du_3;// from kzeros 3 body decays
        TH1F *kzeros_anumu_du_3;// from kzeros 3 body decays
        TH1F *kzeros_nue_du_3;// from kzeros 3 body decays
        TH1F *kzeros_anue_du_3;// from kzeros 3 body decays

        TH1F *kzeros_nue_du_mu;// from muons coming from the kzeros
        TH1F *kzeros_anumu_du_mu;// from muons coming from the kzeros
        TH1F *kzeros_anue_du_mu;// from muons coming from the kzeros
        TH1F *kzeros_numu_du_mu;// from muons coming from the kzeros

        //geometry decay................
        TH1F *h_piminus_DIF_r;            //distribution of pi- vs distance to axis r 
        TH1F *h_piminus_DIF_z;            //distribution of pi- vs distance to axis z 
        TH2F *h_piminus_DIF_xy;           //distribution of pi- longitudinal space xy
        TH2F *h_piminus_DIF_zr;           //distribution of pi- longitudinal space zr
        TH2F *h_piminus_DIF_thp;          //impulsion versus polar angle theta of pi-
        TH2F *h_piminus_DIF_thp_EnuW;     //impulsion versus polar angle theta of pi- with energy neutrino as a weight
        TH2F *h_piminus_DIF_thpinu;       //polar angle theta of neutrino versus polar angle theta of pi-
        TH2F *h_piminus_DIF_Ethnu;        //polar angle of neutrino versus energy of neutrino
        TH1F *h_piminus_DIF_EnuForw;      //distribution of pi- vs energy of neutrino

        TH1F *h_Kminus_DIF_r;
        TH1F *h_Kminus_DIF_z;
        TH2F *h_Kminus_DIF_xy;
        TH2F *h_Kminus_DIF_zr;
        TH2F *h_Kminus_DIF_thp;
        TH2F *h_Kminus_DIF_thp_EnuW;
        TH2F *h_Kminus_DIF_thpinu;
        TH2F *h_Kminus_DIF_Ethnu;
        TH1F *h_Kminus_DIF_EnuForw;

        TH1F *h_piplus_DIF_r;
        TH1F *h_piplus_DIF_z;
        TH2F *h_piplus_DIF_xy;
        TH2F *h_piplus_DIF_zr;
        TH2F *h_piplus_DIF_thp;
        TH2F *h_piplus_DIF_thp_EnuW;  
        TH2F *h_piplus_DIF_thpinu;
        TH2F *h_piplus_DIF_Ethnu;
        TH1F *h_piplus_DIF_EnuForw;

        TH1F *h_Kplus_DIF_r;
        TH1F *h_Kplus_DIF_z;
        TH2F *h_Kplus_DIF_xy;
        TH2F *h_Kplus_DIF_zr;
        TH2F *h_Kplus_DIF_thp;
        TH2F *h_Kplus_DIF_thp_EnuW;
        TH2F *h_Kplus_DIF_thpinu;
        TH2F *h_Kplus_DIF_Ethnu;
        TH1F *h_Kplus_DIF_EnuForw;

        void FillHisto(Double_t EnerTarg, Double_t LengthTarg) ;
        void FillHisto2(Double_t p_muplus_PI,
		        Double_t p_muminus_PI,
		        Double_t p_muplus_K,
		        Double_t p_muminus_K,
		        Double_t E_numu_PI,
		        Double_t E_anumu_PI,
		        Double_t E_numu_K,
		        Double_t E_anumu_K) ;

        void FillGeometryDecay(Double_t x,
			       Double_t y,
			       Double_t z,
			       Double_t px,
			       Double_t py,
			       Double_t pz,
			       Double_t pxnu,
			       Double_t pynu,
			       Double_t pznu,
			       Int_t moth,
			       Int_t q) ;

        void FillFinalOps(Double_t ProcEvtsScale, Double_t K_REP) ;

        void FillProba2(Int_t q,
                        G4String partyp,
                        Int_t parfl,
                        Double_t enuLAB,
                        Double_t weight,
                        Double_t weightIDEAL) ;

        void FillProba3K(Int_t decflag,
                         Double_t e_nu,
                         Double_t weight) ;

        void FillprobaMu(Int_t ipart,
                         Int_t parflag,
                         Double_t e_nu,
                         Double_t weight,
                         Double_t weightnue) ;

        void PrepareGLoBESFlux(G4String,G4double);

	void FillSumFluxesFlavours();

};

#endif
