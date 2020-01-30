#ifndef SBEventAction_h
#define SBEventAction_h 1

#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class SBRunAction;
class SBProba;
class SBEventActionMessenger;

class SBEventAction : public G4UserEventAction
{
public:
  SBEventAction(SBRunAction*,SBProba*);
  virtual ~SBEventAction();

  void  BeginOfEventAction(const G4Event*);
  void    EndOfEventAction(const G4Event*);
    
  void AddTarg(G4double de, G4double dl) {EnergyTarg += de; TrackLTarg += dl;};
  void AddGap(G4double de, G4double dl) {EnergyGap += de; TrackLGap += dl;};
                     
  void SetPrintModulo(G4int val){printModulo = val;};
  void SetEKinProt(G4double val){EKinProt = val;};
  void SetPower(G4double val){PowerMW = val;};
  void SetSBVerbosity(G4int val){SBVerbosity = val;};
  G4int GetSBVerbosity(){return SBVerbosity;};

  void DumpChannelsInfo();

  G4double EKinProt;
  G4double PowerMW;

  G4int PIPLUS_TRACKID;
  G4int PIMINUS_TRACKID;
  G4int KPLUS_TRACKID;
  G4int KMINUS_TRACKID;
  G4int KZEROL_TRACKID;
  G4int KZEROS_TRACKID;

  G4double p_muplus_PI;
  G4double p_muminus_PI;
  G4double p_muplus_K;
  G4double p_muminus_K;

  G4double E_numu_PI;
  G4double E_anumu_PI;
  G4double E_numu_K;
  G4double E_anumu_K;

  G4double EForw_numu_PI;
  G4double EForw_anumu_PI;
  G4double EForw_numu_K;
  G4double EForw_anumu_K;

  G4double W_numu_PI;
  G4double W_anumu_PI;

  G4double W_numu_K;
  G4double W_anumu_K;

  G4double NDIF_piplus;
  G4double NDIF_piminus;
  G4double NDIF_Kplus;
  G4double NDIF_Kminus;
  G4double NDIF_KzeroL;
  G4double NDIF_KzeroS;

  G4ThreeVector x_DIF_K_kplus;
  G4ThreeVector p_DIF_K_kplus;
  G4ThreeVector x_DIF_K_kminus;
  G4ThreeVector p_DIF_K_kminus;

  G4ThreeVector x_DIF_K_kzeroL;
  G4ThreeVector p_DIF_K_kzeroL;
  G4ThreeVector x_DIF_K_kzeroS;
  G4ThreeVector p_DIF_K_kzeroS;

  G4ThreeVector x_DIF_PI_piplus;
  G4ThreeVector p_DIF_PI_piplus;
  G4ThreeVector x_DIF_PI_piminus;
  G4ThreeVector p_DIF_PI_piminus;

  G4ThreeVector x_DIF_PI_muplus;
  G4ThreeVector p_DIF_PI_muplus;
  G4ThreeVector x_DIF_PI_muminus;
  G4ThreeVector p_DIF_PI_muminus;
  G4ThreeVector x_DIF_PI_numu;
  G4ThreeVector p_DIF_PI_numu;
  G4ThreeVector x_DIF_PI_anumu;
  G4ThreeVector p_DIF_PI_anumu;

  G4ThreeVector x_DIF_K_muplus;
  G4ThreeVector p_DIF_K_muplus;
  G4ThreeVector x_DIF_K_muminus;
  G4ThreeVector p_DIF_K_muminus;
  G4ThreeVector x_DIF_K_numu;
  G4ThreeVector p_DIF_K_numu;
  G4ThreeVector x_DIF_K_anumu;
  G4ThreeVector p_DIF_K_anumu;

  G4int NBody;
  G4bool pich2B;
  G4bool kch2B;
  G4bool kch3B;
  G4bool kzero3B;
  G4int decflag;

  G4int n1;
  G4int n2;
  G4int n3;
  G4int n4;
  G4int n5;
  G4int n6;
  G4int n7;
  G4int n8;
  G4int n9;
  G4int n10;
  G4int n11;
  G4int n12;
  G4int n13;
  G4int n14;
  G4int n15;
  G4int n16;
  G4int n17;
  G4int n18;
  G4int n19;
  G4int n20;
  G4int n30;
  G4int n31;

  G4ThreeVector p_nu1[20];
  G4ThreeVector p_nu2[20];
  G4ThreeVector p_nu3[20];
  G4ThreeVector p_nu4[20];
  G4ThreeVector p_nu5[20];
  G4ThreeVector p_nu6[20];
  G4ThreeVector p_nu7[20];
  G4ThreeVector p_nu8[20];
  G4ThreeVector p_nu9[20];
  G4ThreeVector p_nu10[20];
  G4ThreeVector p_nu11[20];
  G4ThreeVector p_nu12[20];
  G4ThreeVector p_nu13[20];
  G4ThreeVector p_nu14[20];
  G4ThreeVector p_nu15[20];
  G4ThreeVector p_nu16[20];
  G4ThreeVector p_nu17[20];
  G4ThreeVector p_nu18[20];
  G4ThreeVector p_nu19[20];
  G4ThreeVector p_nu20[20];  
  G4ThreeVector p_nu30[20];
  G4ThreeVector p_nu31[20];

  G4ThreeVector p_mu1[20];
  G4ThreeVector p_mu2[20];
  G4ThreeVector p_mu3[20];
  G4ThreeVector p_mu4[20];
  G4ThreeVector p_mu5[20];
  G4ThreeVector p_mu6[20];
  G4ThreeVector p_mu7[20];
  G4ThreeVector p_mu8[20];
  G4ThreeVector p_mu9[20];
  G4ThreeVector p_mu10[20];
  G4ThreeVector p_mu11[20];
  G4ThreeVector p_mu12[20];
  G4ThreeVector p_mu13[20];
  G4ThreeVector p_mu14[20];
  G4ThreeVector p_mu15[20];
  G4ThreeVector p_mu16[20];
  G4ThreeVector p_mu17[20];
  G4ThreeVector p_mu18[20];
  G4ThreeVector p_mu19[20];
  G4ThreeVector p_mu20[20];  
  G4ThreeVector p_mu30[20];
  G4ThreeVector p_mu31[20];

  G4ThreeVector x_mu1[20];
  G4ThreeVector x_mu2[20];
  G4ThreeVector x_mu3[20];
  G4ThreeVector x_mu4[20];
  G4ThreeVector x_mu5[20];
  G4ThreeVector x_mu6[20];
  G4ThreeVector x_mu7[20];
  G4ThreeVector x_mu8[20];
  G4ThreeVector x_mu9[20];
  G4ThreeVector x_mu10[20];
  G4ThreeVector x_mu11[20];
  G4ThreeVector x_mu12[20];
  G4ThreeVector x_mu13[20];
  G4ThreeVector x_mu14[20];
  G4ThreeVector x_mu15[20];
  G4ThreeVector x_mu16[20];
  G4ThreeVector x_mu17[20];
  G4ThreeVector x_mu18[20];
  G4ThreeVector x_mu19[20];
  G4ThreeVector x_mu20[20];  
  G4ThreeVector x_mu30[20];
  G4ThreeVector x_mu31[20];
  
  G4ThreeVector p_par1[20];
  G4ThreeVector p_par2[20];
  G4ThreeVector p_par3[20];
  G4ThreeVector p_par4[20];
  G4ThreeVector p_par5[20];
  G4ThreeVector p_par6[20];
  G4ThreeVector p_par7[20];
  G4ThreeVector p_par8[20];
  G4ThreeVector p_par9[20];
  G4ThreeVector p_par10[20];
  G4ThreeVector p_par11[20];
  G4ThreeVector p_par12[20];
  G4ThreeVector p_par13[20];
  G4ThreeVector p_par14[20];
  G4ThreeVector p_par15[20];
  G4ThreeVector p_par16[20];
  G4ThreeVector p_par17[20];
  G4ThreeVector p_par18[20];
  G4ThreeVector p_par19[20];
  G4ThreeVector p_par20[20];  
  G4ThreeVector p_par30[20];
  G4ThreeVector p_par31[20];

  G4int id_par1[20];
  G4int id_par2[20];
  G4int id_par3[20];
  G4int id_par4[20];
  G4int id_par5[20];
  G4int id_par6[20];
  G4int id_par7[20];
  G4int id_par8[20];
  G4int id_par9[20];
  G4int id_par10[20];
  G4int id_par11[20];
  G4int id_par12[20];
  G4int id_par13[20];
  G4int id_par14[20];
  G4int id_par15[20];
  G4int id_par16[20];
  G4int id_par17[20];
  G4int id_par18[20];
  G4int id_par19[20];
  G4int id_par20[20];  
  G4int id_par30[20];
  G4int id_par31[20];

  G4int id1[20];
  G4int id2[20];
  G4int id3[20];
  G4int id4[20];
  G4int id5[20];
  G4int id6[20];
  G4int id7[20];
  G4int id8[20];
  G4int id9[20];
  G4int id10[20];
  G4int id11[20];
  G4int id12[20];
  G4int id13[20];
  G4int id14[20];
  G4int id15[20];
  G4int id16[20];
  G4int id17[20];
  G4int id18[20];
  G4int id19[20];
  G4int id20[20];  
  G4int id30[20];
  G4int id31[20];

  G4int id_sec_pip1_1[20];
  G4int id_sec_pip1_2[20];
  G4int id_sec_pip1_3[20];
  G4int id_sec_pip1_4[20];
  G4int id_sec_pip1_5[20];
  G4int id_sec_pip1_6[20];
  G4int id_sec_pip1_7[20];
  G4int id_sec_pip1_8[20];
  G4int id_sec_pip1_9[20];
  G4int id_sec_pip1_10[20];
  G4int id_sec_pip1_11[20];
  G4int id_sec_pip1_12[20];
  G4int id_sec_pip1_13[20];
  G4int id_sec_pip1_14[20];
  G4int id_sec_pip1_15[20];
  G4int id_sec_pip1_16[20];
  G4int id_sec_pip1_17[20];
  G4int id_sec_pip1_18[20];
  G4int id_sec_pip1_19[20];
  G4int id_sec_pip1_20[20];  
  G4int id_sec_pip1_30[20];
  G4int id_sec_pip1_31[20];

  G4int id_sec_pim1_1[20];
  G4int id_sec_pim1_2[20];
  G4int id_sec_pim1_3[20];
  G4int id_sec_pim1_4[20];
  G4int id_sec_pim1_5[20];
  G4int id_sec_pim1_6[20];
  G4int id_sec_pim1_7[20];
  G4int id_sec_pim1_8[20];
  G4int id_sec_pim1_9[20];
  G4int id_sec_pim1_10[20];
  G4int id_sec_pim1_11[20];
  G4int id_sec_pim1_12[20];
  G4int id_sec_pim1_13[20];
  G4int id_sec_pim1_14[20];
  G4int id_sec_pim1_15[20];
  G4int id_sec_pim1_16[20];
  G4int id_sec_pim1_17[20];
  G4int id_sec_pim1_18[20];
  G4int id_sec_pim1_19[20];
  G4int id_sec_pim1_20[20];  
  G4int id_sec_pim1_30[20];
  G4int id_sec_pim1_31[20];

  G4int id_sec_pip2_1[20];
  G4int id_sec_pip2_2[20];
  G4int id_sec_pip2_3[20];
  G4int id_sec_pip2_4[20];
  G4int id_sec_pip2_5[20];
  G4int id_sec_pip2_6[20];
  G4int id_sec_pip2_7[20];
  G4int id_sec_pip2_8[20];
  G4int id_sec_pip2_9[20];
  G4int id_sec_pip2_10[20];
  G4int id_sec_pip2_11[20];
  G4int id_sec_pip2_12[20];
  G4int id_sec_pip2_13[20];
  G4int id_sec_pip2_14[20];
  G4int id_sec_pip2_15[20];
  G4int id_sec_pip2_16[20];
  G4int id_sec_pip2_17[20];
  G4int id_sec_pip2_18[20];
  G4int id_sec_pip2_19[20];
  G4int id_sec_pip2_20[20];  
  G4int id_sec_pip2_30[20];
  G4int id_sec_pip2_31[20];

  G4int id_sec_pim2_1[20];
  G4int id_sec_pim2_2[20];
  G4int id_sec_pim2_3[20];
  G4int id_sec_pim2_4[20];
  G4int id_sec_pim2_5[20];
  G4int id_sec_pim2_6[20];
  G4int id_sec_pim2_7[20];
  G4int id_sec_pim2_8[20];
  G4int id_sec_pim2_9[20];
  G4int id_sec_pim2_10[20];
  G4int id_sec_pim2_11[20];
  G4int id_sec_pim2_12[20];
  G4int id_sec_pim2_13[20];
  G4int id_sec_pim2_14[20];
  G4int id_sec_pim2_15[20];
  G4int id_sec_pim2_16[20];
  G4int id_sec_pim2_17[20];
  G4int id_sec_pim2_18[20];
  G4int id_sec_pim2_19[20];
  G4int id_sec_pim2_20[20];  
  G4int id_sec_pim2_30[20];
  G4int id_sec_pim2_31[20];


  G4int class_par30[20];
  G4int class_par31[20];

  G4ThreeVector p_K_3B;
  G4double enuLAB_K_3B;

  // primary, secondary, tertiary particle flag

  G4int ParLevel_piplus;
  G4int ParLevel_piminus;
  G4int ParLevel_muplus;
  G4int ParLevel_muminus;

  G4int MAXIND;
  G4int vPDG_ID[10000];
  G4int vPDG_PAR_ID[10000];
  G4int vPDG_PARPAR_ID[10000];

  G4String vName_ID[10000];
  G4String vName_PAR_ID[10000];
  G4String vName_PARPAR_ID[10000];

  // Entering tunnel
  // --
  G4ThreeVector x_EXI_piplus;
  G4ThreeVector p_EXI_piplus;
  G4ThreeVector x_EXI_piminus;
  G4ThreeVector p_EXI_piminus;
  G4ThreeVector x_EXI_kplus;
  G4ThreeVector p_EXI_kplus;
  G4ThreeVector x_EXI_kminus;
  G4ThreeVector p_EXI_kminus;
  G4ThreeVector x_EXI_muplus;
  G4ThreeVector p_EXI_muplus;
  G4ThreeVector x_EXI_muminus;
  G4ThreeVector p_EXI_muminus;
  G4ThreeVector x_EXI_k0L;
  G4ThreeVector p_EXI_k0L;

  G4double angnue0OLD;
  G4double angnumu0OLD;
  G4double anganue0OLD;
  G4double anganumu0OLD;

private:
  SBRunAction*  runAct;
  SBProba* Probability;

  G4double  EnergyTarg, EnergyGap;
  G4double  TrackLTarg, TrackLGap;

  
  G4int     printModulo;
  G4int     SBVerbosity;
  SBEventActionMessenger*  eventMessenger;

};

#endif

    
