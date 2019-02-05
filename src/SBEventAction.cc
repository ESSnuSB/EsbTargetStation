#include "SBEventAction.hh"
#include "SBRunAction.hh"
#include "SBProba.hh"
#include "SBEventActionMessenger.hh"

#include "G4Event.hh"
#include "G4TrajectoryContainer.hh"
#include "G4VTrajectory.hh"
#include "G4VVisManager.hh"
#include "G4UnitsTable.hh"

#include "G4KaonPlus.hh"
#include "G4PionPlus.hh"

#include "Randomize.hh"
#include <iomanip>

SBEventAction::SBEventAction(SBRunAction* run, SBProba* proba) :runAct(run),Probability(proba),printModulo(1),eventMessenger(0)
{
  eventMessenger = new SBEventActionMessenger(this);

  // Proton Beam Parameters 
  // -- 
  // Parameters can be overridden with macro file
  // --
  PowerMW = 4.;//MW
  EKinProt = 4.5*CLHEP::GeV;

  G4cout<<"CALL SBEventAction::SBEventAction "<<G4endl;
}

SBEventAction::~SBEventAction()
{
  G4cout<<"CALL SBEventAction::~SBEventAction "<<G4endl;
  delete eventMessenger;
}

void SBEventAction::BeginOfEventAction(const G4Event* evt)
{  
  G4cout<<"CALL SBEventAction::BeginOfEventAction "<<G4endl;
  G4int evtNb = evt->GetEventID();

  // Check if input file end has been reachedin that case end run.
  // --
  runAct->CheckEOF();

  if (evtNb%printModulo == 0) { 
    G4cout << "\n---> Begin of event: " << evtNb << G4endl;
    CLHEP::HepRandom::showEngineStatus();
  }
 
  // Initialisation per event
  // --
  EnergyTarg = EnergyGap = 0.;
  TrackLTarg = TrackLGap = 0.;

  PIPLUS_TRACKID 	= 0;
  PIMINUS_TRACKID 	= 0;
  KPLUS_TRACKID 	= 0;
  KMINUS_TRACKID 	= 0;
  KZEROS_TRACKID 	= 0;
  KZEROL_TRACKID 	= 0;

  p_muplus_PI		= -99999.;
  p_muminus_PI		= -99999.;
  p_muplus_K		= -99999.;
  p_muminus_K		= -99999.;

  E_numu_PI		= -99999.;
  E_anumu_PI		= -99999.;
  E_numu_K		= -99999.;
  E_anumu_K		= -99999.;

  EForw_numu_PI		= -99999.;
  EForw_anumu_PI	= -99999.;
  EForw_numu_K		= -99999.;
  EForw_anumu_K		= -99999.;

  W_numu_PI		= -99999.;
  W_anumu_PI		= -99999.;

  W_numu_K		= -99999.;
  W_anumu_K		= -99999.;

  NDIF_piplus		= 0;
  NDIF_piminus		= 0;
  NDIF_Kplus		= 0;
  NDIF_Kminus		= 0;
  NDIF_KzeroS		= 0;
  NDIF_KzeroL		= 0;

  p_DIF_PI_piplus	= G4ThreeVector(0,0,0);
  p_DIF_PI_piminus	= G4ThreeVector(0,0,0);
  p_DIF_K_kplus		= G4ThreeVector(0,0,0);
  p_DIF_K_kminus	= G4ThreeVector(0,0,0);

  p_DIF_PI_numu		= G4ThreeVector(0,0,0);
  p_DIF_PI_anumu	= G4ThreeVector(0,0,0);
  p_DIF_K_numu		= G4ThreeVector(0,0,0);
  p_DIF_K_anumu		= G4ThreeVector(0,0,0);

  p_DIF_PI_muplus	= G4ThreeVector(0,0,0);
  p_DIF_PI_muminus	= G4ThreeVector(0,0,0);
  p_DIF_K_muplus	= G4ThreeVector(0,0,0);
  p_DIF_K_muminus	= G4ThreeVector(0,0,0);

  p_DIF_K_kzeroL	= G4ThreeVector(0,0,0);
  p_DIF_K_kzeroS	= G4ThreeVector(0,0,0);

  x_DIF_PI_piplus	= G4ThreeVector(0,0,0);
  x_DIF_PI_piminus	= G4ThreeVector(0,0,0);
  x_DIF_K_kplus		= G4ThreeVector(0,0,0);
  x_DIF_K_kminus	= G4ThreeVector(0,0,0);

  x_DIF_K_kzeroL	= G4ThreeVector(0,0,0);
  x_DIF_K_kzeroS	= G4ThreeVector(0,0,0);

  x_DIF_PI_numu		= G4ThreeVector(0,0,0);
  x_DIF_PI_anumu	= G4ThreeVector(0,0,0);
  x_DIF_K_numu		= G4ThreeVector(0,0,0);
  x_DIF_K_anumu		= G4ThreeVector(0,0,0);

  x_DIF_PI_muplus	= G4ThreeVector(0,0,0);
  x_DIF_PI_muminus	= G4ThreeVector(0,0,0);
  x_DIF_K_muplus	= G4ThreeVector(0,0,0);
  x_DIF_K_muminus	= G4ThreeVector(0,0,0);
  
  ParLevel_piplus	= -1;
  ParLevel_piminus	= -1;
  ParLevel_muplus	= -1;
  ParLevel_muminus	= -1;

  MAXIND=10000;// must be equal to the one in SBSteppingAction.hh ...
  for(int i=0;i<MAXIND;i++){
    vPDG_ID[i]		= 0;
    vPDG_PAR_ID[i]	= 0;
    vPDG_PARPAR_ID[i]	= 0;
  }

  NBody		= 0;
  pich2B	= false;
  kch2B		= false;
  kch3B		= false;
  kzero3B	= false;
  decflag	= 0;
  p_K_3B	= G4ThreeVector(0,0,0);
  enuLAB_K_3B	= 0;

  n1		= 0;
  n2		= 0;
  n3		= 0;
  n4		= 0;
  n5		= 0;
  n6		= 0;
  n7		= 0;
  n8		= 0;
  n9		= 0;
  n10		= 0;
  n11		= 0;
  n12		= 0;
  n13		= 0;
  n14		= 0;
  n15		= 0;
  n16		= 0;
  n17		= 0;
  n18		= 0;
  n19		= 0;
  n20		= 0;
  n30		= 0;
  n31		= 0;
  for(int i=0;i<20;i++){
    p_nu1[i]	= G4ThreeVector(0,0,0);
    p_nu2[i]	= G4ThreeVector(0,0,0);
    p_nu3[i]	= G4ThreeVector(0,0,0);
    p_nu4[i]	= G4ThreeVector(0,0,0);
    p_nu5[i]	= G4ThreeVector(0,0,0);
    p_nu6[i]	= G4ThreeVector(0,0,0);
    p_nu7[i]	= G4ThreeVector(0,0,0);
    p_nu8[i]	= G4ThreeVector(0,0,0);
    p_nu9[i]	= G4ThreeVector(0,0,0);
    p_nu10[i]	= G4ThreeVector(0,0,0);
    p_nu11[i]	= G4ThreeVector(0,0,0);
    p_nu12[i]	= G4ThreeVector(0,0,0);
    p_nu13[i]   = G4ThreeVector(0,0,0);
    p_nu14[i]	= G4ThreeVector(0,0,0);
    p_nu15[i]	= G4ThreeVector(0,0,0);
    p_nu16[i]	= G4ThreeVector(0,0,0);
    p_nu17[i]	= G4ThreeVector(0,0,0);
    p_nu18[i]	= G4ThreeVector(0,0,0);
    p_nu19[i]	= G4ThreeVector(0,0,0);
    p_nu20[i]	= G4ThreeVector(0,0,0);
    p_nu30[i]	= G4ThreeVector(0,0,0);
    p_nu31[i]	= G4ThreeVector(0,0,0);

    p_mu1[i]	= G4ThreeVector(0,0,0);
    p_mu2[i]	= G4ThreeVector(0,0,0);
    p_mu3[i]	= G4ThreeVector(0,0,0);
    p_mu4[i]	= G4ThreeVector(0,0,0);
    p_mu5[i]	= G4ThreeVector(0,0,0);
    p_mu6[i]	= G4ThreeVector(0,0,0);
    p_mu7[i]	= G4ThreeVector(0,0,0);
    p_mu8[i]	= G4ThreeVector(0,0,0);
    p_mu9[i]	= G4ThreeVector(0,0,0);
    p_mu10[i]	= G4ThreeVector(0,0,0);
    p_mu11[i]	= G4ThreeVector(0,0,0);
    p_mu12[i]	= G4ThreeVector(0,0,0);
    p_mu13[i]	= G4ThreeVector(0,0,0);
    p_mu14[i]	= G4ThreeVector(0,0,0);
    p_mu15[i]	= G4ThreeVector(0,0,0);
    p_mu16[i]	= G4ThreeVector(0,0,0);
    p_mu17[i]	= G4ThreeVector(0,0,0);
    p_mu18[i]	= G4ThreeVector(0,0,0);
    p_mu19[i]	= G4ThreeVector(0,0,0);
    p_mu20[i]	= G4ThreeVector(0,0,0);
    p_mu30[i]	= G4ThreeVector(0,0,0);
    p_mu31[i]	= G4ThreeVector(0,0,0);

    x_mu1[i]	= G4ThreeVector(0,0,0);
    x_mu2[i]	= G4ThreeVector(0,0,0);
    x_mu3[i]	= G4ThreeVector(0,0,0);
    x_mu4[i]	= G4ThreeVector(0,0,0);
    x_mu5[i]	= G4ThreeVector(0,0,0);
    x_mu6[i]	= G4ThreeVector(0,0,0);
    x_mu7[i]	= G4ThreeVector(0,0,0);
    x_mu8[i]	= G4ThreeVector(0,0,0);
    x_mu9[i]	= G4ThreeVector(0,0,0);
    x_mu10[i]	= G4ThreeVector(0,0,0);
    x_mu11[i]	= G4ThreeVector(0,0,0);
    x_mu12[i]	= G4ThreeVector(0,0,0);
    x_mu13[i]	= G4ThreeVector(0,0,0);
    x_mu14[i]	= G4ThreeVector(0,0,0);
    x_mu15[i]	= G4ThreeVector(0,0,0);
    x_mu16[i]	= G4ThreeVector(0,0,0);
    x_mu17[i]	= G4ThreeVector(0,0,0);
    x_mu18[i]	= G4ThreeVector(0,0,0);
    x_mu19[i]	= G4ThreeVector(0,0,0);
    x_mu20[i]	= G4ThreeVector(0,0,0);
    x_mu30[i]	= G4ThreeVector(0,0,0);
    x_mu31[i]	= G4ThreeVector(0,0,0);
    
    p_par1[i]	= G4ThreeVector(0,0,0);
    p_par2[i]	= G4ThreeVector(0,0,0);
    p_par3[i]	= G4ThreeVector(0,0,0);
    p_par4[i]	= G4ThreeVector(0,0,0);
    p_par5[i]	= G4ThreeVector(0,0,0);
    p_par6[i]	= G4ThreeVector(0,0,0);
    p_par7[i]	= G4ThreeVector(0,0,0);
    p_par8[i]	= G4ThreeVector(0,0,0);
    p_par9[i]	= G4ThreeVector(0,0,0);
    p_par10[i]	= G4ThreeVector(0,0,0);
    p_par11[i]	= G4ThreeVector(0,0,0);
    p_par12[i]	= G4ThreeVector(0,0,0);
    p_par13[i]	= G4ThreeVector(0,0,0);
    p_par14[i]	= G4ThreeVector(0,0,0);
    p_par15[i]	= G4ThreeVector(0,0,0);
    p_par16[i]	= G4ThreeVector(0,0,0);
    p_par17[i]	= G4ThreeVector(0,0,0);
    p_par18[i]	= G4ThreeVector(0,0,0);
    p_par19[i]	= G4ThreeVector(0,0,0);
    p_par20[i]	= G4ThreeVector(0,0,0);
    p_par30[i]	= G4ThreeVector(0,0,0);
    p_par31[i]	= G4ThreeVector(0,0,0);

    id1[i]	= 0;
    id2[i]	= 0;
    id3[i]	= 0;
    id4[i]	= 0;
    id5[i]	= 0;
    id6[i]	= 0;
    id7[i]	= 0;
    id8[i]	= 0;
    id9[i]	= 0;
    id10[i]	= 0;
    id11[i]	= 0;
    id12[i]	= 0;
    id13[i]	= 0;
    id14[i]	= 0;
    id15[i]	= 0;
    id16[i]	= 0;
    id17[i]	= 0;
    id18[i]	= 0;
    id19[i]	= 0;
    id20[i]	= 0;  
    id30[i]	= 0;
    id31[i]	= 0;

    id_par1[i]	= 0;
    id_par2[i]	= 0;
    id_par3[i]	= 0;
    id_par4[i]	= 0;
    id_par5[i]	= 0;
    id_par6[i]	= 0; 
    id_par7[i]	= 0;
    id_par8[i]	= 0;
    id_par9[i]	= 0;
    id_par10[i]	= 0;
    id_par11[i]	= 0;
    id_par12[i]	= 0;
    id_par13[i]	= 0;
    id_par14[i]	= 0;
    id_par15[i]	= 0;
    id_par16[i]	= 0;
    id_par17[i]	= 0;
    id_par18[i]	= 0;
    id_par19[i]	= 0;
    id_par20[i]	= 0;  
    id_par30[i]	= 0;
    id_par31[i]	= 0;

    class_par30[i]=0;
    class_par31[i]=0;
    
  }

  x_EXI_piplus	= G4ThreeVector(0,0,0);
  p_EXI_piplus	= G4ThreeVector(0,0,0);

  x_EXI_piminus	= G4ThreeVector(0,0,0);
  p_EXI_piminus	= G4ThreeVector(0,0,0);

  x_EXI_kplus	= G4ThreeVector(0,0,0);
  p_EXI_kplus	= G4ThreeVector(0,0,0);

  x_EXI_kminus	= G4ThreeVector(0,0,0);
  p_EXI_kminus	= G4ThreeVector(0,0,0);

  x_EXI_muplus	= G4ThreeVector(0,0,0);
  p_EXI_muplus	= G4ThreeVector(0,0,0);

  x_EXI_muminus	= G4ThreeVector(0,0,0);
  p_EXI_muminus	= G4ThreeVector(0,0,0);

  x_EXI_k0L	= G4ThreeVector(0,0,0);
  p_EXI_k0L	= G4ThreeVector(0,0,0);

  angnue0OLD 	= -999;
  angnumu0OLD 	= -999;
  anganue0OLD 	= -999;
  anganumu0OLD 	= -999;
}

void SBEventAction::EndOfEventAction(const G4Event* evt)
{
  G4cout<<"CALL SBEventAction::EndOfEventAction "<<G4endl;

  // Accumulates statistic
  // --
  runAct->fillPerEvent(EnergyTarg, EnergyGap, TrackLTarg, TrackLGap);

  //G4cout <<"PROVA "<< p_muplus_PI/CLHEP::GeV << G4endl=0;

  E_numu_PI	= p_DIF_PI_numu.mag();
  E_anumu_PI	= p_DIF_PI_anumu.mag();
  E_numu_K	= p_DIF_K_numu.mag();
  E_anumu_K	= p_DIF_K_anumu.mag();

  p_muplus_PI=p_DIF_PI_muplus.mag();
  p_muminus_PI=p_DIF_PI_muminus.mag();
  p_muplus_K=p_DIF_K_muplus.mag();
  p_muminus_K=p_DIF_K_muminus.mag();

  runAct->fillPerEvent1(p_muplus_PI,
			p_muminus_PI,
			p_muplus_K,
			p_muminus_K,
			E_numu_PI,
			E_anumu_PI,
			E_numu_K,
			E_anumu_K);

  G4double mpi 	= G4PionPlus::PionPlus()->GetPDGMass();
  G4double mK 	= G4KaonPlus::KaonPlus()->GetPDGMass();

  EForw_numu_PI=Probability->EnuForward(p_DIF_PI_piplus,mpi);
  EForw_anumu_PI=Probability->EnuForward(p_DIF_PI_piminus,mpi);
  EForw_numu_K=Probability->EnuForward(p_DIF_K_kplus,mK);
  EForw_anumu_K=Probability->EnuForward(p_DIF_K_kminus,mK);

  G4double PWN = Probability->PowerNorm(EKinProt,PowerMW);       

  if(evt->GetEventID()==0) G4cout <<"SBEventAction: Ek(p) "<< EKinProt << " GeV. Power " << PowerMW << " (MW). Norm fact.: "<< PWN << G4endl;

  if(SBVerbosity>0)G4cout <<"SBEventAction: PROBABS "<< W_numu_PI << " " << W_anumu_PI <<" "<< W_numu_K << " " << W_anumu_K << G4endl;

  if((ParLevel_piminus!=-1)||
     (ParLevel_piplus!=-1)||
     (ParLevel_muminus!=-1)||
     (ParLevel_muplus!=-1)
     ) G4cout<<"SBEventAction hierarchy ";
  if(ParLevel_piminus!=-1)G4cout<<" pi- "<<ParLevel_piminus;
  if(ParLevel_piplus!=-1)G4cout<<" pi+ "<<ParLevel_piplus;
  if(ParLevel_muminus!=-1)G4cout<<" mu- "<<ParLevel_muminus;
  if(ParLevel_muplus!=-1)G4cout<<" mu+ "<<ParLevel_muplus;
  G4cout<<G4endl;

  DumpChannelsInfo();

  //======================================

  for(int j=0;j<n1;j++){
    Probability->proba2(p_par1[j],"kaon",0,1,PWN);
  }
  for(int j=0;j<n2;j++){
    Probability->proba2(p_par2[j],"kaon",0,-1,PWN);
  }
  //-----------------------

  for(int j=0;j<n30;j++){
    Probability->proba2(p_par30[j],"pion",class_par30[j],1,PWN);
  }
  for(int j=0;j<n31;j++){
    Probability->proba2(p_par31[j],"pion",class_par31[j],-1,PWN);
  }

  //-----------------------
  for(int j=0;j<n7;j++){
    Probability->proba3K(p_par7[j],p_nu7[j].mag(),PWN,7);
  }
  for(int j=0;j<n8;j++){
    Probability->proba3K(p_par8[j],p_nu8[j].mag(),PWN,8);
  }
  for(int j=0;j<n9;j++){
    Probability->proba3K(p_par9[j],p_nu9[j].mag(),PWN,9);
  }
  for(int j=0;j<n10;j++){
    Probability->proba3K(p_par10[j],p_nu10[j].mag(),PWN,10);
  }
  for(int j=0;j<n13;j++){
    Probability->proba3K(p_par13[j],p_nu13[j].mag(),PWN,13);
  }
  for(int j=0;j<n14;j++){
    Probability->proba3K(p_par14[j],p_nu14[j].mag(),PWN,14);
  }
  for(int j=0;j<n15;j++){
    Probability->proba3K(p_par15[j],p_nu15[j].mag(),PWN,15);
  }
  for(int j=0;j<n16;j++){
    Probability->proba3K(p_par16[j],p_nu16[j].mag(),PWN,16);
  }
  //-----------------------
  // decay channels producing muons
  for(int j=0;j<n1;j++){
    Probability->probaMu(x_mu1[j],p_mu1[j],p_par1[j],mK,1,2,PWN);
  }
  for(int j=0;j<n2;j++){
    Probability->probaMu(x_mu2[j],p_mu2[j],p_par2[j],mK,-1,2,PWN);
  }
  for(int j=0;j<n9;j++){
    Probability->probaMu(x_mu9[j],p_mu9[j],p_par9[j],mK,1,2,PWN);
  }
  for(int j=0;j<n10;j++){
    Probability->probaMu(x_mu10[j],p_mu10[j],p_par10[j],mK,-1,2,PWN);
  }
  for(int j=0;j<n15;j++){
    Probability->probaMu(x_mu15[j],p_mu15[j],p_par15[j],mK,1,3,PWN);
  }
  for(int j=0;j<n16;j++){
    Probability->probaMu(x_mu16[j],p_mu16[j],p_par16[j],mK,-1,3,PWN);
  }
  //----
  for(int j=0;j<n30;j++){
    if(class_par30[j]==0)Probability->probaMu(x_mu30[j],p_mu30[j],p_par30[j],mpi,1,1,PWN);
    if(class_par30[j]==1)Probability->probaMu(x_mu30[j],p_mu30[j],p_par30[j],mpi,1,2,PWN);
    if(class_par30[j]==2)Probability->probaMu(x_mu30[j],p_mu30[j],p_par30[j],mpi,1,3,PWN);
  }
  for(int j=0;j<n31;j++){
    if(class_par31[j]==0)Probability->probaMu(x_mu31[j],p_mu31[j],p_par31[j],mpi,-1,1,PWN);
    if(class_par31[j]==1)Probability->probaMu(x_mu31[j],p_mu31[j],p_par31[j],mpi,-1,2,PWN);
    if(class_par31[j]==2)Probability->probaMu(x_mu31[j],p_mu31[j],p_par31[j],mpi,-1,3,PWN);
  }
  //-----------------------

  G4int moth=0,chargeDIF=0;
  G4double xx=0,yy=0,zz=0;
  G4double pxx=0,pyy=0,pzz=0;
  G4double pxnu=0,pynu=0,pznu=0;

  if(x_DIF_PI_piplus!=G4ThreeVector(0,0,0)){
    moth=1;
    chargeDIF=1;
    xx = x_DIF_PI_piplus.x();
    yy = x_DIF_PI_piplus.y();
    zz = x_DIF_PI_piplus.z();
    pxx = p_DIF_PI_piplus.x();
    pyy = p_DIF_PI_piplus.y();
    pzz = p_DIF_PI_piplus.z();
    pxnu = p_DIF_PI_numu.x();
    pynu = p_DIF_PI_numu.y();
    pznu = p_DIF_PI_numu.z();
  }
  if(x_DIF_K_kplus!=G4ThreeVector(0,0,0)){
    moth=2;
    chargeDIF=1;
    xx = x_DIF_K_kplus.x();
    yy = x_DIF_K_kplus.y();
    zz = x_DIF_K_kplus.z();
    pxx = p_DIF_K_kplus.x();
    pyy = p_DIF_K_kplus.y();
    pzz = p_DIF_K_kplus.z();
    pxnu = p_DIF_K_numu.x();
    pynu = p_DIF_K_numu.y();
    pznu = p_DIF_K_numu.z();
  }

  if(x_DIF_PI_piminus!=G4ThreeVector(0,0,0)){
    moth=1;
    chargeDIF=-1;
    xx = x_DIF_PI_piminus.x();
    yy = x_DIF_PI_piminus.y();
    zz = x_DIF_PI_piminus.z();
    pxx = p_DIF_PI_piminus.x();
    pyy = p_DIF_PI_piminus.y();
    pzz = p_DIF_PI_piminus.z();
    pxnu = p_DIF_PI_anumu.x();
    pynu = p_DIF_PI_anumu.y();
    pznu = p_DIF_PI_anumu.z();
  }
  if(x_DIF_K_kminus!=G4ThreeVector(0,0,0)){
    moth=2;
    chargeDIF=-1;
    xx = x_DIF_K_kminus.x();
    yy = x_DIF_K_kminus.y();
    zz = x_DIF_K_kminus.z();
    pxx = p_DIF_K_kminus.x();
    pyy = p_DIF_K_kminus.y();
    pzz = p_DIF_K_kminus.z();
    pxnu = p_DIF_K_anumu.x();
    pynu = p_DIF_K_anumu.y();
    pznu = p_DIF_K_anumu.z();
  }

  if((NDIF_piplus+NDIF_piminus+NDIF_Kplus+NDIF_Kminus+NDIF_KzeroS+NDIF_KzeroL)>0)
    G4cout<<"SBEventAction: decays in flight observed for this input track: " <<
      " pi+ "<<NDIF_piplus<<
      " pi- "<<NDIF_piminus<<
      " K+ "<<NDIF_Kplus<<
      " K- "<<NDIF_Kminus<<
      " K0L "<<NDIF_KzeroL<<
      " K0S "<<NDIF_KzeroS<<G4endl;

  runAct->DIFStat(NDIF_piplus,
		  NDIF_piminus,
		  NDIF_Kplus,
		  NDIF_Kminus,
		  NDIF_KzeroS,
		  NDIF_KzeroL);

  if(moth)runAct->fillDIF(xx,yy,zz,
			  pxx,pyy,pzz,
			  pxnu,pynu,pznu,
			  moth,chargeDIF);

  /***********/
  xx = x_EXI_piplus.x();
  yy = x_EXI_piplus.y();
  zz = x_EXI_piplus.z();
  pxx = p_EXI_piplus.x();
  pyy = p_EXI_piplus.y();
  pzz = p_EXI_piplus.z();
  runAct->fillEXIT(xx,yy,zz,pxx,pyy,pzz,1);

  xx = x_EXI_piminus.x();
  yy = x_EXI_piminus.y();
  zz = x_EXI_piminus.z();
  pxx = p_EXI_piminus.x();
  pyy = p_EXI_piminus.y();
  pzz = p_EXI_piminus.z();
  runAct->fillEXIT(xx,yy,zz,pxx,pyy,pzz,2);

  // Fill histo with tracks' parameters at target exit
  // --
  runAct->fillTARG(E_numu_PI);

  // Print per event (modulo n)
  // --
  G4int evtNb = evt->GetEventID();
  if (evtNb%printModulo == 0) {
    //G4cout<<"---> End of event: "<<evtNb<<G4endl;	
    if(SBVerbosity>0)
    G4cout<<"SBEventAction: Target: total energy: " <<std::setw(7)<<G4BestUnit(EnergyTarg,"Energy")<<" total track length: "<<std::setw(7)<<G4BestUnit(TrackLTarg,"Length")<<G4endl;
       //<< "        Gap: total energy: " << std::setw(7)
       //                                 << G4BestUnit(EnergyGap,"Energy")
       //<< "       total track length: " << std::setw(7)
       //                                 << G4BestUnit(TrackLGap,"Length")
       //<< G4endl;  
  }
  //m_hh->theta->Fill(EnergyTarg);

  //check if input file end has been reached
  //in that case end run.
  //runAct->CheckEOF();
}  


void SBEventAction::DumpChannelsInfo(){
  //======================================
  if(n1)G4cout << "Decay 1: K+ -> mu+ numu "<<n1<<G4endl;
  for(int j=0;j<n1;j++){
    if(n1)G4cout << j <<" p_par1 "<< p_par1[j]/CLHEP::GeV <<" p_mu1 "<<p_mu1[j]/CLHEP::GeV<<" x_mu1 "<<x_mu1[j]/CLHEP::m<<" p_nu1 "<<p_nu1[j]/CLHEP::GeV<<G4endl;
  }
  if(n2)G4cout << "Decay 2: "<<n2<<G4endl;
  for(int j=0;j<n2;j++){
    if(n2)G4cout << j <<" p_par2 "<< p_par2[j]/CLHEP::GeV <<" p_mu2 "<<p_mu2[j]/CLHEP::GeV<<" x_mu2 "<<x_mu2[j]/CLHEP::m<<" p_nu2 "<<p_nu2[j]/CLHEP::GeV<<G4endl;
  }
  if(n3)G4cout << "Decay 3: "<<n3<<G4endl;
  for(int j=0;j<n3;j++){
    if(n3)G4cout << j <<" p_par3 "<< p_par3[j]/CLHEP::GeV <<" p_mu3 "<<p_mu3[j]/CLHEP::GeV<<" x_mu3 "<<x_mu3[j]/CLHEP::m<<" p_nu3 "<<p_nu3[j]/CLHEP::GeV<<G4endl;
  }
  if(n4)G4cout << "Decay 4: "<<n4<<G4endl;
  for(int j=0;j<n4;j++){
    if(n4)G4cout << j <<" p_par4 "<< p_par4[j]/CLHEP::GeV <<" p_mu4 "<<p_mu4[j]/CLHEP::GeV<<" x_mu4 "<<x_mu4[j]/CLHEP::m<<" p_nu4 "<<p_nu4[j]/CLHEP::GeV<<G4endl;
  }
  if(n5)G4cout << "Decay 5: "<<n5<<G4endl;
  for(int j=0;j<n5;j++){
    if(n5)G4cout << j <<" p_par5 "<< p_par5[j]/CLHEP::GeV <<" p_mu5 "<<p_mu5[j]/CLHEP::GeV<<" x_mu5 "<<x_mu5[j]/CLHEP::m<<" p_nu5 "<<p_nu5[j]/CLHEP::GeV<<G4endl;
  }
  if(n6)G4cout << "Decay 6: "<<n6<<G4endl;
  for(int j=0;j<n6;j++){
    if(n6)G4cout << j <<" p_par6 "<< p_par6[j]/CLHEP::GeV <<" p_mu6 "<<p_mu6[j]/CLHEP::GeV<<" x_mu6 "<<x_mu6[j]/CLHEP::m<<" p_nu6 "<<p_nu6[j]/CLHEP::GeV<<G4endl;
  }
  if(n7)G4cout << "Decay 7: "<<n7<<G4endl;
  for(int j=0;j<n7;j++){
    if(n7)G4cout << j <<" p_par7 "<< p_par7[j]/CLHEP::GeV <<" p_mu7 "<<p_mu7[j]/CLHEP::GeV<<" x_mu7 "<<x_mu7[j]/CLHEP::m<<" p_nu7 "<<p_nu7[j]/CLHEP::GeV<<G4endl;
  }
  if(n8)G4cout << "Decay 8: "<<n8<<G4endl;
  for(int j=0;j<n8;j++){
    if(n8)G4cout << j <<" p_par8 "<< p_par8[j]/CLHEP::GeV <<" p_mu8 "<<p_mu8[j]/CLHEP::GeV<<" x_mu8 "<<x_mu8[j]/CLHEP::m<<" p_nu8 "<<p_nu8[j]/CLHEP::GeV<<G4endl;
  }
  if(n9)G4cout << "Decay 9: "<<n9<<G4endl;
  for(int j=0;j<n9;j++){
    if(n9)G4cout << j <<" p_par9 "<< p_par9[j]/CLHEP::GeV <<" p_mu9 "<<p_mu9[j]/CLHEP::GeV<<" x_mu9 "<<x_mu9[j]/CLHEP::m<<" p_nu9 "<<p_nu9[j]/CLHEP::GeV<<G4endl;
  }
  if(n10)G4cout << "Decay 10: "<<n10<<G4endl;
  for(int j=0;j<n10;j++){
    if(n10)G4cout << j <<" p_par10 "<< p_par10[j]/CLHEP::GeV <<" p_mu10 "<<p_mu10[j]/CLHEP::GeV<<" x_mu10 "<<x_mu10[j]/CLHEP::m<<" p_nu10 "<<p_nu10[j]/CLHEP::GeV<<G4endl;
  }
  if(n11)G4cout << "Decay 11: "<<n11<<G4endl;
  for(int j=0;j<n11;j++){
    if(n11)G4cout << j <<" p_par11 "<< p_par11[j]/CLHEP::GeV <<" p_mu11 "<<p_mu11[j]/CLHEP::GeV<<" x_mu11 "<<x_mu11[j]/CLHEP::m<<" p_nu11 "<<p_nu11[j]/CLHEP::GeV<<G4endl;
  }
  if(n12)G4cout << "Decay 12: "<<n12<<G4endl;
  for(int j=0;j<n12;j++){
    if(n12)G4cout << j <<" p_par12 "<< p_par12[j]/CLHEP::GeV <<" p_mu12 "<<p_mu12[j]/CLHEP::GeV<<" x_mu12 "<<x_mu12[j]/CLHEP::m<<" p_nu12 "<<p_nu12[j]/CLHEP::GeV<<G4endl;
  }
  if(n13)G4cout << "Decay 13: "<<n13<<G4endl;
  for(int j=0;j<n13;j++){
    if(n13)G4cout << j <<" p_par13 "<< p_par13[j]/CLHEP::GeV <<" p_mu13 "<<p_mu13[j]/CLHEP::GeV<<" x_mu13 "<<x_mu13[j]/CLHEP::m<<" p_nu13 "<<p_nu13[j]/CLHEP::GeV<<G4endl;
  }
  if(n14)G4cout << "Decay 14: "<<n14<<G4endl;
  for(int j=0;j<n14;j++){
    if(n14)G4cout << j <<" p_par14 "<< p_par14[j]/CLHEP::GeV <<" p_mu14 "<<p_mu14[j]/CLHEP::GeV<<" x_mu14 "<<x_mu14[j]/CLHEP::m<<" p_nu14 "<<p_nu14[j]/CLHEP::GeV<<G4endl;
  }
  if(n15)G4cout << "Decay 15: "<<n15<<G4endl;
  for(int j=0;j<n15;j++){
    if(n15)G4cout << j <<" p_par15 "<< p_par15[j]/CLHEP::GeV <<" p_mu15 "<<p_mu15[j]/CLHEP::GeV<<" x_mu15 "<<x_mu15[j]/CLHEP::m<<" p_nu15 "<<p_nu15[j]/CLHEP::GeV<<G4endl;
  }
  if(n16)G4cout << "Decay 16: "<<n16<<G4endl;
  for(int j=0;j<n16;j++){
    if(n16)G4cout << j <<" p_par16 "<< p_par16[j]/CLHEP::GeV <<" p_mu16 "<<p_mu16[j]/CLHEP::GeV<<" x_mu16 "<<x_mu16[j]/CLHEP::m<<" p_nu16 "<<p_nu16[j]/CLHEP::GeV<<G4endl;
  }
  if(n17)G4cout << "Decay 17: "<<n17<<G4endl;
  for(int j=0;j<n17;j++){
    if(n17)G4cout << j <<" p_par17 "<< p_par17[j]/CLHEP::GeV <<" p_mu17 "<<p_mu17[j]/CLHEP::GeV<<" x_mu17 "<<x_mu17[j]/CLHEP::m<<" p_nu17 "<<p_nu17[j]/CLHEP::GeV<<G4endl;
  }
  if(n18)G4cout << "Decay 18: "<<n18<<G4endl;
  for(int j=0;j<n18;j++){
    if(n18)G4cout << j <<" p_par18 "<< p_par18[j]/CLHEP::GeV <<" p_mu18 "<<p_mu18[j]/CLHEP::GeV<<" x_mu18 "<<x_mu18[j]/CLHEP::m<<" p_nu18 "<<p_nu18[j]/CLHEP::GeV<<G4endl;
  }
  if(n19)G4cout << "Decay 19: K0S -> pi+ pi- "<<n19<<G4endl;
  for(int j=0;j<n19;j++){
    if(n19)G4cout << j <<" p_par19 "<< p_par19[j]/CLHEP::GeV <<" p_mu19 "<<p_mu19[j]/CLHEP::GeV<<" x_mu19 "<<x_mu19[j]/CLHEP::m<<" p_nu19 "<<p_nu19[j]/CLHEP::GeV<<
	     " id_sec_pip1 "<<id_sec_pip1_19[j]<<" id_sec_pim1 "<<id_sec_pim1_19[j]<<G4endl;
  }
  if(n20)G4cout << "Decay 20: "<<n20<<G4endl;
  for(int j=0;j<n20;j++){
    if(n20)G4cout << j <<" p_par20 "<< p_par20[j]/CLHEP::GeV <<" p_mu20 "<<p_mu20[j]/CLHEP::GeV<<" x_mu20 "<<x_mu20[j]/CLHEP::m<<" p_nu20 "<<p_nu20[j]/CLHEP::GeV<<G4endl;
  }


  //======================================

  // add the others
  if(n30)G4cout << "Decay 30: "<<n30<<G4endl;
  for(int j=0;j<n30;j++){
    if(n30)G4cout << j <<" p_par30 "<< p_par30[j]/CLHEP::GeV <<" p_mu30 "<<p_mu30[j]/CLHEP::GeV<<" x_mu30 "<<x_mu30[j]/CLHEP::m<<" p_nu30 "<<p_nu30[j]/CLHEP::GeV<<
	     " class "<<class_par30[j]<<" id_par30 "<<id_par30[j]<<G4endl;
  }
  if(n31)G4cout << "Decay 31: "<<n31<<G4endl;
  for(int j=0;j<n31;j++){
    if(n31)G4cout << j <<" p_par31 "<< p_par31[j]/CLHEP::GeV <<" p_mu31 "<<p_mu31[j]/CLHEP::GeV<<" x_mu31 "<<x_mu31[j]/CLHEP::m<<" p_nu31 "<<p_nu31[j]/CLHEP::GeV<<
	     " class "<<class_par31[j]<<" id_par31 "<<id_par31[j]<<G4endl;
  }
}
