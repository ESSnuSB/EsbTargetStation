#include "SBProba.hh"
#include "SBDetectorConstruction.hh"
#include "SBRunAction.hh"
#include "SBEventAction.hh"
#include "SBAnalysisManager.hh"

#include "TF1.h"
#include "G4PionPlus.hh"
#include "G4PionZero.hh"
#include "G4KaonZero.hh"
#include "G4MuonPlus.hh"
#include "G4KaonPlus.hh"
#include "G4Electron.hh"

#include "Randomize.hh"

//to make it a ROOT class
//ClassImp(SBProba)

SBProba::SBProba(SBDetectorConstruction* SBDC, SBRunAction* SBRA/*,SBEventAction* SBEA*/):SBDetector(SBDC),runAct(SBRA)/*,evtAct(SBEA)*/
{
  NBINE=100;
  surf = 100.*CLHEP::m*CLHEP::m; 	// 100m^2
  dist = 100.*CLHEP::km;  		// 100 10^3 m
  multiplicity=1;
  
  mmu  = G4MuonPlus::MuonPlus()->GetPDGMass();
  mK0  = G4KaonZero::KaonZero()->GetPDGMass();
  mK   = G4KaonPlus::KaonPlus()->GetPDGMass();
  me   = G4Electron::Electron()->GetPDGMass();
  mpi  = G4PionPlus::PionPlus()->GetPDGMass();
  mpi0 = G4PionZero::PionZero()->GetPDGMass();

  ctmu = 658.654*CLHEP::m;//
  G4cout << "SBProba::SBProba " << surf/(CLHEP::m*CLHEP::m) << " m2 " <<dist/CLHEP::km << " km " << multiplicity <<" "<< G4endl;
}

SBProba::~SBProba(){;}
// PWN=probability

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//   2 Body decay:  pi -> mu numu
//                   K -> mu numu

void SBProba::proba2(G4ThreeVector vP, 
		     G4String partyp,
		     G4int parfl,
		     G4int q,
		     G4double PWN)
{
  G4double M=0;
  G4double norm=0, weight=0, weightIDEAL=0, E=0;
  G4double beta=0, gamma=0;
  G4double P = vP.mag();
  G4double cth = 0;

  G4double OffAxisAngle=SBDetector->GetOffAxisAngle();
  G4double OffAxisPhi=SBDetector->GetOffAxisPhi();

  G4double eOME[3]={0.,0.,0.};
  eOME[0]=sin(OffAxisAngle)*cos(OffAxisPhi);
  eOME[1]=sin(OffAxisAngle)*sin(OffAxisPhi);
  eOME[2]=cos(OffAxisAngle);

  if(P)cth=((eOME[0]*vP.x())+(eOME[1]*vP.y())+(eOME[2]*vP.z()))/P;
  
  if(partyp=="kaon")M=mK;
  if(partyp=="pion")M=mpi;

  if(P){
    E=sqrt(P*P+M*M);
    beta = P/E;
    gamma = E/M;
    norm = 1.e17*surf/(dist*dist)/(4.*acos(-1.))/multiplicity;
    weight = norm*PWN*(1.-beta*beta)/pow(1.-beta*cth,2);
    weightIDEAL = norm*PWN*(1.-beta*beta)/pow(1.-beta*1.,2);

    G4double enuCM = (M*M-mmu*mmu)/(2.*M);
    G4double enuLAB = 1.+beta*(beta-cth)/(beta*cth-1.);
    enuLAB = enuCM*gamma*enuLAB;
    
    SBAnalysisManager* analysis = SBAnalysisManager::getInstance();
    analysis->FillProba2(q,partyp,parfl,enuLAB,weight,weightIDEAL);
    
    G4cout << "SBProba::proba2 "
	   << surf/(CLHEP::m*CLHEP::m)<<" m2 " 
	   << dist/CLHEP::km<<" km " 
	   << multiplicity<<" "
	   << weight<<" m(GeV) "
	   << M/CLHEP::GeV<<" beta "
	   << beta<<" gamma "
	   << gamma
	   << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double SBProba::PowerNorm(G4double Ek,G4double PowerMW){
  G4double year =1.e7;//s convention
  G4double eV_J = 1.602176*1.e-19;
  G4double pot_s = PowerMW*1.e6/(Ek*1.e9*eV_J);
  G4double pot_y = pot_s*year;
  G4double Pw = pot_y/1.e23;
  return Pw;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//////////////////////////////////////////////////////////////////////////////////////////////inutile
G4double SBProba::EnuForward(G4ThreeVector vP, G4double M){
  //   2 Body decay:  pi -> mu numu
  //                   K -> mu numu
  G4double E=0, beta=0, gamma=0, enuCM=0, enuLAB=0;
  G4double P = vP.mag();
  G4double cth = 0;

  //if(P)cth = vP.z()/P;// check ...

  G4double OffAxisAngle=SBDetector->GetOffAxisAngle();
  G4double OffAxisPhi=SBDetector->GetOffAxisPhi();
  G4double eOME[3]={0.,0.,0.};
  eOME[0]=sin(OffAxisAngle)*cos(OffAxisPhi);
  eOME[1]=sin(OffAxisAngle)*sin(OffAxisPhi);
  eOME[2]=cos(OffAxisAngle);
  if(P)cth=((eOME[0]*vP.x())+(eOME[1]*vP.y())+(eOME[2]*vP.z()))/P;


  //if(P>0.000001){
  if(P>0.000001){
    E=sqrt(P*P+M*M);
    beta = P/E;
    gamma = E/M;
    enuCM = (M*M-mmu*mmu)/(2.*M);
    enuLAB = 1.+beta*(beta-cth)/(beta*cth-1.);
    enuLAB = enuCM*gamma*enuLAB;
  }
  //      G4cout << "exit proba2" << ipart << G4endl;  
  return enuLAB;
}
///////////////////////////////////////////////////////////////////////////////////////////
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SBProba::proba3K(G4ThreeVector pK, 
		      G4double enuLAB,
		      G4double PWN,
		      G4int decflag=7){

  // 3 body kaon decays

  G4cout << "SBProba::proba3K: pK(GeV) "<<pK/CLHEP::GeV<<" enuLAB(GeV) "<<enuLAB/CLHEP::GeV<<" PWN "<<PWN<<" decflag "<<decflag<<G4endl;
  G4cout << "SBProba::proba3K: surf " << surf/(CLHEP::m*CLHEP::m) << " m2 dist " <<dist/CLHEP::km << " km " << multiplicity <<" "<< G4endl;

  G4double cthK=0,leptM=0,piM=0;

  if((decflag==7)||(decflag==8)){leptM=me; piM=mpi0;}
  if((decflag==9)||(decflag==10)){leptM=mmu; piM=mpi0;}  
  if((decflag==13)||(decflag==14)){leptM=me; piM=mpi;}
  if((decflag==15)||(decflag==16)){leptM=mmu; piM=mpi;}

  G4cout << "SBProba::proba3K: lepton M: "<<leptM/CLHEP::MeV<<
    " (MeV). pion M: "<<piM/CLHEP::MeV<<" (MeV)."<<G4endl;

  G4double PK = pK.mag();
  //if(PK)cthK=pK.z()/PK;    

  G4double OffAxisAngle=SBDetector->GetOffAxisAngle();
  G4double OffAxisPhi=SBDetector->GetOffAxisPhi();
  G4double eOME[3]={0.,0.,0.};
  eOME[0]=sin(OffAxisAngle)*cos(OffAxisPhi);
  eOME[1]=sin(OffAxisAngle)*sin(OffAxisPhi);
  eOME[2]=cos(OffAxisAngle);
  if(PK)cthK=((eOME[0]*pK.x())+(eOME[1]*pK.y())+(eOME[2]*pK.z()))/PK;

  G4double EK = sqrt(PK*PK+mK0*mK0);
  G4double betaK = PK/EK;
  G4double gammaK = EK/mK0;
  
  G4double norm = (1e17*surf)/(pow(dist,2)*4.*acos(-1.)*multiplicity);
  G4double cthstr = (betaK-cthK)/((betaK*cthK)-1.);//cos(nu) in K r.f.
  G4double fact = ((1.-(betaK*betaK))/pow(1.-betaK*cthK,2));
  G4double fact1=(1./(gammaK*(1.+(betaK*cthstr))));

  G4cout << "SBProba::proba3K: norm "<<norm << " fact "<<fact<<" cthstr "<<cthstr<<G4endl;

  G4double weight = (norm*fact*fact1)*2./(mK0-leptM-piM);//? mumble mumble
  weight*=PWN;

  G4double weight2 = norm*fact;
  weight2*=PWN;

  G4cout << "SBProba::proba3K:"<<
    " EK (GeV): "<<EK/CLHEP::GeV<<
    " pK (GeV): "<<pK/CLHEP::GeV<<
    " betaK: "<<betaK<<
    " gammaK: "<<gammaK<<
    " cthstr: "<<cthstr<<
    " mult: "<<multiplicity<<
    " Enulab (GeV): "<<enuLAB/CLHEP::GeV<<
    " norm: "<<norm<<
    " fact: "<<fact<<
    " fact1: "<<fact1<<
    " weight: "<<weight<<
    " weight2: "<<weight2<<G4endl;
 
  //G4double de=1;//MeV
  de=1;//MeV
  //G4double enuzero = G4UniformRand();
  //G4double E0=EK-leptM-piM;
  //G4double E0=mK-leptM-piM;
  G4double E0=200;//MeV
  G4int nsampl = int(E0/de);
  G4double e_nu=0;
  G4cout <<"EK-mpi-mlept "<<EK-leptM-piM<<" de "<<de <<" NSAMPL " << nsampl <<G4endl;

  TF1 *f = new TF1("f","2.47459e-5-7.44997e-06*x+5.43311e-07*x*x-2.51867e-09*x*x*x",0,200);
  
  //double ww=0;
  double enustar=0;
  weight=weight*de*gammaK*(1+betaK*cthstr);
  for(int j=0;j<nsampl;j++){
    //enustar = E0 - (j+1-enuzero)*de;
    /*  if(enustar<200){
      ww=2.47459e-5
	-7.44997e-06*enustar
	+5.43311e-07*pow(enustar,2)
	-2.51867e-09*pow(enustar,3);
    */

    enustar=f->GetRandom();
    e_nu = enustar*gammaK*(1+betaK*cthstr);//transform to lab

    //G4cout <<"check "<<e_nu<<" "<<ww<<G4endl;
    

    //////////////////////////////////////////////////////////////
    SBAnalysisManager* analysis = SBAnalysisManager::getInstance();
    analysis->FillProba3K(decflag,e_nu,weight);

    /*
    if(decflag==7)runAct->kaons_nue_du_3->Fill(e_nu/CLHEP::GeV,weight);
    if(decflag==8)runAct->kaons_anue_du_3->Fill(e_nu/CLHEP::GeV,weight);
    if(decflag==9)runAct->kaons_numu_du_3->Fill(e_nu/CLHEP::GeV,weight);
    if(decflag==10)runAct->kaons_anumu_du_3->Fill(e_nu/CLHEP::GeV,weight);
    
    if(decflag==13)runAct->kzeros_nue_du_3->Fill(e_nu/CLHEP::GeV,weight);
    if(decflag==14)runAct->kzeros_anue_du_3->Fill(e_nu/CLHEP::GeV,weight);
    if(decflag==15)runAct->kzeros_numu_du_3->Fill(e_nu/CLHEP::GeV,weight);
    if(decflag==16)runAct->kzeros_anumu_du_3->Fill(e_nu/CLHEP::GeV,weight);
    */
  }

  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

/***************************************/
/***************************************/
/***************************************/

void SBProba::probaMu(
		      G4ThreeVector xmu, 
		      G4ThreeVector pmu,
		      G4ThreeVector ppar, 
		      G4double mpar = 139.,
		      G4int iipart=-1,
		      G4int parflag=0,
		      G4double PWN=1.
		      ){
  // muon decay
  ipart=iipart;
  G4double Emu=0, Epar=0, betaMu=0;
  betaPar=0;
  G4double gammaMu=0, norm=0, weight=0, weight0=0;
  gammaPar=0;
  G4double polaT=0, polaL=0, func=0, fact=0, weightnue=0, weightnue0=0;
  G4double threshold=0, e_nu=0, a=0, probaDKmu=0,cthMu=0,cthMu0=0;
  G4int nsampl=0;
  
  //G4cout <<"enter probaMu " << ipart << G4endl;

  G4double Pmu=pmu.mag();
  G4double Ppar=ppar.mag();

  G4double OffAxisAngle=SBDetector->GetOffAxisAngle();
  G4double OffAxisPhi=SBDetector->GetOffAxisPhi();
  G4double eOME[3]={0.,0.,0.};
  eOME[0]=sin(OffAxisAngle)*cos(OffAxisPhi);
  eOME[1]=sin(OffAxisAngle)*sin(OffAxisPhi);
  eOME[2]=cos(OffAxisAngle);
  if(Pmu)cthMu=((eOME[0]*pmu.x())+(eOME[1]*pmu.y())+(eOME[2]*pmu.z()))/Pmu;
  if(Pmu)cthMu0=pmu.z()/Pmu;

  //G4cout << "brius Ppar "<< Ppar << " Pmu "<< Pmu << G4endl;

  //if(Pmu<0.000001){
  //G4cout << "ProbaMu: P tooo small. return. Pmu = "<< Pmu << G4endl;
  //return;
  //}

  Emu = sqrt(Pmu*Pmu+mmu*mmu);
  betaMu = Pmu/Emu;
  gammaMu = Emu/mmu;

  Epar = sqrt(Ppar*Ppar+mpar*mpar);
  betaPar = Ppar/Epar;
  gammaPar = Epar/mpar;

  runAct->h_MUDEC_pmupar->Fill(Pmu/CLHEP::GeV,Ppar/CLHEP::GeV);

  if((betaMu>=1)&&(cthMu>=1)){
    G4cout <<"SBproba::probaMu ERROR, costheta = "<<cthMu<<" beta = " << betaMu << G4endl;
    return;
  }
  
  // temp!
  //G4int itra=0,itrackSav=0,ievent=0,ieventSav=0;

  G4double pstar = (mpar*mpar-mmu*mmu)/(2.*mpar);
  //G4cout <<"pstar "<< pstar << G4endl;

  //G4double thmu=acos(pmu.z()/Pmu);
  //G4double thpar=acos(ppar.z()/Ppar);
  //thpar=0;
  G4double thpimu = acos((pmu.x()*ppar.x()+pmu.y()*ppar.y()+pmu.z()*ppar.z())/(Pmu*Ppar));

  runAct->h_MUDEC_thpimu_lab->Fill(thpimu);

  //G4cout << "zero "<<Pmu<<" "<<Emu<<" "<<gammaPar<<" "<<betaPar<<" "<<gammaPar*Epar-betaPar*gammaPar*Ppar*cos(thpar)<<" thpimu "<<thpimu<<G4endl;
  //cthstarMu=sin(thpar)*(gammaPar/pstar)*(P*cos(thmu-thpar)-betaPar*Emu);
  //cthstarMu=gammaPar*(Pmu*cos(thmu)-betaPar*Emu)/sqrt(pow(gammaPar*(Pmu*cos(thmu)-betaPar*Emu),2)+pow(Pmu*sin(thmu),2));

  cthstarMu=(gammaPar*(Pmu*cos(thpimu)-betaPar*Emu))/pstar;

  //G4cout <<"cthstarMu "<< cthstarMu << G4endl;

  if(pow(cthstarMu,2)<1.){
    polaT = gammaPar*betaPar*sqrt(1.-pow(cthstarMu,2))/(gammaMu*betaMu);
  }else{
    G4cout << "SBProba::probaMu error: abs(cthstarMu)>1 " << G4endl;
    return;
  }

  if(polaT>1){
    G4cout <<"SBProba::probaMu ERROR, polaT = " << polaT << G4endl;
    G4cout <<"SBProba::probaMu ERROR, pi info = " << gammaPar << " "  << betaPar << G4endl;
    G4cout <<"SBProba::probaMu ERROR, mu info = " << gammaMu <<" "<<betaMu<<" "<<cthstarMu<<G4endl;
    return;
  }

  polaL = sqrt(1.-pow(polaT,2));

  mpi = G4PionPlus::PionPlus()->GetPDGMass();

  threshold = ((mpi*mpi)-(mmu*mmu))/((mpi*mpi)+(mmu*mmu));
  threshold = -threshold/betaPar;

  if(ipart==(-1)) polaL = -polaL;//muminus
  if(cthstarMu>threshold) polaL = -polaL;
      
  if(ipart==1)   runAct->h_MUDEC_polaL_muplus->Fill(polaL);
  if(ipart==(-1))runAct->h_MUDEC_polaL_muminus->Fill(polaL);

  //if((ieventSav!=ievent)||(itrackSav!=itra)){ 
  // no information about the polarization...
  //polaL = 0;
  //}
  
  //G4cout <<"TESTPROBA1 "<< pmu << " " << xmu <<G4endl;

  probaDKmu=dkproba(pmu,xmu);

  //G4cout <<"TESTPROBA "<< probaDKmu <<G4endl;

  //100  FORMAT('probaMu',E14.6,E14.6,E14.6,E14.6,E14.6,E14.6,E14.6,E14.6);

  de = 1;//MeV (G4 standard unit)
  //de = 200;//MeV (G4 standard unit)

  //norm = (de/1000.)*1.e17*surf/pow(dist,2)/(2.*acos(-1.))*probaDKmu/multiplicity;// is dividing by 1/2pi correct ? check factor 2 ins in "fact"
  norm = (de)*1.e17*surf/pow(dist,2)/(4.*acos(-1.))*probaDKmu/multiplicity;

  //cthstarNu = (betaMu-cthMu0)/(betaMu*cthMu0-1.);//of neutrino BUG ?
  cthstarNu = (betaMu-cthMu)/(betaMu*cthMu-1.);//of neutrino BUG ?
  //cthstarNu = (betaMu-cos(thpimu))/(betaMu*cos(thpimu)-1.);//of neutrino
 
  G4cout <<"costh* nu "<<cthstarNu<<G4endl;

  func = (1.-betaMu*betaMu)/pow(1.-betaMu*cthMu,2);
  fact = 2./(Emu*(1.+betaMu*cthstarNu));

  runAct->h_MUDEC_cthstarMu->Fill(cthstarMu);
  runAct->h_MUDEC_cthstarNu->Fill(cthstarNu);

  //     call a random number to avoid systematic effect
  //CALL GRNDM(enuzero,1);
  G4double enuzero = G4UniformRand();// temp
  nsampl = int(Emu/de);

  G4cout <<"Emu "<<Emu<<" de "<<de <<" NSAMPL " << nsampl <<" weight "<<norm*func*fact*PWN<<G4endl;
 
  //G4cout << "CICCIO "<<cthstarMu <<" "<<cthstarNu << G4endl;
  //G4cout << "probaMu power normalization "<<PWN<< G4endl;

  G4double succ_samples =0;
  for(int j=0;j<nsampl;j++){
    // sample neutrino energy in the lab frame from its max value (~Emu) to zero
    e_nu = Emu - (j+1-enuzero)*de;

    //Enu* = Enu/ [gammamu ( 1 + betamu Costhstar)]
    // a is 2Enu*/mmu
    a = fact*e_nu;

    //G4cout << "CICCIO1 "<<Emu <<" "<<e_nu <<" "<<a<<" "<<fact<<" "<<betaMu<<" "<<cthstarNu<<G4endl;
    if((a>=0)&&(a<=1)){
      succ_samples++;

      G4double f0m = 2.*pow(a,2)*(3.-2.*a);
      G4double f1m = 2.*pow(a,2)*(1.-2.*a);

      G4double f0e = 12.*pow(a,2)*(1.-a);
      G4double f1e = f0e;

      //polaL=1;//temp

      //weight = pow(a,2)*( (3.-2.*a) - polaL*cthstarNu*(1.-2.*a) );
      weight0 = f0m - polaL*cthstarNu*f1m;
      weight = (norm*func*fact*PWN)*weight0;

      //weightnue = (pow(a,2)-pow(a,3)) - polaL*cthstarNu*(pow(a,2)-pow(a,3));
      //weightnue = 6.*norm*func*fact*PWN*weightnue; 
      weightnue0 = f0e - polaL*cthstarNu*f1e;
      weightnue = (norm*func*fact*PWN)*weightnue0;
      /*     
	     if((weight<0)||(weightnue<0))
	     G4cout << "probaMu brius "
	     <<ipart<<" "
	     <<e_nu<<" "
	     <<weight<<" "
	     <<weightnue<<" "
	     <<weight0<<" "
	     <<weightnue0<<" "
	     <<norm<<" "
	     <<probaDKmu<<" "
	     <<func<<" "
	     <<fact<<" "
	     <<G4endl;
      */
      G4double Enustar = mmu*a/2.;

      if(ipart==(-1)){//muminus
	runAct->h_MUDEC_Enumustar->Fill(Enustar,weight0);
	runAct->h_MUDEC_Eanuestar->Fill(Enustar,weightnue0);
	//runAct->h_MUDEC_Enumustar->Fill(Enustar);
	//runAct->h_MUDEC_Eanuestar->Fill(Enustar);
      } else if (ipart==1){// muplus
	runAct->h_MUDEC_Eanumustar->Fill(Enustar,weight0);
	runAct->h_MUDEC_Enuestar->Fill(Enustar,weightnue0);
	//runAct->h_MUDEC_Eanumustar->Fill(Enustar);
	//runAct->h_MUDEC_Enuestar->Fill(Enustar);
      }

      //G4cout << "brius " << Enustar << " " << weight <<" "<<weightnue<< G4endl;
      ///////////////////////////////////////////////////////////////////////////////////////////////

      SBAnalysisManager* analysis = SBAnalysisManager::getInstance();
      analysis->FillprobaMu(ipart,parflag,e_nu,weight,weightnue);
      
      ////////////////////////////////////////////////////////////////////////////////////////////////////
      
    } else {// a in [0,1]
      //G4cout <<"SBProba::probaMu WARNING (sampling ineff.), a = " << a << G4endl;
    } // a in [0,1]
  }// loop
  

  G4cout <<"SBProba::probaMu successful samplings: "<< succ_samples<<" over "<<nsampl<<G4endl;
  //      if(ipart==6) G4cout <<"exit probaMu " << Probae(5)<<" " <<Probnm(5) << G4endl;
  //      G4cout <<" exit proba " << ipart << G4endl


  G4double rmu =sqrt(pow(xmu.x(),2)+pow(xmu.y(),2));

  runAct->h_MUDEC_probaTunnel->Fill(probaDKmu);

  runAct->h_MUDEC_probaTunnel_pmu->Fill(Pmu/CLHEP::GeV,probaDKmu);
  runAct->h_MUDEC_probaTunnel_thmu->Fill(acos(cthMu0),probaDKmu);
  runAct->h_MUDEC_probaTunnel_rmu->Fill(rmu/CLHEP::m,probaDKmu);

  runAct->h_MUDEC_muxy->Fill(xmu.x()/CLHEP::m,xmu.y()/CLHEP::m);
  runAct->h_MUDEC_muzr->Fill(xmu.z()/CLHEP::m,rmu/CLHEP::m);

  return;

}

/***************************************/
/***************************************/
/***************************************/

G4double SBProba::dkproba(G4ThreeVector pmu, G4ThreeVector xmuon){


  //G4int vbl = evtAct->GetSBVerbosity();
  G4int vbl = 1;
  G4ThreeVector xmu,umu;
  G4double r0=SBDetector->GetTunnelRadius();
  //G4double z0=0.5*SBDetector->GetHallSizeZ()-SBDetector->GetTunnelLength();
  G4double z0=SBDetector->GetHallSizeZ();
  
  if(vbl>0)
    G4cout << "dkproba pmu "<< pmu/CLHEP::GeV << " (GeV) r0 "<< r0/CLHEP::m <<" (m) z0 "<<z0/CLHEP::m<<" (m) xmuon "<< xmuon/CLHEP::m<<" (m)"<< G4endl;

  // put the origin at the beginning of the tunnel (correctly of all the volume indeed!)
  xmu.setX(xmuon.x());
  xmu.setY(xmuon.y());
  //xmu.setZ(xmuon.z()-z0);
  xmu.setZ(xmuon.z()+0.5*z0);
  
  G4double pathmu=0,probmu=0,param=0;
  G4double apmu=pmu.mag();
  if(!apmu){
    G4cout << " dkproba ERROR pmu = 0 "<< G4endl;
    probmu = 0;
    return probmu;
  }
  G4double path0mu = ctmu*apmu/mmu;
  if(vbl>0)
    G4cout << "dkproba xmu "<< xmu/CLHEP::m <<  " (m) ctmu "<< ctmu/CLHEP::m << " (m) pmu " << apmu/CLHEP::GeV <<" (GeV) mmu "<< mmu/CLHEP::GeV <<" (GeV) path0mu "<<path0mu/CLHEP::m<<" (m) "<< G4endl;

  umu.setX(pmu.x()/pmu.mag());
  umu.setY(pmu.y()/pmu.mag());
  umu.setZ(pmu.z()/pmu.mag());

  //dvdot(umu,umu,2,ap);
  G4double ap=(umu.x()*umu.x())+(umu.y()*umu.y());

  if(ap==0.){
    G4cout <<" dkproba:: ap = 0. axis parallel muon!"<<G4endl;
    param=z0/umu.z();
    pathmu=z0;
    probmu=1.-exp(-pathmu/path0mu);
    return probmu;
  }
  
  //dvdot(umu,xmu,2,bp);
  //dvdot(xmu,xmu,2,cp);
  G4double bp = (umu.x()*xmu.x())+(umu.y()*xmu.y());
  G4double cp = (xmu.x()*xmu.x())+(xmu.y()*xmu.y())-pow(r0,2);

  if((bp*bp-ap*cp)<0.){
    probmu = 0.;
    G4cout <<" dkproba ERROR sqrt < 0 ! a "<<ap<<" b "<<bp<<" c "<<cp<<G4endl;
    return probmu;
  }

  param = (-bp+sqrt(bp*bp-ap*cp))/ap;
  G4ThreeVector posfin = G4ThreeVector(0,0,0);

  // Pion on the wall
  if (param*umu.z()<=(z0-xmu.z())) {
    if(vbl>0)
      G4cout << "dkproba: mu exit on the wall" <<G4endl;
    posfin.setX(umu.x()*param+xmu.x());
    posfin.setY(umu.y()*param+xmu.y());
    posfin.setZ(umu.z()*param+xmu.z());
  } else {
    // Pion on the end plate
    if(vbl>0)
      G4cout << "dkproba: mu exit on the end-plate" <<G4endl;
    param = (z0-xmu.z())/umu.z();
    posfin.setX(umu.x()*param+xmu.x());
    posfin.setY(umu.y()*param+xmu.y());
    posfin.setZ(z0);
  }
  
  // the path available for decay
  G4ThreeVector vdiff = posfin-xmu;
  pathmu = vdiff.mag();

  if(pathmu<=0.){
    probmu = 0.;
    G4cout <<" dkproba: pathmu <= 0 ! posfin "<< posfin <<
      " param "<<param<<
      " a "<<ap<<
      " b "<<bp<<
      " c "<<cp<<
      G4endl;
  }

  runAct->h_MUDEC_pathTunnel->Fill(pathmu/CLHEP::m);

  // the probability is

  probmu=1.-exp(-pathmu/path0mu);
  if(vbl>0)
    G4cout << " dkproba: probmu " << probmu << " pathmu " <<pathmu/CLHEP::m <<" (m) path0mu "<<path0mu/CLHEP::m<<" (m) "<< G4endl;
  
  return probmu;
}
