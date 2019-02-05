#include "SBSteppingAction.hh"
#include "SBDetectorConstruction.hh"
#include "SBPrimaryGeneratorAction.hh"
#include "SBEventAction.hh"
#include "SBRunAction.hh"

#include "G4Step.hh"
#include "G4VProcess.hh"

#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4KaonZeroLong.hh"
#include "G4KaonZeroShort.hh"
#include "G4KaonZero.hh"
#include "G4AntiKaonZero.hh"
#include "G4Lambda.hh"
#include "G4AntiLambda.hh"

SBSteppingAction::SBSteppingAction(SBDetectorConstruction* det,SBEventAction* evt,SBPrimaryGeneratorAction* primGen, SBRunAction* run)
  :detector(det), evtAct(evt), primaryGen(primGen), runAct(run)
{ 
  pOutFileExitTunnel.open("ExitTunnel.data");
  pOutFileExitTarget.open("ExitTarget.data");
}

SBSteppingAction::~SBSteppingAction()
{
  pOutFileExitTunnel.close();
  pOutFileExitTarget.close();
}

void SBSteppingAction::UserSteppingAction(const G4Step* aStep)
{

  G4bool DEB = false;
  if(evtAct->GetSBVerbosity()>0)DEB=true;
  
  G4String INUMU 	= G4NeutrinoMu::NeutrinoMu()->GetParticleName();
  G4String IANUMU 	= G4AntiNeutrinoMu::AntiNeutrinoMu()->GetParticleName();
  G4String INUE 	= G4NeutrinoE::NeutrinoE()->GetParticleName();
  G4String IANUE 	= G4AntiNeutrinoE::AntiNeutrinoE()->GetParticleName();
  G4String IPIMINUS 	= G4PionMinus::PionMinus()->GetParticleName();
  G4String IPIPLUS 	= G4PionPlus::PionPlus()->GetParticleName();
  G4String IMUMINUS 	= G4MuonMinus::MuonMinus()->GetParticleName();
  G4String IMUPLUS 	= G4MuonPlus::MuonPlus()->GetParticleName();
  G4String IKMINUS 	= G4KaonMinus::KaonMinus()->GetParticleName();
  G4String IKPLUS 	= G4KaonPlus::KaonPlus()->GetParticleName();

  //G4int PDGpiplus =  G4PionPlus::PionPlus()->GetPDGEncoding();
  //G4int PDGpiminus =  G4PionMinus::PionMinus()->GetPDGEncoding();
  G4int PDGkplus 	=  G4KaonPlus::KaonPlus()->GetPDGEncoding();
  G4int PDGkminus 	=  G4KaonMinus::KaonMinus()->GetPDGEncoding();
  G4int PDGk0L 		=  G4KaonZeroLong::KaonZeroLong()->GetPDGEncoding();
  G4int PDGk0S 		=  G4KaonZeroShort::KaonZeroShort()->GetPDGEncoding();
  G4int PDGk0 		=  G4KaonZero::KaonZero()->GetPDGEncoding();
  G4int PDGak0 		=  G4AntiKaonZero::AntiKaonZero()->GetPDGEncoding();
  G4int PDGL 		=  G4Lambda::Lambda()->GetPDGEncoding();
  G4int PDGaL 		=  G4AntiLambda::AntiLambda()->GetPDGEncoding();

  //G4cout<<"CALL SBSteppingAction::UserSteppingAction"<<G4endl;

  G4VPhysicalVolume* volume0 = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();  

  G4double edep 	= aStep->GetTotalEnergyDeposit();
  G4double stepl 	= 0.;
  if(aStep->GetTrack()->GetDefinition()->GetPDGCharge() != 0.)
    stepl = aStep->GetStepLength();
  
  if (volume0 == detector->GetTarget()) evtAct->AddTarg(edep,stepl);

  G4StepPoint* point1 	= aStep->GetPreStepPoint();
  G4StepPoint* point2 	= aStep->GetPostStepPoint();

  G4ThreeVector pos1 	= point1->GetPosition();
  G4ThreeVector pos2 	= point2->GetPosition();  
  G4ThreeVector mom1 	= point1->GetMomentum();
  G4ThreeVector mom2 	= point2->GetMomentum();  

  G4TouchableHandle touch1 = point1->GetTouchableHandle();
  G4VPhysicalVolume* volume = touch1->GetVolume();//To get the current volume:
  G4String name = volume->GetName();

  G4TouchableHandle touch2 = point2->GetTouchableHandle();

  G4Track* track 	= aStep->GetTrack();
  G4int PID 		= track->GetDynamicParticle()->GetPDGcode();

  G4ThreeVector momentum = track->GetMomentum();
  G4ThreeVector pos 	 = track->GetPosition();
  //G4double kinEnergy     = track->GetKineticEnergy();
  //G4double globalTime    = track->GetGlobalTime();  

  G4double xprim=0,yprim=0,zprim=0,pxprim=0,pyprim=0,pzprim=0;
  xprim 	= primaryGen->GetPrimaryX();
  yprim 	= primaryGen->GetPrimaryY();
  zprim 	= primaryGen->GetPrimaryZ();
  pxprim 	= primaryGen->GetPrimaryPX();
  pyprim 	= primaryGen->GetPrimaryPY();
  pzprim 	= primaryGen->GetPrimaryPZ();

  G4double pTOTprim 	= sqrt(pxprim*pxprim+pyprim*pyprim+pzprim*pzprim);
  G4double pTOT 	= sqrt(momentum.x()*momentum.x()+momentum.y()*momentum.y()+momentum.z()*momentum.z());
  G4double pTOT2 	= sqrt(mom2.x()*mom2.x()+mom2.y()*mom2.y()+mom2.z()*mom2.z());
  G4double polarangle 	= 0;
  if(pTOT2!=0.)polarangle = acos(mom2.z()/pTOT2);
  G4ThreeVector mdir 	= track->GetMomentumDirection();

  if(runAct->h_geom_g4->GetEntries()<100){
    if ((point1->GetStepStatus()==fGeomBoundary)||(point1->GetStepStatus()==fWorldBoundary)){
      runAct->h_geom_g4->Fill(pos1.z()/CLHEP::m,sqrt(pow(pos1.x()/CLHEP::m,2)+pow(pos1.y()/CLHEP::m,2)));    
    }
  }

  G4int NMA=10;
  bool EXITFILEDISABLE=false;
  if(EXITFILEDISABLE==false){
    G4int pid_out = 13;
    G4StepPoint* thePrePoint = aStep->GetPreStepPoint();
    G4VPhysicalVolume* thePrePV = thePrePoint->GetPhysicalVolume();
    G4StepPoint* thePostPoint = aStep->GetPostStepPoint();
    G4VPhysicalVolume* thePostPV = thePostPoint->GetPhysicalVolume();
    if(thePostPoint->GetStepStatus()==fGeomBoundary ){
      if(thePostPV){
	if((thePostPV->GetName()=="World")&& 
	   (thePrePV->GetName()=="TARG")) {      
	  if(track->GetDefinition()->GetParticleName()=="pi+")pid_out=13;
	  if(track->GetDefinition()->GetParticleName()=="pi-")pid_out=14;
	  if(track->GetDefinition()->GetParticleName()=="kaon+")pid_out=15;
	  if(track->GetDefinition()->GetParticleName()=="kaon-")pid_out=16;
	  if(track->GetDefinition()->GetParticleName()=="kaon0")pid_out=24;
	  if(track->GetDefinition()->GetParticleName()=="anti_kaon0")pid_out=25;
	  if(track->GetDefinition()->GetParticleName()=="kaon0S")pid_out=19;
	  if(track->GetDefinition()->GetParticleName()=="kaon0L")pid_out=12;

	  if(track->GetDefinition()->GetParticleName()=="kaon+"){
	    if(primaryGen->n_kplus_exit++<NMA){
	      primaryGen->x_kplus_exit[primaryGen->n_kplus_exit]=pos2.x();
	      primaryGen->y_kplus_exit[primaryGen->n_kplus_exit]=pos2.y();
	      primaryGen->z_kplus_exit[primaryGen->n_kplus_exit]=pos2.z();
	      primaryGen->px_kplus_exit[primaryGen->n_kplus_exit]=momentum.x();
	      primaryGen->py_kplus_exit[primaryGen->n_kplus_exit]=momentum.y();
	      primaryGen->pz_kplus_exit[primaryGen->n_kplus_exit]=momentum.z();
	      primaryGen->n_kplus_exit++;
	    }
	  }
	  if(track->GetDefinition()->GetParticleName()=="kaon-"){
	    if(primaryGen->n_kminus_exit++<NMA){
	      primaryGen->x_kminus_exit[primaryGen->n_kminus_exit]=pos2.x();
	      primaryGen->y_kminus_exit[primaryGen->n_kminus_exit]=pos2.y();
	      primaryGen->z_kminus_exit[primaryGen->n_kminus_exit]=pos2.z();
	      primaryGen->px_kminus_exit[primaryGen->n_kminus_exit]=momentum.x();
	      primaryGen->py_kminus_exit[primaryGen->n_kminus_exit]=momentum.y();
	      primaryGen->pz_kminus_exit[primaryGen->n_kminus_exit]=momentum.z();
	      primaryGen->n_kminus_exit++;
	    }
	  }
	  if(track->GetDefinition()->GetParticleName()=="kaon0S"){
	    if(primaryGen->n_k0s_exit++<NMA){
	      primaryGen->x_k0s_exit[primaryGen->n_k0s_exit]=pos2.x();
	      primaryGen->y_k0s_exit[primaryGen->n_k0s_exit]=pos2.y();
	      primaryGen->z_k0s_exit[primaryGen->n_k0s_exit]=pos2.z();
	      primaryGen->px_k0s_exit[primaryGen->n_k0s_exit]=momentum.x();
	      primaryGen->py_k0s_exit[primaryGen->n_k0s_exit]=momentum.y();
	      primaryGen->pz_k0s_exit[primaryGen->n_k0s_exit]=momentum.z();
	      primaryGen->n_k0s_exit++;
	    }
	  }
	  if(track->GetDefinition()->GetParticleName()=="kaon0L"){
	    if(primaryGen->n_k0l_exit++<NMA){
	      primaryGen->x_k0l_exit[primaryGen->n_k0l_exit]=pos2.x();
	      primaryGen->y_k0l_exit[primaryGen->n_k0l_exit]=pos2.y();
	      primaryGen->z_k0l_exit[primaryGen->n_k0l_exit]=pos2.z();
	      primaryGen->px_k0l_exit[primaryGen->n_k0l_exit]=momentum.x();
	      primaryGen->py_k0l_exit[primaryGen->n_k0l_exit]=momentum.y();
	      primaryGen->pz_k0l_exit[primaryGen->n_k0l_exit]=momentum.z();
	      primaryGen->n_k0l_exit++;
	    }
	  }

	  if((track->GetDefinition()->GetParticleName()=="pi+")||
	     (track->GetDefinition()->GetParticleName()=="pi-")||
	     (track->GetDefinition()->GetParticleName()=="kaon+")||
	     (track->GetDefinition()->GetParticleName()=="kaon-")||
	     (track->GetDefinition()->GetParticleName()=="kaon0S")||
	     (track->GetDefinition()->GetParticleName()=="kaon0L")||
	     (track->GetDefinition()->GetParticleName()=="kaon0")||
	     (track->GetDefinition()->GetParticleName()=="anti_kaon0")){
	    pOutFileExitTarget <<primaryGen->GetFLUKAPOTs()<<" "
			       <<pid_out<<" "
			       <<pos2.x()/CLHEP::cm<<" "
			       <<pos2.y()/CLHEP::cm<<" "
			       <<(pos2.z()-detector->GetTargetZetaPos()+0.5*detector->GetTargetThickness())/CLHEP::cm<<" "
			       <<momentum.x()/CLHEP::GeV<<" "
			       <<momentum.y()/CLHEP::GeV<<" "
			       <<momentum.z()/CLHEP::GeV<<" 1. 1. 1."
			       <<G4endl;
	  }
	}
      }
    }
    
    // Identify piplus, piminus, mu+, mu- at the end of the TUNL and store the position and momentum 
    // at the event level. 
    // --
    if ((point1->GetStepStatus() == fGeomBoundary)&&(name=="TUNL")){
      pOutFileExitTunnel <<pos1.x()/CLHEP::cm<<" "<<pos1.y()/CLHEP::cm<<" "<<pos1.z()/CLHEP::cm<<" "
			 << momentum.x()/CLHEP::GeV<<" "<<momentum.y()/CLHEP::GeV<<" "<<momentum.z()/CLHEP::GeV<<" 1. " << PID << " "
			 << xprim<<" "<<yprim<<" "<<zprim<<" "<<pxprim<<" "<<pyprim<<" "<<pzprim<<G4endl;
      
      if(DEB)G4cout <<"SBSteppingAction: "<< track->GetDefinition()->GetParticleName() << " enters TUNNEL. pos = "<<pos1<<" p = "<< pTOT << G4endl;
      
      if(track->GetDefinition()->GetParticleName()=="pi+"){
	evtAct->x_EXI_piplus=pos1;
	evtAct->p_EXI_piplus=momentum;
      }

      if(track->GetDefinition()->GetParticleName()=="pi-"){
	evtAct->x_EXI_piminus=pos1;
	evtAct->p_EXI_piminus=momentum;
      }
      
      if(track->GetDefinition()->GetParticleName()=="mu+"){
	evtAct->x_EXI_muplus=pos1;
	evtAct->p_EXI_muplus=momentum;
      }
      if(track->GetDefinition()->GetParticleName()=="mu-"){
	evtAct->x_EXI_muminus=pos1;
	evtAct->p_EXI_muminus=momentum;
      }
    }
 
  }
  
  const G4VProcess* pProcess = point2->GetProcessDefinedStep();
  G4String theProcessName=" ";
  if(pProcess)theProcessName = pProcess->GetProcessName();
  const G4String theParticleName = track->GetDefinition()->GetParticleName();

  G4int index = track->GetTrackID();
  G4int indexp = track->GetParentID();
  G4int IMAX=evtAct->MAXIND;
  if((index<IMAX)&&(indexp<IMAX)){
    evtAct->vPDG_ID[index]=PID;
    evtAct->vPDG_PAR_ID[index]=evtAct->vPDG_ID[indexp];
    evtAct->vPDG_PARPAR_ID[index]=evtAct->vPDG_PAR_ID[indexp];//?
    // names
    evtAct->vName_ID[index]=theParticleName;
    evtAct->vName_PAR_ID[index]=evtAct->vName_ID[indexp];
    evtAct->vName_PARPAR_ID[index]=evtAct->vName_PAR_ID[indexp];//?
  }else{}

  G4int PARID=0;// = track->GetParentID();
  G4int parPDG=0;// = evtAct->vPDG_ID[PARID];
  G4int parparPDG=0;// = evtAct->vPDG_PAR_ID[PARID];
  G4int parparparPDG=0;// = evtAct->vPDG_PARPAR_ID[PARID];
  
  G4String parName="";// = evtAct->vName_ID[PARID];
  G4String parparName="";// = evtAct->vName_PAR_ID[PARID];
  G4String parparparName="";// = evtAct->vName_PARPAR_ID[PARID];
  
  PARID = track->GetParentID();
  if(PARID<IMAX){
    parPDG = evtAct->vPDG_ID[PARID];
    parparPDG = evtAct->vPDG_PAR_ID[PARID];
    parparparPDG = evtAct->vPDG_PARPAR_ID[PARID];
    
    parName = evtAct->vName_ID[PARID];
    parparName = evtAct->vName_PAR_ID[PARID];
    parparparName = evtAct->vName_PARPAR_ID[PARID];
  } else {}
  
  G4bool kaonfath = false, kaonGfath = false, kaonGGfath = false;
  G4bool k0fath = false, k0Gfath = false, k0GGfath = false;

  kaonfath =((parPDG==PDGkplus)||(parPDG==PDGkminus));
  kaonGfath =((parparPDG==PDGkplus)||(parparPDG==PDGkminus));
  kaonGGfath = ((parparparPDG==PDGkplus)||(parparparPDG==PDGkminus));

  k0fath =
    ((parPDG==PDGk0)||
     (parPDG==PDGak0)||
     (parPDG==PDGL)||// added lambdas in the definition
     (parPDG==PDGaL)||
     (parPDG==PDGk0S)||
     (parPDG==PDGk0L));
  k0Gfath =
    ((parparPDG==PDGk0)||
     (parparPDG==PDGak0)||
     (parparPDG==PDGL)||// added lambdas in the definition
     (parparPDG==PDGaL)||
     (parparPDG==PDGk0S)||
     (parparPDG==PDGk0L));
  k0GGfath =
    ((parparparPDG==PDGk0)||
     (parparparPDG==PDGak0)||
     (parparparPDG==PDGL)||// added lambdas in the definition
     (parparparPDG==PDGaL)||
     (parparparPDG==PDGk0S)||
     (parparparPDG==PDGk0L));


  // Kill electron muons gamma pi0
  // --
  G4bool KILLMUENABLE = true;
  if(KILLMUENABLE){
    if(theParticleName=="mu-"){
      aStep->GetTrack()->SetTrackStatus(fStopAndKill);
      aStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
      //G4cout <<"SBSteppingAction: mu- killed" << G4endl;
    }
    if(theParticleName=="mu+"){
      aStep->GetTrack()->SetTrackStatus(fStopAndKill);
      aStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
      //G4cout <<"SBSteppingAction: mu+ killed" << G4endl;
    }
  }
  if(theParticleName=="e+"){
    aStep->GetTrack()->SetTrackStatus(fStopAndKill);
    aStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
    //G4cout <<"SBSteppingAction: e+ killed" << G4endl;
  }
  if(theParticleName=="e-"){
    aStep->GetTrack()->SetTrackStatus(fStopAndKill);
    aStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
    //G4cout <<"SBSteppingAction: e- killed" << G4endl;
  }
  if(theParticleName=="gamma"){
    aStep->GetTrack()->SetTrackStatus(fStopAndKill);
    aStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
    //G4cout <<"SBSteppingAction: gamma killed" << G4endl;
  }
  if(theParticleName=="pi0"){
    aStep->GetTrack()->SetTrackStatus(fStopAndKill);
    aStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
    //G4cout <<"SBSteppingAction: pi0 killed" << G4endl;
  }

  //***************************************************
  //***************************************************
  //***************************************************
  // check nue
  // direct estimate without weighting
  // meaningful only if kaon replication is off and muons/etc are not killed

  bool NUECHECK_ENABLE=false;
  if(NUECHECK_ENABLE==true){
      
    if((PID==12)||(PID==14)||(PID==-12)||(PID==-14)){//nu
      G4ThreeVector nup = track->GetMomentum();
      //G4double nupx = nup.x();
      //G4double nupy = nup.y();
      //G4double nupz = nup.z();
      //G4double nuE = sqrt((nupx*nupx)+(nupy*nupy)+(nupz*nupz));
	
      G4double OffAxisAngle=11.8/810.;//NOVA
      G4double OffAxisPhi=0;
      G4double eOME[3]={0.,0.,0.};
      eOME[0]=sin(OffAxisAngle)*cos(OffAxisPhi);
      eOME[1]=sin(OffAxisAngle)*sin(OffAxisPhi);
      eOME[2]=cos(OffAxisAngle);
	
      //G4double angnue=0, angnumu=0;
      //G4double angnue0=0, angnumu0=0;	
      //G4double anganue=0, anganumu=0;
      //G4double anganue0=0, anganumu0=0;	

      G4cout<< "neutrino "<<PID<<" tagged"<<G4endl;

      //G4double cone[8]={0.1,0.05,0.025,0.0125,0.00625,0.003125,0.0015625,0.00078125};

/*
      if(nuE!=0){
	if(PID==12){
	  angnue=acos(((eOME[0]*nupx)+(eOME[1]*nupy)+(eOME[2]*nupz))/nuE);	  
	  angnue0=acos(nupz/nuE);
	    
	  if(angnue0!=evtAct->angnue0OLD){// to avoid double counting
	    if(angnue0<cone[0])runAct->h_G4nue_1->Fill(nuE/CLHEP::GeV);
	    if(angnue0<cone[1])runAct->h_G4nue_2->Fill(nuE/CLHEP::GeV);
	    if(angnue0<cone[2])runAct->h_G4nue_3->Fill(nuE/CLHEP::GeV);
	    if(angnue0<cone[3])runAct->h_G4nue_4->Fill(nuE/CLHEP::GeV);
	    if(angnue0<cone[4])runAct->h_G4nue_5->Fill(nuE/CLHEP::GeV);
	    if(angnue0<cone[5])runAct->h_G4nue_6->Fill(nuE/CLHEP::GeV);
	    if(angnue0<cone[6])runAct->h_G4nue_7->Fill(nuE/CLHEP::GeV);
	    if(angnue0<cone[7])runAct->h_G4nue_8->Fill(nuE/CLHEP::GeV);
	      
	    if(angnue<cone[0])runAct->h_G4nueOA_1->Fill(nuE/CLHEP::GeV);
	    if(angnue<cone[1])runAct->h_G4nueOA_2->Fill(nuE/CLHEP::GeV);
	    if(angnue<cone[2])runAct->h_G4nueOA_3->Fill(nuE/CLHEP::GeV);
	    if(angnue<cone[3])runAct->h_G4nueOA_4->Fill(nuE/CLHEP::GeV);
	    if(angnue<cone[4])runAct->h_G4nueOA_5->Fill(nuE/CLHEP::GeV);	      
	    if(angnue<cone[5])runAct->h_G4nueOA_6->Fill(nuE/CLHEP::GeV);	      
	    if(angnue<cone[6])runAct->h_G4nueOA_7->Fill(nuE/CLHEP::GeV);	      
	    if(angnue<cone[7])runAct->h_G4nueOA_8->Fill(nuE/CLHEP::GeV);	      
            G4cout <<"SBSteppingAction: nue_OFFAXIS FIRST ONE "<< nuE/CLHEP::GeV<< G4endl;
	    evtAct->angnue0OLD=angnue0;
	  }
	}//nue
	  	  
	if(PID==14){
	  angnumu=acos(((eOME[0]*nupx)+(eOME[1]*nupy)+(eOME[2]*nupz))/nuE);	  
	  angnumu0=acos(nupz/nuE);
	    
	  if(angnumu0!=evtAct->angnumu0OLD){// to avoid double counting
	    if(angnumu0<cone[0])runAct->h_G4numu_1->Fill(nuE/CLHEP::GeV);
	    if(angnumu0<cone[1])runAct->h_G4numu_2->Fill(nuE/CLHEP::GeV);
	    if(angnumu0<cone[2])runAct->h_G4numu_3->Fill(nuE/CLHEP::GeV);
	    if(angnumu0<cone[3])runAct->h_G4numu_4->Fill(nuE/CLHEP::GeV);
	    if(angnumu0<cone[4])runAct->h_G4numu_5->Fill(nuE/CLHEP::GeV);
	    if(angnumu0<cone[5])runAct->h_G4numu_6->Fill(nuE/CLHEP::GeV);
	    if(angnumu0<cone[6])runAct->h_G4numu_7->Fill(nuE/CLHEP::GeV);
	    if(angnumu0<cone[7])runAct->h_G4numu_8->Fill(nuE/CLHEP::GeV);
	      
	    if(angnumu<cone[0])runAct->h_G4numuOA_1->Fill(nuE/CLHEP::GeV);
	    if(angnumu<cone[1])runAct->h_G4numuOA_2->Fill(nuE/CLHEP::GeV);
	    if(angnumu<cone[2])runAct->h_G4numuOA_3->Fill(nuE/CLHEP::GeV);
	    if(angnumu<cone[3])runAct->h_G4numuOA_4->Fill(nuE/CLHEP::GeV);
	    if(angnumu<cone[4])runAct->h_G4numuOA_5->Fill(nuE/CLHEP::GeV);		
	    if(angnumu<cone[5])runAct->h_G4numuOA_6->Fill(nuE/CLHEP::GeV);		
	    if(angnumu<cone[6])runAct->h_G4numuOA_7->Fill(nuE/CLHEP::GeV);		
	    if(angnumu<cone[7])runAct->h_G4numuOA_8->Fill(nuE/CLHEP::GeV);		
	    G4cout <<"SBSteppingAction: numu_OFFAXIS FIRST ONE "<< nuE/CLHEP::GeV<< G4endl;
	    evtAct->angnumu0OLD=angnumu0;
	  }	    
	}//numu	  

	if(PID==-12){
	  anganue=acos(((eOME[0]*nupx)+(eOME[1]*nupy)+(eOME[2]*nupz))/nuE);	  
	  anganue0=acos(nupz/nuE);
	    
	  if(anganue0!=evtAct->anganue0OLD){// to avoid double counting
	    if(anganue0<cone[0])runAct->h_G4anue_1->Fill(nuE/CLHEP::GeV);
	    if(anganue0<cone[1])runAct->h_G4anue_2->Fill(nuE/CLHEP::GeV);
	    if(anganue0<cone[2])runAct->h_G4anue_3->Fill(nuE/CLHEP::GeV);
	    if(anganue0<cone[3])runAct->h_G4anue_4->Fill(nuE/CLHEP::GeV);
	    if(anganue0<cone[4])runAct->h_G4anue_5->Fill(nuE/CLHEP::GeV);
	    if(anganue0<cone[5])runAct->h_G4anue_6->Fill(nuE/CLHEP::GeV);
	    if(anganue0<cone[6])runAct->h_G4anue_7->Fill(nuE/CLHEP::GeV);
	    if(anganue0<cone[7])runAct->h_G4anue_8->Fill(nuE/CLHEP::GeV);
	      
	    if(anganue<cone[0])runAct->h_G4anueOA_1->Fill(nuE/CLHEP::GeV);
	    if(anganue<cone[1])runAct->h_G4anueOA_2->Fill(nuE/CLHEP::GeV);
	    if(anganue<cone[2])runAct->h_G4anueOA_3->Fill(nuE/CLHEP::GeV);
	    if(anganue<cone[3])runAct->h_G4anueOA_4->Fill(nuE/CLHEP::GeV);
	    if(anganue<cone[4])runAct->h_G4anueOA_5->Fill(nuE/CLHEP::GeV);	      
	    if(anganue<cone[5])runAct->h_G4anueOA_6->Fill(nuE/CLHEP::GeV);	      
	    if(anganue<cone[6])runAct->h_G4anueOA_7->Fill(nuE/CLHEP::GeV);	      
	    if(anganue<cone[7])runAct->h_G4anueOA_8->Fill(nuE/CLHEP::GeV);	      
	    G4cout <<"SBSteppingAction: anue_OFFAXIS FIRST ONE "<< nuE/CLHEP::GeV<< G4endl;
	    evtAct->anganue0OLD=anganue0;
	  }
	}//anue
	  	  
	if(PID==-14){
	  anganumu=acos(((eOME[0]*nupx)+(eOME[1]*nupy)+(eOME[2]*nupz))/nuE);	  
	  anganumu0=acos(nupz/nuE);
	    
	  if(anganumu0!=evtAct->anganumu0OLD){// to avoid double counting
	    if(anganumu0<cone[0])runAct->h_G4anumu_1->Fill(nuE/CLHEP::GeV);
	    if(anganumu0<cone[1])runAct->h_G4anumu_2->Fill(nuE/CLHEP::GeV);
	    if(anganumu0<cone[2])runAct->h_G4anumu_3->Fill(nuE/CLHEP::GeV);
	    if(anganumu0<cone[3])runAct->h_G4anumu_4->Fill(nuE/CLHEP::GeV);
	    if(anganumu0<cone[4])runAct->h_G4anumu_5->Fill(nuE/CLHEP::GeV);
	    if(anganumu0<cone[5])runAct->h_G4anumu_6->Fill(nuE/CLHEP::GeV);
	    if(anganumu0<cone[6])runAct->h_G4anumu_7->Fill(nuE/CLHEP::GeV);
	    if(anganumu0<cone[7])runAct->h_G4anumu_8->Fill(nuE/CLHEP::GeV);
	      
	    if(anganumu<cone[0])runAct->h_G4anumuOA_1->Fill(nuE/CLHEP::GeV);
	    if(anganumu<cone[1])runAct->h_G4anumuOA_2->Fill(nuE/CLHEP::GeV);
	    if(anganumu<cone[2])runAct->h_G4anumuOA_3->Fill(nuE/CLHEP::GeV);
	    if(anganumu<cone[3])runAct->h_G4anumuOA_4->Fill(nuE/CLHEP::GeV);
	    if(anganumu<cone[4])runAct->h_G4anumuOA_5->Fill(nuE/CLHEP::GeV);		
	    if(anganumu<cone[5])runAct->h_G4anumuOA_6->Fill(nuE/CLHEP::GeV);		
	    if(anganumu<cone[6])runAct->h_G4anumuOA_7->Fill(nuE/CLHEP::GeV);		
	    if(anganumu<cone[7])runAct->h_G4anumuOA_8->Fill(nuE/CLHEP::GeV);		
	    G4cout <<"SBSteppingAction: anumu_OFFAXIS FIRST ONE "<< nuE/CLHEP::GeV<< G4endl;
	    evtAct->anganumu0OLD=anganumu0;
	  }	    
	}//anumu  

      }//E!=0
*/

    }//nuecheck
  }// enablea

  //***************************************************
  //***************************************************
  //***************************************************

  if(theProcessName=="Decay"){
    //but then ask additionally for:
    G4StepStatus stepStatus = point2->GetStepStatus();

    //pProcess->DumpInfo();
    
    G4cout << 
      "SBSteppingAction: Decay of PARTICLE "<<theParticleName<<
      " PDGcode: "<<PID<<
      " track ID: "<<track->GetTrackID()<<
      " parent ID: "<< PARID <<
      " step status "<<stepStatus<< G4endl;

    if (stepStatus!=fAtRestDoItProc) { // decay in flight;
      //if (stepStatus==fAtRestDoItProc) { // at rest

      //if(DEB)G4cout << "SBSteppingAction: Decay IN FLIGHT of PARTICLE " << PID <<" "<<theParticleName<<" track ID: "<<track->GetTrackID()<<". parent ID: "<< track->GetParentID()<<". Step status: "<<stepStatus<<" "<<G4endl;

      if(PID==211){
	//if(DEB)
	G4cout << "SBSteppingAction: PI+ Decay IN FLIGHT " << PID <<" "<<theParticleName<<" "<<track->GetTrackID()<< G4endl;
	evtAct->PIPLUS_TRACKID = track->GetTrackID();
	evtAct->x_DIF_PI_piplus=track->GetPosition();
	evtAct->p_DIF_PI_piplus=track->GetMomentum();
	evtAct->NDIF_piplus++;

	//G4cout << "PI+ PARENT IDENTITY "<<parPDG<<" from a "<<parparPDG<<G4endl;
	G4cout << "PI+ <~ "<<parName<<" <~ "<<parparName<<G4endl;

	if(kaonfath||kaonGfath||kaonGGfath){// fails if pion does twice PionMinusInelastic ... 
	  evtAct->ParLevel_piplus = 1;
	  G4cout << "K+/- => PI+ "<<G4endl;
	}else if(k0fath||k0Gfath||k0GGfath){
	  evtAct->ParLevel_piplus = 2;
	  G4cout << "K0 => PI+ "<<G4endl;
	}else{
	  evtAct->ParLevel_piplus = 0;
	  G4cout << "direct PI+"<<G4endl;
	}

	/*
	  if(PARID==evtAct->KPLUS_TRACKID){
	  evtAct->ParLevel_piplus = 1;
	  G4cout << "SBSteppingAction: pi+ from a K+ tagged. pi+ ParLevel "<<evtAct->ParLevel_piplus<<G4endl;
	  }
	  if(PARID==evtAct->KMINUS_TRACKID){
	  evtAct->ParLevel_piplus = 1;
	  G4cout << "SBSteppingAction: pi+ from a K- tagged. pi+ ParLevel "<<evtAct->ParLevel_piplus<<G4endl;
	  }
	  if(PARID==evtAct->KZEROS_TRACKID){
	  evtAct->ParLevel_piplus = 1;
	  G4cout << "SBSteppingAction: pi+ from a K0S tagged. pi+ ParLevel "<<evtAct->ParLevel_piplus<<G4endl;
	  }
	  if(PARID==evtAct->KZEROL_TRACKID){
	  evtAct->ParLevel_piplus = 1;
	  G4cout << "SBSteppingAction: pi+ from a K0L tagged. pi+ ParLevel "<<evtAct->ParLevel_piplus<<G4endl;
	  }
	*/
      } else if(PID==-211){
	//if(DEB)
	G4cout << "SBSteppingAction: PI- Decay IN FLIGHT " << PID <<" "<<theParticleName<<" "<<track->GetTrackID()<< G4endl;
	evtAct->PIMINUS_TRACKID = track->GetTrackID();
	evtAct->x_DIF_PI_piminus=track->GetPosition();
	evtAct->p_DIF_PI_piminus=track->GetMomentum();
	evtAct->NDIF_piminus++;

	//G4cout << "PI- PARENT IDENTITY "<<parPDG <<" from a "<<parparPDG<<G4endl;
	G4cout << "PI- <~ "<<parName<<" <~ "<<parparName<<G4endl;

	if(kaonfath||kaonGfath||kaonGGfath){// fails if pion does twice PionMinusInelastic ... 
	  evtAct->ParLevel_piminus = 1;
	  G4cout << "K+/- => PI- "<<G4endl;
	}else if(k0fath||k0Gfath||k0GGfath){
	  evtAct->ParLevel_piminus = 2;
	  G4cout << "K0 => PI- "<<G4endl;
	}else{
	  evtAct->ParLevel_piminus = 0;
	  G4cout << "direct PI-"<<G4endl;
	}

	/*
	  if(PARID==evtAct->KMINUS_TRACKID){
	  evtAct->ParLevel_piminus = 1;
	  G4cout << "SBSteppingAction: pi- from a K- tagged. pi- ParLevel " <<evtAct->ParLevel_piminus<<G4endl;
	  }
	  if(PARID==evtAct->KPLUS_TRACKID){
	  evtAct->ParLevel_piminus = 1;
	  G4cout << "SBSteppingAction: pi- from a K+ tagged. pi- ParLevel " <<evtAct->ParLevel_piminus<<G4endl;
	  }
	  if(PARID==evtAct->KZEROS_TRACKID){
	  evtAct->ParLevel_piminus = 1;
	  G4cout << "SBSteppingAction: pi- from a K0S tagged. pi- ParLevel " <<evtAct->ParLevel_piminus<<G4endl;
	  }
	  if(PARID==evtAct->KZEROL_TRACKID){
	  evtAct->ParLevel_piminus = 1;
	  G4cout << "SBSteppingAction: pi- from a K0L tagged. pi- ParLevel " <<evtAct->ParLevel_piminus<<G4endl;
	  }
	*/
      } else if(PID==321){
	//if(DEB)
	G4cout << "SBSteppingAction: K+ Decay IN FLIGHT " << PID <<" "<<theParticleName<<" "<<track->GetTrackID()<< G4endl;
	evtAct->KPLUS_TRACKID = track->GetTrackID();
	evtAct->x_DIF_K_kplus=track->GetPosition();
	evtAct->p_DIF_K_kplus=track->GetMomentum();
	evtAct->NDIF_Kplus++;
      } else if(PID==-321){
	//if(DEB)
	G4cout << "SBSteppingAction: K- Decay IN FLIGHT " << PID <<" "<<theParticleName<<" "<<track->GetTrackID()<< G4endl;
	evtAct->KMINUS_TRACKID = track->GetTrackID();
	evtAct->x_DIF_K_kminus=track->GetPosition();
	evtAct->p_DIF_K_kminus=track->GetMomentum();
	evtAct->NDIF_Kminus++;
      } else if(PID==310){
	//if(DEB)
	G4cout << "SBSteppingAction: K0S Decay IN FLIGHT " << PID <<" "<<theParticleName<<" "<<track->GetTrackID()<< G4endl;
	evtAct->KZEROS_TRACKID = track->GetTrackID();
	evtAct->x_DIF_K_kzeroS=track->GetPosition();
	evtAct->p_DIF_K_kzeroS=track->GetMomentum();
	evtAct->NDIF_KzeroS++;
      } else if(PID==130){
	//if(DEB)
	G4cout << "SBSteppingAction: K0L Decay IN FLIGHT " << PID <<" "<<theParticleName<<" "<<track->GetTrackID()<< G4endl;
	evtAct->KZEROL_TRACKID = track->GetTrackID();
	evtAct->x_DIF_K_kzeroL=track->GetPosition();
	evtAct->p_DIF_K_kzeroL=track->GetMomentum();
	evtAct->NDIF_KzeroL++;


      } else if(PID==13){
	//if(DEB)
	G4cout << "SBSteppingAction: MU- Decay IN FLIGHT " << PID <<" "<<theParticleName<<" "<<track->GetTrackID()<<" x= " << track->GetPosition()<<" p= "<<track->GetMomentum()<<G4endl;

      } else if(PID==-13){
	//if(DEB)
	G4cout << "SBSteppingAction: MU+ Decay IN FLIGHT " << PID <<" "<<theParticleName<<" "<<track->GetTrackID()<<" x= " << track->GetPosition()<<" p= "<<track->GetMomentum()<<G4endl;

      }

    }//DIF
    
    //888888888888888888888888888
    //
    // check branching ratios
    //    
    if(theParticleName=="pi+"){runAct->h_decay_NUM->Fill(30.);}
    if(theParticleName=="pi-"){runAct->h_decay_NUM->Fill(31.);}
    if(theParticleName=="kaon+"){
      runAct->h_decay_NUM->Fill(1.);
      runAct->h_decay_NUM->Fill(3.);
      runAct->h_decay_NUM->Fill(5.);
      runAct->h_decay_NUM->Fill(7.);
      runAct->h_decay_NUM->Fill(9.);
      runAct->h_decay_NUM->Fill(11.);
    }
    if(theParticleName=="kaon-"){
      runAct->h_decay_NUM->Fill(2.);
      runAct->h_decay_NUM->Fill(4.);
      runAct->h_decay_NUM->Fill(6.);
      runAct->h_decay_NUM->Fill(8.);
      runAct->h_decay_NUM->Fill(10.);
      runAct->h_decay_NUM->Fill(12.);
    }
    if(theParticleName=="kaon0L"){
      runAct->h_decay_NUM->Fill(13.);
      runAct->h_decay_NUM->Fill(14.);
      runAct->h_decay_NUM->Fill(15.);
      runAct->h_decay_NUM->Fill(16.);
      runAct->h_decay_NUM->Fill(17.);
      runAct->h_decay_NUM->Fill(18.);
    }
    if(theParticleName=="kaon0S"){
      runAct->h_decay_NUM->Fill(19.);
      runAct->h_decay_NUM->Fill(20.);
    }
    //888888888888888888888888888
    //
    // Detect the final state of the decay and categorize
    //

    const G4TrackVector* trkList = aStep->GetSecondary();
    G4TrackVector::const_iterator ite;
    G4ThreeVector p_nutrack=G4ThreeVector(0,0,0);
    G4ThreeVector p_parent=G4ThreeVector(0,0,0);
    G4ThreeVector p_sec_mu=G4ThreeVector(0,0,0);
    G4ThreeVector x_sec_mu=G4ThreeVector(0,0,0);
    G4int PAID = 0, ID=0; 
    G4int id_sec_pip1 = 0; 
    G4int id_sec_pim1 = 0; 
    G4int id_sec_pip2 = 0; 
    G4int id_sec_pim2 = 0; 
    //G4String name_parent[5]={""};
    bool firstpip=true;
    bool firstpim=true;
    G4int Nsec = 0;
    G4String SecNames[5]={""};
    for( ite = trkList->begin(); ite != trkList->end(); ite++){
      /*if(DEB)*/
      G4String CrePro =(*ite)->GetCreatorProcess()->GetProcessName();
      G4String SecParName=(*ite)->GetDefinition()->GetParticleName(); 
      if(CrePro=="Decay"){  
	PAID = track->GetParentID(); 
	ID = track->GetTrackID(); 
	p_parent=track->GetMomentum();

	G4cout <<"SBSteppingAction: DIF SECONDARY "<<SecParName<<" created by "<<CrePro<<G4endl;
	G4cout <<"SBSteppingAction: PARENT ID "<<PAID<<" TRACK ID "<<ID<<G4endl;

	//name_parent=track->GetDefinition()->GetParticleName(); 

	//G4cout<<"PROVA "<<(G4int)(*ite)->GetTrackID()<<" "<<(*ite)->GetParentID()<<G4endl;

	if((SecParName=="nu_e")||(SecParName=="anti_nu_e")||
	   (SecParName=="nu_mu")||(SecParName=="anti_nu_mu"))
	  p_nutrack = (*ite)->GetMomentum();
	
	if((SecParName=="mu+")||(SecParName=="mu-")){
	  p_sec_mu = (*ite)->GetMomentum();
	  x_sec_mu = (*ite)->GetPosition();
	}

	if((SecParName=="pi+")&&firstpip){
	  id_sec_pip1 = (*ite)->GetTrackID();
	  G4cout << " ID first decay pi+ "<<id_sec_pip1<< G4endl;
	  firstpip=false;
	}
	if((SecParName=="pi-")&&firstpim){
	  id_sec_pim1 = (*ite)->GetTrackID();
	  G4cout << " ID first decay pi- "<<id_sec_pim1<< G4endl;
	  firstpim=false;
	}
	if((SecParName=="pi+")&&(!firstpip)){
	  id_sec_pip2 = (*ite)->GetTrackID();
	  G4cout << " ID second decay pi- "<<id_sec_pip2<< G4endl;
	}
	if((SecParName=="pi-")&&(!firstpim)){
	  id_sec_pim2 = (*ite)->GetTrackID();
	  G4cout << " ID second decay pi- "<<id_sec_pim2<< G4endl;
	}
	
	SecNames[Nsec]=SecParName;
	Nsec++;
      }
    }
    //if((*ite)->GetDefinition()==G4Electron::Electron()) 
    G4cout <<"SBSteppingAction: "<< Nsec << " secondaries found."<<G4endl;
    
    evtAct->NBody=Nsec;
    
    G4int N_eplus=0,N_eminus=0,N_muplus=0,N_muminus=0,N_nue=0,N_anue=0,N_numu=0,N_anumu=0,N_piplus=0,N_piminus=0,N_pi0=0;
    for(int i=0;i<Nsec;i++){
      if(SecNames[i]=="e+")N_eplus++;
      if(SecNames[i]=="e-")N_eminus++;
      if(SecNames[i]=="mu+")N_muplus++;
      if(SecNames[i]=="mu-")N_muminus++;
      if(SecNames[i]=="nu_e")N_nue++;
      if(SecNames[i]=="anti_nu_e")N_anue++;
      if(SecNames[i]=="nu_mu")N_numu++;
      if(SecNames[i]=="anti_nu_mu")N_anumu++;
      if(SecNames[i]=="pi+")N_piplus++;
      if(SecNames[i]=="pi-")N_piminus++;
      if(SecNames[i]=="pi0")N_pi0++;
    }
    
    G4double parP=0;
    
    if(Nsec==2){
      if(theParticleName=="pi+"){
	if((N_muplus==1)&&(N_numu==1)){
	  evtAct->decflag=30;
	  evtAct->p_nu30[evtAct->n30]=p_nutrack;
	  evtAct->p_par30[evtAct->n30]=p_parent;
	  evtAct->id_par30[evtAct->n30]=PAID;
	  evtAct->id30[evtAct->n30]=ID;
	  evtAct->p_mu30[evtAct->n30]=p_sec_mu;
	  evtAct->x_mu30[evtAct->n30]=x_sec_mu;
	  //runAct->h_prim_match->Fill(evtAct->p_par30[evtAct->n30].mag()/CLHEP::GeV,pTOTprim/CLHEP::GeV);
	  parP = evtAct->p_par30[evtAct->n30].mag();
	  runAct->h_prim_match_pi->Fill((parP-pTOTprim)/CLHEP::GeV,parP/pTOTprim);
	  runAct->h_decay_flags->Fill(evtAct->decflag);

	  //G4cout << "birbo "<<primaryGen->particleName<<" "<<evtAct->class_par30[evtAct->n30]<<" "<<p_parent<<G4endl;
	  // works only for input from file ... :(
	  /*
	    if(primaryGen->particleName=="pi+")evtAct->class_par30[evtAct->n30]=0;
	    if(primaryGen->particleName=="pi-")evtAct->class_par30[evtAct->n30]=0;
	    if(primaryGen->particleName=="kaon+")evtAct->class_par30[evtAct->n30]=1;
	    if(primaryGen->particleName=="kaon-")evtAct->class_par30[evtAct->n30]=1;
	    if(primaryGen->particleName=="kaon0")evtAct->class_par30[evtAct->n30]=2;
	    if(primaryGen->particleName=="anti_kaon0")evtAct->class_par30[evtAct->n30]=2;
	    if(primaryGen->particleName=="kaon0S")evtAct->class_par30[evtAct->n30]=2;
	    if(primaryGen->particleName=="kaon0L")evtAct->class_par30[evtAct->n30]=2;
	  */
	    
	  // check if primary pion or from K charged or from K neutral

	  G4int mypaid=evtAct->id_par30[evtAct->n30];
	  for(int k=0;k<evtAct->n3;k++){//K+ -> pi+ pi0
	    if (mypaid==evtAct->id3[k])evtAct->class_par30[evtAct->n30]=1;
	  }
	  for(int k=0;k<evtAct->n4;k++){//K- -> pi- pi0
	    if (mypaid==evtAct->id4[k])evtAct->class_par30[evtAct->n30]=1;
	  }
	  for(int k=0;k<evtAct->n5;k++){//K+ -> pi+ pi+ pi-
	    if (mypaid==evtAct->id5[k])evtAct->class_par30[evtAct->n30]=1;
	  }
	  for(int k=0;k<evtAct->n6;k++){//K- -> pi- pi- pi+
	    if (mypaid==evtAct->id6[k])evtAct->class_par30[evtAct->n30]=1;
	  }
	  for(int k=0;k<evtAct->n11;k++){//K+ -> pi+ pi0 pi0
	    if (mypaid==evtAct->id11[k])evtAct->class_par30[evtAct->n30]=1;
	  }
	  for(int k=0;k<evtAct->n12;k++){//K- -> pi- pi0 pi0
	    if (mypaid==evtAct->id12[k])evtAct->class_par30[evtAct->n30]=1;
	  }
	  for(int k=0;k<evtAct->n13;k++){//K0L -> e+ nue pi-
	    if (mypaid==evtAct->id13[k])evtAct->class_par30[evtAct->n30]=2;
	  }
	  for(int k=0;k<evtAct->n14;k++){//K0L -> e- anue pi+
	    if (mypaid==evtAct->id14[k])evtAct->class_par30[evtAct->n30]=2;
	  }
	  for(int k=0;k<evtAct->n15;k++){//K0L -> mu+ numu pi-
	    if (mypaid==evtAct->id15[k])evtAct->class_par30[evtAct->n30]=2;
	  }
	  for(int k=0;k<evtAct->n16;k++){//K0L -> mu- anumu pi+
	    if (mypaid==evtAct->id16[k])evtAct->class_par30[evtAct->n30]=2;
	  }
	  for(int k=0;k<evtAct->n18;k++){//K0L -> pi+ pi- pi0
	    if (mypaid==evtAct->id18[k])evtAct->class_par30[evtAct->n30]=2;
	  }
	  for(int k=0;k<evtAct->n19;k++){//K0S -> pi+ pi-
	    if (mypaid==evtAct->id19[k])evtAct->class_par30[evtAct->n30]=2;
	  }

	  /*
	    G4int mypaid=evtAct->id_par30[evtAct->n30];

	    for(int k=0;k<evtAct->n3;k++){//K+ -> pi+ pi0
	    if (mypaid==evtAct->id_par3[k])evtAct->class_par30[evtAct->n30]=1;
	    }

	    for(int k=0;k<evtAct->n3;k++){//K+ -> pi+ pi0
	    if (mypaid==evtAct->id_sec_pip1_3[k])evtAct->class_par30[evtAct->n30]=1;
	    }
	    for(int k=0;k<evtAct->n4;k++){//K- -> pi- pi0
	    if (mypaid==evtAct->id_sec_pim1_4[k])evtAct->class_par30[evtAct->n30]=1;
	    }
	    for(int k=0;k<evtAct->n5;k++){//K+ -> pi+ pi+ pi-
	    if (mypaid==evtAct->id_sec_pip1_5[k])evtAct->class_par30[evtAct->n30]=1;
	    if (mypaid==evtAct->id_sec_pip2_5[k])evtAct->class_par30[evtAct->n30]=1;
	    if (mypaid==evtAct->id_sec_pim1_5[k])evtAct->class_par30[evtAct->n30]=1;
	    }
	    for(int k=0;k<evtAct->n6;k++){//K- -> pi- pi- pi+
	    if (mypaid==evtAct->id_sec_pim1_6[k])evtAct->class_par30[evtAct->n30]=1;
	    if (mypaid==evtAct->id_sec_pim2_6[k])evtAct->class_par30[evtAct->n30]=1;
	    if (mypaid==evtAct->id_sec_pip1_6[k])evtAct->class_par30[evtAct->n30]=1;
	    }
	    for(int k=0;k<evtAct->n11;k++){//K+ -> pi+ pi0 pi0
	    if (mypaid==evtAct->id_sec_pip1_11[k])evtAct->class_par30[evtAct->n30]=1;
	    }
	    for(int k=0;k<evtAct->n12;k++){//K- -> pi- pi0 pi0
	    if (mypaid==evtAct->id_sec_pim1_12[k])evtAct->class_par30[evtAct->n30]=1;
	    }
	    for(int k=0;k<evtAct->n13;k++){//K0L -> e+ nue pi-
	    if (mypaid==evtAct->id_sec_pim1_13[k])evtAct->class_par30[evtAct->n30]=2;
	    }
	    for(int k=0;k<evtAct->n14;k++){//K0L -> e- anue pi+
	    if (mypaid==evtAct->id_sec_pip1_14[k])evtAct->class_par30[evtAct->n30]=2;
	    }
	    for(int k=0;k<evtAct->n15;k++){//K0L -> mu+ numu pi-
	    if (mypaid==evtAct->id_sec_pim1_15[k])evtAct->class_par30[evtAct->n30]=2;
	    }
	    for(int k=0;k<evtAct->n16;k++){//K0L -> mu- anumu pi+
	    if (mypaid==evtAct->id_sec_pip1_16[k])evtAct->class_par30[evtAct->n30]=2;
	    }
	    for(int k=0;k<evtAct->n18;k++){//K0L -> pi+ pi- pi0
	    if (mypaid==evtAct->id_sec_pip1_18[k])evtAct->class_par30[evtAct->n30]=2;
	    if (mypaid==evtAct->id_sec_pim1_18[k])evtAct->class_par30[evtAct->n30]=2;
	    }
	    for(int k=0;k<evtAct->n19;k++){//K0S -> pi+ pi-
	    if (mypaid==evtAct->id_sec_pip1_19[k])evtAct->class_par30[evtAct->n30]=2;
	    if (mypaid==evtAct->id_sec_pim1_19[k])evtAct->class_par30[evtAct->n30]=2;
	    }

	  */

	  evtAct->n30++;
	}

	if((N_eplus==1)&&(N_nue==1)){G4cout<<"decadimento pi+ -> e+ nue"<<G4endl;}

      }
      if(theParticleName=="pi-"){
	if((N_muminus==1)&&(N_anumu==1)){
	  evtAct->decflag=31;
	  evtAct->p_nu31[evtAct->n31]=p_nutrack;
	  evtAct->p_par31[evtAct->n31]=p_parent;
	  evtAct->id_par31[evtAct->n31]=PAID;
	  evtAct->id30[evtAct->n31]=ID;
	  evtAct->p_mu31[evtAct->n31]=p_sec_mu;
	  evtAct->x_mu31[evtAct->n31]=x_sec_mu;
	  //runAct->h_prim_match->Fill(evtAct->p_par31[evtAct->n31].mag()/CLHEP::GeV,pTOTprim/CLHEP::GeV);
	  parP = evtAct->p_par31[evtAct->n31].mag();
	  runAct->h_prim_match_pi->Fill((parP-pTOTprim)/CLHEP::GeV,parP/pTOTprim);
	  runAct->h_decay_flags->Fill(evtAct->decflag);
	  /*
	    if(primaryGen->particleName=="pi+")evtAct->class_par31[evtAct->n31]=0;
	    if(primaryGen->particleName=="pi-")evtAct->class_par31[evtAct->n31]=0;
	    if(primaryGen->particleName=="kaon+")evtAct->class_par31[evtAct->n31]=1;
	    if(primaryGen->particleName=="kaon-")evtAct->class_par31[evtAct->n31]=1;
	    if(primaryGen->particleName=="kaon0")evtAct->class_par31[evtAct->n31]=2;
	    if(primaryGen->particleName=="anti_kaon0")evtAct->class_par31[evtAct->n31]=2;
	    if(primaryGen->particleName=="kaon0S")evtAct->class_par31[evtAct->n31]=2;
	    if(primaryGen->particleName=="kaon0L")evtAct->class_par31[evtAct->n31]=2;
	  */

	  // check if primary pion or from K charged or from K neutral

	  G4int mypaid=evtAct->id_par31[evtAct->n31];
	  for(int k=0;k<evtAct->n3;k++){//K+ -> pi+ pi0
	    if (mypaid==evtAct->id3[k])evtAct->class_par31[evtAct->n31]=1;
	  }
	  for(int k=0;k<evtAct->n4;k++){//K- -> pi- pi0
	    if (mypaid==evtAct->id4[k])evtAct->class_par31[evtAct->n31]=1;
	  }
	  for(int k=0;k<evtAct->n5;k++){//K+ -> pi+ pi+ pi-
	    if (mypaid==evtAct->id5[k])evtAct->class_par31[evtAct->n31]=1;
	  }
	  for(int k=0;k<evtAct->n6;k++){//K- -> pi- pi- pi+
	    if (mypaid==evtAct->id6[k])evtAct->class_par31[evtAct->n31]=1;
	  }
	  for(int k=0;k<evtAct->n11;k++){//K+ -> pi+ pi0 pi0
	    if (mypaid==evtAct->id11[k])evtAct->class_par31[evtAct->n31]=1;
	  }
	  for(int k=0;k<evtAct->n12;k++){//K- -> pi- pi0 pi0
	    if (mypaid==evtAct->id12[k])evtAct->class_par31[evtAct->n31]=1;
	  }
	  for(int k=0;k<evtAct->n13;k++){//K0L -> e+ nue pi-
	    if (mypaid==evtAct->id13[k])evtAct->class_par31[evtAct->n31]=2;
	  }
	  for(int k=0;k<evtAct->n14;k++){//K0L -> e- anue pi+
	    if (mypaid==evtAct->id14[k])evtAct->class_par31[evtAct->n31]=2;
	  }
	  for(int k=0;k<evtAct->n15;k++){//K0L -> mu+ numu pi-
	    if (mypaid==evtAct->id15[k])evtAct->class_par31[evtAct->n31]=2;
	  }
	  for(int k=0;k<evtAct->n16;k++){//K0L -> mu- anumu pi+
	    if (mypaid==evtAct->id16[k])evtAct->class_par31[evtAct->n31]=2;
	  }
	  for(int k=0;k<evtAct->n18;k++){//K0L -> pi+ pi- pi0
	    if (mypaid==evtAct->id18[k])evtAct->class_par31[evtAct->n31]=2;
	  }
	  for(int k=0;k<evtAct->n19;k++){//K0S -> pi+ pi-
	    if (mypaid==evtAct->id19[k])evtAct->class_par31[evtAct->n31]=2;
	  }

	  /*
	    G4int mypaid=evtAct->id_par31[evtAct->n31];
	    for(int k=0;k<evtAct->n3;k++){//K+ -> pi+ pi0
	    if (mypaid==evtAct->id_sec_pip1_3[k])evtAct->class_par31[evtAct->n31]=1;
	    }
	    for(int k=0;k<evtAct->n4;k++){//K- -> pi- pi0
	    if (mypaid==evtAct->id_sec_pim1_4[k])evtAct->class_par31[evtAct->n31]=1;
	    }
	    for(int k=0;k<evtAct->n5;k++){//K+ -> pi+ pi+ pi-
	    if (mypaid==evtAct->id_sec_pip1_5[k])evtAct->class_par31[evtAct->n31]=1;
	    if (mypaid==evtAct->id_sec_pip2_5[k])evtAct->class_par31[evtAct->n31]=1;
	    if (mypaid==evtAct->id_sec_pim1_5[k])evtAct->class_par31[evtAct->n31]=1;
	    }
	    for(int k=0;k<evtAct->n6;k++){//K- -> pi- pi- pi+
	    if (mypaid==evtAct->id_sec_pim1_6[k])evtAct->class_par31[evtAct->n31]=1;
	    if (mypaid==evtAct->id_sec_pim2_6[k])evtAct->class_par31[evtAct->n31]=1;
	    if (mypaid==evtAct->id_sec_pip1_6[k])evtAct->class_par31[evtAct->n31]=1;
	    }
	    for(int k=0;k<evtAct->n11;k++){//K+ -> pi+ pi0 pi0
	    if (mypaid==evtAct->id_sec_pip1_11[k])evtAct->class_par31[evtAct->n31]=1;
	    }
	    for(int k=0;k<evtAct->n12;k++){//K- -> pi- pi0 pi0
	    if (mypaid==evtAct->id_sec_pim1_12[k])evtAct->class_par31[evtAct->n31]=1;
	    }
	    for(int k=0;k<evtAct->n13;k++){//K0L -> e+ nue pi-
	    if (mypaid==evtAct->id_sec_pim1_13[k])evtAct->class_par31[evtAct->n31]=2;
	    }
	    for(int k=0;k<evtAct->n14;k++){//K0L -> e- anue pi+
	    if (mypaid==evtAct->id_sec_pip1_14[k])evtAct->class_par31[evtAct->n31]=2;
	    }
	    for(int k=0;k<evtAct->n15;k++){//K0L -> mu+ numu pi-
	    if (mypaid==evtAct->id_sec_pim1_15[k])evtAct->class_par31[evtAct->n31]=2;
	    }
	    for(int k=0;k<evtAct->n16;k++){//K0L -> mu- anumu pi+
	    if (mypaid==evtAct->id_sec_pip1_16[k])evtAct->class_par31[evtAct->n31]=2;
	    }
	    for(int k=0;k<evtAct->n18;k++){//K0L -> pi+ pi- pi0
	    if (mypaid==evtAct->id_sec_pip1_18[k])evtAct->class_par31[evtAct->n31]=2;
	    if (mypaid==evtAct->id_sec_pim1_18[k])evtAct->class_par31[evtAct->n31]=2;
	    }
	    for(int k=0;k<evtAct->n19;k++){//K0S -> pi+ pi-
	    if (mypaid==evtAct->id_sec_pip1_19[k])evtAct->class_par31[evtAct->n31]=2;
	    if (mypaid==evtAct->id_sec_pim1_19[k])evtAct->class_par31[evtAct->n31]=2;
	    }
	  */

	  evtAct->n31++;
	}

	if((N_eminus==1)&&(N_anue==1)){G4cout<<"decadimento pi- -> e- antinue"<<G4endl;}

      }
      if(theParticleName=="kaon+"){
	if((N_muplus==1)&&(N_numu==1)){
	  evtAct->decflag=1;
	  evtAct->p_nu1[evtAct->n1]=p_nutrack;
	  evtAct->p_par1[evtAct->n1]=p_parent;
	  evtAct->id_par1[evtAct->n1]=PAID;
	  evtAct->id1[evtAct->n1]=ID;
	  evtAct->p_mu1[evtAct->n1]=p_sec_mu;
	  evtAct->x_mu1[evtAct->n1]=x_sec_mu;
	  //runAct->h_prim_match->Fill(evtAct->p_par1[evtAct->n1].mag()/CLHEP::GeV,pTOTprim/CLHEP::GeV);
	  parP = evtAct->p_par1[evtAct->n1].mag();
	  runAct->h_prim_match_ka->Fill((parP-pTOTprim)/CLHEP::GeV,parP/pTOTprim);
	  runAct->h_decay_flags->Fill(evtAct->decflag);
	  evtAct->n1++;
	}
	if((N_piplus==1)&&(N_pi0==1)){
	  evtAct->decflag=3;
	  evtAct->p_nu3[evtAct->n3]=p_nutrack;
	  evtAct->p_par3[evtAct->n3]=p_parent;
	  evtAct->id_par3[evtAct->n3]=PAID;
	  evtAct->id3[evtAct->n3]=ID;
	  evtAct->p_mu3[evtAct->n3]=p_sec_mu;
	  evtAct->x_mu3[evtAct->n3]=x_sec_mu;

	  evtAct->id_sec_pip1_3[evtAct->n3]=id_sec_pip1;

	  //runAct->h_prim_match->Fill(evtAct->p_par3[evtAct->n3].mag()/CLHEP::GeV,pTOTprim/CLHEP::GeV);
	  parP = evtAct->p_par3[evtAct->n3].mag();
	  runAct->h_prim_match_ka->Fill((parP-pTOTprim)/CLHEP::GeV,parP/pTOTprim);
	  runAct->h_decay_flags->Fill(evtAct->decflag);
	  evtAct->n3++;
	}
      }
      if(theParticleName=="kaon-"){
	if((N_muminus==1)&&(N_anumu==1)){
	  evtAct->decflag=2;
	  evtAct->p_nu2[evtAct->n2]=p_nutrack;
	  evtAct->p_par2[evtAct->n2]=p_parent;
	  evtAct->id_par2[evtAct->n2]=PAID;
	  evtAct->id2[evtAct->n2]=ID;
	  evtAct->p_mu2[evtAct->n2]=p_sec_mu;
	  evtAct->x_mu2[evtAct->n2]=x_sec_mu;
	  //runAct->h_prim_match->Fill(evtAct->p_par2[evtAct->n2].mag()/CLHEP::GeV,pTOTprim/CLHEP::GeV);
	  parP = evtAct->p_par2[evtAct->n2].mag();
	  runAct->h_prim_match_ka->Fill((parP-pTOTprim)/CLHEP::GeV,parP/pTOTprim);
	  runAct->h_decay_flags->Fill(evtAct->decflag);
	  evtAct->n2++;
	}
	if((N_piminus==1)&&(N_pi0==1)){
	  evtAct->decflag=4;
	  evtAct->p_nu4[evtAct->n4]=p_nutrack;
	  evtAct->p_par4[evtAct->n4]=p_parent;
	  evtAct->id_par4[evtAct->n4]=PAID;
	  evtAct->id4[evtAct->n4]=ID;
	  evtAct->p_mu4[evtAct->n4]=p_sec_mu;
	  evtAct->x_mu4[evtAct->n4]=x_sec_mu;

	  evtAct->id_sec_pim1_4[evtAct->n4]=id_sec_pim1;

	  //runAct->h_prim_match->Fill(evtAct->p_par4[evtAct->n4].mag()/CLHEP::GeV,pTOTprim/CLHEP::GeV);
	  parP = evtAct->p_par4[evtAct->n4].mag();
	  runAct->h_prim_match_ka->Fill((parP-pTOTprim)/CLHEP::GeV,parP/pTOTprim);
	  runAct->h_decay_flags->Fill(evtAct->decflag);
	  evtAct->n4++;
	}
      }
      if(theParticleName=="kaon0S"){
	if((N_piplus==1)&&(N_piminus==1)){
	  evtAct->decflag=19;
	  evtAct->p_nu19[evtAct->n19]=p_nutrack;
	  evtAct->p_par19[evtAct->n19]=p_parent;
	  evtAct->id_par19[evtAct->n19]=PAID;
	  evtAct->id19[evtAct->n19]=ID;
	  evtAct->p_mu19[evtAct->n19]=p_sec_mu;
	  evtAct->x_mu19[evtAct->n19]=x_sec_mu;

	  evtAct->id_sec_pip1_19[evtAct->n19]=id_sec_pip1;
	  evtAct->id_sec_pim1_19[evtAct->n19]=id_sec_pim1;

	  //runAct->h_prim_match->Fill(evtAct->p_par19[evtAct->n19].mag()/CLHEP::GeV,pTOTprim/CLHEP::GeV);
	  parP = evtAct->p_par19[evtAct->n19].mag();
	  runAct->h_prim_match_k0->Fill((parP-pTOTprim)/CLHEP::GeV,parP/pTOTprim);
	  runAct->h_decay_flags->Fill(evtAct->decflag);
	  evtAct->n19++;
	}
	if(N_pi0==2){
	  evtAct->decflag=20;
	  evtAct->p_nu20[evtAct->n20]=p_nutrack;
	  evtAct->p_par20[evtAct->n20]=p_parent;
	  evtAct->id_par20[evtAct->n20]=PAID;
	  evtAct->id20[evtAct->n20]=ID;
	  evtAct->p_mu20[evtAct->n20]=p_sec_mu;
	  evtAct->x_mu20[evtAct->n20]=x_sec_mu;
	  //runAct->h_prim_match->Fill(evtAct->p_par20[evtAct->n20].mag()/CLHEP::GeV,pTOTprim/CLHEP::GeV);
	  parP = evtAct->p_par20[evtAct->n20].mag();
	  runAct->h_prim_match_k0->Fill((parP-pTOTprim)/CLHEP::GeV,parP/pTOTprim);
	  runAct->h_decay_flags->Fill(evtAct->decflag);
	  evtAct->n20++;
	}
      }
    }//Nsec 2
    
    if(Nsec==3){
      if(theParticleName=="kaon+"){
	if((N_piplus==2)&&(N_piminus==1)){
	  evtAct->decflag=5;
	  evtAct->p_nu5[evtAct->n5]=p_nutrack;
	  evtAct->p_par5[evtAct->n5]=p_parent;
	  evtAct->id_par5[evtAct->n5]=PAID;
	  evtAct->id5[evtAct->n5]=ID;
	  evtAct->p_mu5[evtAct->n5]=p_sec_mu;
	  evtAct->x_mu5[evtAct->n5]=x_sec_mu;

	  evtAct->id_sec_pip1_5[evtAct->n5]=id_sec_pip1;
	  evtAct->id_sec_pip2_5[evtAct->n5]=id_sec_pip2;
	  evtAct->id_sec_pim1_5[evtAct->n5]=id_sec_pim1;

	  //runAct->h_prim_match->Fill(evtAct->p_par5[evtAct->n5].mag()/CLHEP::GeV,pTOTprim/CLHEP::GeV);
	  parP = evtAct->p_par5[evtAct->n5].mag();
	  runAct->h_prim_match_ka->Fill((parP-pTOTprim)/CLHEP::GeV,parP/pTOTprim);
	  runAct->h_decay_flags->Fill(evtAct->decflag);
	  evtAct->n5++;
	}
	if((N_eplus==1)&&(N_nue==1)&&(N_pi0==1)){
	  evtAct->decflag=7;
	  evtAct->p_nu7[evtAct->n7]=p_nutrack;
	  evtAct->p_par7[evtAct->n7]=p_parent;
	  evtAct->id_par7[evtAct->n7]=PAID;
	  evtAct->id7[evtAct->n7]=ID;
	  evtAct->p_mu7[evtAct->n7]=p_sec_mu;
	  evtAct->x_mu7[evtAct->n7]=x_sec_mu;
	  //runAct->h_prim_match->Fill(evtAct->p_par7[evtAct->n7].mag()/CLHEP::GeV,pTOTprim/CLHEP::GeV);
	  parP = evtAct->p_par7[evtAct->n7].mag();
	  runAct->h_prim_match_ka->Fill((parP-pTOTprim)/CLHEP::GeV,parP/pTOTprim);
	  runAct->h_decay_flags->Fill(evtAct->decflag);
	  evtAct->n7++;
	  G4cout <<"SBSteppingAction: K3 "<<p_parent<<" "<<p_nutrack<<" "<<evtAct->decflag<<G4endl;
	}
	if((N_muplus==1)&&(N_numu==1)&&(N_pi0==1)){
	  evtAct->decflag=9;
	  evtAct->p_nu9[evtAct->n9]=p_nutrack;
	  evtAct->p_par9[evtAct->n9]=p_parent;
	  evtAct->id_par9[evtAct->n9]=PAID;
	  evtAct->id9[evtAct->n9]=ID;
	  evtAct->p_mu9[evtAct->n9]=p_sec_mu;
	  evtAct->x_mu9[evtAct->n9]=x_sec_mu;
	  //runAct->h_prim_match->Fill(evtAct->p_par9[evtAct->n9].mag()/CLHEP::GeV,pTOTprim/CLHEP::GeV);
	  parP = evtAct->p_par9[evtAct->n9].mag();
	  runAct->h_prim_match_ka->Fill((parP-pTOTprim)/CLHEP::GeV,parP/pTOTprim);
	  runAct->h_decay_flags->Fill(evtAct->decflag);
	  evtAct->n9++;
	  G4cout <<"SBSteppingAction: K3 "<<p_parent<<" "<<p_nutrack<<" "<<evtAct->decflag<<G4endl;
	}
	if((N_piplus==1)&&(N_pi0==2)){
	  evtAct->decflag=11;
	  evtAct->p_nu11[evtAct->n11]=p_nutrack;
	  evtAct->p_par11[evtAct->n11]=p_parent;
	  evtAct->id_par11[evtAct->n11]=PAID;
	  evtAct->id11[evtAct->n11]=ID;
	  evtAct->p_mu11[evtAct->n11]=p_sec_mu;
	  evtAct->x_mu11[evtAct->n11]=x_sec_mu;

	  evtAct->id_sec_pip1_11[evtAct->n11]=id_sec_pip1;

	  //runAct->h_prim_match->Fill(evtAct->p_par11[evtAct->n11].mag()/CLHEP::GeV,pTOTprim/CLHEP::GeV);
	  parP = evtAct->p_par11[evtAct->n11].mag();
	  runAct->h_prim_match_ka->Fill((parP-pTOTprim)/CLHEP::GeV,parP/pTOTprim);
	  runAct->h_decay_flags->Fill(evtAct->decflag);
	  evtAct->n11++;
	}
      }//k+
      if(theParticleName=="kaon-"){
	if((N_piminus==2)&&(N_piplus==1)){
	  evtAct->decflag=6;
	  evtAct->p_nu6[evtAct->n6]=p_nutrack;
	  evtAct->p_par6[evtAct->n6]=p_parent;
	  evtAct->id_par6[evtAct->n6]=PAID;
	  evtAct->id6[evtAct->n6]=ID;
	  evtAct->p_mu6[evtAct->n6]=p_sec_mu;
	  evtAct->x_mu6[evtAct->n6]=x_sec_mu;

	  evtAct->id_sec_pip1_6[evtAct->n6]=id_sec_pip1;
	  evtAct->id_sec_pim1_6[evtAct->n6]=id_sec_pim1;
	  evtAct->id_sec_pim2_6[evtAct->n6]=id_sec_pim2;

	  //runAct->h_prim_match->Fill(evtAct->p_par6[evtAct->n6].mag()/CLHEP::GeV,pTOTprim/CLHEP::GeV);
	  parP = evtAct->p_par6[evtAct->n6].mag();
	  runAct->h_prim_match_ka->Fill((parP-pTOTprim)/CLHEP::GeV,parP/pTOTprim);
	  runAct->h_decay_flags->Fill(evtAct->decflag);
	  evtAct->n6++;
	}
	if((N_eminus==1)&&(N_anue==1)&&(N_pi0==1)){
	  evtAct->decflag=8;
	  evtAct->p_nu8[evtAct->n8]=p_nutrack;
	  evtAct->p_par8[evtAct->n8]=p_parent;
	  evtAct->id_par8[evtAct->n8]=PAID;
	  evtAct->id8[evtAct->n8]=ID;
	  evtAct->p_mu8[evtAct->n8]=p_sec_mu;
	  evtAct->x_mu8[evtAct->n8]=x_sec_mu;
	  //runAct->h_prim_match->Fill(evtAct->p_par8[evtAct->n8].mag()/CLHEP::GeV,pTOTprim/CLHEP::GeV);
	  parP = evtAct->p_par8[evtAct->n8].mag();
	  runAct->h_prim_match_ka->Fill((parP-pTOTprim)/CLHEP::GeV,parP/pTOTprim);
	  runAct->h_decay_flags->Fill(evtAct->decflag);
	  evtAct->n8++;
	  G4cout <<"SBSteppingAction: K3 "<<p_parent<<" "<<p_nutrack<<" "<<evtAct->decflag<<G4endl;
	}
	if((N_muminus==1)&&(N_anumu==1)&&(N_pi0==1)){
	  evtAct->decflag=10;
	  evtAct->p_nu10[evtAct->n10]=p_nutrack;
	  evtAct->p_par10[evtAct->n10]=p_parent;
	  evtAct->id_par10[evtAct->n10]=PAID;
	  evtAct->id10[evtAct->n10]=ID;
	  evtAct->p_mu10[evtAct->n10]=p_sec_mu;
	  evtAct->x_mu10[evtAct->n10]=x_sec_mu;
	  //runAct->h_prim_match->Fill(evtAct->p_par10[evtAct->n10].mag()/CLHEP::GeV,pTOTprim/CLHEP::GeV);
	  parP = evtAct->p_par10[evtAct->n10].mag();
	  runAct->h_prim_match_ka->Fill((parP-pTOTprim)/CLHEP::GeV,parP/pTOTprim);
	  runAct->h_decay_flags->Fill(evtAct->decflag);
	  evtAct->n10++;
	  G4cout <<"SBSteppingAction: K3 "<<p_parent<<" "<<p_nutrack<<" "<<evtAct->decflag<<G4endl;
	}
	if((N_piminus==1)&&(N_pi0==2)){
	  evtAct->decflag=12;
	  evtAct->p_nu12[evtAct->n12]=p_nutrack;
	  evtAct->p_par12[evtAct->n12]=p_parent;
	  evtAct->id_par12[evtAct->n12]=PAID;
	  evtAct->id12[evtAct->n12]=ID;
	  evtAct->p_mu12[evtAct->n12]=p_sec_mu;
	  evtAct->x_mu12[evtAct->n12]=x_sec_mu;

	  evtAct->id_sec_pim1_12[evtAct->n12]=id_sec_pim1;

	  //runAct->h_prim_match->Fill(evtAct->p_par12[evtAct->n12].mag()/CLHEP::GeV,pTOTprim/CLHEP::GeV);
	  parP = evtAct->p_par12[evtAct->n12].mag();
	  runAct->h_prim_match_ka->Fill((parP-pTOTprim)/CLHEP::GeV,parP/pTOTprim);
	  runAct->h_decay_flags->Fill(evtAct->decflag);
	  evtAct->n12++;
	}
      }//k-
      if(theParticleName=="kaon0L"){
	if((N_eplus==1)&&(N_nue==1)&&(N_piminus==1)){

	  evtAct->decflag=13;
	  evtAct->p_nu13[evtAct->n13]=p_nutrack;
	  evtAct->p_par13[evtAct->n13]=p_parent;
	  evtAct->id_par13[evtAct->n13]=PAID;
	  evtAct->id13[evtAct->n13]=ID;
	  evtAct->p_mu13[evtAct->n13]=p_sec_mu;
	  evtAct->x_mu13[evtAct->n13]=x_sec_mu;

	  evtAct->id_sec_pim1_13[evtAct->n13]=id_sec_pim1;

	  //runAct->h_prim_match->Fill(evtAct->p_par13[evtAct->n13].mag()/CLHEP::GeV,pTOTprim/CLHEP::GeV);
	  parP = evtAct->p_par13[evtAct->n13].mag();
	  runAct->h_prim_match_k0->Fill((parP-pTOTprim)/CLHEP::GeV,parP/pTOTprim);
	  runAct->h_decay_flags->Fill(evtAct->decflag);
	  evtAct->n13++;
	  G4cout <<"SBSteppingAction: K3 "<<p_parent<<" "<<p_nutrack<<" "<<evtAct->decflag<<G4endl;

	}
	if((N_eminus==1)&&(N_anue==1)&&(N_piplus==1)){

	  evtAct->decflag=14;
	  evtAct->p_nu14[evtAct->n14]=p_nutrack;
	  evtAct->p_par14[evtAct->n14]=p_parent;
	  evtAct->id_par14[evtAct->n14]=PAID;
	  evtAct->id14[evtAct->n14]=ID;
	  evtAct->p_mu14[evtAct->n14]=p_sec_mu;
	  evtAct->x_mu14[evtAct->n14]=x_sec_mu;

	  evtAct->id_sec_pip1_14[evtAct->n14]=id_sec_pip1;

	  //runAct->h_prim_match->Fill(evtAct->p_par14[evtAct->n14].mag()/CLHEP::GeV,pTOTprim/CLHEP::GeV);
	  parP = evtAct->p_par14[evtAct->n14].mag();
	  runAct->h_prim_match_k0->Fill((parP-pTOTprim)/CLHEP::GeV,parP/pTOTprim);
	  runAct->h_decay_flags->Fill(evtAct->decflag);
	  evtAct->n14++;
	  G4cout <<"SBSteppingAction: K3 "<<p_parent<<" "<<p_nutrack<<" "<<evtAct->decflag<<G4endl;

	}
	if((N_muplus==1)&&(N_numu==1)&&(N_piminus==1)){
	  evtAct->decflag=15;
	  evtAct->p_nu15[evtAct->n15]=p_nutrack;
	  evtAct->p_par15[evtAct->n15]=p_parent;
	  evtAct->id_par15[evtAct->n15]=PAID;
	  evtAct->id15[evtAct->n15]=ID;
	  evtAct->p_mu15[evtAct->n15]=p_sec_mu;
	  evtAct->x_mu15[evtAct->n15]=x_sec_mu;

	  evtAct->id_sec_pim1_15[evtAct->n15]=id_sec_pim1;

	  //runAct->h_prim_match->Fill(evtAct->p_par15[evtAct->n15].mag()/CLHEP::GeV,pTOTprim/CLHEP::GeV);
	  parP = evtAct->p_par15[evtAct->n15].mag();
	  runAct->h_prim_match_k0->Fill((parP-pTOTprim)/CLHEP::GeV,parP/pTOTprim);
	  runAct->h_decay_flags->Fill(evtAct->decflag);
	  evtAct->n15++;
	  G4cout <<"SBSteppingAction: K3 "<<p_parent<<" "<<p_nutrack<<" "<<evtAct->decflag<<G4endl;
	}
	if((N_muminus==1)&&(N_anumu==1)&&(N_piplus==1)){
	  evtAct->decflag=16;
	  evtAct->p_nu16[evtAct->n16]=p_nutrack;
	  evtAct->p_par16[evtAct->n16]=p_parent;
	  evtAct->id_par16[evtAct->n16]=PAID;
	  evtAct->id16[evtAct->n16]=ID;
	  evtAct->p_mu16[evtAct->n16]=p_sec_mu;
	  evtAct->x_mu16[evtAct->n16]=x_sec_mu;

	  evtAct->id_sec_pip1_16[evtAct->n16]=id_sec_pip1;

	  //runAct->h_prim_match->Fill(evtAct->p_par16[evtAct->n16].mag()/CLHEP::GeV,pTOTprim/CLHEP::GeV);
	  parP = evtAct->p_par16[evtAct->n16].mag();
	  runAct->h_prim_match_k0->Fill((parP-pTOTprim)/CLHEP::GeV,parP/pTOTprim);
	  runAct->h_decay_flags->Fill(evtAct->decflag);
	  evtAct->n16++;
	  G4cout <<"SBSteppingAction: K3 "<<p_parent<<" "<<p_nutrack<<" "<<evtAct->decflag<<G4endl;
	}
	if(N_pi0==3){
	  evtAct->decflag=17;
	  evtAct->p_nu17[evtAct->n17]=p_nutrack;
	  evtAct->p_par17[evtAct->n17]=p_parent;
	  evtAct->id_par17[evtAct->n17]=PAID;
	  evtAct->id17[evtAct->n17]=ID;
	  evtAct->p_mu17[evtAct->n17]=p_sec_mu;
	  evtAct->x_mu17[evtAct->n17]=x_sec_mu;
	  //runAct->h_prim_match->Fill(evtAct->p_par17[evtAct->n17].mag()/CLHEP::GeV,pTOTprim/CLHEP::GeV);
	  parP = evtAct->p_par17[evtAct->n17].mag();
	  runAct->h_prim_match_k0->Fill((parP-pTOTprim)/CLHEP::GeV,parP/pTOTprim);
	  runAct->h_decay_flags->Fill(evtAct->decflag);
	  evtAct->n17++;
	}
	if((N_piplus==1)&&(N_piminus==1)&&(N_pi0==1)){
	  evtAct->decflag=18;
	  evtAct->p_nu18[evtAct->n18]=p_nutrack;
	  evtAct->p_par18[evtAct->n18]=p_parent;
	  evtAct->id_par18[evtAct->n18]=PAID;
	  evtAct->id18[evtAct->n18]=ID;
	  evtAct->p_mu18[evtAct->n18]=p_sec_mu;
	  evtAct->x_mu18[evtAct->n18]=x_sec_mu;

	  evtAct->id_sec_pip1_18[evtAct->n18]=id_sec_pip1;
	  evtAct->id_sec_pim1_18[evtAct->n18]=id_sec_pim1;

	  //runAct->h_prim_match->Fill(evtAct->p_par18[evtAct->n18].mag()/CLHEP::GeV,pTOTprim/CLHEP::GeV);
	  parP = evtAct->p_par18[evtAct->n18].mag();
	  runAct->h_prim_match_k0->Fill((parP-pTOTprim)/CLHEP::GeV,parP/pTOTprim);
	  runAct->h_decay_flags->Fill(evtAct->decflag);
	  evtAct->n18++;
	}
      }//k0L

      if((theParticleName=="kaon+")||
	 (theParticleName=="kaon-")||
	 (theParticleName=="kaon0L")){
	evtAct->p_K_3B=track->GetMomentum();
	evtAct->enuLAB_K_3B=p_nutrack.mag();
      }
      
    }//3 BODY

    if((evtAct->decflag==13)||// 3 body K0 with direct neutrinos
       (evtAct->decflag==14)||
       (evtAct->decflag==15)||
       (evtAct->decflag==16)
       ) evtAct->kzero3B=true;
    
    if((evtAct->decflag==7)||// 3 body K+- with direct neutrinos
       (evtAct->decflag==8)||
       (evtAct->decflag==9)||
       (evtAct->decflag==10)
       ) evtAct->kch3B=true;

    if((evtAct->decflag==1)||// 2 body K+- with direct neutrinos
       (evtAct->decflag==2)
       ) evtAct->kch2B=true;

    if((evtAct->decflag==31)||// 2 body pi+- with direct neutrinos
       (evtAct->decflag==30)
       ) evtAct->pich2B=true;

    if(evtAct->pich2B) G4cout << "SBSteppingAction: CLASS. 2 body pi+/- with direct nus " <<G4endl;
    if(evtAct->kch2B) G4cout << "SBSteppingAction: CLASS. 2 body K+/- with direct nus " <<G4endl;
    if(evtAct->kch3B) G4cout << "SBSteppingAction: CLASS. 3 body K+/- with direct nus " <<G4endl;
    if(evtAct->kzero3B) G4cout << "SBSteppingAction: CLASS. 3 body K0L with direct nus " <<G4endl;

    if(evtAct->decflag!=0){
      G4cout << "SBSteppingAction: Decay flag "<< evtAct->decflag << G4endl;
      //runAct->h_decay_flags->Fill(double(evtAct->decflag));
      /*
	runAct->h_decay_flags->Fill(1,double(evtAct->n1));
	runAct->h_decay_flags->Fill(2,double(evtAct->n2));
	runAct->h_decay_flags->Fill(3,double(evtAct->n3));
	runAct->h_decay_flags->Fill(4,double(evtAct->n4));
	runAct->h_decay_flags->Fill(5,double(evtAct->n5));
	runAct->h_decay_flags->Fill(6,double(evtAct->n6));
	runAct->h_decay_flags->Fill(7,double(evtAct->n7));
	runAct->h_decay_flags->Fill(8,double(evtAct->n8));
	runAct->h_decay_flags->Fill(9,double(evtAct->n9));
	runAct->h_decay_flags->Fill(10,double(evtAct->n10));
	runAct->h_decay_flags->Fill(11,double(evtAct->n11));
	runAct->h_decay_flags->Fill(12,double(evtAct->n12));
	runAct->h_decay_flags->Fill(13,double(evtAct->n13));
	runAct->h_decay_flags->Fill(14,double(evtAct->n14));
	runAct->h_decay_flags->Fill(15,double(evtAct->n15));
	runAct->h_decay_flags->Fill(16,double(evtAct->n16));
	runAct->h_decay_flags->Fill(17,double(evtAct->n17));
	runAct->h_decay_flags->Fill(18,double(evtAct->n18));
	runAct->h_decay_flags->Fill(19,double(evtAct->n19));
	runAct->h_decay_flags->Fill(20,double(evtAct->n20));
	runAct->h_decay_flags->Fill(30,double(evtAct->n30));
	runAct->h_decay_flags->Fill(31,double(evtAct->n31));
      */
    }

    //888888888888888888888888888

  } // decay;
  

    /*********************************/
    // necessary ? skip or keep just as debug
    //  
  if(track->GetCreatorProcess() != 0){
    G4String ParticName = track->GetDefinition()->GetParticleName();
    G4String ProcName = track->GetCreatorProcess()->GetProcessName();
    G4ThreeVector imp = track->GetMomentum();
    //G4ThreeVector pos = track->GetPosition();
    pos = track->GetPosition();
    //if(DEB)G4cout << ParticName << " from " << ProcName << " tagged!" << G4endl;
    G4int PARENTID = track->GetParentID(); 
    if(ProcName=="Decay"){
      if((ParticName==IMUPLUS)){
	//G4cout << "MU+ PARENT IDENTITY "<<evtAct->vPDG_ID[PARENTID]<<G4endl;
	if(PARENTID==evtAct->PIPLUS_TRACKID){//pi+ ~> mu+ 
	  if(DEB)G4cout << "       "<<ParticName << " from " << ProcName << " tagged! parent ID "<<PARENTID<<" "<<evtAct->PIPLUS_TRACKID<<" "<<imp<<" "<<pos<<G4endl;
	  if(evtAct->p_muplus_PI==-99999.){// only at first step
	    if(ParticName=="mu+"){
	      evtAct->x_DIF_PI_muplus = pos; 
	      evtAct->p_DIF_PI_muplus = imp;
	      
	      //G4cout << "MU+ PARENT IDENTITY "<<parPDG<<" from a "<<parparPDG<<G4endl;
	      G4cout << "MU+ <~ "<<parName<<" <~ "<<parparName<<" <~ "<<parparparName<<G4endl;
	      /*if((parPDG==PDGpiplus)||(parPDG==PDGpiminus)){*/
	      if((kaonGfath)||(kaonGGfath)){
		evtAct->ParLevel_muplus = 2;
		G4cout << "cascade   K+/- => PI => MU+ "<<G4endl;
	      }else if((k0Gfath)||(k0GGfath)){
		evtAct->ParLevel_muplus = 3;
		G4cout << "cascade   K0 => PI => MU+ "<<G4endl;
	      }else{
		evtAct->ParLevel_muplus = 1;
		G4cout << "direct   PI => MU+ "<<G4endl;
	      }
	      /*	      }
			      else if(kaonfath){
			      evtAct->ParLevel_muplus = 4;
			      G4cout << "direct    K => MU+ "<<G4endl;
			      }
	      */

	      //ooooooooooooooooooooo
	      /*
		if(PARENTID==evtAct->PIPLUS_TRACKID){
		evtAct->ParLevel_muplus = evtAct->ParLevel_piplus + 1;
		G4cout << "SBSteppingAction: mu+ from a pi+ tagged. mu+ ParLevel: "<< evtAct->ParLevel_muplus <<G4endl;
		}
	      */
	      //ooooooooooooooooooooo
	    }
	    if(DEB)G4cout <<" FIRST mu+ "<<evtAct->x_DIF_PI_muplus<<" "<<evtAct->p_DIF_PI_muplus<< G4endl;
	    evtAct->p_muplus_PI=evtAct->p_DIF_PI_muplus.mag();
	  }
	}
	if(PARENTID==evtAct->KPLUS_TRACKID){//K+ ~> mu+ 
	  if(DEB)G4cout << "       "<<ParticName << " from " << ProcName << " tagged! parent ID "<<PARENTID<<" "<<evtAct->KPLUS_TRACKID<<" "<<imp<<" "<<pos<<G4endl;
	  if(evtAct->p_muplus_K==-99999.){// only at first step
	    if(ParticName=="mu+"){
	      evtAct->x_DIF_K_muplus = pos; 
	      evtAct->p_DIF_K_muplus = imp;
	    }
	    if(DEB)G4cout <<" FIRST mu+ "<<evtAct->x_DIF_K_muplus<<" "<<evtAct->p_DIF_K_muplus<< G4endl;
	    evtAct->p_muplus_K=evtAct->p_DIF_K_muplus.mag();
	  }	  
	}
      }//mu+
      if((ParticName==IMUMINUS)){
	//G4cout << "MU- PARENT IDENTITY "<<evtAct->vPDG_ID[PARENTID]<<G4endl;
	if(PARENTID==evtAct->PIMINUS_TRACKID){//pi- ~> mu- 
	  if(DEB)G4cout << "       "<<ParticName << " from " << ProcName << " tagged! parent ID "<<PARENTID<<" "<<evtAct->PIMINUS_TRACKID<<" "<<imp<<" "<<pos<<G4endl;
	  if(evtAct->p_muminus_PI==-99999.){// only at first step
	    if(ParticName=="mu-"){
	      evtAct->x_DIF_PI_muminus = pos; 
	      evtAct->p_DIF_PI_muminus = imp;

	      //G4cout << "MU- PARENT IDENTITY "<<parPDG<<" from a "<<parparPDG<<G4endl;
	      G4cout << "MU- <~ "<<parName<<" <~ "<<parparName<<" <~ "<<parparparName<<G4endl;
	      /*if((parPDG==PDGpiplus)||(parPDG==PDGpiminus)){*/
	      if((kaonGfath)||(kaonGGfath)){
		evtAct->ParLevel_muminus = 2;
		G4cout << "cascade   K+/- => PI => MU- "<<G4endl;
	      }else if((k0Gfath)||(k0GGfath)){
		evtAct->ParLevel_muminus = 3;
		G4cout << "cascade   K0 => PI => MU- "<<G4endl;
	      }else{
		evtAct->ParLevel_muminus = 1;
		G4cout << "direct   PI => MU- "<<G4endl;
	      }
	      /*	      }
			      else if(kaonfath){
			      evtAct->ParLevel_muminus = 4;
			      G4cout << "direct   K => MU- "<<G4endl;
			      }
	      */

	      //ooooooooooooooooooooo  
	      /*
		if(PARENTID==evtAct->PIMINUS_TRACKID){
		evtAct->ParLevel_muminus = evtAct->ParLevel_piminus + 1;
		G4cout << "SBSteppingAction: mu- from a pi- tagged. mu- ParLevel: "<< evtAct->ParLevel_muminus <<G4endl;
		}
	      */
	      //ooooooooooooooooooooo
	    }
	    if(DEB)G4cout <<" FIRST mu- "<<evtAct->x_DIF_PI_muminus<<" "<<evtAct->p_DIF_PI_muminus<< G4endl;
	    evtAct->p_muminus_PI=evtAct->p_DIF_PI_muminus.mag();
	  }
	}
	if(PARENTID==evtAct->KMINUS_TRACKID){//K- ~> mu- 
	  if(DEB)G4cout << "       "<<ParticName << " from " << ProcName << " tagged! parent ID "<<PARENTID<<" "<<evtAct->KMINUS_TRACKID<<" "<<imp<<" "<<pos<<G4endl;
	  if(evtAct->p_muminus_K==-99999.){// only at first step
	    if(ParticName=="mu-"){
	      evtAct->x_DIF_K_muminus = pos; 
	      evtAct->p_DIF_K_muminus = imp;
	    }
	    if(DEB)G4cout <<" FIRST mu- "<<evtAct->x_DIF_K_muminus<<" "<<evtAct->p_DIF_K_muminus<< G4endl;
	    evtAct->p_muminus_K=evtAct->p_DIF_K_muminus.mag();
	  }	  
	}
      }//mu-


      if((ParticName==INUMU)){
	if(PARENTID==evtAct->PIPLUS_TRACKID){//pi+ ~> numu 
	  if(DEB)G4cout << "       "<<ParticName << " from " << ProcName << " tagged! parent ID "<<PARENTID<<" "<<evtAct->PIPLUS_TRACKID<<" "<<imp<<" "<<pos<<G4endl;
	  if(evtAct->E_numu_PI==-99999.){// only at first step
	    if(ParticName=="nu_mu"){
	      evtAct->x_DIF_PI_numu = pos; 
	      evtAct->p_DIF_PI_numu = imp;
	    }
	    if(DEB)G4cout <<"PION FIRST nu_mu "<<evtAct->x_DIF_PI_numu<<" "<<evtAct->p_DIF_PI_numu<< G4endl;
	    evtAct->E_numu_PI=evtAct->p_DIF_PI_numu.mag();
	  }
	}
	if(PARENTID==evtAct->KPLUS_TRACKID){//K+ ~> numu 
	  if(DEB)G4cout << "       "<<ParticName << " from " << ProcName << " tagged! parent ID "<<PARENTID<<" "<<evtAct->KPLUS_TRACKID<<" "<<imp<<" "<<pos<<G4endl;
	  if(evtAct->E_numu_K==-99999.){// only at first step
	    if(ParticName=="nu_mu"){
	      evtAct->x_DIF_K_numu = pos; 
	      evtAct->p_DIF_K_numu = imp;
	    }
	    if(DEB)G4cout <<"KAPPA FIRST nu_mu "<<evtAct->x_DIF_K_numu<<" "<<evtAct->p_DIF_K_numu<< G4endl;
	    evtAct->E_numu_K=evtAct->p_DIF_K_numu.mag();
	  }
	}
      }//numu


      if((ParticName==IANUMU)){
	if(PARENTID==evtAct->PIMINUS_TRACKID){//pi- ~> antinumu 
	  if(DEB)G4cout << "       "<<ParticName << " from " << ProcName << " tagged! parent ID "<<PARENTID<<" "<<evtAct->PIPLUS_TRACKID<<" "<<imp<<" "<<pos<<G4endl;	  
	  if(evtAct->E_anumu_PI==-99999.){// only at first step
	    if(ParticName=="anti_nu_mu"){
	      evtAct->x_DIF_PI_anumu = pos; 
	      evtAct->p_DIF_PI_anumu = imp;
	    }
	    if(DEB)G4cout <<" FIRST anti_nu_mu "<<evtAct->x_DIF_PI_anumu<<" "<<evtAct->p_DIF_PI_anumu<< G4endl;
	    evtAct->E_anumu_PI=evtAct->p_DIF_PI_anumu.mag();
	  }
	}
	if(PARENTID==evtAct->KMINUS_TRACKID){//K- ~> antinumu 
	  if(DEB)G4cout << "       "<<ParticName << " from " << ProcName << " tagged! parent ID "<<PARENTID<<" "<<evtAct->KPLUS_TRACKID<<" "<<imp<<" "<<pos<<G4endl;
	  if(evtAct->E_anumu_K==-99999.){// only at first step
	    if(ParticName=="anti_nu_mu"){
	      evtAct->x_DIF_K_anumu = pos; 
	      evtAct->p_DIF_K_anumu = imp;
	    }
	    if(DEB)G4cout <<" FIRST anti_nu_mu "<<evtAct->x_DIF_K_anumu<<" "<<evtAct->p_DIF_K_anumu<< G4endl;
	    evtAct->E_anumu_K=evtAct->p_DIF_K_anumu.mag();
	  }
	}
      }//antinumu
      
    }//Decay
  }//CreatorProcess

  //if(track->GetParentID()==0 && track->GetUserInformation()==0) {
  //if(theProcessName == "Decay" && theParticleName=="pi+") {
  //if(track->GetKineticEnergy()>0.){
  //TrackInformation* anInfo = new TrackInformation(track);
  //track->SetUserInformation(anInfo);
  //}
  // if(DEB)G4cout << "pi plus decay"<< G4endl;
  //}
  //}
  

  // from gdecay.f modified
  /*
    G4double  mom1, mom, ener;
    G4bool twoBodiesProba = false;
    G4bool Kaon3bodies = false;
    // geant3
    G4int IGAM = 1;
    G4int IPOS = 2;
    G4int IELE = 3;
    G4int INEU = 4;
    G4int IMUPLUS = 5;
    G4int IMINUS = 6;
    G4int IPIPLUS = 8;
    G4int IPIMINUS = 9;
    G4int IKPLUS = 11;
    G4int IKMINUS = 12;
    G4int IK0L = 10;
    // manca MUMINUS nel check ?

    if(ID==1&&(ipart==IPIPLUS||ipart==IPIMINUS||ipart==IKPLUS||ipart==IKMINUS)) {
    mom1 = sqrt(pcm(1,1)**2+pcm(2,1)**2+pcm(3,1)**2);
    cthstar = PCM(3,1)/mom1;
    mom = VECT(7);
    ener = SQRT(mom*mom+AMASS*AMASS);
    betaPar  = mom/ener;
    gammaPar  = ener/AMASS;
    parentMass = AMASS;
    ieventSav = ievent;
    itrackSav = itra;
    twoBodiesProba = true;
    }
    if((ipart==IK0L)&&(ID==IPOS||ID==IELE||ID==INEU||ID==IMUPLUS)){
    mom1 = sqrt(pcm(1,2)**2+pcm(2,2)**2+pcm(3,2)**2);
    cthstar = PCM(3,2)/mom1;
    mom = VECT(7);
    ener = SQRT(mom*mom+AMASS*AMASS);
    betaPar  = mom/ener;
    gammaPar  = ener/AMASS;
    Kaon3bodies = true;
    }
    if((ipart==IKPLUS||ipart==IKMINUS)&&(ID==INEU||ID==IMUPLUS)){
    mom1 = sqrt(pcm(1,1)**2+pcm(2,1)**2+pcm(3,1)**2);
    cthstar = PCM(3,1)/mom1;
    mom = VECT(7);
    ener = SQRT(mom*mom+AMASS*AMASS);
    betaPar  = mom/ener;
    gammaPar  = ener/AMASS;
    Kaon3bodies = true;
    }


    if(twoBodiesProba)proba2(mom,cth);
    if(Kaon3bodies)proba3K(mom,cth);
 
  */

}

/***********************************************/
G4int SBSteppingAction::FLUKA_G4_PID(G4int pid_in=0){
  
  const G4int N = 30;
  G4int i_FLU[N]={0};
  G4int i_G4[N]={0};
  
  i_FLU[0]=13;//piplus
  i_FLU[1]=14;//piminus
  i_FLU[2]=23;//pizero
  i_FLU[3]=4;//eplus
  i_FLU[4]=3;//eminus
  i_FLU[5]=1;//proton
  i_FLU[6]=2;//antiproton
  i_FLU[7]=8;//neutron
  i_FLU[8]=9;//antineutron
  i_FLU[9]=7;//gamma
  i_FLU[10]=15;//kplus
  i_FLU[11]=16;//kminus
  i_FLU[12]=24;//kzero
  i_FLU[13]=25;//antikzero
  i_FLU[14]=17;//lambda
  i_FLU[15]=18;//antilambda
  i_FLU[16]=5;//nue
  i_FLU[17]=6;//antinue
  i_FLU[18]=27;//numu
  i_FLU[19]=28;//antinumu
  i_FLU[20]=43;//nutau
  i_FLU[21]=44;//antinutau
  i_FLU[22]=-6;//alpha
  i_FLU[23]=-5;//he3
  i_FLU[24]=-3;//deu
  i_FLU[25]=-4;//tri
  i_FLU[26]=10;//muplus
  i_FLU[27]=11;//muminus
  i_FLU[28]=21;//sigplus
  i_FLU[29]=20;//sigminus

  i_G4[0]=211;//piplus
  i_G4[1]=-211;//piminus
  i_G4[2]=111;//pizero
  i_G4[3]=-11;//eplus
  i_G4[4]=11;//eminus
  i_G4[5]=2212;//proton
  i_G4[6]=-2212;//antiproton
  i_G4[7]=2112;//neutron
  i_G4[8]=-2112;//antineutron
  i_G4[9]=22;//gamma
  i_G4[10]=321;//kplus
  i_G4[11]=-321;//kminus
  i_G4[12]=311;//kzero
  i_G4[13]=-311;//antikzero
  i_G4[14]=3122;//lambda
  i_G4[15]=-3122;//antilambda
  i_G4[16]=12;//nue
  i_G4[17]=-12;//antinue
  i_G4[18]=14;//numu
  i_G4[19]=-14;//antinumu
  i_G4[20]=16;//nutau
  i_G4[21]=-16;//antinutau
  i_G4[22]=0;//alpha NOT PRESENT
  i_G4[23]=0;//he3 NOT PRESENT
  i_G4[24]=0;//deu NOT PRESENT
  i_G4[25]=0;//tri NOT PRESENT
  i_G4[26]=-13;//muplus
  i_G4[27]=13;//muminus
  i_G4[28]=3222;//sigplus
  i_G4[29]=3112;//sigminus
 
  G4int pid_out=0; 

  for (G4int i=0;i<N;i++){
    if(i_G4[i]==pid_in) pid_out = i_FLU[i];
  }

  return pid_out;
}
