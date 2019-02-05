#include "SBPrimaryGeneratorAction.hh"
#include "SBDetectorConstruction.hh"
#include "SBPrimaryGeneratorMessenger.hh"
#include "SBRunAction.hh"

#include "G4Event.hh"
#include "G4HEPEvtInterface.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4KaonZero.hh"
#include "Randomize.hh"
#include "TMath.h"

using namespace TMath;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo

SBPrimaryGeneratorAction::SBPrimaryGeneratorAction(SBDetectorConstruction* SBDC):SBDetector(SBDC),rndmFlag("off")
{
  // Create a messenger for this class
  // -- 
  gunMessenger = new SBPrimaryGeneratorMessenger(this);

  procEvts 	= 0;
  procTracks 	= 0;
  OLDidev 	= 0;

  K_REPL=1;//kaon replication

  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  
  // Particule Gun generation
  // --
  G4ParticleTable* 	particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName="proton");
  particleGun->SetParticleDefinition(particle);
  
  x0 = 100.*CLHEP::cm;                                               
  y0 = 100.*CLHEP::cm;                                               
  G4double z0 = SBDetector->GetTargetZetaPos()-0.5*SBDetector->GetTargetThickness()-500;
  gunEKin = 4.5*CLHEP::GeV;

  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));  
  particleGun->SetParticleEnergy(gunEKin);        
    
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  particleGun->SetParticlePosition(G4ThreeVector(-x0,y0,z0));  
  particleGun->SetParticleEnergy(gunEKin); 
     
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  particleGun->SetParticlePosition(G4ThreeVector(x0,-y0,z0));  
  particleGun->SetParticleEnergy(gunEKin);
  
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  particleGun->SetParticlePosition(G4ThreeVector(-x0,-y0,z0));  
  particleGun->SetParticleEnergy(gunEKin);
  
  G4cout <<"SBPrimaryGenerationAction::SBPrimaryGeneratorAction gun kinetic energy is "<< gunEKin/CLHEP::GeV <<" GeV. "<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo

SBPrimaryGeneratorAction::~SBPrimaryGeneratorAction()
{
  G4cout<<"CALL SBPrimaryGeneratorAction::~SBPrimaryGeneratorAction "<<G4endl;
  fin.close();

  G4cout<<"SBPrimaryGenerationAction: TOTAL processed events: "<<procEvts<<" (may be different from FLUKA pots if preselection is applied in the input file)"<<G4endl;
  G4cout<<"SBPrimaryGenerationAction: TOTAL processed single tracks: "<<procTracks<<G4endl;

  delete particleGun;
  delete gunMessenger;
  delete HEPEvt;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo

void SBPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4cout<<"========================================================"<<G4endl;
  G4cout<<"===============   S  T  A  R  T  ======================="<<G4endl;
  G4cout<<"========================================================"<<G4endl;

  if(ExtGen==1)externalGEN=true;
  if(ExtGen==0)externalGEN=false;
  
  myEOF=false;
  
  G4int evid =anEvent->GetEventID();
  
  if(evid==0){
    G4cout << ExtGen <<" GENERATOR"<<G4endl;
    if(externalGEN)G4cout<<"SBPrimaryGeneratorAction USING EXTERNAL GENERATOR"<<G4endl;
    if(!externalGEN)G4cout<<"SBPrimaryGeneratorAction USING GEANT4 GENERATOR"<<G4endl;  
  }

  if(externalGEN){
    if(evid==0){
      InputFileOK=true;
      G4cout <<"SBPrimaryGenerationAction: Opening input ASCII file" << filenamein <<" " <<fin.fail()<<G4endl;
      fin.open(filenamein);
      if(fin.fail()){
	G4cout <<"SBPrimaryGenerationAction: Input ASCII file "<<filenamein <<" NOT found. Check name in control card (.mac). Exiting." << G4endl;
	InputFileOK=false;
	return;
      } else {
	G4cout <<"SBPrimaryGenerationAction: Opened input ASCII file" << filenamein <<" "<<fin.fail()<< G4endl;
      }
    }

    if(!InputFileOK){
      G4cout << "SBPrimaryGenerationAction: problem with input ASCII file "<< filenamein << G4endl;
      return;
    }

  }

  G4double z0 = SBDetector->GetTargetZetaPos()-0.5*SBDetector->GetTargetThickness()-500;

  if (rndmFlag == "on"){
    y0 = (SBDetector->GetTargetDiameter())*(G4UniformRand()-0.5);
    z0 = (SBDetector->GetTargetDiameter())*(G4UniformRand()-0.5);
  }  

  G4cout<<"0 "<<G4endl;

  if(!externalGEN){

    if(evid==0)G4cout << "SBPrimaryGenerationAction: Using ParticleGun for primary generations"<< G4endl;

    G4ParticleTable* partTab = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* mypart;
    mypart = partTab->FindParticle("proton");
    particleGun->SetParticleDefinition(mypart);
   
    particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
    particleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
    particleGun->SetParticleEnergy(gunEKin);
    G4cout <<"SBPrimaryGenerationAction: gun kinetic energy is "<< gunEKin/CLHEP::GeV <<" GeV. "<<G4endl;
    particleGun->GeneratePrimaryVertex(anEvent);

    particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
    particleGun->SetParticlePosition(G4ThreeVector(-x0,y0,z0));  
    G4cout <<"SBPrimaryGenerationAction: gun kinetic energy is "<< gunEKin/CLHEP::GeV <<" GeV. "<<G4endl;
    particleGun->SetParticleEnergy(gunEKin); 
    particleGun->GeneratePrimaryVertex(anEvent);
   
    particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
    particleGun->SetParticlePosition(G4ThreeVector(x0,-y0,z0));  
    G4cout <<"SBPrimaryGenerationAction: gun kinetic energy is "<< gunEKin/CLHEP::GeV <<" GeV. "<<G4endl;
    particleGun->SetParticleEnergy(gunEKin);
    particleGun->GeneratePrimaryVertex(anEvent);

    particleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
    particleGun->SetParticlePosition(G4ThreeVector(-x0,-y0,z0));  
    G4cout <<"SBPrimaryGenerationAction: gun kinetic energy is "<< gunEKin/CLHEP::GeV <<" GeV. "<<G4endl;
    particleGun->SetParticleEnergy(gunEKin);
    particleGun->GeneratePrimaryVertex(anEvent);   
 
    G4cout << "SBPrimaryGenerationAction: GEANT4 primary vtx =( "<<x0/CLHEP::cm<<", "<<y0/CLHEP::cm<<", "<<z0/CLHEP::cm<<" ) (cm)"<< G4endl;
    FLUKApots=FLUKApots+1.;

    G4cout << "SBPrimaryGenerationAction: nK+ "<<n_kplus_exit<<G4endl;
    G4cout << "SBPrimaryGenerationAction: nK- "<<n_kminus_exit<<G4endl;
    G4cout << "SBPrimaryGenerationAction: nK0S "<<n_k0s_exit<<G4endl;
    G4cout << "SBPrimaryGenerationAction: nK0L "<<n_k0l_exit<<G4endl;
    G4cout << "SBPrimaryGenerationAction: REPLICATE kaons "<<K_REPL-1<<" times."<<G4endl;

    G4double kekin=0;
    
    mypart = partTab->FindParticle("kaon+");
    particleGun->SetParticleDefinition(mypart);
    for(int ii=0;ii<n_kplus_exit;ii++){
      G4cout<<"k+ "<<x_kplus_exit[ii]<<" "<<y_kplus_exit[ii]<<" "<<z_kplus_exit[ii]<<" "<<px_kplus_exit[ii]<<" "<<py_kplus_exit[ii]<<" "<<pz_kplus_exit[ii]<<G4endl;
      particleGun->SetParticleMomentumDirection(G4ThreeVector(px_kplus_exit[ii],py_kplus_exit[ii],pz_kplus_exit[ii]));
      particleGun->SetParticlePosition(G4ThreeVector(x_kplus_exit[ii],y_kplus_exit[ii],z_kplus_exit[ii]));
      kekin = sqrt(pow(px_kplus_exit[ii],2)+pow(py_kplus_exit[ii],2)+pow(pz_kplus_exit[ii],2))-mypart->GetPDGMass();
      particleGun->SetParticleEnergy(kekin);
      for(int ll=0;ll<K_REPL-1;ll++){   
	particleGun->GeneratePrimaryVertex(anEvent);
      }
    }
    mypart = partTab->FindParticle("kaon-");
    particleGun->SetParticleDefinition(mypart);    
    for(int ii=0;ii<n_kminus_exit;ii++){
      G4cout<<"k- "<<x_kminus_exit[ii]<<" "<<y_kminus_exit[ii]<<" "<<z_kminus_exit[ii]<<" "<<px_kminus_exit[ii]<<" "<<py_kminus_exit[ii]<<" "<<pz_kminus_exit[ii]<<G4endl;
      particleGun->SetParticleMomentumDirection(G4ThreeVector(px_kminus_exit[ii],py_kminus_exit[ii],pz_kminus_exit[ii]));
      particleGun->SetParticlePosition(G4ThreeVector(x_kminus_exit[ii],y_kminus_exit[ii],z_kminus_exit[ii]));
      kekin = sqrt(pow(px_kminus_exit[ii],2)+pow(py_kminus_exit[ii],2)+pow(pz_kminus_exit[ii],2))-mypart->GetPDGMass();
      particleGun->SetParticleEnergy(kekin);
      for(int ll=0;ll<K_REPL-1;ll++){   
	particleGun->GeneratePrimaryVertex(anEvent);
      }
    }
    mypart = partTab->FindParticle("kaon0S");
    particleGun->SetParticleDefinition(mypart);    
    for(int ii=0;ii<n_k0s_exit;ii++){
      G4cout<<"k0s "<<x_k0s_exit[ii]<<" "<<y_k0s_exit[ii]<<" "<<z_k0s_exit[ii]<<" "<<px_k0s_exit[ii]<<" "<<py_k0s_exit[ii]<<" "<<pz_k0s_exit[ii]<<G4endl;
      particleGun->SetParticleMomentumDirection(G4ThreeVector(px_k0s_exit[ii],py_k0s_exit[ii],pz_k0s_exit[ii]));
      particleGun->SetParticlePosition(G4ThreeVector(x_k0s_exit[ii],y_k0s_exit[ii],z_k0s_exit[ii]));
      kekin = sqrt(pow(px_k0s_exit[ii],2)+pow(py_k0s_exit[ii],2)+pow(pz_k0s_exit[ii],2))-mypart->GetPDGMass();
      particleGun->SetParticleEnergy(kekin);
      for(int ll=0;ll<K_REPL-1;ll++){   
	particleGun->GeneratePrimaryVertex(anEvent);
      }
    }
    mypart = partTab->FindParticle("kaon0L");
    particleGun->SetParticleDefinition(mypart);    
    for(int ii=0;ii<n_k0l_exit;ii++){
      G4cout<<"k0l "<<x_k0l_exit[ii]<<" "<<y_k0l_exit[ii]<<" "<<z_k0l_exit[ii]<<" "<<px_k0l_exit[ii]<<" "<<py_k0l_exit[ii]<<" "<<pz_k0l_exit[ii]<<G4endl;
      particleGun->SetParticleMomentumDirection(G4ThreeVector(px_k0l_exit[ii],py_k0l_exit[ii],pz_k0l_exit[ii]));
      particleGun->SetParticlePosition(G4ThreeVector(x_k0l_exit[ii],y_k0l_exit[ii],z_k0l_exit[ii]));
      kekin = sqrt(pow(px_k0l_exit[ii],2)+pow(py_k0l_exit[ii],2)+pow(pz_k0l_exit[ii],2))-mypart->GetPDGMass();
      particleGun->SetParticleEnergy(kekin);
      for(int ll=0;ll<K_REPL-1;ll++){   
	particleGun->GeneratePrimaryVertex(anEvent);
      }
    }

    n_kplus_exit=0;
    n_kminus_exit=0;
    n_k0s_exit=0;
    n_k0l_exit=0;

    //////////////////////////////////////////////////////////////////////////////////////////////
  } else {
    if(evid==0)G4cout << "SBPrimaryGenerationAction: using ASCII file input for primary generations"<< G4endl;

    G4cout<<"1 "<<G4endl;

    double idev=0,pid=0,x=0*CLHEP::cm,y=0*CLHEP::cm,z=0*CLHEP::cm,px=0*CLHEP::GeV,py=0*CLHEP::GeV,pz=0*CLHEP::GeV,dum0=0,dum1=0,dum2=0;

    if(fin.eof())myEOF=true;
    fin >> idev >> pid >> x >> y >> z >> px >> py >> pz >> dum0 >> dum1 >> dum2;
    G4cout<<"2 "<<G4endl;

    if(idev!=0)FLUKApots = idev;// NB. here FLUKa evt ID is assumed to start from zero!!!
      
    if(idev!=OLDidev){
      procEvts++;
      OLDidev=idev;
    }
    procTracks++;
      
    G4cout<<"3 "<<G4endl;
    particleName="pi-";
    G4double mass = 0.139*CLHEP::GeV;
      
    bool skip = false;
      
    G4bool Antoine = false;
    if(Antoine){
      if(pid==-3)particleName="deuteron";
      else if(pid==-4)particleName="triton";
      else if(pid==-6)particleName="alpha";
      else if(pid==1)particleName="proton";
      else if(pid==2)particleName="anti_proton";//check name
      else if(pid==3)particleName="e-";
      else if(pid==4)particleName="e+";
      else if(pid==7)particleName="gamma";
      else if(pid==8)particleName="neutron";
      else if(pid==9)particleName="anti_neutron";//check name
      else if(pid==10)particleName="mu+";
      else if(pid==11)particleName="mu-";
      else if(pid==12)particleName="kaon0L";
      else if(pid==13)particleName="pi+";
      else if(pid==14)particleName="pi-";
      else if(pid==15)particleName="kaon+";
      else if(pid==16)particleName="kaon-";
      else if(pid==17)particleName="lambda";//check name
      else if(pid==18)particleName="kaon0S";//? 18 in converison by Antoine gukine.F
      //else if(pid==19)particleName="kaon0S";//
      else if(pid==20)particleName="sigma-";//check name
      else if(pid==21)particleName="sigma+";//check name
      else if(pid==22)particleName="sigma0";//check name
      else if(pid==23)particleName="pi0";
      else if(pid==24)particleName="kaon0";
      //else if(pid==25)particleName="anti_kaon0";
      else {
	G4cout << "SBPrimaryGenerationAction: FLUKA track with pid " << pid <<" will not be not tracked by GEANT4 " << G4endl;
	skip=true;
      }
    }else{
      if(pid==-3)particleName="deuteron";
      else if(pid==-4)particleName="triton";
      else if(pid==-6)particleName="alpha";
      else if(pid==1)particleName="proton";
      else if(pid==2)particleName="anti_proton";//check name
      else if(pid==3)particleName="e-";
      else if(pid==4)particleName="e+";
      else if(pid==7)particleName="gamma";
      else if(pid==8)particleName="neutron";
      else if(pid==9)particleName="anti_neutron";//check name
      else if(pid==10)particleName="mu+";
      else if(pid==11)particleName="mu-";
      else if(pid==12)particleName="kaon0L";
      else if(pid==13)particleName="pi+";
      else if(pid==14)particleName="pi-";
      else if(pid==15)particleName="kaon+";
      else if(pid==16)particleName="kaon-";
      else if(pid==17)particleName="lambda";//check name
      else if(pid==18)particleName="anti_lambda";//check name
      else if(pid==19)particleName="kaon0S";//? check 18 in converison by Antoine gukine.F
      else if(pid==20)particleName="sigma-";//check name
      else if(pid==21)particleName="sigma+";//check name
      else if(pid==22)particleName="sigma0";//check name
      else if(pid==23)particleName="pi0";
      else if(pid==24)particleName="kaon0";
      else if(pid==25)particleName="anti_kaon0";
      else {
	G4cout << "SBPrimaryGenerationAction: FLUKA track with pid " << pid <<" will not be not tracked by GEANT4 " << G4endl;
	skip=true;
      }
    }
 
    G4cout<<"4 "<<particleName<<G4endl;

    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle = particleTable->FindParticle(particleName);

    G4cout<<"4a "<<particle<<G4endl;
    if(particleName=="kaon0")G4cout<<"kaon0 shoot"<<G4endl;

    particleGun->SetParticleDefinition(particle);

    mass=particle->GetPDGMass();
    G4cout<<"5 "<<G4endl;

    //tmp! shoot only muons =========
    //G4cout <<"MUON SHOOTING"<<G4endl;
    //particle = particleTable->FindParticle("mu+");
    //particleGun->SetParticleDefinition(particle);
    //mass=0.105*CLHEP::GeV;
    //px=0;
    //py=0;
    //pz=0.001*CLHEP::GeV*G4UniformRand();
    //===============================      

    px*=1000;
    py*=1000;
    pz*=1000;
      
    x*=10;
    y*=10;
    z*=10;
      
    XPrimary=x;
    YPrimary=y;
    ZPrimary=z;
    PXPrimary=px;
    PYPrimary=py;
    PZPrimary=pz;
      
    ZPrimary_G4frame=z+z0;
      
    G4double En = sqrt((mass*mass)+(px*px)+(py*py)+(pz*pz));
    G4double Ekin = En - mass;
      
    G4cout << "FLUKA ID evt "<< idev <<" "<< particleName << G4endl;
    G4cout << "ASCII: x "<< x/CLHEP::cm <<
      " y "<< y/CLHEP::cm <<
      " z "<< z/CLHEP::cm <<" cm"<<G4endl; 
      
    G4cout << "DETECTOR COORD: x "<< x/CLHEP::cm <<
      " y "<< y/CLHEP::cm <<
      " z "<< z/CLHEP::cm <<" cm"<<G4endl; 

    G4cout << "z detector "<< z/CLHEP::cm <<" cm"<< G4endl;
    G4cout << "m "<< mass/CLHEP::GeV <<" GeV"<< G4endl;
    G4cout << "px "<< px/CLHEP::GeV <<" GeV"<< G4endl;
    G4cout << "py "<< py/CLHEP::GeV <<" GeV"<< G4endl;
    G4cout << "pz "<< pz/CLHEP::GeV <<" GeV"<< G4endl;
    G4cout << "p "<< sqrt((px*px)+(py*py)+(pz*pz))/CLHEP::GeV <<" GeV"<< G4endl;
    G4cout << "Etot "<< En/CLHEP::GeV <<" GeV"<< G4endl;
    G4cout << "Ekin "<< Ekin/CLHEP::GeV <<" GeV"<< G4endl;
      
    if(evid==0)G4cout << "SBPrimaryGenerationAction: z start of target at generation level "<< z0/CLHEP::m <<" m"<< G4endl;
      
    if(!skip){
      particleGun->SetParticleMomentumDirection(G4ThreeVector(px,py,pz));
      particleGun->SetParticlePosition(G4ThreeVector(x,y,z));
      particleGun->SetParticleEnergy(Ekin);
      particleGun->GeneratePrimaryVertex(anEvent);

      //
      // kaon duplication
      //
	
      //G4int K_REPL = 1;
      
      if((particleName=="kaon0")||
	 (particleName=="anti_kaon0")||
	 (particleName=="kaon+")||
	 (particleName=="kaon-")	   
	 ){
	G4cout <<"SBPrimaryGeneratorAction: "<<particleName<<" will be replicated "<<K_REPL-1<<" times."<<G4endl;
	for(int ll=0;ll<K_REPL-1;ll++){
	  G4cout<<"ciao "<<ll<<G4endl;	    
	  particleGun->GeneratePrimaryVertex(anEvent);
	}	  
      }
    }
  }
}
