#ifndef SBProba_h
#define SBProba_h 1
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "SBDetectorConstruction.hh"
#include "SBRunAction.hh"
//#include "SBEventAction.hh"

class SBDetectorConstruction;
class SBRunAction;
class SBProba
{
public:

  SBProba(SBDetectorConstruction*,SBRunAction*/*,SBEventAction**/);
  ~SBProba();

  void SetSurf(G4double val){surf=val;};
  void SetDist(G4double val){dist=val;};

  G4double EnuForward(G4ThreeVector,G4double);
  G4double PowerNorm(G4double,G4double);
  void proba2(G4ThreeVector,G4String,G4int,G4int,G4double);
  void proba3K(G4ThreeVector,G4double,G4double,G4int);
  void probaMu(G4ThreeVector,G4ThreeVector,G4ThreeVector,G4double,G4int,G4int,G4double);
  G4double dkproba(G4ThreeVector,G4ThreeVector);

 private: 

  G4double dist;
  G4double surf;
  // in MeV/c2
  G4double me;//e mass
  G4double mmu;//mu mass
  G4double mpi;//pi+- mass
  G4double mpi0;//pi0 mass
  G4double mK0;//K0 mass
  G4double mK;//K+- mass

  SBDetectorConstruction* SBDetector;
  SBRunAction* runAct;
  //SBEventAction* evtAct;

  G4int ipart;

  G4double cthstarMu;
  G4double cthstarNu;
  G4double betaPar;
  G4double gammaPar;
  //G4double ammu;
  //G4double ammu2;
  //G4double ampi;
  //G4double ampi2;
  G4double ctmu; 
  G4int NBINE;// Number of bins in neutrino's energy
  G4double de;                   // energy bin size
  
  //G4double ProbnmP[100];       //   nu_mu : pi+ -> mu+ nu_mu
  //G4double ProbamP[100];       //  anu_mu : pi- -> mu- anu_mu
  //G4double Probne[100];        //   nu_e  : mu+ -> e+ nu_e anu_mu
  //G4double Probnm[100];        //  anu_mu : mu- -> e- anu_e nu_mu
  //G4double Probae[100];        //  anu_e  : mu- -> e- anu_e nu_mu
  //G4double Probam[100];        //   nu_mu : mu+ -> e+ anu_e anu_mu  
  //
  //     Decay products of charged kaons
  //  
  //G4double ProbnmK[100];       //   nu_mu : K+ -> mu+ nu_mu
  //G4double ProbamK[100];       //  anu_mu : K- -> mu- anu_mu
  
  //G4double ProbCKnue[100];      // 
  //G4double ProbCKanue[100];     //   3 bodies decays of K+ and K- 
  //G4double ProbCKnumu[100];     //
  //G4double ProbCKanumu[100];    //
  
  //     neutrino from the decay of pion form a charged kaon decay chain.
  
  //G4double ProbnmCKpi[100];    //   nu_mu : pi+ -> mu+ nu_mu
  //G4double ProbamCKpi[100];    //  anu_mu : pi- -> mu- anu_mu
  
  //     neutrino from the decay of muons from a charged kaon decay chain.
  
  //G4double ProbneCKmu[100];    //   nu_e  : mu+ -> e+ nu_e anu_mu
  //G4double ProbnmCKmu[100];    //  anu_mu : mu- -> e- anu_e nu_mu
  //G4double ProbaeCKmu[100];    //  anu_e  : mu- -> e- anu_e nu_mu
  //G4double ProbamCKmu[100];    //   nu_mu : mu+ -> e+ anu_e anu_mu
  
  //
  //     Decay products of neutral kaons
  //
  
  //G4double ProbK0nue[100];       // 
  //G4double ProbK0anue[100];      //   3 bodies decays of K0 long
  //G4double ProbK0numu[100];      //
  //G4double ProbK0anumu[100];     //
  
  //     neutrino from the decay of pion form a neutral kaon decay chain.
  
  //G4double ProbnmK0pi[100];    //   nu_mu : pi+ -> mu+ nu_mu
  //G4double ProbamK0pi[100];    //  anu_mu : pi- -> mu- anu_mu
  
  //     neutrino from the decay of muons from a neutral kaon decay chain.
  
  //G4double ProbneK0mu[100];    //   nu_e  : mu+ -> e+ nu_e anu_mu
  //G4double ProbnmK0mu[100];    //  anu_mu : mu- -> e- anu_e nu_mu
  //G4double ProbaeK0mu[100];    //  anu_e  : mu- -> e- anu_e nu_mu
  //G4double ProbamK0mu[100];    //   nu_mu : mu+ -> e+ anu_e anu_mu
  
  //G4int isFromKaon;
  G4int multiplicity;
  //to make it a ROOT class
  //ClassDef(SBDetectorConstruction,2)//SBDetectorConstruction
};
#endif
