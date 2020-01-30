#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "Randomize.hh"
#include "SBProba.hh"
#include "SBDetectorConstruction.hh"
#include "QGSP_BERT.hh"
#include "SBPrimaryGeneratorAction.hh"
#include "SBRunAction.hh"
#include "SBEventAction.hh"
#include "SBSteppingAction.hh"
#include "SBSteppingVerbose.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#if defined(G4UI_USE_TCSH)
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#elif defined(G4UI_USE_XM)
#include "G4UIXm.hh"
#elif defined(G4UI_USE_WIN32)
#include "G4UIWin32.hh"
#elif defined(G4UI_USE_QT)
#include "G4UIQt.hh"
#include "G4Qt.hh"
#else
#include "G4UIterminal.hh"
#endif

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

namespace {
  void PrintUsage() {
    G4cerr << " Usage: " << G4endl;
    G4cerr << " exampleB4c [-m macro ] [-u UIsession] [-t nThreads]" << G4endl;
    G4cerr << "   note: -t option is available only for multi-threaded mode."
           << G4endl;
  }
}

int main(int argc,char** argv)
{

  // Evaluate arguments
  // --
  if ( argc > 7 ) {
    PrintUsage();
    return 1;
  }

  G4String macro;
  G4String session;

  // Detect interactive mode (if no macro provided) and define UI session
  // -- 
  G4UIExecutive* ui = 0;
  if ( ! macro.size() ) {
    ui = new G4UIExecutive(argc, argv, session);
  }

  // Choose the Random engine
  // --
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  
  // User Verbose output class
  // --
  G4VSteppingVerbose::SetInstance(new SBSteppingVerbose);

  // Define Manager
  // --
  G4RunManager * runManager = new G4RunManager;

  // Set mandatory initialization classes
  // --
  SBDetectorConstruction* detector = new SBDetectorConstruction;
  runManager->SetUserInitialization(detector);

  // Set physics list classes
  // --
  //G4VUserPhysicsList* physics = new SBPhysicsList;
  G4VUserPhysicsList* physics = new QGSP_BERT;
  runManager->SetUserInitialization(physics);

  // Set user action classes
  // --
  //G4VUserPrimaryGeneratorAction* gen_action = new SBPrimaryGeneratorAction(detector);
  SBPrimaryGeneratorAction* gen_action = new SBPrimaryGeneratorAction(detector);
  runManager->SetUserAction(gen_action);

  // Run Action
  // --
  SBRunAction* run_action = new SBRunAction(gen_action);  
  runManager->SetUserAction(run_action);

  // Probability
  // --
  SBProba *probab = new SBProba(detector,run_action);
  SBEventAction* event_action = new SBEventAction(run_action,probab);
  runManager->SetUserAction(event_action);

  // Stepping Action
  // --
  G4UserSteppingAction* stepping_action = new SBSteppingAction(detector, event_action, gen_action, run_action);
  runManager->SetUserAction(stepping_action);

  //G4BlineTracer* theBlineTool = new G4BlineTracer();    //track mag_field mache pas ( need G4MageticField)
  
  // Initialize G4 kernel
  // --
  runManager->Initialize();
  
//#ifdef G4VIS_USE
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
//#endif
  
  // Get the pointer to the User Interface manager
  //
  G4UImanager* UI = G4UImanager::GetUIpointer();      
     
  int dum=0;
  G4cout << Form("\\rm -f g4_??.eps") << G4endl;
  dum=system(Form("\\rm -f g4_??.eps"));
  G4cout << Form("\\rm -f g4_??.prim") << G4endl;
  dum=system(Form("\\rm -f g4_??.prim"));

  G4String command = "/control/execute ";
  G4String fileName = "";
  if(argc>1)fileName=argv[1];

  if (argc!=1)   // batch mode
    {
      UI->ApplyCommand(command+fileName);    
    }
  else           // interactive mode : define visualization UI terminal
    {

     G4UIsession* session = 0;

    // Get the pointer to the User Interface manager
    // --
    auto UImanager = G4UImanager::GetUIpointer();

    UImanager->ApplyCommand("/control/execute init_vis.mac");
//    if (ui->IsGUI()) {
      UImanager->ApplyCommand("/control/execute gui.mac");
//    }
    ui->SessionStart();
    delete ui;


/*
#if defined(G4UI_USE_TCSH)
      G4cout << " ****************** TCSH *************************" << G4endl;
      session = new G4UIterminal(new G4UItcsh);
#elif defined(G4UI_USE_XM)
      G4cout << " ******************  XM  *************************" << G4endl;
      session = new G4UIXm(argc,argv);
      UI->ApplyCommand("/control/execute visTutor/gui.mac");      
#elif defined(G4UI_USE_WIN32)
      G4cout << " ****************** WIN32 ************************" << G4endl;
      session = new G4UIWin32();
      UI->ApplyCommand("/control/execute visTutor/gui.mac");      
#elif defined(G4UI_USE_QT)
      G4cout << " ******************  QT  *************************" << G4endl;
      session = new G4UIQt(argc,argv);
      UI->ApplyCommand("/control/execute visTutor/gui.mac");      
#else
*/ 
     session = new G4UIterminal();
/*
#endif  
#ifdef G4VIS_USE
      //G4cout << " ****************** /control/execute myrun.mac *************************" << G4endl;
      //UI->ApplyCommand("/control/execute myrun.mac");     
#endif
*/
      session->SessionStart();
      delete session;
    }

  // Job termination
  // --
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !

  //G4cout <<"VERBOSITY "<<event_action->GetSBVerbosity()<<G4endl;

  //G4cout << Form("cp SB_simul_v0.root %s.root",fileName.c_str()) << G4endl;
  //dum=system(Form("cp SB_simul_v0.root %s.root",fileName.c_str()));

  //G4cout << Form("cp nufluxes_GLOBESformat.txt %s_GLOBESformat.txt",fileName.c_str()) << G4endl;
  //dum=system(Form("cp nufluxes_GLOBESformat.txt %s_GLOBESformat.txt",fileName.c_str()));

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

