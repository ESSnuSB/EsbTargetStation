class SBSteppingVerbose;

#ifndef SBSteppingVerbose_h
#define SBSteppingVerbose_h 1
#include "G4SteppingVerbose.hh"

class SBSteppingVerbose : public G4SteppingVerbose
{
 public:   

   SBSteppingVerbose();
  ~SBSteppingVerbose();

   void StepInfo();
   void TrackingStarted();

};

#endif
