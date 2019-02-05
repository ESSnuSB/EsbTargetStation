#ifndef SBField_H
#define SBField_H 1

#include "globals.hh"
#include "G4MagneticField.hh"
#include "SBDetectorConstruction.hh"

class SBDetectorConstruction;

class SBField : public G4MagneticField
{
public:
  SBField(G4double,G4int,SBDetectorConstruction*);
  ~SBField();
  
  void GetFieldValue( const G4double point[3], G4double *Bfield ) const;
  
  void SetCurrent(G4double vI){I0=vI;};
  //void SetRadius(G4double vR){R0=vR;};
  
  G4double GetCurrent(){return I0;};
  //G4double GetRadius(){return R0;};
  
private:
  //SBDetectorConstruction* SBDetector;
  G4double I0;
  G4double x0;
  G4double y0;
  G4double K0;
  //G4double rmax_sq;
  //G4double zmax;
};

#endif
