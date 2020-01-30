#include "G4UnitsTable.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include "SBField.hh"
#include "SBDetectorConstruction.hh"

//SBField::SBField(G4double I_OneHornIntensity, G4int hID, SBDetectorConstruction* SBDC):SBDetector(SBDC)
SBField::SBField(G4double I_OneHornIntensity, G4int hID, SBDetectorConstruction* SBDC)
{
  //I0 = 300000.*ampere;
  //R0 = 4.*cm;//InCondRad[0] + InCondDepth[1];
  I0=I_OneHornIntensity;
  G4double MU04PI = 1.e-7*(CLHEP::tesla*CLHEP::m)/CLHEP::ampere;      ///(µ0/4Pi)
  K0=MU04PI *2.*I0;

  G4cout << "I0 = " <<I0/CLHEP::ampere <<" A "<<G4endl; 

  x0 = SBDC->Horn_xpos;                                               ////(modif)
  y0 = SBDC->Horn_ypos;                                               ////(modif)

  G4cout << " Horn Id    = " << hID << G4endl ; 
  
  if(hID==0){x0*=1.;y0*=1.;}
  if(hID==1){x0*=-1.;y0*=1.;}                                                        //
  if(hID==2){x0*=1.;y0*=-1.;}                                                        //
  if(hID==3){x0*=-1.;y0*=-1.;}                                                //

  //G4cout << " Horn_xpos = " <<  x0 << " cm ; Horn_ypos = " << y0 << G4endl;
  G4cout <<"SBField: Horn position (x,y) cm "<< x0/CLHEP::cm <<" "<<y0/CLHEP::cm<<G4endl;
}

SBField::~SBField(){;}

void SBField::GetFieldValue(const G4double point[3],G4double *Bfield) const
{
  G4bool DEB=false;
  //G4bool DEB=true;
  G4double r = 0*CLHEP::cm;
  // COrrect value for r
  r=sqrt(((point[0]-x0)*(point[0]-x0))+((point[1]-y0)*(point[1]-y0)));
  
  // Test Value
  //r=sqrt((point[0]*point[0])+(point[1]*point[1]));
  
  //G4cout << " Get Field Value ; Horn_xpos = " <<  x0 << " cm ; Horn_ypos = " << y0 << G4endl;
  //G4cout << " Rayon r = " << r << G4endl;
  
  if(r>1*CLHEP::cm) //if((sqr(point[0])+sqr(point[1]))<rmax_sq)
  {
    Bfield[0] = -(K0/r)*((point[1]-y0)/r);
    Bfield[1] =  (K0/r)*((point[0]-x0)/r);
    Bfield[2] = 0;
        
  } else { 
    Bfield[0] = 0.;
    Bfield[1] = 0.;
    Bfield[2] = 0.;
     
  }

  if(DEB)G4cout <<"x "<<point[0]/CLHEP::cm<<" (cm) y "<<point[1]/CLHEP::cm<<" (cm) z "<<point[2]/CLHEP::cm<<" (cm) Bx " << Bfield[0]/CLHEP::tesla <<" (T) By "<<Bfield[1]/CLHEP::tesla<<" (T) Bz "<<Bfield[2]/CLHEP::tesla<<" (T)"<<G4endl;

}
