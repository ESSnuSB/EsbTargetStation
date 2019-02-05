//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo
//
// SBDetectorconstructionNew ->  Detector construction.    
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo

#ifndef SBDetectorConstruction_h
#define SBDetectorConstruction_h 1

#include "SBDetectorConstructionMessenger.hh" 
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include <Rtypes.h>

class G4Box;
class G4Tubs;
class G4Polycone;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4SubtractionSolid;
class G4Material;
class SBDetectorConstructionMessenger;

class SBDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  
  SBDetectorConstruction();
  ~SBDetectorConstruction();
     
private:
   
  G4double OffAxisAngle;
  G4double OffAxisPhi;  
 
  G4Material* TargetMaterial;
  G4double TargetThickness;
  G4double TargetHalfThickness;
  G4double TargetDiameter;
  G4double TargetZetaPos;
     
  G4Material* defaultMaterial;

  G4Material* HornMaterial;
  G4Material* TunnelMaterial;

  G4Material* FillingMaterial;

  G4double rO_Hall;
  G4double L_Hall;

//....oooOO0OOooo....  Inner conductor  ....oooOO0OOooo....
  G4double InCondHalfLength;      

//....oooOO0OOooo....  Outer Conductor  ....oooOO0OOooo.... 
  G4double OutCondHalfLength;

//....oooOO0OOooo....  Double Skin (Water surround the horn and the double skin surround the water)  ....oooOO0OOooo.... ////////////////////////////
  G4double WaterThickness; 
  G4double SkinThickness;

//....oooOO0OOooo....  Target  ....oooOO0OOooo....
  G4double TargetRadius; 
  G4double TargetLength;

// ....oooOO0OOooo....  Tunnel  ....oooOO0OOooo....
  G4double TunnelRadius; 
  G4double TunnelLength;

//....oooOO0OOooo....  Horn  ....oooOO0OOooo....
  G4int HornDesign;
  G4int DoSkin;
  G4int DoRefl;
  G4int InjSchema;
  G4int Horn_N;
  G4double HornDisplacRad;
  G4double L1;               // = 58.9*cm;
  G4double L2;               // = 46.8*cm;
  G4double L3;               // = 60.3*cm;
  G4double L4;               // = 47.5*cm;
  G4double L5;               // = 10.8*cm;
  G4double t1;               // = 3*cm;
  G4double t2;               // = 3*cm;
  G4double t3;               // = 3*cm;
  G4double t4;               // = 10*cm;
  G4double R;                // = 12*cm;  //radial extension  (target stick to the horn)
  G4double r3;               // = 50.8*cm;//curvature radius 
  G4double R4;               // = 27.2*cm //curvature radius
  G4double R2;               // = 19.1*cm;//radial extension  
  G4double R3;               // = 35.9*cm;//radial extension

  G4double Horn_I1;          //300000A                                                   ////////////////////////////////////////////////
  G4double Horn_I2;          //600000A

  G4double OutCondLnkTotLength;
  G4double OutCondLnk_zpos;
  
  G4double HalfMotherVoluLength;
  G4double HalfTunnelLength; 

  G4double Tunnel_zpos;
  G4double totL;
  G4double Clearance;
  G4double HalfHornLength;
  G4double HornRadius;
  
  G4double current1;
  G4double current2;
  
  G4int numZPlanes;
  G4int numZPlanesOut;
  G4double zPlane[24];
  G4double rInner[24];
  G4double rOuter[24];
  G4double zPlaneOut[24];
  G4double rInnerOut[24];
  G4double rOuterOut[24];

  G4double zPlaneB[24];
  G4double rInnerB[24];
  G4double rOuterB[24];
  G4double zPlaneBOut[24];
  G4double rInnerBOut[24];
  G4double rOuterBOut[24];
  /*
  G4double rInnerOut1[24];
  G4double rOuterOut1[24];
  G4double rInner1[24];
  G4double rOuter1[24];
  */
  G4double yInnerOut[24];
  G4double yOuter[24];
  G4double xOut[24];
  G4double xIn[24];
  G4double yInner[24];
  G4double yOuterOut[24];
  //-----------------------------------//////////////////////////////////////modif
  // MINOS like
  //-----------------------------------
  G4double SkinDepth_MS1;
  G4double SkinDepth_MS2;
  G4int NMS1;
  G4double Z0_MS1;
  G4double zMS1[100],riMS1[100],roMS1[100];
  G4double zMS2[100],riMS2[100],roMS2[100];
  G4int NMS2;
  G4double Z0_MS2;
  G4double zMS1_skin[100],riMS1_skin[100],roMS1_skin[100];
  G4double zMS2_skin[100],riMS2_skin[100],roMS2_skin[100];
  //-----------------------------------             //////////////////////////////////////////////////


//....oooOO0OOooo....  Solids Construction  ....oooOO0OOooo....
  
  SBDetectorConstructionMessenger* detectorConstructionMessenger;

  G4VPhysicalVolume* Construct();

  G4Tubs*            solidWorld;
  G4LogicalVolume*   logicWorld;
  G4VPhysicalVolume* physiWorld;

  G4Tubs*            solidTunnel;
  G4LogicalVolume*   logicTunnel;
  G4VPhysicalVolume* physiTunnel;

  G4Polycone* solidInnerConductor;
  G4LogicalVolume* logicInnerConductor;
  G4VPhysicalVolume* physiInnerConductor[4];
  
  G4Polycone* solidVirtualOuterConductor;
  G4SubtractionSolid* solidOuterConductor;
  G4LogicalVolume* logicOuterConductor;
  G4VPhysicalVolume* physiOuterConductor[4];

  G4Polycone* solidBField1;
  G4LogicalVolume* logicBField1[4];
  G4VPhysicalVolume* physiBField1[4];

  G4Polycone* solidVirtualBField2;
  G4SubtractionSolid* solidBField2;
  G4LogicalVolume* logicBField2[4];
  G4VPhysicalVolume* physiBField2[4];
   
  G4Tubs* solidTarget;
  G4LogicalVolume* logicTarget;
  G4VPhysicalVolume* physiTarget[4];

  void DefineMaterials();
  G4VPhysicalVolume* ConstructSystem();    

  void SetVolumesVisibility();
  void DoMagField1(bool);
  void DoMagField2(bool);
  void DoInnerConductor(bool);
  void DoOuterConductor(bool);
    
  void DoTunnel(bool);
  void DoTarget(bool);

public:

  void SetHornParameters();
  
  void SetOffAxisAngle(G4double);
  void SetOffAxisPhi(G4double);
 
  G4double GetOffAxisAngle(){return OffAxisAngle;};
  G4double GetOffAxisPhi(){return OffAxisPhi;};
  
  G4double Horn_xpos;
  G4double Horn_ypos;
  G4double Horn_zpos;
  G4double Horn_rpos;

  void UpdateGeometry();  
  void PrintHornParameters(); 
                    
  G4double GetHallSizeX(){return 2.*rO_Hall;}; 
  G4double GetHallSizeY(){return 2.*rO_Hall;};
  G4double GetHallSizeZ(){return L_Hall;};

  G4double GetTunnelLength(){return TunnelLength;};
  G4double GetTunnelRadius(){return TunnelRadius;};

  void SetTunnelLength(G4double val){TunnelLength=val;};
  void SetTunnelRadius(G4double val){TunnelRadius=val;};
  void SetHornDistanceToXAxis(G4double val){Horn_xpos=val; G4cout << " Horn_xpos = " << Horn_xpos << G4endl; };
  void SetHornDistanceToYAxis(G4double val){Horn_ypos=val; G4cout << " Horn_ypos = " << Horn_ypos << G4endl; };
  void SetHornDistanceToZAxis(G4double val){Horn_zpos=val; G4cout << " Horn_zpos = " << Horn_zpos << G4endl; };
  void SetHornDistanceToRAxis(G4double val){Horn_rpos=val; G4cout << " Horn_rpos = " << Horn_rpos << G4endl; };
  void SetHorncurrent1(G4double val){current1=val;         G4cout << " current1 = " << current1 << G4endl;   };
  void SetHorncurrent2(G4double val){current2=val;         G4cout << " current2 = " << current2 << G4endl;   };
  G4double GetInnerConductorLength(){return 2.*InCondHalfLength;};
     
  G4double GetTargetThickness(){return TargetThickness;}; 
  G4double GetTargetDiameter(){return TargetDiameter;}; 
  G4Material* GetTargetMaterial(){return TargetMaterial;};
  G4double GetTargetZetaPos(){return TargetZetaPos;};
          
  const G4VPhysicalVolume* GetphysiWorld(){return physiWorld;}           
  const G4VPhysicalVolume* GetTarget(){return physiTarget[0];}
//to make it a ROOT class
  //ClassDef(SBDetectorConstruction,2)//SBDetectorConstruction
};
#endif
   
  
