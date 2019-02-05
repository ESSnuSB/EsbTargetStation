//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo
//
// SBDetectorconstructionNew ->  Detector construction.    
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo

#include "G4NistManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "SBDetectorConstruction.hh"
#include "SBDetectorConstructionMessenger.hh"
#include "QGSP_BERT.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Polycone.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4SubtractionSolid.hh"  

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"

#include "G4PropagatorInField.hh"
#include "G4ClassicalRK4.hh"
#include "SBField.hh"


#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ChordFinder.hh"

#include "G4Mag_UsualEqRhs.hh"

#include "TMath.h"

#include "Riostream.h"
#include "Rtypes.h"
#include "TROOT.h"
#include "TLine.h"
#include "TVirtualPad.h"
#include "TClass.h"
#include "TVirtualX.h"

#include "G4RunManager.hh"

//to make it a ROOT class
//ClassImp(SBDetectorConstruction)

void SBDetectorConstruction::SetOffAxisAngle(G4double val){OffAxisAngle=val;}
void SBDetectorConstruction::SetOffAxisPhi(G4double val){OffAxisPhi=val;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo

SBDetectorConstruction::SBDetectorConstruction():TargetMaterial(0),defaultMaterial(0),solidWorld(0),logicWorld(0),physiWorld(0)
{

  // Default parameter values     (target +tunnel)
  // --
  TargetThickness=780.*mm;
  TargetDiameter=11*mm;
  TargetHalfThickness=(TargetThickness*0.5); 

  TunnelLength=40000.*mm;
  HalfTunnelLength=(TunnelLength*0.5);
  TunnelRadius=2000.*mm;

  // Default parameter for the horn position.
  // ----------------------------------------
  Horn_xpos = 1000.0*mm;
  Horn_ypos = 1000.0*mm;
  Horn_zpos = -0.5*L_Hall+2.*HalfMotherVoluLength-HalfHornLength;
  Horn_rpos = sqrt((Horn_xpos*Horn_xpos)+(Horn_ypos*Horn_ypos));

  current1 = -300000.*ampere;
  current2 = -300000.*ampere;
  
  // Materials
  // ---------
  DefineMaterials();
  
  // create commands for interactive definition
  // ------------------------------------------
  detectorConstructionMessenger=new SBDetectorConstructionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo

SBDetectorConstruction::~SBDetectorConstruction()
{
  delete detectorConstructionMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo

G4VPhysicalVolume* SBDetectorConstruction::Construct()
{
  return ConstructSystem();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo

G4VPhysicalVolume* SBDetectorConstruction::ConstructSystem()
{
  
  bool placeInnerConductor=true;
  bool placeBfield1=true;
  bool placeOuterConductor=true;
  bool placeBfield2=true;
  bool placeTarget=true;
  bool placeTunnel=true;
  
  placeInnerConductor=true;
  placeBfield1=true;
  placeOuterConductor=true;
  placeBfield2=true;
  placeTarget=true;
  placeTunnel=true;
  
  // Clean old geometry, if any
  // --
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // Define Horn Parmeter
  // --
  SetHornParameters();

  // WorldVolume
  // --
  solidWorld=new G4Tubs("World",0.*cm,rO_Hall,L_Hall/2.,0.*deg,360.*deg);               
  logicWorld=new G4LogicalVolume(solidWorld,defaultMaterial,"World");
  physiWorld=new G4PVPlacement(0,G4ThreeVector(),logicWorld,"World",0,false,0);
  G4cout << physiWorld->GetName() <<" created."<<G4endl;

  // Contruct Tunnel
  // --
  DoTunnel(placeTunnel);

  // Construct Horn
  // --
  DoInnerConductor(placeInnerConductor);
  DoOuterConductor(placeOuterConductor);
  DoMagField1(placeBfield1);
  DoMagField2(placeBfield2);
  DoTarget(placeTarget);

  // Set Volume Visibility
  // --
  SetVolumesVisibility();

  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo

void SBDetectorConstruction::UpdateGeometry()
{
  G4cout << "SBDetectorConstruction::UpdateGeometry "<< G4endl;  
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructSystem());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo

void SBDetectorConstruction::SetHornParameters()
{
  G4cout << "SBDetectorConstruction::SetHornParameters "<< G4endl; 
 
  // Target
  // --
  TargetRadius = TargetDiameter/2.;

  // Variables for the Horn
  // --
  L1 = 589.0*mm;
  L2 = 468.0*mm;
  L3 = 603.0*mm;
  L4 = 475.0*mm;
  L5 = 10.8*mm;
  t1 = 3.0*mm;
  t2 = 3.0*mm;
  t3 = 3.0*mm;
  t4 = 10.0*mm;
   R = 12.0*mm;  
  r3 = 50.8*mm; 
  R4 = 272.0*mm;          
  R2 = 191.0*mm;  
  R3 = 359.0*mm;
  totL = t2+r3+L1+L2+L3+L4+L5+R4+t3;      // 2474,6*mm
  OutCondHalfLength = totL*0.5;           // 1237,3*mm 
  InCondHalfLength = (totL-t2-t3)*0.5;    // 1234,3*mm
 
  Clearance = 3000.*mm;
  HalfHornLength = OutCondHalfLength;
  G4cout <<"HalfHornLength= "<<HalfHornLength<<" mm"<<G4endl;
  rO_Hall=4000.*mm;    
 
  HalfMotherVoluLength = HalfHornLength+Clearance/2.;
  G4cout <<"HalfMotherVoluLength= "<<HalfMotherVoluLength<<" mm"<<G4endl;

  L_Hall = TunnelLength+2.*HalfMotherVoluLength;
  G4cout <<"L_Hall= "<<L_Hall<<" mm"<<G4endl; 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo
void SBDetectorConstruction::DefineMaterials()
{
  G4cout << "SBDetectorConstruction::DefineMaterials "<< G4endl;  
  G4String symbol;
  G4double a, z, density;
  G4int ncomponents, natoms;
  G4double fractionmass;
  
  G4Element* H  = new G4Element("Hydrogen",symbol="H" , z= 1., a= 1.01*g/mole);
  G4Element* C  = new G4Element("Carbon"  ,symbol="C" , z= 6., a= 12.01*g/mole);
  G4Element* N  = new G4Element("Nitrogen",symbol="N" , z= 7., a= 14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"  ,symbol="O" , z= 8., a= 16.00*g/mole);
  G4Element* Si = new G4Element("Silicon",symbol="Si" , z= 14., a= 28.09*g/mole);
  G4Element* Al = new G4Element("Aluminum",symbol="Al" , z= 13., a= 26.9815386*g/mole);
  G4Element* Be = new G4Element("Berillium",symbol="Be" , z= 4., a= 9.0121823*g/mole);

  G4Material* Aluminum=new G4Material("Aluminum",z=13.,a=26.98*g/mole,density=2.700*g/cm3);

  G4Material* myGraphite=new G4Material("myGraphite",z=6.,a=12.01*g/mole,density=1.85*g/cm3);
  //G4Material* myTantalum=new G4Material("myTantalum",z=73.,a=180.94788*g/mole,density=16.69*g/cm3);

  G4Material* H2O=new G4Material("Water", density= 1.000*g/cm3, ncomponents=2);
  H2O->AddElement(H,natoms=2);
  H2O->AddElement(O,natoms=1); 
  H2O->GetIonisation()->SetMeanExcitationEnergy(75.0*eV);
  
  G4Material* Sci=new G4Material("Scintillator", density= 1.032*g/cm3, ncomponents=2);
  Sci->AddElement(C,natoms=9);
  Sci->AddElement(H,natoms=10);
  
  G4Material* Myl=new G4Material("Mylar", density= 1.397*g/cm3, ncomponents=3);
  Myl->AddElement(C,natoms=10);
  Myl->AddElement(H,natoms=8);
  Myl->AddElement(O,natoms=4);
  
  G4Material* SiO2=new G4Material("quartz",density= 2.200*g/cm3, ncomponents=2);
  SiO2->AddElement(Si,natoms=1);
  SiO2->AddElement(O,natoms=2);
  
  G4Material* Air=new G4Material("Air",density= 1.290*mg/cm3,ncomponents=2);
  Air->AddElement(N,fractionmass=0.7);
  Air->AddElement(O,fractionmass=0.3);
 
  G4Material* AlBeMet=new G4Material("AlBeMet",density=0.2071*g/cm3,ncomponents=2);
  AlBeMet->AddElement(Al,fractionmass=38.*perCent);
  AlBeMet->AddElement(Be,fractionmass=62.*perCent);
  
  G4Material* CO2=new G4Material("CarbonicGas", density= 1.842*mg/cm3, ncomponents=2,kStateGas, 325.*kelvin, 50.*atmosphere);
  CO2->AddElement(C,natoms=1);
  CO2->AddElement(O,natoms=2);
 
  G4Material* steam=new G4Material("WaterSteam",density=0.3*mg/cm3,ncomponents=1,kStateGas,500.*kelvin,2.*atmosphere);
  steam->AddMaterial(H2O,fractionmass=1.);

  G4Material* Vacuum=new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,kStateGas, 2.73*kelvin, 3.e-18*pascal);

  G4Material* beam=new G4Material("Beam", density= 1.e-5*g/cm3, ncomponents=1,kStateGas, STP_Temperature, 2.e-2*bar);

  beam->AddMaterial(Air, fractionmass=1.);
  G4Material* AirTunnel=new G4Material("AirTunnel",density= 0.00154*g/cm3,ncomponents=1,kStateGas,STP_Temperature,0.001238*atmosphere); 
  AirTunnel->AddMaterial(Air,fractionmass=1.);

  if(0)G4cout << *(G4Material::GetMaterialTable()) << G4endl;
 
  defaultMaterial=Air;
  HornMaterial=Aluminum;
  TargetMaterial=myGraphite;
 
  //HornMaterial=Vacuum; 
  //TargetMaterial=Vacuum;
  
  FillingMaterial=Vacuum;
  TunnelMaterial=AirTunnel;
  //TunnelMaterial=Vacuum;
  //TunnelMaterial=Air;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo

void SBDetectorConstruction::DoInnerConductor(bool place)
{
  G4cout << "SBDetectorConstruction::DoInnerConductor "<< G4endl; 
  // Warning:l'origine du systeme de coordonnées globales est au centre du world volume 
  // -- 
  numZPlanes=23;

 // Initialization
 // --
 for(int i=0;i<24;i++){
   zPlane[i]=0;
   rInner[i]=0;
   rOuter[i]=0;

   xIn[i]=0;
   yOuter[i]=0;
  }

 // Z_plane
 // --
 zPlane[0]=-InCondHalfLength;
 for(int i=1;i<=8;i++){
  zPlane[i]=zPlane[i-1]+(r3/8);
  }
 zPlane[9]=-(InCondHalfLength-r3-L1);
 zPlane[10]=-(InCondHalfLength-r3-L1-L2);  ////modifié
 zPlane[11]=(InCondHalfLength-L4-L5-R4);
 zPlane[12]=(InCondHalfLength-L5-R4);
 zPlane[13]=(InCondHalfLength-R4);

 for(int i=14;i<=23;i++){
  zPlane[i]=zPlane[i-1]+(R4)/10;
 }

 G4cout <<"zPlane[0]= "<<zPlane[0]<<G4endl;
 G4cout <<"zPlane[1]= "<<zPlane[1]<<G4endl;
 G4cout <<"zPlane[2]= "<<zPlane[2]<<G4endl;
 G4cout <<"zPlane[3]= "<<zPlane[3]<<G4endl;
 G4cout <<"zPlane[4]= "<<zPlane[4]<<G4endl;
 G4cout <<"zPlane[5]= "<<zPlane[5]<<G4endl;
 G4cout <<"zPlane[6]= "<<zPlane[6]<<G4endl;
 G4cout <<"zPlane[7]= "<<zPlane[7]<<G4endl;
 G4cout <<"zPlane[8]= "<<zPlane[8]<<G4endl;
 G4cout <<"zPlane[9]= "<<zPlane[9]<<G4endl;
 G4cout <<"zPlane[10]= "<<zPlane[10]<<G4endl;
 G4cout <<"zPlane[11]= "<<zPlane[11]<<G4endl;
 G4cout <<"zPlane[12]= "<<zPlane[12]<<G4endl;
 G4cout <<"zPlane[13]= "<<zPlane[13]<<G4endl;
 G4cout <<"zPlane[14]= "<<zPlane[14]<<G4endl;
 G4cout <<"zPlane[15]= "<<zPlane[15]<<G4endl;
 G4cout <<"zPlane[16]= "<<zPlane[16]<<G4endl;
 G4cout <<"zPlane[17]= "<<zPlane[17]<<G4endl;
 G4cout <<"zPlane[18]= "<<zPlane[18]<<G4endl;
 G4cout <<"zPlane[19]= "<<zPlane[19]<<G4endl;
 G4cout <<"zPlane[20]= "<<zPlane[20]<<G4endl;
 G4cout <<"zPlane[21]= "<<zPlane[21]<<G4endl;
 G4cout <<"zPlane[22]= "<<zPlane[22]<<G4endl;
 G4cout <<"zPlane[23]= "<<zPlane[23]<<G4endl;

 G4cout <<"/////////////////////////////////////////////////////////////"<<G4endl;
 
 // R_inner 
 // --
 rInner[0]=t1+r3+R;
 for(int i=1;i<=8;i++){
  rInner[i]=rInner[i-1]-(r3)/8;
 }
 rInner[9]=R+t1;
 rInner[10]=R+R2+t1;
 rInner[11]=R+R2+t1;
 rInner[12]=R+t1;
 rInner[13]=R+t1;

 //for(int i=14;i<=23;i++){
 //rInner[i]=rInner[i-1]+(R2-t3)/10;
 //}
 
 /////////////////////////////////////////////////////////////////////////////////
 //rInner1[13]=R+t1;
 
 xIn[13]=(R4)/10;
 
 for(int i=14;i<=23;i++){
 xIn[i]=xIn[i-1]+(R4)/10;
 }
 for(int i=13;i<=23;i++){
 yInner[i]=sqrt((1-(pow(xIn[i],2)/pow(R4,2)))*pow(R2-t1,2));
 }

 for(int i=14;i<=23;i++){
 rInner[i]=R+t1+(R2-t1-yInner[i-1])-0.5;
 }
 /////////////////////////////////////////////////////////////////////////////////
 /////////////////////////////////////////////////////////////////////////////////
 //rOuter1[13]=R+R2+R3;
  
 for(int i=13;i<=23;i++){
 yOuter[i]=sqrt((1-(pow(xIn[i],2)/pow(R4,2)))*pow(R3,2));
 }

 for(int i=14;i<=23;i++){
 rOuter[i]=R+R2+R3-(R3-yOuter[i-1]);
 }
 /////////////////////////////////////////////////////////////////////////////////

 G4cout <<"rInner[0]= "<<rInner[0]<<G4endl;
 G4cout <<"rInner[1]= "<<rInner[1]<<G4endl;
 G4cout <<"rInner[2]= "<<rInner[2]<<G4endl;
 G4cout <<"rInner[3]= "<<rInner[3]<<G4endl;
 G4cout <<"rInner[4]= "<<rInner[4]<<G4endl;
 G4cout <<"rInner[5]= "<<rInner[5]<<G4endl;
 G4cout <<"rInner[6]= "<<rInner[6]<<G4endl;
 G4cout <<"rInner[7]= "<<rInner[7]<<G4endl;
 G4cout <<"rInner[8]= "<<rInner[8]<<G4endl;
 G4cout <<"rInner[9]= "<<rInner[9]<<G4endl;
 G4cout <<"rInner[10]= "<<rInner[10]<<G4endl;
 G4cout <<"rInner[11]= "<<rInner[11]<<G4endl;
 G4cout <<"rInner[12]= "<<rInner[12]<<G4endl;
 G4cout <<"rInner[13]= "<<rInner[13]<<G4endl;
 G4cout <<"rInner[14]= "<<rInner[14]<<G4endl;
 G4cout <<"rInner[15]= "<<rInner[15]<<G4endl;
 G4cout <<"rInner[16]= "<<rInner[16]<<G4endl;
 G4cout <<"rInner[17]= "<<rInner[17]<<G4endl;
 G4cout <<"rInner[18]= "<<rInner[18]<<G4endl;
 G4cout <<"rInner[19]= "<<rInner[19]<<G4endl;
 G4cout <<"rInner[20]= "<<rInner[20]<<G4endl;
 G4cout <<"rInner[21]= "<<rInner[21]<<G4endl;
 G4cout <<"rInner[22]= "<<rInner[22]<<G4endl;
 G4cout <<"rInner[23]= "<<rInner[23]<<G4endl;
 
 G4cout <<"/////////////////////////////////////////////////////////////"<<G4endl;

 // R_outer
 // --
 rOuter[0]=R+R2+R3-r3;  ///rOuterOut[0]=R+R2+R3-r3;
 for(int i=1;i<=8;i++){
  rOuter[i]=rOuter[i-1]+(r3)/8;
  }
 for(int i=9;i<=13;i++){
  rOuter[i]=R+R2+R3;   ///
  }
 //for(int i=14;i<=23;i++){
  //rOuter[i]=rOuter[i-1]-(R3)/10;   ///
  //}

/////////////////////////////////////////////////////////////////////////////////
 //rOuter1[13]=R+R2+R3;
  
 for(int i=13;i<=23;i++){
 yOuter[i]=sqrt((1-(pow(xIn[i],2)/pow(R4,2)))*pow(R3,2));
 }

 for(int i=14;i<=23;i++){
 rOuter[i]=R+R2+R3-(R3-yOuter[i-1]);
 }

 /////////////////////////////////////////////////////////////////////////////////
 G4cout <<"rOuter[0]= "<<rOuter[0]<<G4endl;
 G4cout <<"rOuter[1]= "<<rOuter[1]<<G4endl;
 G4cout <<"rOuter[2]= "<<rOuter[2]<<G4endl;
 G4cout <<"rOuter[3]= "<<rOuter[3]<<G4endl;
 G4cout <<"rOuter[4]= "<<rOuter[4]<<G4endl;
 G4cout <<"rOuter[5]= "<<rOuter[5]<<G4endl;
 G4cout <<"rOuter[6]= "<<rOuter[6]<<G4endl;
 G4cout <<"rOuter[7]= "<<rOuter[7]<<G4endl;
 G4cout <<"rOuter[8]= "<<rOuter[8]<<G4endl;
 G4cout <<"rOuter[9]= "<<rOuter[9]<<G4endl;
 G4cout <<"rOuter[10]= "<<rOuter[10]<<G4endl;
 G4cout <<"rOuter[11]= "<<rOuter[11]<<G4endl;
 G4cout <<"rOuter[12]= "<<rOuter[12]<<G4endl;
 G4cout <<"rOuter[13]= "<<rOuter[13]<<G4endl;
 G4cout <<"rOuter[14]= "<<rOuter[14]<<G4endl;
 G4cout <<"rOuter[15]= "<<rOuter[15]<<G4endl;
 G4cout <<"rOuter[16]= "<<rOuter[16]<<G4endl;
 G4cout <<"rOuter[17]= "<<rOuter[17]<<G4endl;
 G4cout <<"rOuter[18]= "<<rOuter[18]<<G4endl;
 G4cout <<"rOuter[19]= "<<rOuter[19]<<G4endl;
 G4cout <<"rOuter[20]= "<<rOuter[20]<<G4endl;
 G4cout <<"rOuter[21]= "<<rOuter[21]<<G4endl;
 G4cout <<"rOuter[22]= "<<rOuter[22]<<G4endl;
 G4cout <<"rOuter[23]= "<<rOuter[23]<<G4endl;
 G4cout <<"/////////////////////////////////////////////////////////////"<<G4endl;

 G4cout << " Before " << G4endl;
 solidInnerConductor=new G4Polycone("INCO",0.*deg,360.*deg,numZPlanes=23,zPlane,rInner,rOuter);
 G4cout << " After  " << G4endl;
 logicInnerConductor=new G4LogicalVolume(solidInnerConductor,HornMaterial,"INCO");

 if(place){
 physiInnerConductor[0]=new G4PVPlacement(0,G4ThreeVector(Horn_xpos,Horn_ypos,Horn_zpos),logicInnerConductor,"INCO",logicWorld,false,0); 
 if(place)G4cout << " created.physiInnerConductor[0]"<<G4endl;
 physiInnerConductor[1]=new G4PVPlacement(0,G4ThreeVector(-Horn_xpos,Horn_ypos,Horn_zpos),logicInnerConductor,"INCO",logicWorld,false,0); 
 if(place)G4cout << " created.physiInnerConductor[1]"<<G4endl;
 physiInnerConductor[2]=new G4PVPlacement(0,G4ThreeVector(Horn_xpos,-Horn_ypos,Horn_zpos),logicInnerConductor,"INCO",logicWorld,false,0); 
 if(place)G4cout << " created.physiInnerConductor[2]"<<G4endl;
 physiInnerConductor[3]=new G4PVPlacement(0,G4ThreeVector(-Horn_xpos,-Horn_ypos,Horn_zpos),logicInnerConductor,"INCO",logicWorld,false,0); 
 if(place)G4cout << " created.physiInnerConductor[3]"<<G4endl;
 }
 G4cout <<"/////////////////////////////////////////////////////////////"<<G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo

void SBDetectorConstruction::DoOuterConductor(bool place)
{
  G4cout << "SBDetectorConstruction::DoOuterConductor "<< G4endl; 
/////warning:l'origine du systeme de coordonnées globales est au centre du world volume 
  numZPlanesOut=23;
//initialization
  for(int i=0;i<23;i++){
    zPlaneOut[i]=0;
    rInnerOut[i]=0;
    rOuterOut[i]=0;
    
    //rInnerOut1[i]=0;
    //rOuterOut1[i]=0;
    yInnerOut[i]=0;
    xOut[i]=0;
    yOuterOut[i]=0;
    }
//Z_plane
 zPlaneOut[0]=-OutCondHalfLength;
 for(int i=1;i<=8;i++){
  zPlaneOut[i]=zPlaneOut[i-1]+(t2+r3)/8;
  }
 zPlaneOut[9]=-(OutCondHalfLength-L1-r3-t2);
 zPlaneOut[10]=-(OutCondHalfLength-r3-L1-L2-t2);       ////modifié
 zPlaneOut[11]=(OutCondHalfLength-L4-L5-R4-t3); //(InCondHalfLength-L4-L5-R4)
 zPlaneOut[12]=(OutCondHalfLength-L5-R4-t3);
 zPlaneOut[13]=(OutCondHalfLength-R4-t3);
 for(int i=14;i<=23;i++){
  zPlaneOut[i]=zPlaneOut[i-1]+(R4+t3)/10;
 }

G4cout <<"zPlaneOut[0]= "<<zPlaneOut[0]<<G4endl;
 G4cout <<"zPlaneOut[1]= "<<zPlaneOut[1]<<G4endl;
 G4cout <<"zPlaneOut[2]= "<<zPlaneOut[2]<<G4endl;
 G4cout <<"zPlaneOut[3]= "<<zPlaneOut[3]<<G4endl;
 G4cout <<"zPlaneOut[4]= "<<zPlaneOut[4]<<G4endl;
 G4cout <<"zPlaneOut[5]= "<<zPlaneOut[5]<<G4endl;
 G4cout <<"zPlaneOut[6]= "<<zPlaneOut[6]<<G4endl;
 G4cout <<"zPlaneOut[7]= "<<zPlaneOut[7]<<G4endl;
 G4cout <<"zPlaneOut[8]= "<<zPlaneOut[8]<<G4endl;
 G4cout <<"zPlaneOut[9]= "<<zPlaneOut[9]<<G4endl;
 G4cout <<"zPlaneOut[10]= "<<zPlaneOut[10]<<G4endl;
 G4cout <<"zPlaneOut[11]= "<<zPlaneOut[11]<<G4endl;
 G4cout <<"zPlaneOut[12]= "<<zPlaneOut[12]<<G4endl;
 G4cout <<"zPlaneOut[13]= "<<zPlaneOut[13]<<G4endl;
 G4cout <<"zPlaneOut[14]= "<<zPlaneOut[14]<<G4endl;
 G4cout <<"zPlaneOut[15]= "<<zPlaneOut[15]<<G4endl;
 G4cout <<"zPlaneOut[16]= "<<zPlaneOut[16]<<G4endl;
 G4cout <<"zPlaneOut[17]= "<<zPlaneOut[17]<<G4endl;
 G4cout <<"zPlaneOut[18]= "<<zPlaneOut[18]<<G4endl;
 G4cout <<"zPlaneOut[19]= "<<zPlaneOut[19]<<G4endl;
 G4cout <<"zPlaneOut[20]= "<<zPlaneOut[20]<<G4endl;
 G4cout <<"zPlaneOut[21]= "<<zPlaneOut[21]<<G4endl;
 G4cout <<"zPlaneOut[22]= "<<zPlaneOut[22]<<G4endl;
 G4cout <<"zPlaneOut[23]= "<<zPlaneOut[23]<<G4endl;

G4cout <<"/////////////////////////////////////////////////////////////"<<G4endl;
//R_inner 
 rInnerOut[0]=t1+r3+R;
 for(int i=1;i<=8;i++){
  rInnerOut[i]=rInnerOut[i-1]-(t2+r3)/8;
  }
 rInnerOut[9]=R;
 rInnerOut[10]=R+R2;
 rInnerOut[11]=R+R2;
 rInnerOut[12]=R;
 rInnerOut[13]=R;
 //for(int i=14;i<=23;i++){
  //rInnerOut[i]=rInnerOut[i-1]+(R2)/10;
  //}

 /////////////////////////////////////////////////////////////////////////////////
 //rInnerOut1[13]=R;
 
 xOut[13]=(R4+t3)/10;
 
 for(int i=14;i<=23;i++){
  xOut[i]=xOut[i-1]+(R4+t3)/10;
 }

 for(int i=13;i<=23;i++){
  yInnerOut[i]=sqrt((1-(pow(xOut[i],2)/pow(R4+t3,2)))*pow(R2,2));
 }

 for(int i=14;i<=23;i++){
  rInnerOut[i]=R+(R2-yInnerOut[i-1])-0.5;        //R2-t1 
 }
 /////////////////////////////////////////////////////////////////////////////////
 
 G4cout <<"rInnerOut[0]= "<<rInnerOut[0]<<G4endl;
 G4cout <<"rInnerOut[1]= "<<rInnerOut[1]<<G4endl;
 G4cout <<"rInnerOut[2]= "<<rInnerOut[2]<<G4endl;
 G4cout <<"rInnerOut[3]= "<<rInnerOut[3]<<G4endl;
 G4cout <<"rInnerOut[4]= "<<rInnerOut[4]<<G4endl;
 G4cout <<"rInnerOut[5]= "<<rInnerOut[5]<<G4endl;
 G4cout <<"rInnerOut[6]= "<<rInnerOut[6]<<G4endl;
 G4cout <<"rInnerOut[7]= "<<rInnerOut[7]<<G4endl;
 G4cout <<"rInnerOut[8]= "<<rInnerOut[8]<<G4endl;
 G4cout <<"rInnerOut[9]= "<<rInnerOut[9]<<G4endl;
 G4cout <<"rInnerOut[10]= "<<rInnerOut[10]<<G4endl;
 G4cout <<"rInnerOut[11]= "<<rInnerOut[11]<<G4endl;
 G4cout <<"rInnerOut[12]= "<<rInnerOut[12]<<G4endl;
 G4cout <<"rInnerOut[13]= "<<rInnerOut[13]<<G4endl;
 G4cout <<"rInnerOut[14]= "<<rInnerOut[14]<<G4endl;
 G4cout <<"rInnerOut[15]= "<<rInnerOut[15]<<G4endl;
 G4cout <<"rInnerOut[16]= "<<rInnerOut[16]<<G4endl;
 G4cout <<"rInnerOut[17]= "<<rInnerOut[17]<<G4endl;
 G4cout <<"rInnerOut[18]= "<<rInnerOut[18]<<G4endl;
 G4cout <<"rInnerOut[19]= "<<rInnerOut[19]<<G4endl;
 G4cout <<"rInnerOut[20]= "<<rInnerOut[20]<<G4endl;
 G4cout <<"rInnerOut[21]= "<<rInnerOut[21]<<G4endl;
 G4cout <<"rInnerOut[22]= "<<rInnerOut[22]<<G4endl;
 G4cout <<"rInnerOut[23]= "<<rInnerOut[23]<<G4endl;
 G4cout <<"/////////////////////////////////////////////////////////////" << G4endl;

 rOuterOut[0]=R+R2+R3-r3;        //rOuter[0]=R+R2+R3-r3;  //////////pb =511.2
 for(int i=1;i<=7;i++){
  rOuterOut[i]=rOuterOut[i-1]+(r3+t4)/8;     //rOuter[i]=rOuter[i-1]+(r3)/8;
 }

 rOuterOut[8]=R+R2+R3+t4;
 for(int i=9;i<=13;i++){
  rOuterOut[i]=R+R2+R3+t4;      //rOuter[i]=R+R2+R3-t4;
 }

 for(int i=13;i<=23;i++){
  yOuterOut[i]=sqrt((1-(pow(xOut[i],2)/pow(R4+t3,2)))*pow(R3+t4,2));
 }

 for(int i=14;i<=23;i++){
  rOuterOut[i]=R+R2+R3+t4-(R3+t4-yOuterOut[i-1]);
 }

 /////////////////////////////////////////////////////////////////////////////////
 
 G4cout <<"rOuterOut[0]= "<<rOuterOut[0]<<G4endl;
 G4cout <<"rOuterOut[1]= "<<rOuterOut[1]<<G4endl;
 G4cout <<"rOuterOut[2]= "<<rOuterOut[2]<<G4endl;
 G4cout <<"rOuterOut[3]= "<<rOuterOut[3]<<G4endl;
 G4cout <<"rOuterOut[4]= "<<rOuterOut[4]<<G4endl;
 G4cout <<"rOuterOut[5]= "<<rOuterOut[5]<<G4endl;
 G4cout <<"rOuterOut[6]= "<<rOuterOut[6]<<G4endl;
 G4cout <<"rOuterOut[7]= "<<rOuterOut[7]<<G4endl;
 G4cout <<"rOuterOut[8]= "<<rOuterOut[8]<<G4endl;
 G4cout <<"rOuterOut[9]= "<<rOuterOut[9]<<G4endl;
 G4cout <<"rOuterOut[10]= "<<rOuterOut[10]<<G4endl;
 G4cout <<"rOuterOut[11]= "<<rOuterOut[11]<<G4endl;
 G4cout <<"rOuterOut[12]= "<<rOuterOut[12]<<G4endl;
 G4cout <<"rOuterOut[13]= "<<rOuterOut[13]<<G4endl;
 G4cout <<"rOuterOut[14]= "<<rOuterOut[14]<<G4endl;
 G4cout <<"rOuterOut[15]= "<<rOuterOut[15]<<G4endl;
 G4cout <<"rOuterOut[16]= "<<rOuterOut[16]<<G4endl;
 G4cout <<"rOuterOut[17]= "<<rOuterOut[17]<<G4endl;
 G4cout <<"rOuterOut[18]= "<<rOuterOut[18]<<G4endl;
 G4cout <<"rOuterOut[19]= "<<rOuterOut[19]<<G4endl;
 G4cout <<"rOuterOut[20]= "<<rOuterOut[20]<<G4endl;
 G4cout <<"rOuterOut[21]= "<<rOuterOut[21]<<G4endl;
 G4cout <<"rOuterOut[22]= "<<rOuterOut[22]<<G4endl;
 G4cout <<"rOuterOut[23]= "<<rOuterOut[23]<<G4endl;
 G4cout <<"/////////////////////////////////////////////////////////////"<<G4endl;

 G4cout << " Before " << G4endl;
 solidVirtualOuterConductor=new G4Polycone("OUTC",0.*deg,360.*deg,numZPlanesOut=23,zPlaneOut,rInnerOut,rOuterOut); 
 solidOuterConductor=new G4SubtractionSolid("OUTC-INCO",solidVirtualOuterConductor,solidInnerConductor,0,G4ThreeVector(0,0,0));
 logicOuterConductor=new G4LogicalVolume(solidOuterConductor,HornMaterial,"OUTC-INCO"); //////à verifier
 G4cout << " After  " << G4endl; 

 if(place){
 physiOuterConductor[0]=new G4PVPlacement(0,G4ThreeVector(Horn_xpos,Horn_ypos,Horn_zpos),logicOuterConductor,"OUTC-INCO",logicWorld,false,0); 
 if(place)G4cout <<" created.physiOuterConductor[0]"<<G4endl;
 physiOuterConductor[1]=new G4PVPlacement(0,G4ThreeVector(-Horn_xpos,Horn_ypos,Horn_zpos),logicOuterConductor,"OUTC-INCO",logicWorld,false,0); 
 if(place)G4cout <<" created.physiOuterConductor[1]"<<G4endl;
 physiOuterConductor[2]=new G4PVPlacement(0,G4ThreeVector(Horn_xpos,-Horn_ypos,Horn_zpos),logicOuterConductor,"OUTC-INCO",logicWorld,false,0); 
 if(place)G4cout <<" created.physiOuterConductor[2]"<<G4endl;
 physiOuterConductor[3]=new G4PVPlacement(0,G4ThreeVector(-Horn_xpos,-Horn_ypos,Horn_zpos),logicOuterConductor,"OUTC-INCO",logicWorld,false,0); 
 if(place)G4cout <<" created.physiOuterConductor[0]"<<G4endl;
 }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo

void SBDetectorConstruction::DoTarget(bool place)
{
  G4cout << "SBDetectorConstruction::DoTarget "<< G4endl;           
  G4cout << " Target Thickness " << TargetThickness << G4endl;
  TargetZetaPos=Horn_zpos-OutCondHalfLength+TargetHalfThickness;         /////position cible ok
  solidTarget=new G4Tubs("TARG",0.*cm,TargetDiameter/2.,TargetThickness/2.,0.*deg,360.*deg);
  logicTarget=new G4LogicalVolume(solidTarget,TargetMaterial,"TARG");

  if(place){
     physiTarget[0]=new G4PVPlacement(0,G4ThreeVector(Horn_xpos,Horn_ypos,TargetZetaPos),logicTarget,"TARG",logicWorld,false,0);
     G4cout << physiTarget[0]->GetName() <<"[0] created."<<G4endl;
     physiTarget[1]=new G4PVPlacement(0,G4ThreeVector(Horn_xpos,-Horn_ypos,TargetZetaPos),logicTarget,"TARG",logicWorld,false,0);
     G4cout << physiTarget[1]->GetName() <<"[1] created."<<G4endl;
     physiTarget[2]=new G4PVPlacement(0,G4ThreeVector(-Horn_xpos,Horn_ypos,TargetZetaPos),logicTarget,"TARG",logicWorld,false,0);
     G4cout << physiTarget[2]->GetName() <<"[2] created."<<G4endl;
     physiTarget[3]=new G4PVPlacement(0,G4ThreeVector(-Horn_xpos,-Horn_ypos,TargetZetaPos),logicTarget,"TARG",logicWorld,false,0);
     G4cout << physiTarget[3]->GetName() <<"[3] created."<<G4endl;
     }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo

void SBDetectorConstruction::DoTunnel(bool place)
{
  G4cout << "SBDetectorConstruction::DoTunnel "<< G4endl;  
  Tunnel_zpos=HalfTunnelLength+(Horn_zpos+HalfHornLength)+500;
  G4cout <<"Tunnel_zpos= "<<Tunnel_zpos<<G4endl;
  solidTunnel=new G4Tubs("TUNL",0.*cm,TunnelRadius,0.5*TunnelLength,0.*deg,360.*deg);
  logicTunnel=new G4LogicalVolume(solidTunnel,TunnelMaterial,"TUNL");
  if(place){
    physiTunnel=new G4PVPlacement(0,G4ThreeVector(0,0,Tunnel_zpos),logicTunnel,"TUNL",logicWorld,false,0);
    G4cout << physiTunnel->GetName() <<" created."<<G4endl;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo

void SBDetectorConstruction::DoMagField1(bool place)
{
 //for the inner conductor
G4cout <<"/////////////////////////////////////////////////////////////"<<G4endl;
 G4cout << "SBDetectorConstruction::DoMagField1 "<< G4endl;  
 int HornNumber=4;
 //initialization
 for(int i=0;i<24;i++){
     rInnerB[i]=0;
     rOuterB[i]=0;
     zPlaneB[i]=0;    
     }
//Z_plane,R_inner,R_outer
 for(int i=0;i<24;i++){ 
    zPlaneB[i]=zPlane[i];
    rInnerB[i]=rInner[i];
    rOuterB[i]=rOuter[i];
    }

 G4double stepMinimum1=1.0e-3*mm;
 SBField* myField1[4];
 G4Mag_UsualEqRhs* iEquation1[4];
 G4MagIntegratorStepper* iStepper1[4];
 G4ChordFinder* iChordFinder1[4];
 G4FieldManager* mfieldMgr1[4];

 for(int i=0;i<HornNumber;i++){
    myField1[i]=new SBField(current1,i,this);                 /////////////////SBField redefinir val magB
    iEquation1[i]=new G4Mag_UsualEqRhs(myField1[i]);
    iStepper1[i]=new G4ClassicalRK4(iEquation1[i]);
    iChordFinder1[i]=new G4ChordFinder(myField1[i],stepMinimum1,iStepper1[i]);
    //mfieldMgr1[i]=new G4FieldManager(myField1[i],iChordFinder1[i]);
    mfieldMgr1[i]=new G4FieldManager(myField1[i],iChordFinder1[i],true);
    }
 G4cout << " Before " << G4endl;
 solidBField1=new G4Polycone("BFL1",0.*deg,360.*deg,numZPlanes=23,zPlaneB,rInnerB,rOuterB);
 G4cout << " After " << G4endl;
 for(int i=0;i<HornNumber;i++){
    logicBField1[i]=new G4LogicalVolume(solidBField1,FillingMaterial,Form("BFL1_%d",i),mfieldMgr1[i],0,0);
    }

 if(place){
    physiBField1[0]=new G4PVPlacement(0,G4ThreeVector(Horn_xpos,Horn_ypos,Horn_zpos),logicBField1[0],"BFL1_0",logicWorld,false,0);
    G4cout << physiBField1[0]->GetName() <<"[0] created."<<G4endl;
    physiBField1[1]=new G4PVPlacement(0,G4ThreeVector(-Horn_xpos,Horn_ypos,Horn_zpos),logicBField1[1],"BFL1_1",logicWorld,false,0); //logicBField1[1]
    G4cout << physiBField1[1]->GetName() <<"[1] created."<<G4endl;
    physiBField1[2]=new G4PVPlacement(0,G4ThreeVector(Horn_xpos,-Horn_ypos,Horn_zpos),logicBField1[2],"BFL1_2",logicWorld,false,0);  //logicBField1[2]
    G4cout << physiBField1[2]->GetName() <<"[2] created."<<G4endl;
    physiBField1[3]=new G4PVPlacement(0,G4ThreeVector(-Horn_xpos,-Horn_ypos,Horn_zpos),logicBField1[3],"BFL1_3",logicWorld,false,0); //logicBField1[3]
    G4cout << physiBField1[3]->GetName() <<"[3] created."<<G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo

void SBDetectorConstruction::DoMagField2(bool place)
{
  //for the outer conductor
 G4cout << "SBDetectorConstruction::DoMagField2 "<< G4endl;  
  int HornNumber=4;
//initialization
 for(int i=0;i<24;i++){
    rInnerBOut[i]=0;
    rOuterBOut[i]=0;
    zPlaneBOut[i]=0;    
    }
//Z_plane,R_inner,R_outer
 for(int i=0;i<24;i++){ 
    zPlaneBOut[i]=zPlaneOut[i];
    rInnerBOut[i]=rInnerOut[i];
    rOuterBOut[i]=rOuterOut[i];
    }

 G4double stepMinimum2=1.0e-3*mm;
 SBField* myField2[4];
 G4Mag_UsualEqRhs* iEquation2[4];
 G4MagIntegratorStepper* iStepper2[4];
 G4ChordFinder* iChordFinder2[4];
 G4FieldManager* mfieldMgr2[4];

 for(int i=0;i<HornNumber;i++){
    myField2[i]=new SBField(current2,i,this);
    iEquation2[i]=new G4Mag_UsualEqRhs(myField2[i]);
    iStepper2[i]=new G4ClassicalRK4(iEquation2[i]);
    iChordFinder2[i]=new G4ChordFinder(myField2[i],stepMinimum2,iStepper2[i]);
    mfieldMgr2[i]=new G4FieldManager(myField2[i],iChordFinder2[i],true);
    }
                      
 solidVirtualBField2=new G4Polycone("BFL2",0.*deg,360.*deg,numZPlanesOut=23,zPlaneBOut,rInnerBOut,rOuterBOut); 
 solidBField2=new G4SubtractionSolid("BFL2-BFL1",solidVirtualBField2,solidBField1,0,G4ThreeVector(0,0,0));
 for(int i=0;i<HornNumber;i++){
    logicBField2[i]=new G4LogicalVolume(solidBField2,FillingMaterial,Form("BFL2_%d",i),mfieldMgr2[i],0,0);
    }

 if(place){
    physiBField2[0]=new G4PVPlacement(0,G4ThreeVector(Horn_xpos,Horn_ypos,Horn_zpos),logicBField2[0],"BFL2_0",logicWorld,false,0);
    G4cout << physiBField1[0]->GetName() <<"[0] created."<<G4endl;
    physiBField2[1]=new G4PVPlacement(0,G4ThreeVector(-Horn_xpos,Horn_ypos,Horn_zpos),logicBField2[1],"BFL2_1",logicWorld,false,0);
    G4cout << physiBField1[1]->GetName() <<"[1] created."<<G4endl;
    physiBField2[2]=new G4PVPlacement(0,G4ThreeVector(Horn_xpos,-Horn_ypos,Horn_zpos),logicBField2[2],"BFL2_2",logicWorld,false,0);
    G4cout << physiBField1[2]->GetName() <<"[2] created."<<G4endl;
    physiBField2[3]=new G4PVPlacement(0,G4ThreeVector(-Horn_xpos,-Horn_ypos,Horn_zpos),logicBField2[3],"BFL2_3",logicWorld,false,0);
    G4cout << physiBField1[3]->GetName() <<"[3] created."<<G4endl;
    }
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo

void SBDetectorConstruction::SetVolumesVisibility()   
{
  G4cout << "SBDetectorConstruction::SetVolumesVisibility "<< G4endl;
      
  G4VisAttributes* targetcol= new G4VisAttributes(G4Colour(1.0,0.0,0.0));
  targetcol->SetForceAuxEdgeVisible(true);
  targetcol->SetVisibility(true);
  if(logicTarget)logicTarget->SetVisAttributes(targetcol);

  G4VisAttributes* InnerConductorcol= new G4VisAttributes(G4Colour(0.85,0.5,0.));
  InnerConductorcol->SetForceAuxEdgeVisible(true);
  InnerConductorcol->SetVisibility(true);
  if(logicInnerConductor)logicInnerConductor->SetVisAttributes(InnerConductorcol);

  G4VisAttributes* OuterConductorcol= new G4VisAttributes(G4Colour(0.85,0.5,0.));
  OuterConductorcol->SetForceAuxEdgeVisible(true);
  OuterConductorcol->SetVisibility(true);
  if(logicOuterConductor)logicOuterConductor->SetVisAttributes(OuterConductorcol);

  G4VisAttributes* Tunnelcol= new G4VisAttributes(G4Colour(0.0,1.0,0.0));
  Tunnelcol->SetForceAuxEdgeVisible(true);
  Tunnelcol->SetVisibility(true);
  if(logicTunnel)logicTunnel->SetVisAttributes(Tunnelcol);
  /*
  G4VisAttributes* BField1col= new G4VisAttributes(G4Colour(0.5,0.5,0.5));   //G4Colour(1.0,0.0,0.0)
  //BField1col->SetForceAuxEdgeVisible(true);
  //BField1col->SetVisibility(true);
  for(int i=0;i<HornNumber;i++){
    if(logicBField1[i])logicBField1[i]->SetVisAttributes(BField1col);
    }

  G4VisAttributes* BField2col= new G4VisAttributes(G4Colour(0.5,0.5,0.5));
  //BField2col->SetForceAuxEdgeVisible(true);
  //BField2col->SetVisibility(true);
  for(int i=0;i<HornNumber;i++){
    if(logicBField2[i])logicBField2[i]->SetVisAttributes(BField2col);
    }
  */

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo

