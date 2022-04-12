//
//  MuFilter.cxx
//
//  by A.Buonaura
//

#include "MuFilter.h"
#include "MuFilterPoint.h"

#include "TGeoManager.h"
#include "TString.h"                    // for TString

#include "TClonesArray.h"
#include "TVirtualMC.h"

#include "TGeoBBox.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoBBox.h"
#include "TGeoCompositeShape.h"

//////////
#include "TGeoPara.h"
#include "TGeoArb8.h"
/////////

#include "TParticle.h"
#include "TParticlePDG.h"
#include "TParticleClassPDG.h"
#include "TVirtualMCStack.h"

#include "FairVolume.h"
#include "FairGeoVolume.h"
#include "FairGeoNode.h"
#include "FairRootManager.h"
#include "FairGeoLoader.h"
#include "FairGeoInterface.h"
#include "FairGeoTransform.h"
#include "FairGeoMedia.h"
#include "FairGeoMedium.h"
#include "FairGeoBuilder.h"
#include "FairRun.h"

#include "ShipDetectorList.h"
#include "ShipUnit.h"
#include "ShipStack.h"

#include <stddef.h>                     // for NULL
#include <iostream>                     // for operator<<, basic_ostream,etc
#include <string.h>
#include <cstring>

using std::cout;
using std::endl;
using std::to_string;
using std::string;
using namespace ShipUnit;

MuFilter::MuFilter()
: FairDetector("MuonFilter", "",kTRUE),
  fTrackID(-1),
fVolumeID(-1),
fPos(),
fMom(),
fTime(-1.),
fLength(-1.),
fELoss(-1),
fMuFilterPointCollection(new TClonesArray("MuFilterPoint"))
{
}

MuFilter::MuFilter(const char* name, Bool_t Active,const char* Title)
: FairDetector(name, true, kMuFilter),
  fTrackID(-1),
fVolumeID(-1),
fPos(),
fMom(),
fTime(-1.),
fLength(-1.),
fELoss(-1),
fMuFilterPointCollection(new TClonesArray("MuFilterPoint"))
{
}

MuFilter::~MuFilter()
{
    if (fMuFilterPointCollection) {
        fMuFilterPointCollection->Delete();
        delete fMuFilterPointCollection;
    }
}

void MuFilter::Initialize()
{
	FairDetector::Initialize();
}

// -----  Private method InitMedium
Int_t MuFilter::InitMedium(const char* name)
{
    static FairGeoLoader *geoLoad=FairGeoLoader::Instance();
    static FairGeoInterface *geoFace=geoLoad->getGeoInterface();
    static FairGeoMedia *media=geoFace->getMedia();
    static FairGeoBuilder *geoBuild=geoLoad->getGeoBuilder();

    FairGeoMedium *ShipMedium=media->getMedium(name);

    if (!ShipMedium)
    {
        Fatal("InitMedium","Material %s not defined in media file.", name);
        return -1111;
    }
    TGeoMedium* medium=gGeoManager->GetMedium(name);
    if (medium!=NULL)
        return ShipMedium->getMediumIndex();
    return geoBuild->createMedium(ShipMedium);
}

void MuFilter::ConstructGeometry()
{
	TGeoVolume *top=gGeoManager->FindVolumeFast("Detector");
	if(!top)  LOG(ERROR) << "no Detector volume found " ;

	//Materials 

	InitMedium("iron");
	TGeoMedium *Fe =gGeoManager->GetMedium("iron");
	InitMedium("aluminium");
	TGeoMedium *Al =gGeoManager->GetMedium("aluminium");
	InitMedium("polyvinyltoluene");
	TGeoMedium *Scint =gGeoManager->GetMedium("polyvinyltoluene");
	InitMedium("Concrete");
	TGeoMedium *concrete = gGeoManager->GetMedium("Concrete");

	Float_t nSiPMs[3];             //  number of SiPMs per side
	Float_t nSides[3];             //  number of sides readout
	nSiPMs[0] = conf_ints["MuFilter/nV"];
	nSiPMs[1] = conf_ints["MuFilter/nU"];
	nSiPMs[2] = conf_ints["MuFilter/nD"];
	nSides[0]  = conf_ints["MuFilter/sV"];
	nSides[1]  = conf_ints["MuFilter/sU"];
	nSides[2]  = conf_ints["MuFilter/sD"];

	Int_t fNUpstreamPlanes = conf_ints["MuFilter/NUpstreamPlanes"]; // Number of planes
	Int_t fNUpstreamBars = conf_ints["MuFilter/NUpstreamBars"]; // Number of staggered bars
	Int_t fNDownstreamPlanes =  conf_ints["MuFilter/NDownstreamPlanes"]; // Number of planes
	Int_t fNDownstreamBars =  conf_ints["MuFilter/NDownstreamBars"]; // Number of staggered bars
	Double_t fDownstreamBarX_ver = conf_floats["MuFilter/DownstreamBarX_ver"]; // Staggered bars of upstream section, vertical bars for x measurement
	Double_t fDownstreamBarY_ver = conf_floats["MuFilter/DownstreamBarY_ver"];
	Double_t fDownstreamBarZ_ver = conf_floats["MuFilter/DownstreamBarZ_ver"];
	Double_t fDS4ZGap = conf_floats["MuFilter/DS4ZGap"];

	// position of left bottom edges in survey coordinate system converted to physicist friendly coordinate system
	std::map<int, TVector3 > edge_Veto;
	edge_Veto[1] = TVector3( -conf_floats["MuFilter/Veto1Dx"],conf_floats["MuFilter/Veto1Dz"],conf_floats["MuFilter/Veto1Dy"]) ;
	edge_Veto[2] = TVector3( -conf_floats["MuFilter/Veto2Dx"],conf_floats["MuFilter/Veto2Dz"],conf_floats["MuFilter/Veto2Dy"]) ;
	std::map<int, TVector3 > edge_Iron;
	std::map<int, TVector3 > edge_MuFilter;
	for (int i=1;i<10;i++){
		string si = to_string(i);
		edge_Iron[i] = TVector3( -conf_floats["MuFilter/Iron"+si+"Dx"],conf_floats["MuFilter/Iron"+si+"Dz"],conf_floats["MuFilter/Iron"+si+"Dy"]);
		edge_MuFilter[i]  = TVector3( -conf_floats["MuFilter/Muon"+si+"Dx"],conf_floats["MuFilter/Muon"+si+"Dz"],conf_floats["MuFilter/Muon"+si+"Dy"]);
	}
	// system alignment parameters
	Double_t fVetoShiftX   = conf_floats["MuFilter/VetoShiftX"];
	Double_t fVetoShiftY   = conf_floats["MuFilter/VetoShiftY"];
	Double_t fVetoShiftZ   = conf_floats["MuFilter/VetoShiftZ"];
	Double_t fShiftX     = conf_floats["MuFilter/ShiftX"];
	Double_t fShiftY     = conf_floats["MuFilter/ShiftY"];
	Double_t fShiftZ     = conf_floats["MuFilter/ShiftZ"];

	TVector3 displacement;

	//Definition of the box containing veto planes
	TGeoVolumeAssembly *volVeto = new TGeoVolumeAssembly("volVeto");
	
	//Veto Planes
	Double_t fVetoBarX     = conf_floats["MuFilter/VetoBarX"]; // Veto Bar dimensions
	Double_t fVetoBarY     = conf_floats["MuFilter/VetoBarY"];
	Double_t fVetoBarZ     = conf_floats["MuFilter/VetoBarZ"];
	Double_t fVetoBarGap     = conf_floats["MuFilter/VetoBarGap"];
	Int_t fNVetoPlanes       = conf_ints["MuFilter/NVetoPlanes"];
	Int_t fNVetoBars          = conf_ints["MuFilter/NVetoBars"];
	Double_t fSupportBoxVW = conf_floats["MuFilter/SupportBoxVW"]; // SupportBox dimensions
	// local position of bottom horizontal bar to survey edge
	TVector3 LocBarVeto = TVector3(-conf_floats["MuFilter/VETOLocX"], conf_floats["MuFilter/VETOLocZ"],conf_floats["MuFilter/VETOLocY"]);

	TVector3 VetoBox1 = TVector3(-conf_floats["MuFilter/VETOBoxX1"],conf_floats["MuFilter/VETOBoxZ1"],conf_floats["MuFilter/VETOBoxY1"]); // bottom front left
	TVector3 VetoBox2 = TVector3(-conf_floats["MuFilter/VETOBoxX2"],conf_floats["MuFilter/VETOBoxZ2"],conf_floats["MuFilter/VETOBoxY2"]); // top back right
	TVector3 VetoBoxDim = TVector3( VetoBox1.X()-VetoBox2.X(), VetoBox2.Y()-VetoBox1.Y(), VetoBox2.Z()-VetoBox1.Z() ) ;
	// support box
	TGeoBBox  *supVetoBoxInner  = new TGeoBBox("supVetoBoxI",VetoBoxDim.X()/2,VetoBoxDim.Y()/2,VetoBoxDim.Z()/2);
	TGeoBBox  *supVetoBoxOuter = new TGeoBBox("supVetoBoxO",VetoBoxDim.X()/2+fSupportBoxVW,VetoBoxDim.Y()/2+fSupportBoxVW,VetoBoxDim.Z()/2+fSupportBoxVW);
	TGeoCompositeShape *subVetoBoxShape = new TGeoCompositeShape("subVetoBoxShape", "supVetoBoxO-supVetoBoxI");
	TGeoVolume *subVetoBox = new TGeoVolume("subVetoBox", subVetoBoxShape, Al);     
	subVetoBox->SetLineColor(kGray+1);

	//Veto bars
	TGeoVolume *volVetoBar = gGeoManager->MakeBox("volVetoBar",Scint,fVetoBarX/2., fVetoBarY/2., fVetoBarZ/2.);

	volVetoBar->SetLineColor(kRed-3);
	AddSensitiveVolume(volVetoBar);

	//adding veto planes
	TGeoVolume* volVetoPlane;
	for (int iplane=0; iplane < fNVetoPlanes; iplane++){
	  
	  string name = "volVetoPlane_"+to_string(iplane);
	  volVetoPlane = new TGeoVolumeAssembly(name.c_str());

	  displacement = edge_Veto[iplane+1] + LocBarVeto + TVector3(-fVetoBarX/2, 0, 0);
	  volVeto->AddNode(volVetoPlane,iplane,
				new TGeoTranslation(displacement.X(),displacement.Y(),displacement.Z()));
	 //  VETOBox1 = bottom front left
	  displacement = edge_Veto[iplane+1] +VetoBox1 + TVector3(-VetoBoxDim.X()/2,VetoBoxDim.Y()/2,VetoBoxDim.Z()/2);
	  volVeto->AddNode(subVetoBox,iplane,
		new TGeoTranslation(displacement.X(),displacement.Y(),displacement.Z()));

	  displacement = TVector3(0, 0, 0);
	  for (Int_t ibar = 0; ibar < fNVetoBars; ibar++){
	    Double_t dy_bar =  (fVetoBarY + fVetoBarGap)*ibar; 
	    volVetoPlane->AddNode(volVetoBar, 1e+4+iplane*1e+3+ibar,
				new TGeoTranslation(displacement.X(),displacement.Y()+dy_bar,displacement.Z()));
	  }
	}
	
		//adding to detector volume
	top->AddNode(volVeto, 1,new TGeoTranslation(fVetoShiftX,fVetoShiftY,fVetoShiftZ)) ;

	//*****************************************UPSTREAM SECTION*********************************//

		//Definition of the box containing Fe Blocks + Timing Detector planes 
	TGeoVolumeAssembly *volMuFilter = new TGeoVolumeAssembly("volMuFilter");
	
	//Iron blocks volume definition
	Double_t fFeBlockX = conf_floats["MuFilter/FeX"]; // Passive Iron blocks dimensions
	Double_t fFeBlockY = conf_floats["MuFilter/FeY"];
	Double_t fFeBlockZ = conf_floats["MuFilter/FeZ"];
	Double_t fFeBlockEndX = conf_floats["MuFilter/FeEndX"]; // last Iron block dimensions
	Double_t fFeBlockEndY = conf_floats["MuFilter/FeEndY"];
	Double_t fFeBlockEndZ = conf_floats["MuFilter/FeEndZ"];
	Double_t fFeBlockBotX = conf_floats["MuFilter/FeBotX"]; // bottom Iron block dimensions
	Double_t fFeBlockBotY = conf_floats["MuFilter/FeBotY"];
	Double_t fFeBlockBotZ = conf_floats["MuFilter/FeBotZ"];
	Double_t fSupportBoxW = conf_floats["MuFilter/SupportBoxW"]; // SupportBox dimensions

	TVector3 DSBox1 = TVector3(-conf_floats["MuFilter/DSBoxX1"],conf_floats["MuFilter/DSBoxZ1"],conf_floats["MuFilter/DSBoxY1"]); // bottom front left
	TVector3 DSBox2 = TVector3(-conf_floats["MuFilter/DSBoxX2"],conf_floats["MuFilter/DSBoxZ2"],conf_floats["MuFilter/DSBoxY2"]); // top back right
	TVector3 DSBoxDim = TVector3( DSBox1.X()-DSBox2.X(), DSBox2.Y()-DSBox1.Y(), DSBox2.Z()-DSBox1.Z() ) ;
	TVector3 USBox1 = TVector3(-conf_floats["MuFilter/DSBoxX1"],conf_floats["MuFilter/DSBoxZ1"],conf_floats["MuFilter/USBoxY1"]); // bottom front left
	TVector3 USBox2 = TVector3(-conf_floats["MuFilter/DSBoxX2"],conf_floats["MuFilter/DSBoxZ2"],conf_floats["MuFilter/USBoxY2"]); // top back right
	TVector3 USBoxDim = TVector3( USBox1.X()-USBox2.X(), USBox2.Y()-USBox1.Y(), USBox2.Z()-USBox1.Z() ) ;

	//Iron blocks volume definition
	TGeoVolume *volFeBlock = gGeoManager->MakeBox("volFeBlock",Fe,fFeBlockX/2, fFeBlockY/2, fFeBlockZ/2);
	volFeBlock->SetLineColor(kGreen-4);
	TGeoVolume *volFeBlockEnd = gGeoManager->MakeBox("volFeBlockEnd",Fe,fFeBlockEndX/2, fFeBlockEndY/2, fFeBlockEndZ/2);
	volFeBlockEnd->SetLineColor(kGreen-4);
	TGeoVolume *volBlockBot = gGeoManager->MakeBox("volBlockBot",concrete,fFeBlockBotX/2, fFeBlockBotY/2, fFeBlockBotZ/2);
	volBlockBot->SetLineColor(kGreen-4);

	// support box
	TGeoBBox  *supDSBoxInner  = new TGeoBBox("supDSBoxI",DSBoxDim.X()/2,DSBoxDim.Y()/2,DSBoxDim.Z()/2);
	TGeoBBox  *supDSBoxOuter = new TGeoBBox("supDSBoxO",DSBoxDim.X()/2+fSupportBoxW,DSBoxDim.Y()/2+fSupportBoxW,DSBoxDim.Z()/2+fSupportBoxW);
	TGeoCompositeShape *subDSBoxShape = new TGeoCompositeShape("subDSBoxShape", "supDSBoxO-supDSBoxI");
	TGeoVolume *subDSBox = new TGeoVolume("subDSBox", subDSBoxShape, Al);     
	subDSBox->SetLineColor(kGray+1);
	TGeoBBox  *supUSBoxInner  = new TGeoBBox("supUSBoxI",USBoxDim.X()/2,USBoxDim.Y()/2,USBoxDim.Z()/2);
	TGeoBBox  *supUSBoxOuter = new TGeoBBox("supUSBoxO",USBoxDim.X()/2+fSupportBoxW,USBoxDim.Y()/2+fSupportBoxW,USBoxDim.Z()/2+fSupportBoxW);
	TGeoCompositeShape *subUSBoxShape = new TGeoCompositeShape("subUSBoxShape", "supUSBoxO-supUSBoxI");
	TGeoVolume *subUSBox = new TGeoVolume("subUSBox", subUSBoxShape, Al);     
	subUSBox->SetLineColor(kGray+1);

	top->AddNode(volMuFilter,1,new TGeoTranslation(fShiftX,fShiftY,fShiftZ ));

	Double_t dy = 0;
	Double_t dz = 0;
	//Upstream Detector planes definition
	Double_t fUpstreamDetZ =  conf_floats["MuFilter/UpstreamDetZ"];
	// local position of bottom horizontal bar to survey edge
	TVector3 LocBarUS = TVector3(
		-conf_floats["MuFilter/DSHLocX"],
		conf_floats["MuFilter/DSHLocZ"] - conf_floats["MuFilter/DownstreamBarY"]/2 + conf_floats["MuFilter/UpstreamBarY"]/2,
		conf_floats["MuFilter/DSHLocY"]);

	TGeoVolume* volUpstreamDet;
	Double_t fUpstreamBarX = conf_floats["MuFilter/UpstreamBarX"]; //Staggered bars of upstream section
	Double_t fUpstreamBarY = conf_floats["MuFilter/UpstreamBarY"];
	Double_t fUpstreamBarZ = conf_floats["MuFilter/UpstreamBarZ"];
	Double_t fUpstreamBarGap = conf_floats["MuFilter/UpstreamBarGap"];

	//adding staggered bars, first part, only 11 bars, (single stations, readout on both ends)
	TGeoVolume *volMuUpstreamBar = gGeoManager->MakeBox("volMuUpstreamBar",Scint,fUpstreamBarX/2, fUpstreamBarY/2, fUpstreamBarZ/2);
	volMuUpstreamBar->SetLineColor(kBlue+2);
	AddSensitiveVolume(volMuUpstreamBar);
	for(Int_t l=0; l<fNUpstreamPlanes; l++)
	{
	  string name = "volMuUpstreamDet_"+std::to_string(l);
	  volUpstreamDet = new TGeoVolumeAssembly(name.c_str());

	  displacement = edge_Iron[l+1] - TVector3(fFeBlockX/2,-fFeBlockY/2,-fFeBlockZ/2);
	  volMuFilter->AddNode(volFeBlock,l,
                                    new TGeoTranslation(displacement.X(),displacement.Y(),displacement.Z()));
// place for H8 mockup target 20cm in front of US1
	  if (edge_Iron[9][2] <0.1 && l==0) {
		TGeoVolume *volFeTarget = gGeoManager->MakeBox("volFeTarget",Fe,80./2, 60./2, 29.5/2);
		volFeTarget->SetLineColor(kGreen-4);
		displacement = edge_Iron[l+1] - TVector3(80/2,-60/2,29.5/2+fFeBlockZ );
		volMuFilter->AddNode(volFeTarget,1,
                                    new TGeoTranslation(displacement.X(),displacement.Y(),displacement.Z()));
	}
	  displacement = edge_MuFilter[l+1]+LocBarUS + TVector3(-fUpstreamBarX/2, 0, 0);
	  volMuFilter->AddNode(volUpstreamDet,fNVetoPlanes+l,
                                    new TGeoTranslation(displacement.X(),displacement.Y(),displacement.Z()));

	 //  USBox1 = bottom front left
	  displacement = edge_MuFilter[l+1] +USBox1 + TVector3(-USBoxDim.X()/2,USBoxDim.Y()/2,USBoxDim.Z()/2);
	  volMuFilter->AddNode(subUSBox,l+fNVetoPlanes,
		new TGeoTranslation(displacement.X(),displacement.Y(),displacement.Z()));

	  displacement = TVector3(0, 0, 0);
	  for (Int_t ibar = 0; ibar < fNUpstreamBars; ibar++){
	    Double_t dy_bar =  (fUpstreamBarY + fUpstreamBarGap)*ibar; 
	    volUpstreamDet->AddNode(volMuUpstreamBar,2e+4+l*1e+3+ibar,
				new TGeoTranslation(displacement.X(),displacement.Y()+conf_floats["MuFilter/USOffZ"+to_string(l+1)]+dy_bar,displacement.Z()));
	  }

	}

	//*************************************DOWNSTREAM (high granularity) SECTION*****************//

    // first loop, adding detector main boxes
	TGeoVolume* volDownstreamDet;

	// local position of bottom horizontal bar to survey edge
	TVector3 LocBarH = TVector3(-conf_floats["MuFilter/DSHLocX"],conf_floats["MuFilter/DSHLocZ"],conf_floats["MuFilter/DSHLocY"]);
	// local position of l left vertical bar to survey edge
	TVector3 LocBarV = TVector3(-conf_floats["MuFilter/DSVLocX"],conf_floats["MuFilter/DSVLocZ"],conf_floats["MuFilter/DSVLocY"]);

	Double_t fDownstreamBarX = conf_floats["MuFilter/DownstreamBarX"]; // Staggered bars of upstream section
	Double_t fDownstreamBarY = conf_floats["MuFilter/DownstreamBarY"];
	Double_t fDownstreamBarZ = conf_floats["MuFilter/DownstreamBarZ"];
	Double_t fDownstreamBarGap = conf_floats["MuFilter/DownstreamBarGap"];
	
	TGeoVolume *volMuDownstreamBar_hor = gGeoManager->MakeBox("volMuDownstreamBar_hor",Scint,fDownstreamBarX/2, fDownstreamBarY/2, fDownstreamBarZ/2);
	volMuDownstreamBar_hor->SetLineColor(kAzure+7);
	AddSensitiveVolume(volMuDownstreamBar_hor);

	//vertical bars, for x measurement
	TGeoVolume *volMuDownstreamBar_ver = gGeoManager->MakeBox("volMuDownstreamBar_ver",Scint,fDownstreamBarX_ver/2, fDownstreamBarY_ver/2, fDownstreamBarZ/2);
	volMuDownstreamBar_ver->SetLineColor(kViolet-2);
	AddSensitiveVolume(volMuDownstreamBar_ver);

	for(Int_t l=0; l<fNDownstreamPlanes; l++)
	{
	// add iron blocks
	if (l<fNDownstreamPlanes-1){
		displacement = edge_Iron[l+fNUpstreamPlanes+1] - TVector3(fFeBlockX/2,-fFeBlockY/2,-fFeBlockZ/2);
		volMuFilter->AddNode(volFeBlock,l+fNUpstreamPlanes+fNVetoPlanes,
				new TGeoTranslation(displacement.X(),displacement.Y(),displacement.Z()));
		}else if (edge_Iron[9][2] >0.1) {
// more iron
		displacement = edge_Iron[l+fNUpstreamPlanes+1]  - TVector3(fFeBlockEndX/2,-fFeBlockEndY/2,-fFeBlockEndZ/2);
		volMuFilter->AddNode(volFeBlockEnd,1,
				new TGeoTranslation(displacement.X(),displacement.Y(),displacement.Z()));
		displacement = edge_Iron[l+fNUpstreamPlanes+1]  - TVector3(fFeBlockBotX/2-10.0, fFeBlockBotY/2,fFeBlockBotZ/2-fFeBlockEndZ);
		volMuFilter->AddNode(volBlockBot,1,
				new TGeoTranslation(displacement.X(),displacement.Y(),displacement.Z()));
	}

	// add detectors
	string name = "volMuDownstreamDet_"+std::to_string(l);
	volDownstreamDet = new TGeoVolumeAssembly(name.c_str());
	displacement = edge_MuFilter[l+fNUpstreamPlanes+1] + LocBarH + TVector3(-fDownstreamBarX/2, 0,0);
	
	volMuFilter->AddNode(volDownstreamDet,l+fNUpstreamPlanes+fNVetoPlanes,
				new TGeoTranslation(displacement.X(),displacement.Y(),displacement.Z()));

	//adding bars within each detector box
	if (l!=fNDownstreamPlanes-1) {
		displacement = TVector3(0, 0,0);
		for (Int_t ibar = 0; ibar < fNDownstreamBars; ibar++){
	                 //adding horizontal bars for y
			Double_t dy_bar = (fDownstreamBarY + fDownstreamBarGap)*ibar;
		    	volDownstreamDet->AddNode(volMuDownstreamBar_hor,3e+4+l*1e+3+ibar,
				new TGeoTranslation(displacement.X(),displacement.Y()+dy_bar,displacement.Z()));
			}
		}
	 //  DSBox1 = bottom front left
	displacement = edge_MuFilter[l+fNUpstreamPlanes+1] +DSBox1 +
			TVector3(-DSBoxDim.X()/2,DSBoxDim.Y()/2,DSBoxDim.Z()/2);
		volMuFilter->AddNode(subDSBox,l+fNUpstreamPlanes+fNVetoPlanes,
		new TGeoTranslation(displacement.X(),displacement.Y(),displacement.Z()));
	//adding vertical bars for x
	displacement = LocBarV + TVector3(0, -fDownstreamBarY_ver/2,0) - LocBarH - TVector3(-fDownstreamBarX/2, 0,0);
	for (Int_t i_vbar = 0; i_vbar<fNDownstreamBars; i_vbar++) {
		Double_t dx_bar = (fDownstreamBarX_ver+ fDownstreamBarGap)*i_vbar;
		Int_t i_vbar_rev  = fNDownstreamBars-1-i_vbar;
		volDownstreamDet->AddNode(volMuDownstreamBar_ver,3e+4+l*1e+3+i_vbar_rev+60,
				new TGeoTranslation(displacement.X()+dx_bar,displacement.Y(),displacement.Z())); 
 		// Andrew added a 60 here to make each horizontal + vertical sub-plane contain bars given detIDs as one plane. So the first bar in the vert. sub plane is the 60th etc. 
			}
	}
	// COLDBOX GEOM
	InitMedium("Borated30polyethylene");
	TGeoMedium *Bor30Poly =gGeoManager->GetMedium("Borated30polyethylene");
	InitMedium("PMMA");
	TGeoMedium *Acrylic =gGeoManager->GetMedium("PMMA");

	Float_t Acrylic_width		= 5.0; 		// cm
	Float_t BPoly_width		= 4.0; 		// cm
   	Float_t CBFrontWall_xdim	= 219.;		// cm
	Float_t CBFrontWall_ydim	= 170.72-Acrylic_width; 	// cm
	Float_t CBLatWall_zdim		= 176.0; 	// cm
	Float_t CBTiny_zdim		= 17.0; 	// cm
	Float_t CBExtra_zdim		= 41.0;		// cm
	Float_t CBExtra_xdim		= 67.5;		// cm
	Float_t SlopedWall_zproj	= 110.0;	// cm	
	Float_t CBRearWall_xdim		= CBFrontWall_xdim-SlopedWall_zproj*TMath::Tan(TMath::DegToRad()*15.)-CBExtra_xdim+Acrylic_width;

    	Float_t FeX		= 80.;
	Float_t FeY		= 60.;
	Float_t FeZ		= 20.;


    // ************************ ACRYLIC
    // Shapes definition
    	TGeoBBox *FeBlock_cb 	 	= new TGeoBBox("FeBlock_cb",FeX/2, FeY/2, FeZ/2);
	TGeoBBox *CBFrontWall_a 	= new TGeoBBox("CBFrontWall_a", CBFrontWall_xdim/2., CBFrontWall_ydim/2., Acrylic_width/2.);
	TGeoBBox *CBLateral_a 		= new TGeoBBox("CBLateral_a", Acrylic_width/2., CBFrontWall_ydim/2., (CBLatWall_zdim-2*Acrylic_width)/2.);
	TGeoBBox *CBExtraFront_a	= new TGeoBBox("CBExtraFront_a", CBExtra_xdim/2., CBFrontWall_ydim/2., Acrylic_width/2.);
	TGeoBBox *CBExtraLat_a		= new TGeoBBox("CBExtraLat_a", Acrylic_width/2., CBFrontWall_ydim/2., (CBExtra_zdim-2*Acrylic_width)/2.);
	TGeoBBox *CBTiny1_a		= new TGeoBBox("CBTiny1_a", Acrylic_width/2., CBFrontWall_ydim/2., (CBTiny_zdim-Acrylic_width-BPoly_width)/2.);
	TGeoBBox *CBTiny2_a		= new TGeoBBox("CBTiny2_a", Acrylic_width/2., CBFrontWall_ydim/2., (CBTiny_zdim-Acrylic_width)/2.);
	TGeoBBox *CBRearWall_a		= new TGeoBBox("CBRearWall_a", CBRearWall_xdim/2., CBFrontWall_ydim/2., Acrylic_width/2.);

    	TGeoTranslation *CBWallpos = new TGeoTranslation("CBWallpos", (CBRearWall_xdim-FeX)/2. - 28.5, (FeY-CBFrontWall_ydim)/2., 0);
	CBWallpos->RegisterYourself();

	TGeoCompositeShape *CBWallDownstream = new TGeoCompositeShape("CBWallDownstream", "CBRearWall_a-(FeBlock_cb:CBWallpos)");
    	TGeoPara *CBWallSlope_a = new TGeoPara("CBWallSlope_a", Acrylic_width/2., CBFrontWall_ydim/2., SlopedWall_zproj/2., 0, -15, 0);

    // Acrylic mother shape definition
    TGeoTranslation *FrontWallpos = new TGeoTranslation("FrontWallpos", -CBRearWall_xdim/2.-CBExtra_xdim+Acrylic_width+CBFrontWall_xdim/2., 0, -(SlopedWall_zproj+2*(CBTiny_zdim-Acrylic_width)+Acrylic_width)+BPoly_width);
	FrontWallpos->RegisterYourself();
	TGeoTranslation *Tiny1pos = new TGeoTranslation("Tiny1pos", (CBRearWall_xdim-Acrylic_width)/2., 0, -CBTiny_zdim/2.+BPoly_width/2.);
	Tiny1pos->RegisterYourself();
	TGeoTranslation *SlopeWallpos = new TGeoTranslation("SlopeWallpos", (CBRearWall_xdim+Acrylic_width)/2.+Acrylic_width+SlopedWall_zproj/(2*TMath::Tan(TMath::DegToRad()*85.)),0, -CBTiny_zdim-SlopedWall_zproj/2.+Acrylic_width/2.+BPoly_width);
	SlopeWallpos->RegisterYourself();
	TGeoTranslation *Tiny2pos = new TGeoTranslation("Tiny2pos", 3*Acrylic_width+CBRearWall_xdim/2.+Acrylic_width/2.+SlopedWall_zproj/(TMath::Tan(TMath::DegToRad()*85.)), 0, -(SlopedWall_zproj+2*(CBTiny_zdim-Acrylic_width)+Acrylic_width)+CBTiny_zdim/2.+BPoly_width);
	Tiny2pos->RegisterYourself();
	TGeoTranslation *CBExtraLatpos = new TGeoTranslation("CBExtraLatpos", -CBRearWall_xdim/2.+Acrylic_width/2., 0, Acrylic_width/2.+(CBExtra_zdim-2*Acrylic_width)/2.);
	CBExtraLatpos->RegisterYourself();
	TGeoTranslation *CBExtraFrontpos = new TGeoTranslation("CBExtraFrontpos", -CBRearWall_xdim/2.+Acrylic_width-CBExtra_xdim/2., 0, CBExtra_zdim-Acrylic_width);
	CBExtraFrontpos->RegisterYourself();
	TGeoTranslation *CBLateralpos = new TGeoTranslation("CBLateralpos", -CBRearWall_xdim/2.-CBExtra_xdim+Acrylic_width+Acrylic_width/2., 0, CBExtra_zdim-CBLatWall_zdim/2.-Acrylic_width/2.);
	CBLateralpos->RegisterYourself();
    
    // Acrylic mother volume definition
    	TGeoCompositeShape *COLDBOXA = new TGeoCompositeShape("COLDBOXA", "CBWallDownstream+(CBFrontWall_a:FrontWallpos)+(CBTiny1_a:Tiny1pos)+(CBWallSlope_a:SlopeWallpos)+(CBTiny2_a:Tiny2pos)+(CBExtraLat_a:CBExtraLatpos)+(CBExtraFront_a:CBExtraFrontpos)+(CBLateral_a:CBLateralpos)");
	TGeoVolume *volCOLDBOXA = new TGeoVolume("volCOLDBOXA", COLDBOXA, Acrylic);
    

    // ************************ BORATED POLYETHYLENE
    	Float_t CBFrontWall_xdim_b	= CBFrontWall_xdim-2*Acrylic_width-BPoly_width;		// cm
	Float_t CBFrontWall_ydim_b	= CBFrontWall_ydim-BPoly_width;
	Float_t CBLatWall_zdim_b	= CBLatWall_zdim-2*Acrylic_width; 	// cm
	Float_t CBExtra_xdim_b		= CBExtra_xdim-2*Acrylic_width;		// cm
	Float_t CBRearWall_xdim_b	= CBRearWall_xdim-Acrylic_width;
    // Shapes definition
    TGeoBBox *CBFrontWall_b 	= new TGeoBBox("CBFrontWall_b", CBFrontWall_xdim_b/2.+(Acrylic_width-BPoly_width)/10., CBFrontWall_ydim_b/2., BPoly_width/2.); // (Acrylic_width-BPoly_width)/10. is due to approximations, I guess
	TGeoBBox *CBLateral_b 		= new TGeoBBox("CBLateral_b", BPoly_width/2., CBFrontWall_ydim_b/2., (CBLatWall_zdim-2*BPoly_width)/2.-Acrylic_width);
	TGeoBBox *CBExtraFront_b	= new TGeoBBox("CBExtraFront_b", CBExtra_xdim_b/2., CBFrontWall_ydim_b/2., BPoly_width/2.);
	TGeoBBox *CBExtraLat_b		= new TGeoBBox("CBExtraLat_b", BPoly_width/2., CBFrontWall_ydim_b/2., (CBExtra_zdim-Acrylic_width-BPoly_width)/2.);
	TGeoBBox *CBTiny1_b		= new TGeoBBox("CBTiny1_b", BPoly_width/2., CBFrontWall_ydim_b/2., (CBTiny_zdim-Acrylic_width-BPoly_width)/2.);
	TGeoBBox *CBTiny2_b		= new TGeoBBox("CBTiny2_b", BPoly_width/2., CBFrontWall_ydim_b/2., (CBTiny_zdim-Acrylic_width)/2.);
	TGeoBBox *CBRearWall_b		= new TGeoBBox("CBRearWall_b", CBRearWall_xdim_b/2., CBFrontWall_ydim_b/2., BPoly_width/2.);

    	TGeoPara *CBWallSlope_b = new TGeoPara("CBWallSlope_b", BPoly_width/2., CBFrontWall_ydim_b/2., SlopedWall_zproj/2., 0, -15, 0);

    // Borated Polyethylene mother shape definition
    TGeoTranslation *FrontWallpos_b = new TGeoTranslation("FrontWallpos_b", -CBRearWall_xdim_b/2.-CBExtra_xdim_b+BPoly_width+CBFrontWall_xdim_b/2.+0.1, 0,-SlopedWall_zproj-(CBTiny_zdim-Acrylic_width+BPoly_width)); // +0.1 is due to approximations, I guess 
	FrontWallpos_b->RegisterYourself();
	TGeoTranslation *Tiny1pos_b = new TGeoTranslation("Tiny1pos_b", CBRearWall_xdim_b/2.+BPoly_width/2., 0, -BPoly_width/2.);
	Tiny1pos_b->RegisterYourself();
	TGeoTranslation *SlopeWallpos_b = new TGeoTranslation("SlopeWallpos_b", SlopedWall_zproj/(2*TMath::Tan(TMath::DegToRad()*85.))+CBRearWall_xdim_b/2.+3*BPoly_width,0, -BPoly_width/2.-SlopedWall_zproj/2.-(CBTiny_zdim-Acrylic_width-BPoly_width)/2.);
	SlopeWallpos_b->RegisterYourself();
	TGeoTranslation *Tiny2pos_b = new TGeoTranslation("Tiny2pos_b", 5*BPoly_width+CBRearWall_xdim/2.+SlopedWall_zproj/(TMath::Tan(TMath::DegToRad()*85.))-(Acrylic_width-BPoly_width)/2., 0, -SlopedWall_zproj-(CBTiny_zdim-Acrylic_width));
	Tiny2pos_b->RegisterYourself();
	TGeoTranslation *CBExtraLatpos_b = new TGeoTranslation("CBExtraLatpos_b", -CBRearWall_xdim_b/2.+BPoly_width/2., 0, (CBExtra_zdim-Acrylic_width)/2.);
	CBExtraLatpos_b->RegisterYourself();
	TGeoTranslation *CBExtraFrontpos_b = new TGeoTranslation("CBExtraFrontpos_b", -CBRearWall_xdim_b/2.+BPoly_width-CBExtra_xdim_b/2., 0, CBExtra_zdim-Acrylic_width);
	CBExtraFrontpos_b->RegisterYourself();
	TGeoTranslation *CBLateralpos_b = new TGeoTranslation("CBLateralpos_b", -CBRearWall_xdim_b/2.-CBExtra_xdim_b+BPoly_width+BPoly_width/2., 0, CBExtra_zdim-CBLatWall_zdim_b/2.+BPoly_width/2.-Acrylic_width);
	CBLateralpos_b->RegisterYourself();
    
	// Borated Polyethylene mother volume definition
	TGeoCompositeShape *COLDBOXB = new TGeoCompositeShape("COLDBOXB","CBRearWall_b+(CBTiny1_b:Tiny1pos_b)+(CBExtraLat_b:CBExtraLatpos_b)+(CBWallSlope_b:SlopeWallpos_b)+(CBTiny2_b:Tiny2pos_b)+(CBExtraFront_b:CBExtraFrontpos_b)+(CBLateral_b:CBLateralpos_b)+(CBFrontWall_b:FrontWallpos_b)");
    	TGeoVolume *volCOLDBOXB = new TGeoVolume("volCOLDBOXB", COLDBOXB, Bor30Poly);
    	

	// Acrylic Roof shape definition
	Double_t Roof4_averts [8][2];
	Roof4_averts[0][0] = 0.; Roof4_averts[0][1] = 0.;
	Roof4_averts[1][0] = 0.; Roof4_averts[1][1] = Acrylic_width;
	Roof4_averts[2][0] = SlopedWall_zproj*(TMath::Tan(TMath::DegToRad()*15.)); Roof4_averts[2][1] = Acrylic_width;
	Roof4_averts[3][0] = SlopedWall_zproj*(TMath::Tan(TMath::DegToRad()*15.)); Roof4_averts[3][1] = 0;
	Roof4_averts[4][0] = 0; Roof4_averts[4][1] = 0;
	Roof4_averts[5][0] = 0; Roof4_averts[5][1] = Acrylic_width;
	Roof4_averts[6][0] = 0; Roof4_averts[6][1] = Acrylic_width;
	Roof4_averts[7][0] = 0; Roof4_averts[7][1] = 0;
	
	TGeoBBox *CBRoof1_a 	= new TGeoBBox("CBRoof1_a", CBExtra_xdim/2., Acrylic_width/2., CBLatWall_zdim/2.);
	TGeoBBox *CBRoof2_a	= new TGeoBBox("CBRoof2_a", (CBRearWall_xdim-Acrylic_width)/2., Acrylic_width/2., (CBLatWall_zdim-CBExtra_zdim+Acrylic_width)/2.);
	TGeoBBox *CBRoof3_a	= new TGeoBBox("CBRoof3_a", (SlopedWall_zproj*(TMath::Tan(TMath::DegToRad()*15.)))/2., Acrylic_width/2., CBTiny_zdim/2.);
	TGeoArb8 *CBRoof4_a	= new TGeoArb8("CBRoof4_a", SlopedWall_zproj/2., (Double_t*) Roof4_averts);

	TGeoTranslation *Roof1_apos = new TGeoTranslation("Roof1_apos", -(CBRearWall_xdim-Acrylic_width)/2.-CBExtra_xdim/2., 0, CBExtra_zdim/2.-Acrylic_width/2.);
	Roof1_apos->RegisterYourself();
	TGeoTranslation *Roof3_apos = new TGeoTranslation("Roof3_apos", (CBRearWall_xdim-Acrylic_width)/2.+(SlopedWall_zproj*(TMath::Tan(TMath::DegToRad()*15.)))/2., 0, -(CBLatWall_zdim-CBExtra_zdim+Acrylic_width)/2.+CBTiny_zdim/2.);
	Roof3_apos->RegisterYourself();
	TGeoTranslation *Roof4_apos = new TGeoTranslation("Roof4_apos", (CBRearWall_xdim-Acrylic_width)/2., -Acrylic_width/2., -(CBLatWall_zdim-CBExtra_zdim+Acrylic_width)/2.+CBTiny_zdim+SlopedWall_zproj/2.);
	Roof4_apos->RegisterYourself();

	// Acrylic roof volume definition
	TGeoCompositeShape *CBRoof_a = new TGeoCompositeShape("CBRoof_a", "CBRoof2_a+(CBRoof1_a:Roof1_apos)+(CBRoof3_a:Roof3_apos)+(CBRoof4_a:Roof4_apos)");
	TGeoVolume *volCBRoof_a = new TGeoVolume("volCBRoof_a", CBRoof_a, Acrylic);
	
	
	// Borated Polythylene Roof shape definition
	Double_t Roof4_bverts [8][2];
	Roof4_bverts[0][0] = 0.; Roof4_bverts[0][1] = 0.;
	Roof4_bverts[1][0] = 0.; Roof4_bverts[1][1] = BPoly_width;
	Roof4_bverts[2][0] = SlopedWall_zproj*(TMath::Tan(TMath::DegToRad()*15.)); Roof4_bverts[2][1] = BPoly_width;
	Roof4_bverts[3][0] = SlopedWall_zproj*(TMath::Tan(TMath::DegToRad()*15.)); Roof4_bverts[3][1] = 0;
	Roof4_bverts[4][0] = 0; Roof4_bverts[4][1] = 0;
	Roof4_bverts[5][0] = 0; Roof4_bverts[5][1] = BPoly_width;
	Roof4_bverts[6][0] = 0; Roof4_bverts[6][1] = BPoly_width;
	Roof4_bverts[7][0] = 0; Roof4_bverts[7][1] = 0;

	TGeoBBox *CBRoof1_b 	= new TGeoBBox("CBRoof1_b", CBExtra_xdim_b/2., BPoly_width/2., CBLatWall_zdim_b/2.);
	TGeoBBox *CBRoof2_b	= new TGeoBBox("CBRoof2_b", (CBRearWall_xdim_b-BPoly_width)/2.+BPoly_width/2., BPoly_width/2., (CBLatWall_zdim_b-CBExtra_zdim+Acrylic_width)/2.);
	TGeoBBox *CBRoof3_b	= new TGeoBBox("CBRoof3_b", (SlopedWall_zproj*(TMath::Tan(TMath::DegToRad()*15.)))/2., BPoly_width/2., (CBTiny_zdim-Acrylic_width)/2.);
	TGeoArb8 *CBRoof4_b	= new TGeoArb8("CBRoof4_b", SlopedWall_zproj/2., (Double_t*) Roof4_bverts);

	TGeoTranslation *Roof1_bpos = new TGeoTranslation("Roof1_bpos", -(CBRearWall_xdim-Acrylic_width)/2.-CBExtra_xdim_b/2., 0, CBExtra_zdim/2.-Acrylic_width/2.);
	Roof1_bpos->RegisterYourself();
	TGeoTranslation *Roof3_bpos = new TGeoTranslation("Roof3_bpos", (CBRearWall_xdim-Acrylic_width)/2.+(SlopedWall_zproj*(TMath::Tan(TMath::DegToRad()*15.)))/2., 0, -(CBLatWall_zdim_b-CBExtra_zdim+Acrylic_width)/2.+(CBTiny_zdim-Acrylic_width)/2.);
	Roof3_bpos->RegisterYourself();
	TGeoTranslation *Roof4_bpos = new TGeoTranslation("Roof4_bpos", (CBRearWall_xdim-Acrylic_width)/2., -BPoly_width/2., -(CBLatWall_zdim_b-CBExtra_zdim+Acrylic_width)/2.+CBTiny_zdim-Acrylic_width+SlopedWall_zproj/2.);
	Roof4_bpos->RegisterYourself();

	// Borated Polyethylene volume definition
	TGeoCompositeShape *CBRoof_b = new TGeoCompositeShape("CBRoof_b", "CBRoof2_b+(CBRoof1_b:Roof1_bpos)+(CBRoof3_b:Roof3_bpos)+(CBRoof4_b:Roof4_bpos)");
	TGeoVolume *volCBRoof_b = new TGeoVolume("volCBRoof_b", CBRoof_b, Bor30Poly);
	
	// Volumes positioning
	TGeoVolumeAssembly *volColdBox = new TGeoVolumeAssembly("volColdBox");
	volCOLDBOXA->SetLineColor(kGray-1);
	volCOLDBOXA->SetTransparency(50);
	volCOLDBOXB->SetLineColor(kGray-1);
	volCOLDBOXB->SetTransparency(50);
	volCBRoof_a->SetLineColor(kGray-1);
	volCBRoof_a->SetTransparency(50);
	volCBRoof_b->SetLineColor(kGray-1);
	volCBRoof_b->SetTransparency(50);
	
    	volColdBox->AddNode(volCOLDBOXA, 0, 0);
	volColdBox->AddNode(volCOLDBOXB, 0, new TGeoTranslation(-BPoly_width-Acrylic_width/2., -BPoly_width/2., -Acrylic_width/2.-BPoly_width/2.));
	volColdBox->AddNode(volCBRoof_a, 0, new TGeoTranslation(Acrylic_width/2., CBFrontWall_ydim/2.+Acrylic_width/2., -(CBLatWall_zdim-CBExtra_zdim+Acrylic_width)/2.+Acrylic_width/2.));
	volColdBox->AddNode(volCBRoof_b, 0, new TGeoTranslation(-Acrylic_width/2., CBFrontWall_ydim/2.-BPoly_width/2., -(CBLatWall_zdim-CBExtra_zdim+Acrylic_width)/2.+Acrylic_width/2.));

	displacement = edge_Iron[1] - TVector3(fFeBlockX/2,-fFeBlockY/2,-fFeBlockZ/2);
	volMuFilter->AddNode(volColdBox, 0, new TGeoTranslation(displacement.X()-(CBRearWall_xdim-fFeBlockX)/2.+28.5, displacement.Y()+(CBFrontWall_ydim-fFeBlockY)/2., displacement.Z()+Acrylic_width-fFeBlockZ/2.-BPoly_width+1.));
	//LOG(INFO) << "++++ColBox position (x, y, z): " << displacement.X()-(CBRearWall_xdim-fFeBlockX)/2.+28.5 << "\t" << displacement.Y()+(CBFrontWall_ydim-fFeBlockY)/2. << "\t" << displacement.Z()+Acrylic_width-fFeBlockZ/2.-BPoly_width+1;

}

Bool_t  MuFilter::ProcessHits(FairVolume* vol)
{
	/** This method is called from the MC stepping */
	//Set parameters at entrance of volume. Reset ELoss.
	if ( gMC->IsTrackEntering() ) 
	{
		fELoss  = 0.;
		fTime   = gMC->TrackTime() * 1.0e09;
		fLength = gMC->TrackLength();
		gMC->TrackPosition(fPos);
		gMC->TrackMomentum(fMom);
	}
	// Sum energy loss for all steps in the active volume
	fELoss += gMC->Edep();

	// Create MuFilterPoint at exit of active volume
	if ( gMC->IsTrackExiting()    ||
			gMC->IsTrackStop()       || 
			gMC->IsTrackDisappeared()   ) {
		fTrackID  = gMC->GetStack()->GetCurrentTrackNumber();
		gMC->CurrentVolID(fVolumeID);

		if (fELoss == 0. ) { return kFALSE; }
		TParticle* p=gMC->GetStack()->GetCurrentTrack();
		Int_t pdgCode = p->GetPdgCode();

		TLorentzVector Pos; 
		gMC->TrackPosition(Pos); 
		Double_t xmean = (fPos.X()+Pos.X())/2. ;
		Double_t ymean = (fPos.Y()+Pos.Y())/2. ;
		Double_t zmean = (fPos.Z()+Pos.Z())/2. ;


		AddHit(fTrackID,fVolumeID, TVector3(xmean, ymean,  zmean),
				TVector3(fMom.Px(), fMom.Py(), fMom.Pz()), fTime, fLength,
				fELoss, pdgCode);

		// Increment number of muon det points in TParticle
		ShipStack* stack = (ShipStack*) gMC->GetStack();
		stack->AddPoint(kMuFilter);
	}   

	return kTRUE;
}

void MuFilter::EndOfEvent()
{
    fMuFilterPointCollection->Clear();
}


void MuFilter::Register()
{

    /** This will create a branch in the output tree called
 *      MuFilterPoint, setting the last parameter to kFALSE means:
 *           this collection will not be written to the file, it will exist
 *                only during the simulation.
 *                     */

    FairRootManager::Instance()->Register("MuFilterPoint", "MuFilter",
                                          fMuFilterPointCollection, kTRUE);
}

TClonesArray* MuFilter::GetCollection(Int_t iColl) const
{
    if (iColl == 0) { return fMuFilterPointCollection; }
    else { return NULL; }
}

void MuFilter::Reset()
{
    fMuFilterPointCollection->Clear();
}


MuFilterPoint* MuFilter::AddHit(Int_t trackID,Int_t detID,
                           TVector3 pos, TVector3 mom,
                           Double_t time, Double_t length,
                            Double_t eLoss, Int_t pdgCode)
{
    TClonesArray& clref = *fMuFilterPointCollection;
    Int_t size = clref.GetEntriesFast();
    return new(clref[size]) MuFilterPoint(trackID,detID, pos, mom,
                                        time, length, eLoss, pdgCode);
}

void MuFilter::GetLocalPosition(Int_t fDetectorID, TVector3& vLeft, TVector3& vRight)
{
  GetPosition(fDetectorID, vLeft, vRight);
  TGeoNavigator* nav = gGeoManager->GetCurrentNavigator();
  TString path = nav->GetPath();
  TString loc( path(0, path.Last('/') ) );
  nav->cd(loc);
  Double_t A[3] = {vLeft.X(),vLeft.Y(),vLeft.Z()};
  Double_t B[3] = {vRight.X(),vRight.Y(),vRight.Z()};
  Double_t locA[3] = {0,0,0};
  Double_t locB[3] = {0,0,0};
  nav->MasterToLocal(A, locA);   nav->MasterToLocal(B, locB);
  vLeft.SetXYZ(locA[0],locA[1],locA[2]);
  vRight.SetXYZ(locB[0],locB[1],locB[2]);
}

void MuFilter::GetPosition(Int_t fDetectorID, TVector3& vLeft, TVector3& vRight)
{

  int subsystem     = floor(fDetectorID/10000);
  int plane             = floor(fDetectorID/1000) - 10*subsystem;
  int bar_number   = fDetectorID%1000;

  TString path = "/cave_1/Detector_0/";
  TString barName;

  switch(subsystem) {
  
  case 1: 
      path+="volVeto_1/volVetoPlane_"+std::to_string(plane)+"_"+std::to_string(plane);
      barName = "/volVetoBar_";
      break;
  case 2: 
      path+="volMuFilter_1/volMuUpstreamDet_"+std::to_string(plane)+"_"+std::to_string(plane+2);
      barName = "/volMuUpstreamBar_";
      break;
  case 3: 
      path+="volMuFilter_1/volMuDownstreamDet_"+std::to_string(plane)+"_"+std::to_string(plane+7);
      barName = "/volMuDownstreamBar_";
      if (bar_number<60){barName+="hor_";}
      else{barName+="ver_";}
      break;
  }

  path += barName+std::to_string(fDetectorID);

  TGeoNavigator* nav = gGeoManager->GetCurrentNavigator();
  nav->cd(path);
  LOG(DEBUG) <<path<<" "<<fDetectorID<<" "<<subsystem<<" "<<bar_number;
  TGeoNode* W = nav->GetCurrentNode();
  TGeoBBox* S = dynamic_cast<TGeoBBox*>(W->GetVolume()->GetShape());

  if (subsystem == 3 and bar_number >59){  // vertical bars
      Double_t top[3] = {0,S->GetDY(), 0};
      Double_t bot[3] = {0,-S->GetDY(),0};
      Double_t Gtop[3],Gbot[3];
      nav->LocalToMaster(top, Gtop);   nav->LocalToMaster(bot, Gbot);
      vLeft.SetXYZ(Gtop[0],Gtop[1],Gtop[2]);
      vRight.SetXYZ(Gbot[0],Gbot[1],Gbot[2]);
    }
    else {     // horizontal bars
      Double_t posXend[3] = {S->GetDX(),0,0};
      Double_t negXend[3] = {-S->GetDX(),0,0};
      Double_t GposXend[3],GnegXend[3];
      nav->LocalToMaster(posXend, GposXend);   nav->LocalToMaster(negXend, GnegXend);
      vLeft.SetXYZ(GposXend[0],GposXend[1],GposXend[2]);
      vRight.SetXYZ(GnegXend[0],GnegXend[1],GnegXend[2]);
    }
}

   Int_t MuFilter::GetnSiPMs(Int_t detID){
       int subsystem     = floor(detID/10000)-1;
       if (subsystem==0){return conf_ints["MuFilter/nV"];}
       if (subsystem==1){return conf_ints["MuFilter/nU"];}
       return conf_ints["MuFilter/nD"];
   }
   Int_t MuFilter::GetnSides(Int_t detID){
       int subsystem     = floor(detID/10000)-1;
       if (subsystem==0){return conf_ints["MuFilter/sV"];}
       if (subsystem==1){return conf_ints["MuFilter/sU"];}
       return conf_ints["MuFilter/sD"];
  }

ClassImp(MuFilter)
