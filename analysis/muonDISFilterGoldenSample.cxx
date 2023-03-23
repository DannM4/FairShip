#include <iostream>
#include <vector>

#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TGeoManager.h"
#include "ShipMCTrack.h"

// Cuts
#include "sndBaseCut.h"
#include "sndMinSciFiHitsCut.h"
#include "sndSciFiStationCut.h"
#include "sndVetoCut.h"
#include "sndMinSciFiConsecutivePlanes.h"
#include "sndDSActivityCut.h"
#include "sndUSQDCCut.h"
#include "sndEventDeltat.h"
#include "sndAvgSciFiFiducialCut.h"
#include "sndAvgDSFiducialCut.h"
/* MuonDIS golden sample filter */
/* Usage: 
 * if you're working with muonDIS simulations the arguments are: 1.inputfilte 2.outputfile 3.cutset(0) 4.geometryfile 5.isPMU(0)
 * if you're working with passingMu simulations one can try with: 1.inputfilte 2.outputfile 3.cutset(0) 4.geometryfile(0) 5.isPMU(1)  
 * */

enum Cutset { stage1cuts } ; // one cutset for now

int main(int argc, char ** argv) {

  std::cout << "Starting muon DIS filter" << std::endl;
  
  if (argc < 4 || argc > 6) {
    std::cout << "Three arguments required: input file name (or reg exp), output file name, cut set (0: stage 1 selection) (geofile as 4th arg if MC)" << std::endl;
    exit(-1);
  }

  // Input files
  bool isMC = false;
  TChain * ch = new TChain("rawConv");
  ch->Add(argv[1]);
  if (ch->GetEntries() == 0){
    delete ch;
    ch = new TChain("cbmsim");
    ch->Add(argv[1]);
    if (ch->GetEntries() > 0) {
      isMC = true;
      ch->SetBranchStatus("EventHeader", 0); // To be able to run on old MC.
    } else {
      std::cout << "Didn't find rawConv or cbmsim in input file" << std::endl;
      exit(-1);
    }
  }
  std::cout << "Got input tree" << std::endl;

  // MC truth
  TClonesArray * MCTracks = new TClonesArray("ShipMCTrack", 2*5000);
  if (isMC) ch->SetBranchAddress("MCTrack", &MCTracks);

  // Output file
  TFile * outFile = new TFile(argv[2], "RECREATE");
  std::cout << "Got output file" << std::endl;  

  ch->GetEntry(0);
  ch->GetFile()->Get("BranchList")->Write("BranchList", TObject::kSingleKey);
  ch->GetFile()->Get("TimeBasedBranchList")->Write("TimeBasedBranchList", TObject::kSingleKey);
  ch->GetFile()->Get("FileHeader")->Write();

  // Set up all branches to copy to output TTree.
  TTree * outTree = ch->CloneTree(0);
  std::cout << "Got output tree" << std::endl;  

  // Set up cuts
  std::cout << "Starting cut set up" << std::endl;

  std::vector< sndAnalysis::baseCut * > cutFlow;
  
  int selected_cutset = std::atoi(argv[3]);
  int isPMU = std::atoi(argv[5]);
  TGraph *crsec = NULL;
  TGraph *crsec2 = NULL;
  TFile  *file2 = NULL;
  if (isPMU) std::cout << "Passing muons simulation obtained" << std::endl;

  if (isMC && argc != 6)
  {
    std::cout << "Five arguments are needed for MC" << std::endl;
    exit(-1);
  }
  else if (isMC && !isPMU ) {
	TGeoManager *m1 = TGeoManager::Import(argv[4]); 
	gGeoManager = m1;
	file2 = new TFile("/eos/experiment/sndlhc/MonteCarlo/Pythia6/MuonDIS/muDIScrossSec.root");
      	if(!file2) exit(-1);
    	crsec = (TGraph*) file2->Get("g_13");
	crsec2 = (TGraph*) file2->Get("g_-13");
	outFile->cd();
  }
  

  if (selected_cutset == stage1cuts){ // Stage 1 cuts
    //    cutFlow.push_back( new sndAnalysis::minSciFiHits(1, ch)); // A. Events with SciFi activity
    cutFlow.push_back( new sndAnalysis::vetoCut(ch)); // B. No veto hits
    cutFlow.push_back( new sndAnalysis::sciFiStationCut(0., std::vector<int>(1, 1), ch)); // C. No hits in first SciFi plane
    //cutFlow.push_back( new sndAnalysis::sciFiStationCut(0.05, std::vector<int>(1, 5), ch)); // D. Vertex not in 5th wall
    cutFlow.push_back( new sndAnalysis::avgSciFiFiducialCut(200, 1200, 300, 128*12-200, ch)); // E. Average SciFi hit channel number must be within [200, 1200] (ver) and [300, max-200] (hor)
    cutFlow.push_back( new sndAnalysis::avgDSFiducialCut(70, 105, 10, 50, ch)); // F. Average DS hit bar number must be within [70, 105] (ver) and [10, 50] (hor)
    cutFlow.push_back( new sndAnalysis::minSciFiConsecutivePlanes(ch)); // G. At least two consecutive SciFi planes hit
    cutFlow.push_back( new sndAnalysis::DSActivityCut(ch)); // H. If there is a downstream hit, require hits in all upstream stations.    
    if (not isMC) cutFlow.push_back( new sndAnalysis::eventDeltatCut(-1, 100, ch)); // J. Previous event more than 100 clock cycles away. To avoid deadtime issues.
    //    if (not isMC) cutFlow.push_back( new sndAnalysis::eventDeltatCut(+1, 10, ch)); // K. Next event more than 10 clock cycles away. To avoid event builder issues.
  }
    else {
    std::cout << "Unrecognized cutset. Exitting" << std::endl;
    exit(-1);
  }
  //  cutFlow.push_back( new sndAnalysis::minSciFiHits(35, ch)); // G. At least 35 hits in SciFi
  //  cutFlow.push_back( new sndAnalysis::USQDCCut(600, ch)); // J. Total US QDC > 600

  std::cout << "Done initializing cuts" << std::endl;

  int n_cuts = (int) cutFlow.size();

  // Book histograms
  // Cut-by-cut
  // All cut variables
  std::vector<std::vector<TH1D*> > cut_by_cut_var_histos = std::vector<std::vector<TH1D*> >();
  for (int i_cut = -1; i_cut < n_cuts; i_cut++){
    std::vector<TH1D*> this_cut_by_cut_var_histos = std::vector<TH1D*>();
    for (sndAnalysis::baseCut * cut : cutFlow) {
      for(int i_dim = 0; i_dim < cut->getNbins().size(); i_dim++){
	this_cut_by_cut_var_histos.push_back(new TH1D(("h"+std::to_string(i_cut+1)+"_"+cut->getShortName()+"_"+std::to_string(i_dim)).c_str(), 
						      cut->getShortName().c_str(), 
						      cut->getNbins()[i_dim], cut->getRangeStart()[i_dim], cut->getRangeEnd()[i_dim]));
      }
    }
    cut_by_cut_var_histos.push_back(this_cut_by_cut_var_histos);
  }
  
  std::vector<std::vector<TH1D*> > cut_by_cut_truth_histos = std::vector<std::vector<TH1D*> >();
  if (isMC) {
    for (int i_cut = -1; i_cut < n_cuts; i_cut++){
    std::vector<TH1D*> this_cut_by_cut_truth_histos = std::vector<TH1D*>();
    this_cut_by_cut_truth_histos.push_back(new TH1D(("h"+std::to_string(i_cut+1)+"_Emu").c_str(), "Emu", 300, 0, 3000));
    this_cut_by_cut_truth_histos.push_back(new TH1D(("h"+std::to_string(i_cut+1)+"_muX").c_str(), "muX", 200, -100, 0));
    this_cut_by_cut_truth_histos.push_back(new TH1D(("h"+std::to_string(i_cut+1)+"_muY").c_str(), "muY", 200, 0, 100));
    this_cut_by_cut_truth_histos.push_back(new TH1D(("h"+std::to_string(i_cut+1)+"_muZ").c_str(), "muZ", 200, 280, 380));

    cut_by_cut_truth_histos.push_back(this_cut_by_cut_truth_histos);
    }
  }

  // N-1
  std::vector<TH1D*> n_minus_1_var_histos = std::vector<TH1D*>();
  for (sndAnalysis::baseCut * cut : cutFlow) {
    for(int i_dim = 0; i_dim < cut->getNbins().size(); i_dim++){
      n_minus_1_var_histos.push_back(new TH1D(("n_minus_1_"+cut->getShortName()+"_"+std::to_string(i_dim)).c_str(), 
					      cut->getShortName().c_str(), 
					      cut->getNbins()[i_dim], cut->getRangeStart()[i_dim], cut->getRangeEnd()[i_dim]));
    }
  }

  std::cout << "Done initializing histograms" << std::endl;

  // Cut flow
  TH1D * cutFlowHistogram = new TH1D("cutFlow", "Cut flow;;Number of events passing cut", n_cuts+1, 0, n_cuts+1);
  TH1D * cutFlowHistogram2 = new TH1D("cutFlow2", "Cut flow;;Number of events passing cut / fb^{-1}", n_cuts+1, 0, n_cuts+1);
  for (int i = 2; i <= cutFlowHistogram->GetNbinsX(); i++){
    cutFlowHistogram->GetXaxis()->SetBinLabel(i, cutFlow.at(i-2)->getName().c_str());
    cutFlowHistogram2->GetXaxis()->SetBinLabel(i, cutFlow.at(i-2)->getName().c_str());
  }

  std::cout << "Done initializing cut flow histogram" << std::endl;

  // Get number of entries
  unsigned long int n_entries = ch->GetEntries();

  // Holder for cut results
  std::vector<bool> passes_cut = std::vector<bool>(n_cuts, false);

  std::cout << "Starting event loop" << std::endl;
  
  TVector3 muVTX;
  const char *volName = NULL;

  for (unsigned long int i_entry = 0; i_entry < n_entries; i_entry++){
    if ((i_entry == 316629 || i_entry == n_entries-1)  && isPMU == 1) continue; // exclude problemati events from the PassingMuons simulation
    ch->GetEntry(i_entry);
    float WEIGHT = 1.;
    if (i_entry%10000 == 0) std::cout << "Reading entry " << i_entry << std::endl;
    if (isMC && !isPMU)
    {
      Bool_t intTarget = false;
      Float_t muX = ((ShipMCTrack*) MCTracks->At(0))->GetStartX();
      Float_t muY = ((ShipMCTrack*) MCTracks->At(0))->GetStartY();
      Float_t muZ = ((ShipMCTrack*) MCTracks->At(0))->GetStartZ();
      muVTX = TVector3(muX, muY, muZ);
      if(gGeoManager->FindNode(muVTX.X(), muVTX.Y(), muVTX.Z())->GetMotherVolume() == gGeoManager->FindNode(muVTX.X(), muVTX.Y(), muVTX.Z())->GetMotherVolume())
      {
        if(gGeoManager->FindNode(muVTX.X(), muVTX.Y(), muVTX.Z())->GetMotherVolume() != NULL)
        {
          volName = gGeoManager->FindNode(muVTX.X(), muVTX.Y(), muVTX.Z())->GetMotherVolume()->GetName();
          if(strncmp(volName, "volTarget", 9)==0 || strncmp(volName, "Brick", 5)==0 || strncmp(volName, "Scifi", 5)==0 || strncmp(volName, "Fiber", 5)==0)  intTarget = true;
        }
      }
      if(!intTarget) continue;
      float W = 10.256410*((ShipMCTrack*) MCTracks->At(0))->GetWeight(); // Old fluka normalization
      float wLHC = W/10; // W/nMult
      float sigma = 0;
      float wInter = ((ShipMCTrack*) MCTracks->At(2))->GetWeight();
      if (((ShipMCTrack*) MCTracks->At(0))->GetPdgCode() == 13) sigma = crsec->Eval(((ShipMCTrack*) MCTracks->At(0))->GetEnergy());
      else { sigma = crsec2->Eval(((ShipMCTrack*) MCTracks->At(0))->GetEnergy());}
      float wDIS = 0.6E-3*sigma;
      WEIGHT = wLHC*wInter*wDIS*1E5;
    }
    else if (isMC && isPMU==1) {
	WEIGHT = 8E8/5E7*((ShipMCTrack*) MCTracks->At(0))->GetWeight()*1E5;
    	bool DISinpassing = false;
	int ntracks = MCTracks->GetEntries();
	for (int ntrack = 0; ntrack < ntracks; ntrack++){
		if (((ShipMCTrack*) MCTracks->At(ntrack))->GetProcID() == 23) DISinpassing = true;
		break;
	}
	if (DISinpassing) continue; //exclude DIS interactions in passing muon simulations	
    }

    cutFlowHistogram->Fill(0.);
    cutFlowHistogram2->Fill(0., WEIGHT);

    // Apply cuts
    int n_cuts_passed = 0;
    bool accept_event = true;
    int i_cut = 0;
    for (sndAnalysis::baseCut * cut : cutFlow){
		if (strncmp(cut->getShortName().c_str(), "NoVetoHits", 10) == 0 || strncmp(cut->getShortName().c_str(), "SciFiStation", 12) == 0)
		{ 
			if (!cut->passCut()){
			  if (accept_event) {cutFlowHistogram->Fill(1 + i_cut); cutFlowHistogram2->Fill(1+ i_cut, WEIGHT);}
		  	passes_cut[i_cut] = true;
			  n_cuts_passed++;
			}
			else {
			  accept_event = false;
			  passes_cut[i_cut] = false;
			}
		}
		else{
			if (cut->passCut()){
        if (accept_event) {cutFlowHistogram->Fill(1 + i_cut); cutFlowHistogram2->Fill(1+ i_cut, WEIGHT);}
        passes_cut[i_cut] = true;
        n_cuts_passed++;
      } 
      else {
        accept_event = false;
        passes_cut[i_cut] = false;
      }
		}
    i_cut++;
		}
		if (accept_event) outTree->Fill();


    // Fill histograms
    std::vector<TH1D*>::iterator hist_it;
    // Sequential
    for (int seq_cut = -1; seq_cut < ((int) passes_cut.size()); seq_cut++){
      if (seq_cut >= 0){
	if (not passes_cut[seq_cut]) break;
      }
      hist_it = cut_by_cut_var_histos[seq_cut+1].begin();
      for (sndAnalysis::baseCut * cut : cutFlow) {
	for (int i_dim = 0; i_dim < cut->getPlotVar().size(); i_dim++){
	  (*hist_it)->Fill(cut->getPlotVar()[i_dim]);
	  hist_it++;
	}
      }
      if (isMC) {
	      cut_by_cut_truth_histos[seq_cut+1][0]->Fill(((ShipMCTrack*) MCTracks->At(0))->GetEnergy()); // Emu
	      cut_by_cut_truth_histos[seq_cut+1][1]->Fill(((ShipMCTrack*) MCTracks->At(0))->GetStartX()); // X
	      cut_by_cut_truth_histos[seq_cut+1][2]->Fill(((ShipMCTrack*) MCTracks->At(0))->GetStartY()); // Y
	      cut_by_cut_truth_histos[seq_cut+1][3]->Fill(((ShipMCTrack*) MCTracks->At(0))->GetStartZ()); // Z
      }
    }

    // N-1
    int current_cut = 0;
    hist_it = n_minus_1_var_histos.begin();
    for (sndAnalysis::baseCut * cut : cutFlow) {
      for (int i_dim = 0; i_dim < cut->getPlotVar().size(); i_dim++){
	      if (((not passes_cut[current_cut]) and (n_cuts_passed == (cutFlow.size()-1)))
	        or (n_cuts_passed == cutFlow.size()) ) (*hist_it)->Fill(cut->getPlotVar()[i_dim]);
	      hist_it++;
      }
      current_cut++;
    }
  }
  delete file2;
  delete crsec;
  delete crsec2; 
  
  outFile->Write();

  return 0;
}
