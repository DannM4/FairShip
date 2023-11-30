#include "sndAvgSciFiFiducialCut.h"
#include "sndSciFiTools.h"
#include "Scifi.h"
#include "TChain.h"

namespace sndAnalysis{
  avgSciFiFiducialCut::avgSciFiFiducialCut(double vertical_min_cut, double vertical_max_cut, double horizontal_min_cut, double horizontal_max_cut, TChain * tree, bool reverseCuts, bool setisMC, Scifi *setscifiDet) : sciFiBaseCut(tree){

    reversed = reverseCuts;
    isMC = setisMC;
    scifiDet = setscifiDet;

    vertical_min = vertical_min_cut;
    vertical_max = vertical_max_cut;
    horizontal_min = horizontal_min_cut;
    horizontal_max = horizontal_max_cut;
    
    cutName = "Avg SciFi Ver channel in ["+std::to_string(vertical_min)+","+std::to_string(vertical_max)+"] Hor in ["+std::to_string(horizontal_min)+","+std::to_string(horizontal_max)+"]";

    shortName = "AvgSFChan";
    nbins = std::vector<int>{128*2, 128*2};
    range_start = std::vector<double>{0, 0};
    range_end = std::vector<double>{128*12, 128*12};
    plot_var = std::vector<double>{-1, -1};
  }

  bool avgSciFiFiducialCut::passCut(){
    initializeEvent();
    
    double avg_ver = 0.;
    unsigned int n_ver = 0;
    double avg_hor = 0.;
    unsigned int n_hor = 0;
    //std::cout << "isMC set to " << isMC << std::endl; 
    sndScifiHit * hit;
    TIter hitIterator(scifiDigiHitCollection);
    std::unordered_map<int, std::vector<double>> rangePerStation;
    std::unordered_map<int, double> avg_sf_x;
    std::unordered_map<int, double> avg_sf_y;
    TVector3 a, b;
    if (!isMC){
      sndAnalysis::getAvgScifipos(scifiDigiHitCollection, avg_sf_x, avg_sf_y, scifiDet);
      sndAnalysis::getTimeCorrectedRange(scifiDigiHitCollection, rangePerStation, scifiDet);
    }
    while ( (hit = (sndScifiHit*) hitIterator.Next()) ){
      if (hit->isValid()){
        if(!isMC){
          const Double_t TDC2ns = 1E9/160.316E6;
          scifiDet->GetSiPMPosition(hit->GetDetectorID(), a, b);
	        Double_t L;
          int sfplane = hit->GetStation();
          if (hit->isVertical()) L = b.Y()-avg_sf_y[sfplane];
          else {L = avg_sf_x[sfplane]-a.X();}
          Double_t hit_time = scifiDet->GetCorrectedTime(hit->GetDetectorID(), hit->GetTime()*TDC2ns, L);
          auto it = rangePerStation.find(sfplane);
          bool InTime = sndAnalysis::isInTimeRange(hit_time, it->second[0], it->second[1]);
          if(!InTime) continue;
        }
	int mat = hit->GetMat();
	int sipm = hit->GetSiPM();
	int channel = hit->GetSiPMChan();

	int x = channel + sipm*128 + mat*4*128;

	if (hit->isVertical()){
	  avg_ver += x;
	  n_ver++;
	} else {
	  avg_hor += x;
	  n_hor++;
	}
      }
    }
    
    if ((n_ver+n_hor) == 0) {
      plot_var[0] = -1;
      plot_var[1] = -1;
      return false;
    }
    
    if (n_ver) {
      avg_ver /= n_ver;
      plot_var[0] = avg_ver;
    } else {
      plot_var[0] = -1;
    }

    if (n_hor) {
      avg_hor /= n_hor;
      plot_var[1] = avg_hor;
    } else {
      plot_var[1] = -1;
    }

    if (n_ver == 0) return false;
    if (n_hor == 0) return false;

    if (not reversed) {
      if (avg_hor < horizontal_min) return false;
      if (avg_hor > horizontal_max) return false;
      if (avg_ver < vertical_min) return false;
      if (avg_ver > vertical_max) return false;
    } else {
      if ((avg_hor > horizontal_min) and (avg_hor < horizontal_max)) return false;
      if ((avg_ver > vertical_min) and (avg_ver < vertical_max)) return false;
    }
    return true;
  }
}	     
