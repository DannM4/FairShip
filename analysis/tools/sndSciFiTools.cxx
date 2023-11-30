#include "sndSciFiTools.h"

#include <numeric>
#include <algorithm>
#include "TH1D.h"
#include "TClonesArray.h"
#include "sndScifiHit.h"
#include "Scifi.h"

void sndAnalysis::getSciFiHitsPerStation(TClonesArray * digiHits, std::vector<int> &horizontal_hits, std::vector<int> &vertical_hits){

  // Clear hits per plane vectors
  std::fill(horizontal_hits.begin(), horizontal_hits.end(), 0);
  std::fill(vertical_hits.begin(), vertical_hits.end(), 0);

  // Add valid hits to hits per plane vectors
  sndScifiHit * hit;
  TIter hitIterator(digiHits);

  while ( (hit = (sndScifiHit*) hitIterator.Next()) ){
    if (hit->isValid()){
      int sta = hit->GetStation();
      if (hit->isVertical()){
	vertical_hits[sta-1]++;
      } else {
	horizontal_hits[sta-1]++;
      }
    }
  }
}

int sndAnalysis::getTotalSciFiHits(std::vector<int> &horizontal_hits, std::vector<int> &vertical_hits){
  return std::accumulate(horizontal_hits.begin(),
			 horizontal_hits.end(),
			 std::accumulate(vertical_hits.begin(),
					 vertical_hits.end(), 0));
}

int sndAnalysis::getTotalSciFiHits(TClonesArray * digiHits){
  std::vector<int> horizontal_hits = std::vector<int>(5);
  std::vector<int> vertical_hits = std::vector<int>(5);
    
  getSciFiHitsPerStation(digiHits, horizontal_hits, vertical_hits);

  return getTotalSciFiHits(horizontal_hits, vertical_hits);
}

std::vector<float> sndAnalysis::getFractionalHitsPerPlane(std::vector<int> &horizontal_hits, std::vector<int> &vertical_hits){
    
  int total_hits = getTotalSciFiHits(horizontal_hits, vertical_hits);
    
  std::vector<float> fractional_hits_per_station = std::vector<float>(horizontal_hits.size());
    
  std::transform(horizontal_hits.begin(), horizontal_hits.end(),
		 vertical_hits.begin(),
		 fractional_hits_per_station.begin(),
		 [&total_hits](const auto& hor, const auto& ver){
		   return ((float) hor + ver)/total_hits;
		 });
    
  return fractional_hits_per_station;
    
}

std::vector<float> sndAnalysis::getFractionalHitsPerPlane(TClonesArray * digiHits){
  std::vector<int> horizontal_hits = std::vector<int>(5);
  std::vector<int> vertical_hits = std::vector<int>(5);
    
  getSciFiHitsPerStation(digiHits, horizontal_hits, vertical_hits);
    
  return getFractionalHitsPerPlane(horizontal_hits, vertical_hits);
}
  
  
int sndAnalysis::findStation(std::vector<int> &horizontal_hits, std::vector<int> &vertical_hits, float threshold){

  std::vector<float> frac = getFractionalHitsPerPlane(horizontal_hits, vertical_hits);
    
  std::vector<float> frac_sum = std::vector<float>(frac.size());
    
  std::partial_sum(frac.begin(), frac.end(), frac_sum.begin());

  std::vector<float>::iterator station = std::find_if(frac_sum.begin(), frac_sum.end(),
						      [&threshold](const auto& f) {
							return f > threshold;
						      });
    
  return station - frac_sum.begin() + 1;
    
}

int sndAnalysis::findStation(TClonesArray * digiHits, float threshold){
  std::vector<int> horizontal_hits = std::vector<int>(5);
  std::vector<int> vertical_hits = std::vector<int>(5);
    
  getSciFiHitsPerStation(digiHits, horizontal_hits, vertical_hits);

  return findStation(horizontal_hits, vertical_hits, threshold);
}

void sndAnalysis::getTimeCorrectedRange(TClonesArray *digiHits, std::unordered_map<int, std::vector<double>> &rangePerStation, Scifi *scifiDet){
    std::unordered_map<int, double> avg_sf_x;
    std::unordered_map<int, double> avg_sf_y;
    const Double_t TDC2ns = 1E9/160.316E6;
    sndAnalysis::getAvgScifipos(digiHits, avg_sf_x, avg_sf_y, scifiDet);
    for ( int iplane = 1; iplane < 6; iplane++ ){
        TVector3 a, b;
        TH1D hist_time = TH1D("hist_time", "", 20, 0, 50);
        sndScifiHit * hit;
        TIter hitIterator(digiHits);
        while ( (hit = (sndScifiHit*) hitIterator.Next()) ){
            if (!hit->isValid()) continue;
            if (hit->GetStation()!=iplane) continue;
            //hist_time.Fill(hit->GetTime()*TDC2ns);
            scifiDet->GetSiPMPosition(hit->GetDetectorID(), a, b);
            Double_t L;
            if (hit->isVertical()) L = b.Y()-avg_sf_y[iplane];
            else {L = avg_sf_x[iplane]-a.X();}
            Double_t hit_time = scifiDet->GetCorrectedTime(hit->GetDetectorID(), hit->GetTime()*TDC2ns, L);
            hist_time.Fill(hit_time);
            //cout << "Time " << hit->GetTime()*TDC2ns << endl;
        }
    int ibin = -1;
    ibin = hist_time.GetMaximumBin();
    //cout << "ibin " << ibin << " lowedge " << hist_time.GetBinLowEdge(ibin+1) << endl;
    rangePerStation[iplane].push_back(hist_time.GetBinLowEdge(ibin));
    rangePerStation[iplane].push_back(hist_time.GetBinLowEdge(ibin+3)); // considering 2 bins after the maximum bin
    hist_time.Reset("M");
    }
}

bool sndAnalysis::isInTimeRange(Double_t hit_time, Double_t time_low, Double_t time_up){
    if (hit_time > time_low && hit_time <= time_up) return true;
    else return false;
}

void sndAnalysis::getAvgScifipos(TClonesArray *digiHits, std::unordered_map<int, double> &avg_sf_x, std::unordered_map<int, double> &avg_sf_y, Scifi *scifiDet){
  std::unordered_map<int, int> n_sf_hits_x;
  std::unordered_map<int, int> n_sf_hits_y;
  TVector3 a, b;
  sndScifiHit * hit;
  TIter hitIterator(digiHits);
  while ( (hit = (sndScifiHit*) hitIterator.Next()) ){
    if (!hit->isValid()) continue;
    Int_t plane = hit->GetStation();
    Int_t detID = hit->GetDetectorID();
    scifiDet->GetSiPMPosition(detID, a, b);
    if (hit->isVertical()){
        n_sf_hits_x[plane] +=1;
        avg_sf_x[plane] += (a.X() + b.X())/2.;
        //cout << "hit " << hit->isVertical() << " pos with TV3 " << a.X() << " " << b.X() << endl;
    }
    else{
        n_sf_hits_y[plane] +=1;
        avg_sf_y[plane] += (a.Y() + b.Y())/2.;
        //cout << "hit " << hit->isVertical() << " pos with TV3 " << a.Y() << " " << b.Y() << endl;
    }       
  }
  for (int iplane=1;iplane<6;iplane++){
    if (n_sf_hits_x[iplane]){
        avg_sf_x[iplane] /= n_sf_hits_x[iplane];
    }
    if (n_sf_hits_y[iplane]){
        avg_sf_y[iplane] /= n_sf_hits_y[iplane];
    }
  }
}

