#pragma once

#include "TClonesArray.h"
#include "Scifi.h"

namespace sndAnalysis {

  // Function to get the number of hits per station in SciFi
  void getSciFiHitsPerStation(TClonesArray * digiHits, std::vector<int> &horizontal_hits, std::vector<int> &vertical_hits);

  // Function to get the total number of SciFi hits
  int getTotalSciFiHits(std::vector<int> &horizontal_hits, std::vector<int> &vertical_hits);
  int getTotalSciFiHits(TClonesArray * digiHits);

  // Function to get the SciFi fractional hits per plane
  std::vector<float> getFractionalHitsPerPlane(std::vector<int> &horizontal_hits, std::vector<int> &vertical_hits);
  std::vector<float> getFractionalHitsPerPlane(TClonesArray * digiHits);

  // Function to find the station where the interaction ocurred by checking the station at which the cumulative hit fraction exceeds a threshold
  int findStation(std::vector<int> &horizontal_hits, std::vector<int> &vertical_hits, float threshold);
  int findStation(TClonesArray * digiHits, float threshold);
 
  // Function to get the time interval for the majority of Scifi hits
  void getTimeCorrectedRange(TClonesArray *digiHits, std::unordered_map<int, std::vector<double>> &rangePerStation, Scifi *scifiDet); 
  bool isInTimeRange(Double_t hit_time, Double_t time_low, Double_t time_up);
  void getAvgScifipos(TClonesArray *digiHits, std::unordered_map<int, double> &avg_sf_x, std::unordered_map<int, double> &avg_sf_y, Scifi *scifiDet);
}
