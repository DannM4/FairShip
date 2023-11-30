#pragma once

#include "sndSciFiBaseCut.h"

#include "TChain.h"
#include "sndScifiHit.h"

class Scifi;

namespace sndAnalysis {
  class avgSciFiFiducialCut : public sciFiBaseCut {
  private :
    double vertical_min, vertical_max, horizontal_min, horizontal_max;
    bool reversed, isMC;
    Scifi *scifiDet;
  public :
    avgSciFiFiducialCut(double vertical_min_cut, double vertical_max_cut, double horizontal_min_cut, double horizontal_max_cut, TChain * tree, bool reverseCuts = false, bool isMC = false, Scifi * scifiDet = nullptr);
    ~avgSciFiFiducialCut(){;}

    bool passCut();

  };
};
