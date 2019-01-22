#ifndef REALDATATESTS_FITTINGSHAPE_H
#define REALDATATESTS_FITTINGSHAPE_H

#include "TH2.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TFile.h"

class FittingShape {
public:
  FittingShape(TH2D *shape, Long_t runNum);
  virtual ~FittingShape();
  TGraphErrors * getMeanGraph() { return  meanGraph; };
  TGraphErrors * getSigmaGraph() { return sigmaGraph; };

private:
  TH2D *shape;
  TH1D *currentBin;
  TGraphErrors *meanGraph;
  TGraphErrors *sigmaGraph;
  TGraphErrors *ndfGraph;
  TGraphErrors *chiSq_ndfGraph;
  Long_t runNum;

  void init();
  void runningBins();
};


#endif //REALDATATESTS_FITTINGSHAPE_H
