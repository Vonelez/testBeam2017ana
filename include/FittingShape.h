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

  TGraphErrors *getMeanGraph() { return meanGraph; };

  TGraphErrors *getSigmaGraph() { return sigmaGraph; };

  TGraphErrors *getNdfGraph() { return ndfGraph; };

  TGraphErrors *getChi2_ndfGraph() { return chiSq_ndfGraph; };

  Int_t getLeftEdge() { return leftEdge; };

  Int_t getRightEdge() { return rightEdge; };

private:
  TH2D *shape;
  TH1D *currentBin;
  TGraphErrors *meanGraph;
  TGraphErrors *sigmaGraph;
  TGraphErrors *ndfGraph;
  TGraphErrors *chiSq_ndfGraph;
  Long_t runNum;
  Double_t leftEdge;
  Double_t rightEdge;
  Double_t leftFitEdge;
  Double_t rightFitEdge;

  void init();

  void runningBins();

  void findingEdges();

  void fitRangeEstimation(TH1D *binHist);
};


#endif //REALDATATESTS_FITTINGSHAPE_H
