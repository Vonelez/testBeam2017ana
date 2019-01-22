#ifndef REALDATATESTS_DERIVATIVE_H
#define REALDATATESTS_DERIVATIVE_H

#include "TGraphErrors.h"
#include "TFile.h"
#include "TCanvas.h"

class Derivative {
public:
  Derivative(TGraphErrors *meanGraph, TGraphErrors *sigmaGraph, Int_t leftEdge, Int_t rightEdge);

  virtual ~Derivative();

  TGraphErrors *getDerivGraph() { return derivative; }

  TGraphErrors *getCoordinateResolGraph() { return coordinateResolGraph; }

private:
  TGraphErrors *derivative;
  TGraphErrors *meanGraph;
  TGraphErrors *sigmaGraph;
  TGraphErrors *coordinateResolGraph;

  Int_t leftEdge;
  Int_t rightEdge;

  void init();

  void derivativeCalc();
};


#endif //REALDATATESTS_DERIVATIVE_H
