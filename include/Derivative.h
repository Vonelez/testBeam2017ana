#ifndef REALDATATESTS_DERIVATIVE_H
#define REALDATATESTS_DERIVATIVE_H

#include "TGraphErrors.h"
#include "TFile.h"
#include "TCanvas.h"

class Derivative {
public:
  Derivative(TGraphErrors *meanGraph, TGraphErrors *sigmaGraph);
  virtual ~Derivative();
  TGraphErrors * getDerivGraph() { return derivative; }
  TGraphErrors * getCoordinateResolGraph() { return coordinateResolGraph; }

private:
  TGraphErrors *derivative;
  TGraphErrors *meanGraph;
  TGraphErrors *sigmaGraph;
  TGraphErrors *coordinateResolGraph;

  void init();
  void derivativeCalc();
};


#endif //REALDATATESTS_DERIVATIVE_H
