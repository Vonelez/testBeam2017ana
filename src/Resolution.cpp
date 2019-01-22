#include "../include/Resolution.h"

Resolution::Resolution(TGraphErrors *resolGraph) {
  this->resolGraph = resolGraph;
  init();
  weightedMean();
}

Resolution::~Resolution() = default;

void Resolution::init() {

}

void Resolution::weightedMean() {
  Double_t *y;
  y = resolGraph->GetY();
  Double_t sum = 0;
  Double_t weightsSum = 0;
  Double_t resolution = 0;
  Double_t resolutionError = 0;

  for (int i = 0; i < resolGraph->GetN(); ++i) {
    Double_t error = resolGraph->GetErrorY(i);
    if (y[i] < 0.05) continue;
    sum += y[i] / (error * error);
    weightsSum += 1 / (error * error);
  }

  resolution = sum / weightsSum;
  resolutionError = sqrt(1 / weightsSum);

  std::cout << "---> Resolution (" << resolution << " +/- " << resolutionError << ") um \n";
}