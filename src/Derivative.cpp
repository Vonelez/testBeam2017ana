#include "../include/Derivative.h"

Derivative::Derivative(TGraphErrors *meanGraph, TGraphErrors *sigmaGraph, Int_t leftEdge, Int_t rightEdge) {
  this->meanGraph = meanGraph;
  this->sigmaGraph = sigmaGraph;
  this->leftEdge = leftEdge;
  this->rightEdge = rightEdge;
  derivative = new TGraphErrors();
  coordinateResolGraph = new TGraphErrors();
  derivativeCalc();

}

Derivative::~Derivative() = default;

void Derivative::derivativeCalc() {
  Double_t *x, *y, *ySigma;
  x = meanGraph->GetX();
  y = meanGraph->GetY();
  ySigma = sigmaGraph->GetY();
  Int_t point = 0;
  Double_t dX = 0;
  Double_t derivativeValue = 0;
  for (int i = leftEdge + 2; i < rightEdge - 2; i++) {
    derivativeValue = ((y[i - 1] - y[i - 2]) / (x[i - 1] - x[i - 2]) + (y[i] - y[i - 1]) / (x[i] - x[i - 1]) +
                       (y[i + 1] - y[i]) / (x[i + 1] - x[i]) + (y[i + 2] - y[i + 1]) / (x[i + 2] - x[i + 1])) / 4;
    dX = x[i] - x[i-1];
    derivative->SetPoint(point, x[i], derivativeValue);

    Double_t dY_errorSquare1 = meanGraph->GetErrorY(i - 2) * meanGraph->GetErrorY(i - 2) +
                               meanGraph->GetErrorY(i - 1) * meanGraph->GetErrorY(i - 1);
    Double_t dY_errorSquare2 = meanGraph->GetErrorY(i - 1) * meanGraph->GetErrorY(i - 1) +
                               meanGraph->GetErrorY(i) * meanGraph->GetErrorY(i);
    Double_t dY_errorSquare3 = meanGraph->GetErrorY(i) * meanGraph->GetErrorY(i) +
                               meanGraph->GetErrorY(i + 1) * meanGraph->GetErrorY(i + 1);
    Double_t dY_errorSquare4 = meanGraph->GetErrorY(i + 1) * meanGraph->GetErrorY(i + 1) +
                               meanGraph->GetErrorY(i + 2) * meanGraph->GetErrorY(i + 2);
    Double_t derivativeError = (1. / 4.) * sqrt((1 / dX) * (1 / dX) *
                                                (dY_errorSquare1 + dY_errorSquare2 + dY_errorSquare3 +
                                                 dY_errorSquare4));
    derivative->SetPointError(point, 0., derivativeError);

    Double_t coordResolution = ySigma[i] / abs(derivativeValue);
    Double_t coordResolutionError =
            (ySigma[i] / (derivativeValue * derivativeValue)) * (ySigma[i] / (derivativeValue * derivativeValue)) *
            derivativeError * derivativeError +
            (1 / derivativeValue) * (1 / derivativeValue) * sigmaGraph->GetErrorY(i) * sigmaGraph->GetErrorY(i);
    coordResolutionError = sqrt(coordResolutionError);
    coordinateResolGraph->SetPoint(point, x[i], coordResolution);
    coordinateResolGraph->SetPointError(point, 0., coordResolutionError);

    point++;
  }

}


