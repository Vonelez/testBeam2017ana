#include <iostream>
#include "include/SettingVShape.h"
#include "include/FittingShape.h"
#include "include/Derivative.h"
#include "include/Resolution.h"

int main() {
  TString dir("/Users/andrew/CLionProjects/shipDataAna/experimentalData/");
  TString fileName("Synchronized_data_run_");
  Long_t runNum = 222;
  TString ending(".root");

  TString RootName(dir + fileName + runNum + ending);

  TFile *AnaFile = TFile::Open(RootName, "read");
  TTree *mt;
  TTree *st;
  AnaFile->GetObject("tree", mt);
  AnaFile->GetObject("ADC1", st);

  Int_t binning = 100;

  SettingVShape *setShape = new SettingVShape(st, mt, binning);

  FittingShape *fittingShape = new FittingShape(setShape->GetShape(), runNum);

  Derivative *derivative = new Derivative(fittingShape->getMeanGraph(), fittingShape->getSigmaGraph());

  Resolution(derivative->getCoordinateResolGraph());

  AnaFile->TFile::Close();

  return 0;
}
