#include <iostream>
#include "include/SettingVShape.h"
#include "include/FittingShape.h"
#include "include/Derivative.h"
#include "include/Resolution.h"
#include "include/QA.h"

int main() {
//  TString dir("/Volumes/Andrew's Seagate Drive/SHiP_DATA/");
  TString dir("/eos/experiment/ship/TestBeam/StrawTracker/2017/H2-Tests/synchronized_trees/");
  TString fileName("Synchronized_data_run_");
  Long_t runNum = 164;
  TString ending(".root");

  TString RootName(dir + fileName + runNum + ending);
  TFile *AnaFile;
  AnaFile = TFile::Open(RootName, "read");
  TTree *mt;
  TTree *st;
  AnaFile->GetObject("tree", mt);
  AnaFile->GetObject("ADC1", st);

  Int_t binning = 100;

  auto *setShape = new SettingVShape(st, mt, binning);

  std::cout << "I did /setShape/" << std::endl;

  auto *fittingShape = new FittingShape(setShape->GetShape(), runNum);

  std::cout << "I did /fittingShape/" << std::endl;

  auto *derivative = new Derivative(fittingShape->getMeanGraph(), fittingShape->getSigmaGraph(), fittingShape->getLeftEdge(), fittingShape->getRightEdge());

  std::cout << "I did /derivative/" << std::endl;

  Resolution(derivative->getCoordinateResolGraph());

  std::cout << "I did /Resolution/" << std::endl;

  QA(setShape, fittingShape, derivative, runNum);

  std::cout << "I did /QA/" << std::endl;

//  AnaFile->TFile::Close();

  return 1;
}
