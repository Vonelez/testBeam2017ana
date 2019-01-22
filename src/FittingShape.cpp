#include "../include/FittingShape.h"

FittingShape::FittingShape(TH2D *shape, Long_t runNum) {
  this->shape = shape;
  this->runNum = runNum;
  init();
  runningBins();
  findingEdges();
}

FittingShape::~FittingShape() = default;

void FittingShape::init() {
  currentBin = new TH1D();
  meanGraph = new TGraphErrors();
  sigmaGraph = new TGraphErrors();
  ndfGraph = new TGraphErrors();
  chiSq_ndfGraph = new TGraphErrors();
}

void FittingShape::runningBins() {
  TString imgFolder("img/");
  TString mkd("mkdir -p ");
  system(mkd + imgFolder);
  TString runFolder("run_");
  runFolder += runNum;
  TString createFolder = mkd + imgFolder + runFolder;
  system(createFolder);
  TString binsFileName("/binsPool");
  TString ending(".root");
  TFile binsFile(imgFolder + runFolder + binsFileName + ending, "RECREATE");
  Int_t point = 0;

  for (int i = 1; i < shape->GetNbinsX(); ++i) {
    currentBin = shape->ProjectionY("Current bin", i, i);
    if (currentBin->Integral() < 100) continue;
    fitRangeEstimation(currentBin);
    currentBin->Fit("gaus", "QR", "SAME", leftFitEdge, rightFitEdge);
    TF1 *fitFunc = currentBin->GetFunction("gaus");
    Double_t xAxisPoint = shape->GetXaxis()->GetBinCenter(i);
    meanGraph->SetPoint(point, xAxisPoint, fitFunc->GetParameter(1));
    meanGraph->SetPointError(point, 0., fitFunc->GetParError(1));
    sigmaGraph->SetPoint(point, xAxisPoint, fitFunc->GetParameter(2));
    sigmaGraph->SetPointError(point, 0., fitFunc->GetParError(2));
    ndfGraph->SetPoint(point, xAxisPoint, fitFunc->GetNDF());
    chiSq_ndfGraph->SetPoint(point, xAxisPoint, fitFunc->GetChisquare() / fitFunc->GetNDF());

    TString binName("bin_XAxis_");
    TString num;
    num.Form("%f", xAxisPoint);
    TString fullBinName = binName + num;
    currentBin->GetXaxis()->SetTitle("T (ns)");
    currentBin->GetXaxis()->SetTitleOffset(1.4);
    currentBin->GetYaxis()->SetTitle("N entries");
    currentBin->Write(fullBinName);

    point++;
  }

}

void FittingShape::findingEdges() {
  leftEdge = 0;
  rightEdge = ndfGraph->GetN();
  Double_t *ndfVal;
  ndfVal = ndfGraph->GetY();

  for (int i = 1; i < ndfGraph->GetN() / 2; ++i) {
    if (abs(ndfVal[i] - ndfVal[i - 1]) > 20) {
      leftEdge = i;
      break;
    }
  }

  for (int j = ndfGraph->GetN() - 1; j > ndfGraph->GetN() / 2; --j) {
    if (abs(ndfVal[j] - ndfVal[j + 1]) > 20) {
      rightEdge = j;
      break;
    }
  }
}

void FittingShape::fitRangeEstimation(TH1D *binHist) {
  Double_t maxValue = binHist->GetMaximum();
  Int_t maxValueBin = binHist->GetMaximumBin();
  leftFitEdge = binHist->GetBinCenter(binHist->FindFirstBinAbove());
  rightFitEdge = binHist->GetBinCenter(binHist->FindLastBinAbove());

  for (int i = maxValueBin - 1; i > 0; --i) {
    if (binHist->GetBinContent(i) < maxValue / 10) {
      leftFitEdge = binHist->GetBinCenter(i);
      break;
    }
  }

  for (int j = maxValueBin + 1; j < binHist->FindLastBinAbove(); ++j) {
    if (binHist->GetBinContent(j) < maxValue / 10) {
      rightFitEdge = binHist->GetBinCenter(j);
      break;
    }
  }
}
