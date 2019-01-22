#include "../include/FittingShape.h"

FittingShape::FittingShape(TH2D *shape, Long_t runNum) {
  this->shape = shape;
  this->runNum = runNum;
  init();
  runningBins();
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
  TString runFolder("linarTests_");
  runFolder += runNum;
  TString createFolder = mkd + imgFolder + runFolder;
  system(createFolder);
  TString binsFileName("/binsPool");
  TString ending(".root");
  TFile binsFile(imgFolder + runFolder + binsFileName + ending, "RECREATE");
  Int_t point = 0;

  for (int i = 0; i < shape->GetNbinsX(); ++i) {
    currentBin = shape->ProjectionY("Current bin", i, i);
    if (currentBin->Integral() < 10) continue;
    Double_t leftEdge = currentBin->GetBinCenter(currentBin->FindFirstBinAbove(currentBin->GetMaximum() / 10.));
    Double_t rightEdge = currentBin->GetBinCenter(currentBin->FindLastBinAbove(currentBin->GetMaximum() / 10.));
    currentBin->Fit("gaus", "QR", "SAME", leftEdge, rightEdge);
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

  TFile resultFile(imgFolder + runFolder + "/results" + ending, "RECREATE");
  meanGraph->Write("MeanValue");
  sigmaGraph->Write("Sigma");
  ndfGraph->Write("NDF");
  chiSq_ndfGraph->Write("Chi2/NDF");
}
