#ifndef ANALYSIS_SETTING_VSHAPE_H
#define ANALYSIS_SETTING_VSHAPE_H

#include "../include/MAMBA_presetting.h"
#include "../include/STRAW_presetting.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"

class SettingVShape {
  public:
    SettingVShape(TTree *STRAW_EVENT_tree, TTree *MAMBA_EVENT_tree, Int_t binning);
    virtual ~SettingVShape();
    virtual TH2D * GetShape();
  
  private:
    TTree *STRAW_EVENT_tree;
    TTree *MAMBA_EVENT_tree;
    TH2D *vshape;
    TGraphErrors *testingLinarIdea_1;
    TGraphErrors *testingLinarIdea_2;
    TGraphErrors *testingLinarIdea_3;
    TGraphErrors *derivativePoint1;
    TGraphErrors *derivativePoint2;
    TGraphErrors *derivativePoint3;
    Int_t binning;
    Double_t conversion;
    Int_t mam_nEntries;
    Int_t N;
    Int_t count_processed;
    Double_t derivative;
    Double_t derivativeError;
    Double_t maxPercent;

    virtual void init();
    virtual void merging(TTree *STRAW_EVENT_tree, TTree *MAMBA_EVENT_tree);
    virtual void filling_hists(STRAW_presetting *straw, MAMBA_presetting *mamba);
    virtual void tetstingBinning();
    virtual void derivativeCalc(Int_t interestedBin, Int_t binningInI, Double_t tPoint);
    virtual bool fitIsNice(TF1* fitFunction);

};

#endif