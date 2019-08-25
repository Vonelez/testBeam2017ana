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

  virtual TH2D *GetShape();

private:
  TTree *STRAW_EVENT_tree;
  TTree *MAMBA_EVENT_tree;
  TH2D *vshape;

  Double_t conversion;
  Int_t mam_nEntries;
  Int_t N;
  Int_t count_processed;



  virtual void merging();

  virtual void filling_hists(STRAW_presetting *straw, MAMBA_presetting *mamba);

};

#endif