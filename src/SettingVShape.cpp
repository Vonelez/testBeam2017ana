#include "../include/SettingVShape.h"

SettingVShape::SettingVShape(TTree *STRAW_EVENT_tree, TTree *MAMBA_EVENT_tree, Int_t binning) {
  this->STRAW_EVENT_tree = STRAW_EVENT_tree;
  this->MAMBA_EVENT_tree = MAMBA_EVENT_tree;
  conversion = 25.0 / 1024;
  mam_nEntries = 0;
  N = 0;
  count_processed = 0;
  vshape = new TH2D("InitVShape", "InitVShape", 24 * 1000 / binning, -12.0, 12.0, 550, -100.0, 1000.0);
  vshape->GetXaxis()->SetTitle("U (mm)");
  vshape->GetYaxis()->SetTitle("T (ns)");
  merging();
}

SettingVShape::~SettingVShape() = default;


void SettingVShape::merging() {
  auto *straw = new STRAW_presetting(STRAW_EVENT_tree);
  auto *mamba = new MAMBA_presetting(MAMBA_EVENT_tree);

  mam_nEntries = (Int_t) MAMBA_EVENT_tree->GetEntries(); //? (Int_t) ?
//  N = mam_nEntries;
  N = 10000;
  cout << "requesting " << N << " events..." << endl;

  cout << "Begin analysis..." << endl;
  for (Int_t i = 0; i < N; i++) {
/// progress bar
    int barWidth = 70;
    cout << "\033[0;32m[\033[0m";
    auto pos = (int) (barWidth * ((float) i / (float) N));
    for (int bar = 0; bar < barWidth; ++bar) {
      if (bar < pos) cout << "\033[0;32m=\033[0m";
      else if (bar == pos) cout << "\033[0;32m>\033[0m";
      else cout << " ";
    }
    cout << "\033[0;32m]\033[0m " << "\033[0;31m" << int(((float) i / (float) N) * 100.0 + 1) << "%\r \033[0m";
    cout.flush();
/// progress bar

    STRAW_EVENT_tree->GetEntry(i);
    mamba->CleanArrays();
    MAMBA_EVENT_tree->GetEntry(i);

    filling_hists(straw, mamba);

    count_processed++;
  }
  cout << endl;
  delete straw, mamba;

}

void SettingVShape::filling_hists(STRAW_presetting *straw, MAMBA_presetting *mamba) {
  if (mamba->ntrk == 1) {
    if ((straw->nt00 == 1) && (straw->nt08 == 1)) {
      vshape->Fill(mamba->spUPos[0] * 10, conversion * (straw->tt00 - straw->tt08));
    }
  }
}

TH2D *SettingVShape::GetShape() {
  return vshape;
}