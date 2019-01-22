#include "../include/SettingVShape.h"

SettingVShape::SettingVShape(TTree *STRAW_EVENT_tree, TTree *MAMBA_EVENT_tree, Int_t binning) {
    this->binning = binning;
    this->STRAW_EVENT_tree = STRAW_EVENT_tree;
    this->MAMBA_EVENT_tree = MAMBA_EVENT_tree;
    init();
    merging(STRAW_EVENT_tree, MAMBA_EVENT_tree);
}

SettingVShape::~SettingVShape() = default;

void SettingVShape::init() {
    conversion = 25.0 / 1024;
    mam_nEntries = 0;
    N = 0;
    maxPercent = 0.1;
    count_processed = 0;
    Int_t bins = 22 * 1000 / binning;
    derivative = 0;
    derivativeError = 0;
    vshape = new TH2D("InitVShape", "InitVShape", bins, -11.0, 11.0, 550, -100.0, 1000.0);
    vshape->GetXaxis()->SetTitle("U (mm)");
    vshape->GetYaxis()->SetTitle("T (ns)");
}

void SettingVShape::merging(TTree *STRAW_EVENT_tree, TTree *MAMBA_EVENT_tree) {
    auto *straw = new STRAW_presetting(STRAW_EVENT_tree);
    auto *mamba = new MAMBA_presetting(MAMBA_EVENT_tree);

    mam_nEntries = (Int_t)MAMBA_EVENT_tree->GetEntries();
    N = mam_nEntries;
//    N = 10000;
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
        MAMBA_EVENT_tree->GetEntry(i);

        filling_hists(straw, mamba);

        count_processed++;
    }
    cout << endl;

    //tetstingBinning();

}

void SettingVShape::filling_hists(STRAW_presetting *straw, MAMBA_presetting *mamba) {
    if (mamba->ntrk == 1) {
        if ((straw->nt00 == 1) && (straw->nt08 == 1)) {
            vshape->Fill(mamba->spUPos[0] * 10, conversion * (straw->tt00 - straw->tt08));
        }
    }
}

TH2D * SettingVShape::GetShape() {
    return vshape;
}

void SettingVShape::tetstingBinning() {
    Int_t interestedBins[3] = {1500, 6000, 8500};
    testingLinarIdea_1 = new TGraphErrors();
    derivativePoint1 = new TGraphErrors();
    testingLinarIdea_2 = new TGraphErrors();
    derivativePoint2 = new TGraphErrors();
    testingLinarIdea_3 = new TGraphErrors();
    derivativePoint3 = new TGraphErrors();
    Double_t point1 = vshape->GetXaxis()->GetBinCenter(interestedBins[0]);
    Double_t point2 = vshape->GetXaxis()->GetBinCenter(interestedBins[1]);
    Double_t point3 = vshape->GetXaxis()->GetBinCenter(interestedBins[2]);

    cout << "creating main data folder ..." << endl;
    TString folder("img/");
    TString mkd("mkdir -p ");
    system(mkd + folder);
    TString newFolder("linarTests_222-224");
    TString createFolder = mkd + folder + newFolder;
    system(createFolder);

    cout << "setting TGraphs ..." << endl;
    TString name = {"#sigma_{T} = f(BinSize) for U = "};
    TString derivName = {"Derivative = f(BinSize) for U = "};
    TString num1;
    num1.Form("%f", point1);
    TString full1 = name + num1;
    testingLinarIdea_1->GetXaxis()->SetTitle("BinSize (#mum)");
    testingLinarIdea_1->GetYaxis()->SetTitle("#sigma_{T} (ns)");
    testingLinarIdea_1->SetTitle(full1);
    TString derivFull1 = derivName + num1;
    derivativePoint1->GetXaxis()->SetTitle("BinSize (#mum)");
    derivativePoint1->GetYaxis()->SetTitle("Dervative");
    derivativePoint1->SetTitle(derivFull1);

    TString num2;
    num2.Form("%f", point2);
    TString full2 = name + num2;
    testingLinarIdea_2->GetXaxis()->SetTitle("BinSize (#mum)");
    testingLinarIdea_2->GetYaxis()->SetTitle("#sigma_{T} (ns)");
    testingLinarIdea_2->SetTitle(full2);
    TString derivFull2 = derivName + num2;
    derivativePoint2->GetXaxis()->SetTitle("BinSize (#mum)");
    derivativePoint2->GetYaxis()->SetTitle("Dervative");
    derivativePoint2->SetTitle(derivFull2);

    TString num3;
    num3.Form("%f", point3);
    TString full3 = name + num3;
    testingLinarIdea_3->GetXaxis()->SetTitle("BinSize (#mum)");
    testingLinarIdea_3->GetYaxis()->SetTitle("#sigma_{T} (ns)");
    testingLinarIdea_3->SetTitle(full3);
    TString derivFull3 = derivName + num3;
    derivativePoint3->GetXaxis()->SetTitle("BinSize (#mum)");
    derivativePoint3->GetYaxis()->SetTitle("Dervative");
    derivativePoint3->SetTitle(derivFull3);



    cout << "Starting main algorithm ... " << endl; 
    Int_t counter = 0;

    cout << "---- Bin 1 ----" << endl;

    TString binFolder1("/bin_");
    binFolder1 += interestedBins[0];
    TString createFolder1 = mkd + folder + newFolder + binFolder1;
    system(createFolder1);

    for(Int_t i = 10; i < 40; i++) {
        if (i % 5 == 0) cout << i << " / 40" << endl;
        TH1D *testedBin = new TH1D();
        testedBin = vshape->ProjectionY("Interested Area 1", interestedBins[0] - i, interestedBins[0] + i );
        if (testedBin->Integral() < 10) continue;
        Double_t leftEdge = testedBin->GetBinCenter(testedBin->FindFirstBinAbove(testedBin->GetMaximum() / 10.));
        Double_t rightEdge = testedBin->GetBinCenter(testedBin->FindLastBinAbove(testedBin->GetMaximum() / 10.));
        testedBin->Fit("gaus", "QR", "SAME", leftEdge, rightEdge);
        TF1 *testFunc = testedBin->GetFunction("gaus");
        testingLinarIdea_1->SetPoint(counter, binning * (1 + 2 * i), testFunc->GetParameter(2));
        testingLinarIdea_1->SetPointError(counter, 0., testFunc->GetParError(2));
        Double_t tMeanPoint1 = testFunc->GetParameter(1);
        derivativeCalc(interestedBins[0], i, tMeanPoint1);
        derivativePoint1->SetPoint(counter, binning * (1 + 2 * i), derivative);
        derivativePoint1->SetPointError(counter, 0., derivativeError);
        counter++;

        TString histName("/binning_");
        Long_t num = binning * (1 + 2 * i);
        TString name = folder + newFolder + binFolder1 + histName;

        TString end = {".root"};
        TString full = name + num + end;
        TFile binFile(full, "RECREATE");
        testedBin->Write();
    }

    cout << "---- Bin 2 ----" << endl;
    counter = 0;

    TString binFolder2("/bin_");
    binFolder2 += interestedBins[1];
    TString createFolder2 = mkd + folder + newFolder + binFolder2;
    system(createFolder2);

    for(Int_t i = 10; i < 40; i++) {
        if (i % 5 == 0) cout << i << " / 40" << endl;
        TH1D *testedBin = new TH1D();
        testedBin = vshape->ProjectionY("Interested Area 2", interestedBins[1] - i, interestedBins[1] + i );
        if (testedBin->Integral() < 10) continue;
        Double_t leftEdge = testedBin->GetBinCenter(testedBin->FindFirstBinAbove(testedBin->GetMaximum() / 10.));
        Double_t rightEdge = testedBin->GetBinCenter(testedBin->FindLastBinAbove(testedBin->GetMaximum() / 10.));
        testedBin->Fit("gaus", "QR", "SAME", leftEdge, rightEdge);
        TF1 *testFunc = testedBin->GetFunction("gaus");
        testingLinarIdea_2->SetPoint(counter, binning * (1 + 2 * i), testFunc->GetParameter(2));
        testingLinarIdea_2->SetPointError(counter, 0., testFunc->GetParError(2));
        Double_t tMeanPoint2 = testFunc->GetParameter(1);
        derivativeCalc(interestedBins[1], i, tMeanPoint2);
        derivativePoint2->SetPoint(counter, binning * (1 + 2 * i), derivative);
        derivativePoint2->SetPointError(counter, 0., derivativeError);
        counter++;

        TString histName("/binning_");
        Long_t num = binning * (1 + 2 * i);
        TString name = folder + newFolder + binFolder2 + histName;

        TString end = {".root"};
        TString full = name + num + end;
        TFile binFile(full, "RECREATE");
        testedBin->Write();
    }

    cout << "---- Bin 3 ----" << endl;
    counter = 0;
    
    TString binFolder3("/bin_");
    binFolder3 += interestedBins[2];
    TString createFolder3 = mkd + folder + newFolder + binFolder3;
    system(createFolder3);

    for(Int_t i = 10; i < 40; i++) {
        if (i % 5 == 0) cout << i << " / 40" << endl;
        TH1D *testedBin = new TH1D();
        testedBin = vshape->ProjectionY("Interested Area 3", interestedBins[2] - i, interestedBins[2] + i );
        if (testedBin->Integral() < 10) continue;
        Double_t leftEdge = testedBin->GetBinCenter(testedBin->FindFirstBinAbove(testedBin->GetMaximum() * maxPercent));
        Double_t rightEdge = testedBin->GetBinCenter(testedBin->FindLastBinAbove(testedBin->GetMaximum() * maxPercent));
        testedBin->Fit("gaus", "QR", "SAME", leftEdge, rightEdge);
        TF1 *testFunc = testedBin->GetFunction("gaus");
        while (!fitIsNice(testFunc)) {
            maxPercent += 0.02;
            leftEdge = testedBin->GetBinCenter(testedBin->FindFirstBinAbove(testedBin->GetMaximum() * maxPercent));
            rightEdge = testedBin->GetBinCenter(testedBin->FindLastBinAbove(testedBin->GetMaximum() * maxPercent));
            testedBin->Fit("gaus", "QR", "SAME", leftEdge, rightEdge);
            testFunc = testedBin->GetFunction("gaus");
        }
        testingLinarIdea_3->SetPoint(counter, binning * (1 + 2 * i), testFunc->GetParameter(2));
        testingLinarIdea_3->SetPointError(counter, 0., testFunc->GetParError(2));
        Double_t tMeanPoint3 = testFunc->GetParameter(1);
        derivativeCalc(interestedBins[2], i, tMeanPoint3);
        derivativePoint3->SetPoint(counter, binning * (1 + 2 * i), derivative);
        derivativePoint3->SetPointError(counter, 0., derivativeError);
        counter++;

        TString histName("/binning_");
        Long_t num = binning * (1 + 2 * i);
        TString name = folder + newFolder + binFolder3 + histName;

        TString end = {".root"};
        TString full = name + num + end;
        TFile binFile(full, "RECREATE");
        testedBin->Write();
    }

    TString Name("/RealDataTests");
    TString suf(".root");
    TFile myfile(folder + newFolder + Name + suf, "RECREATE");
    vshape->Write("InitVShape");
    testingLinarIdea_1->Write("InterestedBin1");
    testingLinarIdea_2->Write("InterestedBin2");
    testingLinarIdea_3->Write("InterestedBin3");
    derivativePoint1->Write("DerivativeBin1");
    derivativePoint2->Write("DerivativeBin2");
    derivativePoint3->Write("DerivativeBin3");

    TCanvas *shapes = new TCanvas("Shapes", "Shapes", 1400, 1000);
    shapes->cd();
    vshape->Draw("COLZ");
    TString shape_name("/VShape.pdf");
    shapes->SaveAs(folder + newFolder + shape_name, "Q");
}

void SettingVShape::derivativeCalc(Int_t interestedBin, Int_t binningInI, Double_t tPoint) {
    Int_t leftCentralBin = interestedBin - (2 * binningInI + 1);
    Int_t rightCentralBin = interestedBin + (2 * binningInI + 1);

    Double_t leftBinX_mean = vshape->GetXaxis()->GetBinCenter(leftCentralBin);
    Double_t rightBinX_mean = vshape->GetXaxis()->GetBinCenter(rightCentralBin);
    Double_t centralBinX_mean = vshape->GetXaxis()->GetBinCenter(interestedBin);


    TH1D *leftBin = new TH1D();
    leftBin = vshape->ProjectionY("Left Bin", leftCentralBin - binningInI, leftCentralBin + binningInI);
    Double_t leftEdgeL = leftBin->GetBinCenter(leftBin->FindFirstBinAbove(leftBin->GetMaximum() / 10.));
    Double_t rightEdgeL = leftBin->GetBinCenter(leftBin->FindLastBinAbove(leftBin->GetMaximum() / 10.));
    leftBin->Fit("gaus", "QR", "SAME", leftEdgeL, rightEdgeL);
    TF1 *leftTestFunc = leftBin->GetFunction("gaus");
    Double_t leftBinT_mean = leftTestFunc->GetParameter(1);
    Double_t Y_1_err = leftTestFunc->GetParError(1);


    TH1D *rightBin = new TH1D();
    rightBin = vshape->ProjectionY("Right Bin", rightCentralBin - binningInI, rightCentralBin + binningInI);
    Double_t leftEdgeR = rightBin->GetBinCenter(rightBin->FindFirstBinAbove(rightBin->GetMaximum() / 10.));
    Double_t rightEdgeR = rightBin->GetBinCenter(rightBin->FindLastBinAbove(rightBin->GetMaximum() / 10.));
    rightBin->Fit("gaus", "QR", "SAME", leftEdgeR, rightEdgeR);
    TF1 *rightTestFunc = rightBin->GetFunction("gaus");
    Double_t rightBinT_mean = rightTestFunc->GetParameter(1);
    Double_t Y_2_err = rightTestFunc->GetParError(1);

    derivative = (1./2.) * ((tPoint - leftBinT_mean) / (centralBinX_mean - leftBinX_mean) + (rightBinT_mean - tPoint) / (rightBinX_mean - centralBinX_mean));

    Double_t DY_error = sqrt(Y_1_err*Y_1_err + Y_2_err*Y_2_err);
    Double_t DX = centralBinX_mean - leftBinX_mean;
    derivativeError = (1./2.) * sqrt((1 / DX)*(1 / DX)*DY_error*DY_error);
}

bool SettingVShape::fitIsNice(TF1 *fitFunction) {
    cout << "Bad fit detected!" << endl;
    bool isNice = true;
    if (fitFunction->GetChisquare() / fitFunction->GetNDF() > 2) isNice = false;
    if (fitFunction->GetNDF() < 11.) isNice = true;
    return isNice;
}