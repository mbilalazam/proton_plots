// root -l -q plot_proton_stacked.C

#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TROOT.h"
#include "THStack.h"
#include <iostream>

// ------------------------------------------------------------
// Helper: Rebin preserving errors
// ------------------------------------------------------------
TH1D* RebinPreserveErrors(TH1D* h, int nbins, double xmin, double xmax, const char* name)
{
    TH1D* hR = new TH1D(name, "", nbins, xmin, xmax);
    hR->Sumw2();
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        double x = h->GetBinCenter(i);
        if (x < xmin || x >= xmax) continue;
        int b = hR->FindBin(x);
        hR->SetBinContent(b, hR->GetBinContent(b) + h->GetBinContent(i));
        hR->SetBinError(b, std::sqrt(hR->GetBinError(b)*hR->GetBinError(b) + h->GetBinError(i)*h->GetBinError(i)));
    }
    return hR;
}

// ------------------------------------------------------------
// Draw one canvas: absolute or normalized, for a given variable
// ------------------------------------------------------------
void DrawStack(
    THStack* hs,
    TH1D* hTotal,       // sum of all components (for axis range)
    TLegend* leg,
    const char* xTitle,
    const char* yTitle,
    const char* label,  // e.g. "proton_cosL_intType"
    bool norm,
    double yMax)
{
    TCanvas* c = new TCanvas(label, label, 800, 600);
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(0.15);

    hs->Draw("HIST");
    hs->GetXaxis()->SetTitle(xTitle);
    hs->GetYaxis()->SetTitle(yTitle);
    hs->GetXaxis()->SetTitleOffset(1.2);
    hs->GetYaxis()->SetTitleOffset(1.2);
    hs->SetTitle("Proton Candidates");
    hs->SetMinimum(0);
    hs->SetMaximum(yMax);
    gPad->Update();

    leg->Draw();

    TLatex tL;
    tL.SetNDC();
    tL.DrawLatex(0.20, 0.94, "#bf{DUNE:ND-LAr 2x2}");
    tL.DrawLatex(0.20, 0.84, "#bf{Work in Progress}");

    TString fname = TString(label) + (norm ? "_norm.png" : "_abs.png");
    c->SaveAs(fname);
}

// ------------------------------------------------------------
// MAIN
// ------------------------------------------------------------
void plot_proton_stacked()
{
    gROOT->LoadMacro("protoDUNEStyle.C");
    gROOT->SetStyle("protoDUNEStyle");
    gROOT->ForceStyle();
    gStyle->SetTitleX(0.25);

    // --------------------------------------------------------
    // Binning
    // --------------------------------------------------------
    const int    kNCosL   = 15;
    const double kCosLMin = -1.0, kCosLMax = 1.0;

    const int    kNLen    = 16;
    const double kLenMin  =  3.0, kLenMax  = 83.0;

    // --------------------------------------------------------
    // Open MC file
    // --------------------------------------------------------
    TFile fMC("/global/u1/m/mazam/CAFs_Analysis/9_May_2026/slides/proton_systematics/final_proton_codes/output.root");
    if (fMC.IsZombie()) { std::cerr << "Cannot open MC file\n"; return; }

    // --------------------------------------------------------
    // Load and rebin: interaction type — cosL
    // --------------------------------------------------------
    TH1D* hCosL_QE_raw  = (TH1D*)fMC.Get("reco_matched_proton_cosL_QE");
    TH1D* hCosL_MEC_raw = (TH1D*)fMC.Get("reco_matched_proton_cosL_MEC");
    TH1D* hCosL_RES_raw = (TH1D*)fMC.Get("reco_matched_proton_cosL_RES");
    TH1D* hCosL_DIS_raw = (TH1D*)fMC.Get("reco_matched_proton_cosL_DIS");
    TH1D* hCosL_COH_raw = (TH1D*)fMC.Get("reco_matched_proton_cosL_COH");

    if (!hCosL_QE_raw || !hCosL_MEC_raw || !hCosL_RES_raw || !hCosL_DIS_raw || !hCosL_COH_raw) {
        std::cerr << "Missing interaction-type cosL histograms\n"; return;
    }

    TH1D* hCosL_QE  = RebinPreserveErrors(hCosL_QE_raw,  kNCosL, kCosLMin, kCosLMax, "hCosL_QE");
    TH1D* hCosL_MEC = RebinPreserveErrors(hCosL_MEC_raw, kNCosL, kCosLMin, kCosLMax, "hCosL_MEC");
    TH1D* hCosL_RES = RebinPreserveErrors(hCosL_RES_raw, kNCosL, kCosLMin, kCosLMax, "hCosL_RES");
    TH1D* hCosL_DIS = RebinPreserveErrors(hCosL_DIS_raw, kNCosL, kCosLMin, kCosLMax, "hCosL_DIS");
    TH1D* hCosL_COH = RebinPreserveErrors(hCosL_COH_raw, kNCosL, kCosLMin, kCosLMax, "hCosL_COH");

    // --------------------------------------------------------
    // Load and rebin: interaction type — length
    // --------------------------------------------------------
    TH1D* hLen_QE_raw  = (TH1D*)fMC.Get("reco_matched_proton_length_QE");
    TH1D* hLen_MEC_raw = (TH1D*)fMC.Get("reco_matched_proton_length_MEC");
    TH1D* hLen_RES_raw = (TH1D*)fMC.Get("reco_matched_proton_length_RES");
    TH1D* hLen_DIS_raw = (TH1D*)fMC.Get("reco_matched_proton_length_DIS");
    TH1D* hLen_COH_raw = (TH1D*)fMC.Get("reco_matched_proton_length_COH");

    if (!hLen_QE_raw || !hLen_MEC_raw || !hLen_RES_raw || !hLen_DIS_raw || !hLen_COH_raw) {
        std::cerr << "Missing interaction-type length histograms\n"; return;
    }

    TH1D* hLen_QE  = RebinPreserveErrors(hLen_QE_raw,  kNLen, kLenMin, kLenMax, "hLen_QE");
    TH1D* hLen_MEC = RebinPreserveErrors(hLen_MEC_raw, kNLen, kLenMin, kLenMax, "hLen_MEC");
    TH1D* hLen_RES = RebinPreserveErrors(hLen_RES_raw, kNLen, kLenMin, kLenMax, "hLen_RES");
    TH1D* hLen_DIS = RebinPreserveErrors(hLen_DIS_raw, kNLen, kLenMin, kLenMax, "hLen_DIS");
    TH1D* hLen_COH = RebinPreserveErrors(hLen_COH_raw, kNLen, kLenMin, kLenMax, "hLen_COH");

    // --------------------------------------------------------
    // Load and rebin: particle type — cosL
    // --------------------------------------------------------
    TH1D* hCosL_trueMuon_raw   = (TH1D*)fMC.Get("reco_matched_proton_cosL_trueMuon");
    TH1D* hCosL_trueProton_raw = (TH1D*)fMC.Get("reco_matched_proton_cosL_trueProton");
    TH1D* hCosL_truePion_raw   = (TH1D*)fMC.Get("reco_matched_proton_cosL_truePion");
    TH1D* hCosL_trueOther_raw  = (TH1D*)fMC.Get("reco_matched_proton_cosL_trueOther");

    if (!hCosL_trueMuon_raw || !hCosL_trueProton_raw || !hCosL_truePion_raw || !hCosL_trueOther_raw) {
        std::cerr << "Missing particle-type cosL histograms\n"; return;
    }

    TH1D* hCosL_trueMuon   = RebinPreserveErrors(hCosL_trueMuon_raw,   kNCosL, kCosLMin, kCosLMax, "hCosL_trueMuon");
    TH1D* hCosL_trueProton = RebinPreserveErrors(hCosL_trueProton_raw, kNCosL, kCosLMin, kCosLMax, "hCosL_trueProton");
    TH1D* hCosL_truePion   = RebinPreserveErrors(hCosL_truePion_raw,   kNCosL, kCosLMin, kCosLMax, "hCosL_truePion");
    TH1D* hCosL_trueOther  = RebinPreserveErrors(hCosL_trueOther_raw,  kNCosL, kCosLMin, kCosLMax, "hCosL_trueOther");

    // --------------------------------------------------------
    // Load and rebin: particle type — length
    // --------------------------------------------------------
    TH1D* hLen_trueMuon_raw   = (TH1D*)fMC.Get("reco_matched_proton_length_trueMuon");
    TH1D* hLen_trueProton_raw = (TH1D*)fMC.Get("reco_matched_proton_length_trueProton");
    TH1D* hLen_truePion_raw   = (TH1D*)fMC.Get("reco_matched_proton_length_truePion");
    TH1D* hLen_trueOther_raw  = (TH1D*)fMC.Get("reco_matched_proton_length_trueOther");

    if (!hLen_trueMuon_raw || !hLen_trueProton_raw || !hLen_truePion_raw || !hLen_trueOther_raw) {
        std::cerr << "Missing particle-type length histograms\n"; return;
    }

    TH1D* hLen_trueMuon   = RebinPreserveErrors(hLen_trueMuon_raw,   kNLen, kLenMin, kLenMax, "hLen_trueMuon");
    TH1D* hLen_trueProton = RebinPreserveErrors(hLen_trueProton_raw, kNLen, kLenMin, kLenMax, "hLen_trueProton");
    TH1D* hLen_truePion   = RebinPreserveErrors(hLen_truePion_raw,   kNLen, kLenMin, kLenMax, "hLen_truePion");
    TH1D* hLen_trueOther  = RebinPreserveErrors(hLen_trueOther_raw,  kNLen, kLenMin, kLenMax, "hLen_trueOther");

    // --------------------------------------------------------
    // Colors: interaction type
    // --------------------------------------------------------
    hCosL_QE->SetFillColor(kRed);       hCosL_QE->SetLineColor(kRed);       hCosL_QE->SetLineWidth(0);
    hCosL_MEC->SetFillColor(kBlue);     hCosL_MEC->SetLineColor(kBlue);     hCosL_MEC->SetLineWidth(0);
    hCosL_DIS->SetFillColor(kYellow);   hCosL_DIS->SetLineColor(kYellow);   hCosL_DIS->SetLineWidth(0);
    hCosL_RES->SetFillColor(kGreen+2);  hCosL_RES->SetLineColor(kGreen+2);  hCosL_RES->SetLineWidth(0);
    hCosL_COH->SetFillColor(kCyan);     hCosL_COH->SetLineColor(kCyan);     hCosL_COH->SetLineWidth(0);

    hLen_QE->SetFillColor(kRed);        hLen_QE->SetLineColor(kRed);        hLen_QE->SetLineWidth(0);
    hLen_MEC->SetFillColor(kBlue);      hLen_MEC->SetLineColor(kBlue);      hLen_MEC->SetLineWidth(0);
    hLen_DIS->SetFillColor(kYellow);    hLen_DIS->SetLineColor(kYellow);    hLen_DIS->SetLineWidth(0);
    hLen_RES->SetFillColor(kGreen+2);   hLen_RES->SetLineColor(kGreen+2);   hLen_RES->SetLineWidth(0);
    hLen_COH->SetFillColor(kCyan);      hLen_COH->SetLineColor(kCyan);      hLen_COH->SetLineWidth(0);

    // --------------------------------------------------------
    // Colors: particle type (grayscale)
    // --------------------------------------------------------
    hCosL_trueProton->SetFillColor(kGray+2);  hCosL_trueProton->SetLineColor(kGray+2);  hCosL_trueProton->SetLineWidth(0);
    hCosL_truePion->SetFillColor(kGray+1);    hCosL_truePion->SetLineColor(kGray+1);    hCosL_truePion->SetLineWidth(0);
    hCosL_trueMuon->SetFillColor(kGray);      hCosL_trueMuon->SetLineColor(kGray);      hCosL_trueMuon->SetLineWidth(0);
    hCosL_trueOther->SetFillColor(kWhite);    hCosL_trueOther->SetLineColor(kBlack);    hCosL_trueOther->SetLineWidth(1);

    hLen_trueProton->SetFillColor(kGray+2);   hLen_trueProton->SetLineColor(kGray+2);   hLen_trueProton->SetLineWidth(0);
    hLen_truePion->SetFillColor(kGray+1);     hLen_truePion->SetLineColor(kGray+1);     hLen_truePion->SetLineWidth(0);
    hLen_trueMuon->SetFillColor(kGray);       hLen_trueMuon->SetLineColor(kGray);       hLen_trueMuon->SetLineWidth(0);
    hLen_trueOther->SetFillColor(kWhite);     hLen_trueOther->SetLineColor(kBlack);     hLen_trueOther->SetLineWidth(1);


    // --------------------------------------------------------
    // Colors: particle type (hatched) — test versions
    // --------------------------------------------------------
    TH1D* hCosL_trueProton_h = (TH1D*)hCosL_trueProton->Clone("hCosL_trueProton_h");
    TH1D* hCosL_truePion_h   = (TH1D*)hCosL_truePion->Clone("hCosL_truePion_h");
    TH1D* hCosL_trueMuon_h   = (TH1D*)hCosL_trueMuon->Clone("hCosL_trueMuon_h");
    TH1D* hCosL_trueOther_h  = (TH1D*)hCosL_trueOther->Clone("hCosL_trueOther_h");

    TH1D* hLen_trueProton_h  = (TH1D*)hLen_trueProton->Clone("hLen_trueProton_h");
    TH1D* hLen_truePion_h    = (TH1D*)hLen_truePion->Clone("hLen_truePion_h");
    TH1D* hLen_trueMuon_h    = (TH1D*)hLen_trueMuon->Clone("hLen_trueMuon_h");
    TH1D* hLen_trueOther_h   = (TH1D*)hLen_trueOther->Clone("hLen_trueOther_h");

    hCosL_trueProton_h->SetFillColor(kBlue+1);    hCosL_trueProton_h->SetLineColor(kBlue+1);    hCosL_trueProton_h->SetFillStyle(3004); hCosL_trueProton_h->SetLineWidth(2);
    hCosL_truePion_h->SetFillColor(kRed+1);       hCosL_truePion_h->SetLineColor(kRed+1);       hCosL_truePion_h->SetFillStyle(3005);   hCosL_truePion_h->SetLineWidth(2);
    hCosL_trueMuon_h->SetFillColor(kGreen+2);     hCosL_trueMuon_h->SetLineColor(kGreen+2);     hCosL_trueMuon_h->SetFillStyle(3006);   hCosL_trueMuon_h->SetLineWidth(2);
    hCosL_trueOther_h->SetFillColor(kMagenta+1);  hCosL_trueOther_h->SetLineColor(kMagenta+1);  hCosL_trueOther_h->SetFillStyle(3444);  hCosL_trueOther_h->SetLineWidth(2);

    hLen_trueProton_h->SetFillColor(kBlue+1);     hLen_trueProton_h->SetLineColor(kBlue+1);     hLen_trueProton_h->SetFillStyle(3004);  hLen_trueProton_h->SetLineWidth(2);
    hLen_truePion_h->SetFillColor(kRed+1);        hLen_truePion_h->SetLineColor(kRed+1);        hLen_truePion_h->SetFillStyle(3005);    hLen_truePion_h->SetLineWidth(2);
    hLen_trueMuon_h->SetFillColor(kGreen+2);      hLen_trueMuon_h->SetLineColor(kGreen+2);      hLen_trueMuon_h->SetFillStyle(3006);    hLen_trueMuon_h->SetLineWidth(2);
    hLen_trueOther_h->SetFillColor(kMagenta+1);   hLen_trueOther_h->SetLineColor(kMagenta+1);   hLen_trueOther_h->SetFillStyle(3444);   hLen_trueOther_h->SetLineWidth(2);

    THStack* hsCosL_part_h = new THStack("hsCosL_part_h", "");
    hsCosL_part_h->Add(hCosL_trueOther_h);
    hsCosL_part_h->Add(hCosL_trueMuon_h);
    hsCosL_part_h->Add(hCosL_truePion_h);
    hsCosL_part_h->Add(hCosL_trueProton_h);

    THStack* hsLen_part_h = new THStack("hsLen_part_h", "");
    hsLen_part_h->Add(hLen_trueOther_h);
    hsLen_part_h->Add(hLen_trueMuon_h);
    hsLen_part_h->Add(hLen_truePion_h);
    hsLen_part_h->Add(hLen_trueProton_h);

    // --------------------------------------------------------
    // Print counts
    // --------------------------------------------------------
    std::cout << "\n=== Interaction type — cosL ===" << std::endl;
    std::cout << "QE:  " << hCosL_QE->Integral()  << std::endl;
    std::cout << "MEC: " << hCosL_MEC->Integral() << std::endl;
    std::cout << "RES: " << hCosL_RES->Integral() << std::endl;
    std::cout << "DIS: " << hCosL_DIS->Integral() << std::endl;
    std::cout << "COH: " << hCosL_COH->Integral() << std::endl;
    
    std::cout << "\n=== Particle type — cosL ===" << std::endl;
    std::cout << "True Proton: " << hCosL_trueProton->Integral() << std::endl;
    std::cout << "True Pion:   " << hCosL_truePion->Integral()   << std::endl;
    std::cout << "True Muon:   " << hCosL_trueMuon->Integral()   << std::endl;
    std::cout << "True Other:  " << hCosL_trueOther->Integral()  << std::endl;
        
    // --------------------------------------------------------
    // Build THStacks
    // --------------------------------------------------------
    // Interaction type stacks
    THStack* hsCosL_int = new THStack("hsCosL_int", "");
    hsCosL_int->Add(hCosL_COH);
    hsCosL_int->Add(hCosL_QE);
    hsCosL_int->Add(hCosL_MEC);
    hsCosL_int->Add(hCosL_RES);
    hsCosL_int->Add(hCosL_DIS);

    THStack* hsLen_int = new THStack("hsLen_int", "");
    hsLen_int->Add(hLen_COH);
    hsLen_int->Add(hLen_QE);
    hsLen_int->Add(hLen_MEC);
    hsLen_int->Add(hLen_RES);
    hsLen_int->Add(hLen_DIS);

    // Particle type stacks
    THStack* hsCosL_part = new THStack("hsCosL_part", "");
    hsCosL_part->Add(hCosL_trueOther);
    hsCosL_part->Add(hCosL_trueMuon);
    hsCosL_part->Add(hCosL_truePion);
    hsCosL_part->Add(hCosL_trueProton);

    THStack* hsLen_part = new THStack("hsLen_part", "");
    hsLen_part->Add(hLen_trueOther);
    hsLen_part->Add(hLen_trueMuon);
    hsLen_part->Add(hLen_truePion);
    hsLen_part->Add(hLen_trueProton);

    // --------------------------------------------------------
    // Total histograms (for yMax)
    // --------------------------------------------------------
    TH1D* hCosL_total = (TH1D*)hCosL_QE->Clone("hCosL_total");
    hCosL_total->Add(hCosL_MEC);
    hCosL_total->Add(hCosL_RES);
    hCosL_total->Add(hCosL_DIS);
    hCosL_total->Add(hCosL_COH);

    TH1D* hLen_total = (TH1D*)hLen_QE->Clone("hLen_total");
    hLen_total->Add(hLen_MEC);
    hLen_total->Add(hLen_RES);
    hLen_total->Add(hLen_DIS);
    hLen_total->Add(hLen_COH);

    double yMaxCosL = 1.3 * hCosL_total->GetMaximum();
    double yMaxLen  = 1.3 * hLen_total->GetMaximum();

    // ============================================================
    // ABSOLUTE plots
    // ============================================================
    {
        TLegend* leg = new TLegend(0.54, 0.58, 0.82, 0.90);
        leg->SetTextFont(133); leg->SetTextSize(18); leg->SetFillStyle(0); leg->SetBorderSize(0);
        leg->SetHeader("NuMI RHC MC");
        leg->AddEntry(hCosL_DIS, "DIS", "f"); leg->AddEntry(hCosL_RES, "RES", "f");
        leg->AddEntry(hCosL_MEC, "MEC", "f"); leg->AddEntry(hCosL_QE,  "QE",  "f");
        leg->AddEntry(hCosL_COH, "COH", "f");
        DrawStack(hsCosL_int, hCosL_total, leg, "Reconstructed cos#kern[-0.9]{(}#theta#kern[0.0]{)}", "Counts", "proton_cosL_intType", false, yMaxCosL);
    }
    {
        TLegend* leg = new TLegend(0.54, 0.58, 0.82, 0.90);
        leg->SetTextFont(133); leg->SetTextSize(18); leg->SetFillStyle(0); leg->SetBorderSize(0);
        leg->SetHeader("NuMI RHC MC");
        leg->AddEntry(hLen_DIS, "DIS", "f"); leg->AddEntry(hLen_RES, "RES", "f");
        leg->AddEntry(hLen_MEC, "MEC", "f"); leg->AddEntry(hLen_QE,  "QE",  "f");
        leg->AddEntry(hLen_COH, "COH", "f");
        DrawStack(hsLen_int, hLen_total, leg, "Reconstructed Track Length [cm]", "Counts", "proton_length_intType", false, yMaxLen);
    }
    {
        TLegend* leg = new TLegend(0.54, 0.58, 0.82, 0.90);
        leg->SetTextFont(133); leg->SetTextSize(18); leg->SetFillStyle(0); leg->SetBorderSize(0);
        leg->SetHeader("NuMI RHC MC");
        leg->AddEntry(hCosL_trueProton, "True Proton", "f"); leg->AddEntry(hCosL_truePion,  "True Pion",  "f");
        leg->AddEntry(hCosL_trueMuon,   "True Muon",   "f"); leg->AddEntry(hCosL_trueOther, "True Other", "f");
        DrawStack(hsCosL_part, hCosL_total, leg, "Reconstructed cos#kern[-0.9]{(}#theta#kern[0.0]{)}", "Counts", "proton_cosL_partType", false, yMaxCosL);
    }
    {
        TLegend* leg = new TLegend(0.54, 0.58, 0.82, 0.90);
        leg->SetTextFont(133); leg->SetTextSize(18); leg->SetFillStyle(0); leg->SetBorderSize(0);
        leg->SetHeader("NuMI RHC MC");
        leg->AddEntry(hLen_trueProton, "True Proton", "f"); leg->AddEntry(hLen_truePion,  "True Pion",  "f");
        leg->AddEntry(hLen_trueMuon,   "True Muon",   "f"); leg->AddEntry(hLen_trueOther, "True Other", "f");
        DrawStack(hsLen_part, hLen_total, leg, "Reconstructed Track Length [cm]", "Counts", "proton_length_partType", false, yMaxLen);
    }

    {
        TLegend* leg = new TLegend(0.54, 0.58, 0.82, 0.90);
        leg->SetTextFont(133); leg->SetTextSize(18); leg->SetFillStyle(0); leg->SetBorderSize(0);
        leg->SetHeader("NuMI RHC MC");
        leg->AddEntry(hCosL_trueProton_h, "True Proton", "f"); leg->AddEntry(hCosL_truePion_h,  "True Pion",  "f");
        leg->AddEntry(hCosL_trueMuon_h,   "True Muon",   "f"); leg->AddEntry(hCosL_trueOther_h, "True Other", "f");
        DrawStack(hsCosL_part_h, hCosL_total, leg, "Reconstructed cos#kern[-0.9]{(}#theta#kern[0.0]{)}", "Counts", "proton_cosL_partType_hatch", false, yMaxCosL);
    }
    {
        TLegend* leg = new TLegend(0.54, 0.58, 0.82, 0.90);
        leg->SetTextFont(133); leg->SetTextSize(18); leg->SetFillStyle(0); leg->SetBorderSize(0);
        leg->SetHeader("NuMI RHC MC");
        leg->AddEntry(hLen_trueProton_h, "True Proton", "f"); leg->AddEntry(hLen_truePion_h,  "True Pion",  "f");
        leg->AddEntry(hLen_trueMuon_h,   "True Muon",   "f"); leg->AddEntry(hLen_trueOther_h, "True Other", "f");
        DrawStack(hsLen_part_h, hLen_total, leg, "Reconstructed Track Length [cm]", "Counts", "proton_length_partType_hatch", false, yMaxLen);
    }

    // ============================================================
    // Build normalized clones
    // ============================================================
    double intTotalCosL = hCosL_total->Integral();
    double intTotalLen  = hLen_total->Integral();

    TH1D* hCosL_QE_n  = (TH1D*)hCosL_QE->Clone("hCosL_QE_n");   hCosL_QE_n->Scale(1.0/intTotalCosL);
    TH1D* hCosL_MEC_n = (TH1D*)hCosL_MEC->Clone("hCosL_MEC_n"); hCosL_MEC_n->Scale(1.0/intTotalCosL);
    TH1D* hCosL_RES_n = (TH1D*)hCosL_RES->Clone("hCosL_RES_n"); hCosL_RES_n->Scale(1.0/intTotalCosL);
    TH1D* hCosL_DIS_n = (TH1D*)hCosL_DIS->Clone("hCosL_DIS_n"); hCosL_DIS_n->Scale(1.0/intTotalCosL);
    TH1D* hCosL_COH_n = (TH1D*)hCosL_COH->Clone("hCosL_COH_n"); hCosL_COH_n->Scale(1.0/intTotalCosL);

    TH1D* hLen_QE_n  = (TH1D*)hLen_QE->Clone("hLen_QE_n");   hLen_QE_n->Scale(1.0/intTotalLen);
    TH1D* hLen_MEC_n = (TH1D*)hLen_MEC->Clone("hLen_MEC_n"); hLen_MEC_n->Scale(1.0/intTotalLen);
    TH1D* hLen_RES_n = (TH1D*)hLen_RES->Clone("hLen_RES_n"); hLen_RES_n->Scale(1.0/intTotalLen);
    TH1D* hLen_DIS_n = (TH1D*)hLen_DIS->Clone("hLen_DIS_n"); hLen_DIS_n->Scale(1.0/intTotalLen);
    TH1D* hLen_COH_n = (TH1D*)hLen_COH->Clone("hLen_COH_n"); hLen_COH_n->Scale(1.0/intTotalLen);

    TH1D* hCosL_trueMuon_n   = (TH1D*)hCosL_trueMuon->Clone("hCosL_trueMuon_n");     hCosL_trueMuon_n->Scale(1.0/intTotalCosL);
    TH1D* hCosL_trueProton_n = (TH1D*)hCosL_trueProton->Clone("hCosL_trueProton_n"); hCosL_trueProton_n->Scale(1.0/intTotalCosL);
    TH1D* hCosL_truePion_n   = (TH1D*)hCosL_truePion->Clone("hCosL_truePion_n");     hCosL_truePion_n->Scale(1.0/intTotalCosL);
    TH1D* hCosL_trueOther_n  = (TH1D*)hCosL_trueOther->Clone("hCosL_trueOther_n");   hCosL_trueOther_n->Scale(1.0/intTotalCosL);

    TH1D* hLen_trueMuon_n   = (TH1D*)hLen_trueMuon->Clone("hLen_trueMuon_n");     hLen_trueMuon_n->Scale(1.0/intTotalLen);
    TH1D* hLen_trueProton_n = (TH1D*)hLen_trueProton->Clone("hLen_trueProton_n"); hLen_trueProton_n->Scale(1.0/intTotalLen);
    TH1D* hLen_truePion_n   = (TH1D*)hLen_truePion->Clone("hLen_truePion_n");     hLen_truePion_n->Scale(1.0/intTotalLen);
    TH1D* hLen_trueOther_n  = (TH1D*)hLen_trueOther->Clone("hLen_trueOther_n");   hLen_trueOther_n->Scale(1.0/intTotalLen);

    THStack* hsCosL_int_n = new THStack("hsCosL_int_n", "");
    hsCosL_int_n->Add(hCosL_COH_n); hsCosL_int_n->Add(hCosL_QE_n);
    hsCosL_int_n->Add(hCosL_MEC_n); hsCosL_int_n->Add(hCosL_RES_n); hsCosL_int_n->Add(hCosL_DIS_n);

    THStack* hsLen_int_n = new THStack("hsLen_int_n", "");
    hsLen_int_n->Add(hLen_COH_n); hsLen_int_n->Add(hLen_QE_n);
    hsLen_int_n->Add(hLen_MEC_n); hsLen_int_n->Add(hLen_RES_n); hsLen_int_n->Add(hLen_DIS_n);

    THStack* hsCosL_part_n = new THStack("hsCosL_part_n", "");
    hsCosL_part_n->Add(hCosL_trueOther_n); hsCosL_part_n->Add(hCosL_trueMuon_n);
    hsCosL_part_n->Add(hCosL_truePion_n);  hsCosL_part_n->Add(hCosL_trueProton_n);

    THStack* hsLen_part_n = new THStack("hsLen_part_n", "");
    hsLen_part_n->Add(hLen_trueOther_n); hsLen_part_n->Add(hLen_trueMuon_n);
    hsLen_part_n->Add(hLen_truePion_n);  hsLen_part_n->Add(hLen_trueProton_n);

    TH1D* hCosL_total_n = (TH1D*)hCosL_total->Clone("hCosL_total_n"); hCosL_total_n->Scale(1.0/intTotalCosL);
    TH1D* hLen_total_n  = (TH1D*)hLen_total->Clone("hLen_total_n");   hLen_total_n->Scale(1.0/intTotalLen);

    double yMaxCosL_n = yMaxCosL / intTotalCosL;
    double yMaxLen_n  = yMaxLen  / intTotalLen;
    
    // ============================================================
    // NORMALIZED plots
    // ============================================================
    {
        TLegend* leg = new TLegend(0.54, 0.58, 0.82, 0.90);
        leg->SetTextFont(133); leg->SetTextSize(18); leg->SetFillStyle(0); leg->SetBorderSize(0);
        leg->SetHeader("#splitline{NuMI RHC}{MC is normalized to unity}");
        leg->AddEntry(hCosL_DIS_n, "DIS", "f"); leg->AddEntry(hCosL_RES_n, "RES", "f");
        leg->AddEntry(hCosL_MEC_n, "MEC", "f"); leg->AddEntry(hCosL_QE_n,  "QE",  "f");
        leg->AddEntry(hCosL_COH_n, "COH", "f");
        DrawStack(hsCosL_int_n, hCosL_total_n, leg, "Reconstructed cos#kern[-0.9]{(}#theta#kern[0.0]{)}", "Normalized Counts", "proton_cosL_intType", true, yMaxCosL_n);
    }
    {
        TLegend* leg = new TLegend(0.54, 0.58, 0.82, 0.90);
        leg->SetTextFont(133); leg->SetTextSize(18); leg->SetFillStyle(0); leg->SetBorderSize(0);
        leg->SetHeader("#splitline{NuMI RHC}{MC is normalized to unity}");
        leg->AddEntry(hLen_DIS_n, "DIS", "f"); leg->AddEntry(hLen_RES_n, "RES", "f");
        leg->AddEntry(hLen_MEC_n, "MEC", "f"); leg->AddEntry(hLen_QE_n,  "QE",  "f");
        leg->AddEntry(hLen_COH_n, "COH", "f");
        DrawStack(hsLen_int_n, hLen_total_n, leg, "Reconstructed Track Length [cm]", "Normalized Counts", "proton_length_intType", true, yMaxLen_n);
    }
    {
        TLegend* leg = new TLegend(0.54, 0.58, 0.82, 0.90);
        leg->SetTextFont(133); leg->SetTextSize(18); leg->SetFillStyle(0); leg->SetBorderSize(0);
        leg->SetHeader("#splitline{NuMI RHC}{MC is normalized to unity}");
        leg->AddEntry(hCosL_trueProton_n, "True Proton", "f"); leg->AddEntry(hCosL_truePion_n,  "True Pion",  "f");
        leg->AddEntry(hCosL_trueMuon_n,   "True Muon",   "f"); leg->AddEntry(hCosL_trueOther_n, "True Other", "f");
        DrawStack(hsCosL_part_n, hCosL_total_n, leg, "Reconstructed cos#kern[-0.9]{(}#theta#kern[0.0]{)}", "Normalized Counts", "proton_cosL_partType", true, yMaxCosL_n);
    }
    {
        TLegend* leg = new TLegend(0.54, 0.58, 0.82, 0.90);
        leg->SetTextFont(133); leg->SetTextSize(18); leg->SetFillStyle(0); leg->SetBorderSize(0);
        leg->SetHeader("#splitline{NuMI RHC }{MC is normalized to unity}");
        leg->AddEntry(hLen_trueProton_n, "True Proton", "f"); leg->AddEntry(hLen_truePion_n,  "True Pion",  "f");
        leg->AddEntry(hLen_trueMuon_n,   "True Muon",   "f"); leg->AddEntry(hLen_trueOther_n, "True Other", "f");
        DrawStack(hsLen_part_n, hLen_total_n, leg, "Reconstructed Track Length [cm]", "Normalized Counts", "proton_length_partType", true, yMaxLen_n);
    }

    std::cout << "\nDone. Saved 8 plots." << std::endl;
}
