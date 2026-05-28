///////////////////////////////////////////////////////////////////////////////
// plot_MC_proton_cosTheta_efficiency.C
//
// MC reco proton angular distribution (cos theta) for Case A and Case D.
// Two versions:
//   1. Raw counts
//   2. Shape normalized
// Both with a fractional pad showing Case D / Case A ratio (efficiency).
//
// Compile:
//   g++ -std=c++17 -O2 -g plot_MC_proton_cosTheta_efficiency.C -o plot_MC_proton_cosTheta_efficiency `root-config --cflags --libs`
// Run:
//   ./plot_MC_proton_cosTheta_efficiency
///////////////////////////////////////////////////////////////////////////////

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TLine.h"
#include <iostream>
#include <cmath>
#include <string>

// =============================================================
// Binning
// =============================================================
const int    kNBins = 15;
const double kXMin  = -1.0;
const double kXMax  =  1.0;

// =============================================================
// Paths
// =============================================================
const std::string FILE_A =
    "/global/u1/m/mazam/CAFs_Analysis/9_May_2026/slides/proton_systematics/final_proton_codes/scripts_eff_purity/mc_caseA.root";

const std::string FILE_D =
    "/global/u1/m/mazam/CAFs_Analysis/9_May_2026/slides/proton_systematics/final_proton_codes/scripts_eff_purity/mc_caseD.root";

const std::string OUT_DIR =
    "/global/u1/m/mazam/CAFs_Analysis/9_May_2026/slides/proton_systematics/final_proton_codes/scripts_eff_purity";

// =============================================================
// Rebin preserving errors
// =============================================================
TH1D* RebinPreserveErrors(TH1D* h, int nbins, double xmin, double xmax,
                           const char* name)
{
    TH1D* hR = new TH1D(name, "", nbins, xmin, xmax);
    hR->Sumw2();
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        double x = h->GetBinCenter(i);
        if (x < xmin || x >= xmax) continue;
        int b = hR->FindBin(x);
        hR->SetBinContent(b, hR->GetBinContent(b) + h->GetBinContent(i));
        hR->SetBinError(b, std::sqrt(hR->GetBinError(b)*hR->GetBinError(b)
                                   + h->GetBinError(i)*h->GetBinError(i)));
    }
    return hR;
}

// =============================================================
// Draw split-pad plot
// =============================================================
void DrawSplitPad(TH1D* hA, TH1D* hD, TH1D* hRatio,
                  const char* yTitle,
                  const std::string& outPng,
                  double ymax)
{
    TCanvas* c = new TCanvas("c", "c", 800, 750);

    // Upper pad
    TPad* p1 = new TPad("p1", "", 0, 0.3, 1, 1.0);
    p1->SetLeftMargin(0.18);
    p1->SetRightMargin(0.15);
    p1->SetBottomMargin(0.02);
    p1->SetTopMargin(0.08);
    p1->Draw();
    p1->cd();

    hA->GetXaxis()->SetLabelSize(0);
    hA->GetXaxis()->SetTitleSize(0);
    hA->GetYaxis()->SetTitle(yTitle);
    hA->GetYaxis()->SetTitleOffset(1.3);
    hA->GetYaxis()->SetRangeUser(0, ymax);
    hA->Draw("HIST");
    hD->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.20, 0.68, 0.55, 0.88);
    leg->SetTextFont(133);
    leg->SetTextSize(18);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hA, "Nominal-MINERvA", "l");
    leg->AddEntry(hD, "Final Selection",  "l");
    leg->Draw();

    TLatex tL; tL.SetNDC();
    tL.DrawLatex(0.20, 0.94, "#bf{DUNE:ND-LAr 2x2}");

    // Lower pad — ratio = efficiency as function of cos theta
    c->cd();
    TPad* p2 = new TPad("p2", "", 0, 0.0, 1, 0.3);
    p2->SetLeftMargin(0.18);
    p2->SetRightMargin(0.15);
    p2->SetTopMargin(0.02);
    p2->SetBottomMargin(0.35);
    p2->Draw();
    p2->cd();

    hRatio->SetLineColor(kBlack);
    hRatio->SetLineWidth(2);
    hRatio->SetTitle("");
    hRatio->GetXaxis()->SetTitle("cos(#theta)");
    hRatio->GetYaxis()->SetTitle("Efficiency");
    hRatio->GetYaxis()->SetRangeUser(0, 1.2);
    hRatio->GetYaxis()->SetNdivisions(505);
    hRatio->GetXaxis()->SetLabelSize(0.12);
    hRatio->GetXaxis()->SetTitleSize(0.14);
    hRatio->GetXaxis()->SetTitleOffset(1.0);
    hRatio->GetYaxis()->SetLabelSize(0.10);
    hRatio->GetYaxis()->SetTitleSize(0.12);
    hRatio->GetYaxis()->SetTitleOffset(0.70);
    hRatio->Draw("HIST E");

    // Reference line at 1
    TLine* line = new TLine(kXMin, 1.0, kXMax, 1.0);
    line->SetLineColor(kRed);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw();

    c->SaveAs(outPng.c_str());
    delete c;
}

// =============================================================
// Main
// =============================================================
void plot_MC_proton_cosTheta_efficiency()
{
    gROOT->LoadMacro("protoDUNEStyle.C");
    gROOT->SetStyle("protoDUNEStyle");
    gROOT->ForceStyle();
    gStyle->SetTitleX(0.25);

    // ── Open files ────────────────────────────────────────────
    TFile* fA = TFile::Open(FILE_A.c_str());
    TFile* fD = TFile::Open(FILE_D.c_str());
    if (!fA || fA->IsZombie() || !fD || fD->IsZombie()) {
        std::cerr << "ERROR: Cannot open input files\n"; return;
    }

    TH1D* hA_raw = (TH1D*)fA->Get("reco_matched_proton_cosL");
    TH1D* hD_raw = (TH1D*)fD->Get("reco_matched_proton_cosL");
    if (!hA_raw || !hD_raw) {
        std::cerr << "ERROR: Missing histograms\n"; return;
    }

    // ── Rebin ─────────────────────────────────────────────────
    TH1D* hA = RebinPreserveErrors(hA_raw, kNBins, kXMin, kXMax, "hA");
    TH1D* hD = RebinPreserveErrors(hD_raw, kNBins, kXMin, kXMax, "hD");

    // ── Style ─────────────────────────────────────────────────
    hA->SetLineColor(kBlack);
    hA->SetLineWidth(3);
    hA->SetTitle("Reco Proton Candidates");
    hA->GetXaxis()->SetTitle("cos(#theta)");

    hD->SetLineColor(kRed);
    hD->SetLineWidth(3);

    // ── Ratio (efficiency as function of cos theta) ───────────
    // Always use raw numbers for ratio
    TH1D* hRatio = (TH1D*)hD->Clone("hRatio");
    hRatio->Divide(hA);

    // ── Version 1: Raw counts ─────────────────────────────────
    double ymaxRaw = 1.4 * hA->GetMaximum();

    DrawSplitPad(hA, hD, hRatio,
                 "Counts",
                 OUT_DIR + "/MC_proton_cosTheta_efficiency_raw.png",
                 ymaxRaw);

    // ── Version 2: Shape normalized ───────────────────────────
    TH1D* hA_norm = (TH1D*)hA->Clone("hA_norm");
    TH1D* hD_norm = (TH1D*)hD->Clone("hD_norm");

    if (hA_norm->Integral() > 0) hA_norm->Scale(1.0 / hA_norm->Integral());
    if (hD_norm->Integral() > 0) hD_norm->Scale(1.0 / hD_norm->Integral());

    hA_norm->SetLineColor(kBlack);
    hA_norm->SetLineWidth(3);
    hA_norm->SetTitle("Reco Proton Candidates");
    hA_norm->GetXaxis()->SetTitle("cos(#theta)");

    hD_norm->SetLineColor(kRed);
    hD_norm->SetLineWidth(3);

    double ymaxNorm = 1.4 * std::max(hA_norm->GetMaximum(), hD_norm->GetMaximum());

    DrawSplitPad(hA_norm, hD_norm, hRatio,
                 "Normalized Counts",
                 OUT_DIR + "/MC_proton_cosTheta_efficiency_norm.png",
                 ymaxNorm);

    std::cout << "\nDone. Plots saved in:\n" << OUT_DIR << "\n";
    fA->Close();
    fD->Close();
}

int main()
{
    plot_MC_proton_cosTheta_efficiency();
    return 0;
}
