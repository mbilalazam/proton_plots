///////////////////////////////////////////////////////////////////////////////
// plot_dataMC_proton_detector_band.C
//
// Data vs MC comparison with detector systematic band for reco-matched protons.
// Observables: cosTheta and track length.
//
// Machinery: same as plotSPINE_systematics_MR65_band_v4_dataNorm.C
//   - MR6.1 nominal + 16 variations (8 parameters x up/down)
//   - Per-bin RMS band: sigma_i = sqrt(1/N * sum_v (var_v - nom)^2)
//   - Fractional uncertainty transferred to MR6.5 nominal
//   - MR6.5 normalized to data
//
// Formatting: same as GENIE/G4RW scripts
//   - protoDUNEStyle, orange band, fixed ymax
//
// Two output plots:
//   band_dataMC_cosTheta_detector.png
//   band_dataMC_length_detector.png
//
// Compile and run:
//   g++ -std=c++17 -O2 -g plot_dataMC_proton_detector_band.C -o plot_dataMC_proton_detector_band `root-config --cflags --libs`
//   ./plot_dataMC_proton_detector_band
///////////////////////////////////////////////////////////////////////////////

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TGraphAsymmErrors.h"

#include <vector>
#include <string>
#include <cmath>
#include <iostream>

// =============================================================
// Paths
// =============================================================
const std::string BASE_DIR =
    "/global/u1/m/mazam/CAFs_Analysis/9_May_2026/slides/proton_systematics/detector_systematics";

const std::string DETSYS_DIR =
    "/global/u1/m/mazam/CAFs_Analysis/9_May_2026/slides/proton_systematics/detector_systematics/detector_systematics";

const std::string PHYSYS_DIR =
    "/global/u1/m/mazam/CAFs_Analysis/9_May_2026/slides/proton_systematics/detector_systematics/physics_systematics";

const std::string FILE_NOM   = BASE_DIR + "/output_mr6p1.root";
const std::string FILE_MR65  =
    "/global/u1/m/mazam/CAFs_Analysis/9_May_2026/slides/proton_systematics/output.root";
const std::string FILE_DATA  =
    "/global/u1/m/mazam/sandbox/v11/sandbox_v11_minervaOff_geomContainRecoOnly_Lgt3cm_300cm_endZneg0to5cmVeto_cosZgtNeg0p9.root";

const std::string OUT_DIR =
    "/global/u1/m/mazam/CAFs_Analysis/9_May_2026/slides/proton_systematics/final_proton_codes";

// 16 variation files: 8 parameters x (down, up)
const std::vector<std::string> VAR_FILES = {
    DETSYS_DIR + "/1_MiniRun6_1E18_RHC_E_field_1sig_down_Systematics/output.root",
    DETSYS_DIR + "/2_MiniRun6_1E18_RHC_E_field_1sig_up_Systematics/output.root",
    DETSYS_DIR + "/3_MiniRun6_1E18_RHC_E_lifetime_1sig_down_Systematics/output.root",
    DETSYS_DIR + "/4_MiniRun6_1E18_RHC_E_lifetime_1sig_up_Systematics/output.root",
    DETSYS_DIR + "/5_MiniRun6_1E18_RHC_long_diff_1sig_down_Systematics/output.root",
    DETSYS_DIR + "/6_MiniRun6_1E18_RHC_long_diff_1sig_up_Systematics/output.root",
    DETSYS_DIR + "/7_MiniRun6_1E18_RHC_tran_diff_1sig_down_Systematics/output.root",
    DETSYS_DIR + "/8_MiniRun6_1E18_RHC_tran_diff_1sig_up_Systematics/output.root",
    DETSYS_DIR + "/9_MiniRun6_1E18_RHC_V_drift_1sig_down_Systematics/output.root",
    DETSYS_DIR + "/10_MiniRun6_1E18_RHC_V_drift_1sig_up_Systematics/output.root",
    PHYSYS_DIR + "/11_MiniRun6_1E18_RHC_Birks_A_1sig_down_Systematics/output.root",
    PHYSYS_DIR + "/12_MiniRun6_1E18_RHC_Birks_A_1sig_up_Systematics/output.root",
    PHYSYS_DIR + "/13_MiniRun6_1E18_RHC_Birks_k_1sig_down_Systematics/output.root",
    PHYSYS_DIR + "/14_MiniRun6_1E18_RHC_Birks_k_1sig_up_Systematics/output.root",
    PHYSYS_DIR + "/15_MiniRun6_1E18_RHC_W_ION_1sig_down_Systematics/output.root",
    PHYSYS_DIR + "/16_MiniRun6_1E18_RHC_W_ION_1sig_up_Systematics/output.root",
};

// =============================================================
// Binning
// =============================================================
const int    kNBinsC = 15;
const double kXMinC  = -1.0;
const double kXMaxC  =  1.0;

const int    kNBinsL = 15;
const double kXMinL  =  3.0;
const double kXMaxL  = 83.0;

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
    hR->SetDirectory(0);
    return hR;
}

// =============================================================
// Build detector band:
//   1. Compute fractional RMS from MR6.1 nominal + variations
//   2. Apply fractional uncertainty to MR6.5 scaled MC
// Same machinery as plotSPINE_systematics_MR65_band_v4_dataNorm.C
// =============================================================
TGraphAsymmErrors* BuildDetectorBand(TH1D* h_nom61,
                                      const std::vector<TH1D*>& h_vars,
                                      TH1D* h_mc65_scaled)
{
    int Nvars = h_vars.size();
    TGraphAsymmErrors* gr = new TGraphAsymmErrors();
    int nB = h_mc65_scaled->GetNbinsX();

    for (int i = 1; i <= nB; i++) {
        double x    = h_mc65_scaled->GetBinCenter(i);
        double y65  = h_mc65_scaled->GetBinContent(i);
        double ex   = 0.5 * h_mc65_scaled->GetBinWidth(i);
        double y61  = h_nom61->GetBinContent(i);

        // RMS of variations around MR6.1 nominal
        double sum = 0.0;
        for (int v = 0; v < Nvars; ++v) {
            double delta = h_vars[v]->GetBinContent(i) - y61;
            sum += delta * delta;
        }
        double sigma61 = (Nvars > 0) ? std::sqrt(sum / Nvars) : 0.0;

        // Transfer fractional uncertainty to MR6.5
        double frac   = (y61 > 0.0) ? sigma61 / y61 : 0.0;
        double sigma65 = frac * y65;

        gr->SetPoint(i-1, x, y65);
        gr->SetPointError(i-1, ex, ex, sigma65, sigma65);
    }
    gr->SetFillColorAlpha(kOrange+1, 0.50);
    gr->SetLineColor(kOrange+1);
    gr->SetMarkerSize(0);
    return gr;
}

// =============================================================
// Draw plot: same formatting as GENIE/G4RW scripts
// =============================================================
void DrawPlot(TH1D* hData, TH1D* hMC,
              TGraphAsymmErrors* grBand,
              const char* xTitle,
              const char* title,
              const char* bandLabel,
              const std::string& outPng,
              double ymax = -1)
{
    TCanvas* c = new TCanvas("c", "c", 800, 600);
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(0.15);

    hData->SetLineColor(kBlack);
    hData->SetLineWidth(3);
    hMC->SetLineColor(kRed);
    hMC->SetLineWidth(3);

    hData->SetTitle(title);
    hData->GetXaxis()->SetTitle(xTitle);
    hData->GetYaxis()->SetTitle("Counts");

    double maxY = (ymax > 0) ? ymax : 1.4 * std::max(hData->GetMaximum(), hMC->GetMaximum());
    hData->GetYaxis()->SetRangeUser(0, maxY);

    hData->Draw("HIST E");
    grBand->Draw("2 SAME");
    hMC->Draw("HIST SAME");

	TLegend* leg = new TLegend(0.50, 0.75, 0.80, 0.93);
    leg->SetTextFont(133);
    leg->SetTextSize(18);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hData,   "Data Reco",               "l");
    leg->AddEntry(hMC,     "MC Reco (norm. to Data)", "l");
    leg->AddEntry(grBand,  bandLabel,                 "f");
    leg->Draw();

    TLatex tL; tL.SetNDC();
    tL.DrawLatex(0.20, 0.94, "#bf{DUNE:ND-LAr 2x2}");

    c->SaveAs(outPng.c_str());
    delete c;
}

// =============================================================
// Main
// =============================================================
void plot_dataMC_proton_detector_band()
{
    gROOT->LoadMacro("protoDUNEStyle.C");
    gROOT->SetStyle("protoDUNEStyle");
    gROOT->ForceStyle();
    gStyle->SetTitleX(0.25);
    gSystem->mkdir(OUT_DIR.c_str(), true);

    // ── Open files ────────────────────────────────────────────────────────
    TFile* fNom  = TFile::Open(FILE_NOM.c_str());
    TFile* fMR65 = TFile::Open(FILE_MR65.c_str());
    TFile* fData = TFile::Open(FILE_DATA.c_str());

    if (!fNom  || fNom->IsZombie())  { std::cerr << "ERROR: MR6.1 nominal file\n"; return; }
    if (!fMR65 || fMR65->IsZombie()) { std::cerr << "ERROR: MR6.5 file\n";         return; }
    if (!fData || fData->IsZombie()) { std::cerr << "ERROR: Data file\n";           return; }

    // ── Load MR6.1 nominal histograms ─────────────────────────────────────
    TH1D* h_nom61_cos_raw = (TH1D*)fNom->Get("reco_matched_proton_cosL");
    TH1D* h_nom61_len_raw = (TH1D*)fNom->Get("reco_matched_proton_length");
    if (!h_nom61_cos_raw || !h_nom61_len_raw) {
        std::cerr << "ERROR: Missing MR6.1 nominal histograms\n"; return;
    }
    TH1D* h_nom61_cos = RebinPreserveErrors(h_nom61_cos_raw, kNBinsC, kXMinC, kXMaxC, "h_nom61_cos");
    TH1D* h_nom61_len = RebinPreserveErrors(h_nom61_len_raw, kNBinsL, kXMinL, kXMaxL, "h_nom61_len");

    // ── Load MR6.5 histograms ─────────────────────────────────────────────
    TH1D* h_mr65_cos_raw = (TH1D*)fMR65->Get("reco_matched_proton_cosL");
    TH1D* h_mr65_len_raw = (TH1D*)fMR65->Get("reco_matched_proton_length");
    if (!h_mr65_cos_raw || !h_mr65_len_raw) {
        std::cerr << "ERROR: Missing MR6.5 histograms\n"; return;
    }
    TH1D* h_mr65_cos = RebinPreserveErrors(h_mr65_cos_raw, kNBinsC, kXMinC, kXMaxC, "h_mr65_cos");
    TH1D* h_mr65_len = RebinPreserveErrors(h_mr65_len_raw, kNBinsL, kXMinL, kXMaxL, "h_mr65_len");

    // ── Load data histograms ──────────────────────────────────────────────
    TH1D* hDataCos_raw = (TH1D*)fData->Get("recoProtonCosLAll");
    TH1D* hDataLen_raw = (TH1D*)fData->Get("recoProtonLengthAll");
    if (!hDataCos_raw || !hDataLen_raw) {
        std::cerr << "ERROR: Missing data histograms\n"; return;
    }
    TH1D* hDataCos = RebinPreserveErrors(hDataCos_raw, kNBinsC, kXMinC, kXMaxC, "hDataCos");
    TH1D* hDataLen = RebinPreserveErrors(hDataLen_raw, kNBinsL, kXMinL, kXMaxL, "hDataLen");

    // ── Load 16 variation histograms ──────────────────────────────────────
    std::vector<TH1D*> h_vars_cos, h_vars_len;
    int Nvars = (int)VAR_FILES.size();

    for (int v = 0; v < Nvars; v++) {
        TFile* fVar = TFile::Open(VAR_FILES[v].c_str());
        if (!fVar || fVar->IsZombie()) {
            std::cerr << "ERROR: Cannot open variation file: " << VAR_FILES[v] << "\n";
            return;
        }
        TH1D* hc_raw = (TH1D*)fVar->Get("reco_matched_proton_cosL");
        TH1D* hl_raw = (TH1D*)fVar->Get("reco_matched_proton_length");
        if (!hc_raw || !hl_raw) {
            std::cerr << "ERROR: Missing histograms in: " << VAR_FILES[v] << "\n";
            return;
        }
        // Clone before closing file
        TH1D* hc = RebinPreserveErrors(hc_raw, kNBinsC, kXMinC, kXMaxC,
                                        Form("h_var%d_cos", v));
        TH1D* hl = RebinPreserveErrors(hl_raw, kNBinsL, kXMinL, kXMaxL,
                                        Form("h_var%d_len", v));
        h_vars_cos.push_back(hc);
        h_vars_len.push_back(hl);
        fVar->Close();
        std::cout << "Loaded variation " << v << ": " << VAR_FILES[v] << "\n";
    }

    // ── Scale MR6.5 to data ───────────────────────────────────────────────
    double scaleCos = (h_mr65_cos->Integral("width") > 0.0)
                    ? hDataCos->Integral("width") / h_mr65_cos->Integral("width") : 1.0;
    TH1D* hMC_cos = (TH1D*)h_mr65_cos->Clone("hMC_cos");
    hMC_cos->Scale(scaleCos);

    double scaleLen = (h_mr65_len->Integral("width") > 0.0)
                    ? hDataLen->Integral("width") / h_mr65_len->Integral("width") : 1.0;
    TH1D* hMC_len = (TH1D*)h_mr65_len->Clone("hMC_len");
    hMC_len->Scale(scaleLen);

    // ── Normalize variations to MR6.1 nominal (POT equalisation) ─────────
    double norm_cos = h_nom61_cos->Integral();
    double norm_len = h_nom61_len->Integral();
    for (int v = 0; v < Nvars; v++) {
        double vi_cos = h_vars_cos[v]->Integral();
        double vi_len = h_vars_len[v]->Integral();
        if (vi_cos > 0.0) h_vars_cos[v]->Scale(norm_cos / vi_cos);
        if (vi_len > 0.0) h_vars_len[v]->Scale(norm_len / vi_len);
    }

    // ── Build detector bands ──────────────────────────────────────────────
    auto* gr_cos = BuildDetectorBand(h_nom61_cos, h_vars_cos, hMC_cos);
    auto* gr_len = BuildDetectorBand(h_nom61_len, h_vars_len, hMC_len);

    // ── Draw plots ────────────────────────────────────────────────────────
    DrawPlot(hDataCos, hMC_cos, gr_cos,
             "cos(#theta)", "Protons",
             "Detector Systematics",
             OUT_DIR + "/band_dataMC_cosTheta_detector.png", 50);

    DrawPlot(hDataLen, hMC_len, gr_len,
             "Track Length [cm]", "Protons",
             "Detector Systematics",
             OUT_DIR + "/band_dataMC_length_detector.png", 70);

    std::cout << "\nDone. Plots saved in:\n" << OUT_DIR << "\n";
    fNom->Close();
    fMR65->Close();
    fData->Close();
}

int main()
{
    plot_dataMC_proton_detector_band();
    return 0;
}