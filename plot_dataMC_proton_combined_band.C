///////////////////////////////////////////////////////////////////////////////
// plot_dataMC_proton_combined_band.C
//
// Data vs MC comparison with three systematic bands:
//   - GENIE (sum-of-dials, symmetric)
//   - FSI / G4RW patched (symmetric, same treatment as GENIE)
//   - Total = GENIE ⊕ FSI added in quadrature
//
// Two output plots:
//   band_dataMC_cosTheta_combined.png
//   band_dataMC_length_combined.png
//
// Compile and run:
//   g++ -std=c++17 -O2 -g plot_dataMC_proton_combined_band.C -o plot_dataMC_proton_combined_band `root-config --cflags --libs`
//   ./plot_dataMC_proton_combined_band
///////////////////////////////////////////////////////////////////////////////

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TGraphAsymmErrors.h"

#include <map>
#include <tuple>
#include <array>
#include <iostream>
#include <cmath>
#include <string>
#include <vector>

// =============================================================
// Paths
// =============================================================
const std::string CAF_FILE =
    "/global/u1/m/mazam/CAFs_Analysis/9_May_2026/slides/proton_systematics/output.root";

const std::string DATA_FILE =
    "/global/u1/m/mazam/sandbox/v11/sandbox_v11_minervaOff_geomContainRecoOnly_Lgt3cm_300cm_endZneg0to5cmVeto_cosZgtNeg0p9.root";

const std::string SYST_FILE =
    "/global/u1/m/mazam/Geant4ReWeighting/2x2_G4RW/test/on_CAFS/step2_reco/nusyst_step2_reco_output.root";

const std::string WFILE_PATCHED =
    "/global/cfs/cdirs/dune/users/mazam/G4RW/output_weights/patched_withRock_genieNu/weights_1000_all_LAr_P_cascade_0_4000_5weights.root";

const std::string OUT_DIR =
    "/global/u1/m/mazam/CAFs_Analysis/9_May_2026/slides/proton_systematics/final_proton_codes";

const std::string BASE_DET =
    "/global/u1/m/mazam/CAFs_Analysis/9_May_2026/slides/proton_systematics/detector_systematics";
const std::string DETSYS_DIR = BASE_DET + "/detector_systematics";
const std::string PHYSYS_DIR = BASE_DET + "/physics_systematics";
const std::string FILE_NOM61 = BASE_DET + "/output_mr6p1.root";

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
// Weight map
// =============================================================
using WKey = std::tuple<int,int,int>;
using W2   = std::array<double,2>;   // [0]=w_1p33, [1]=w_0p67

double SanitizeW(double w) {
    // return (w < 0.01 || w > 100.0) ? 1.0 : w;
    return w;
}

void LoadWeights(const std::string& file, std::map<WKey,W2>& wmap)
{
    TFile* f = TFile::Open(file.c_str());
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: Cannot open: " << file << "\n"; return;
    }
    TTree* t = (TTree*)f->Get("g4rwEvt");
    if (!t) {
        std::cerr << "ERROR: No g4rwEvt tree in: " << file << "\n";
        f->Close(); return;
    }
    int file_idx, event_id, track_id;
    double w_1p33, w_0p67, w_1p66, w_0p34, w_2p00;
    t->SetBranchAddress("file_idx",    &file_idx);
    t->SetBranchAddress("event_id",    &event_id);
    t->SetBranchAddress("track_id",    &track_id);
    t->SetBranchAddress("weight_1p33", &w_1p33);
    t->SetBranchAddress("weight_0p67", &w_0p67);
    t->SetBranchAddress("weight_1p66", &w_1p66);
    t->SetBranchAddress("weight_0p34", &w_0p34);
    t->SetBranchAddress("weight_2p00", &w_2p00);
    long N = t->GetEntries();
    for (long i = 0; i < N; i++) {
        t->GetEntry(i);
        wmap[{file_idx, event_id, track_id}] = {SanitizeW(w_1p33), SanitizeW(w_0p67)};
    }
    f->Close();
    std::cout << "Loaded " << N << " entries from: " << file << "\n";
}

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
// Build GENIE band (symmetric) — same as GENIE script
// =============================================================
TGraphAsymmErrors* BuildGENIEBand(TGraphAsymmErrors* grBand, TH1D* hMC)
{
    TGraphAsymmErrors* gr = new TGraphAsymmErrors();
    for (int i = 0; i < grBand->GetN(); ++i) {
        double x, yNom;
        grBand->GetPoint(i, x, yNom);
        int b = hMC->FindBin(x);
        double yMC  = hMC->GetBinContent(b);
        double sigma = grBand->GetErrorYhigh(i);
        double frac  = (yNom > 0.0) ? sigma / yNom : 0.0;
        double exl   = grBand->GetErrorXlow(i);
        double exh   = grBand->GetErrorXhigh(i);
        gr->SetPoint(i, x, yMC);
        gr->SetPointError(i, exl, exh, frac*yMC, frac*yMC);
    }
    gr->SetFillColorAlpha(kOrange+1, 0.50);
    gr->SetLineColor(kOrange+1);
    return gr;
}

// =============================================================
// Build FSI band (symmetric) — same treatment as GENIE
// sigma = 0.5*(eup + edn), applied symmetrically
// =============================================================
TGraphAsymmErrors* BuildFSIBand(TH1D* h_nom, TH1D* h_up, TH1D* h_dn,
                                 TH1D* hMC)
{
    TGraphAsymmErrors* gr = new TGraphAsymmErrors();
    int nB = hMC->GetNbinsX();
    for (int i = 1; i <= nB; i++) {
        double x    = hMC->GetBinCenter(i);
        double yMC  = hMC->GetBinContent(i);
        double ex   = 0.5 * hMC->GetBinWidth(i);
        double yNom = h_nom->GetBinContent(i);
        double eup  = std::max(0.0, h_up->GetBinContent(i) - yNom);
        double edn  = std::max(0.0, yNom - h_dn->GetBinContent(i));
        double sigma = 0.5 * (eup + edn);
        double frac  = (yNom > 0.0) ? sigma / yNom : 0.0;
        gr->SetPoint(i-1, x, yMC);
        gr->SetPointError(i-1, ex, ex, frac*yMC, frac*yMC);
    }
    gr->SetFillColorAlpha(kBlue+1, 0.35);
    gr->SetLineColor(kBlue+1);
    return gr;
}

// =============================================================
// Build detector band
// =============================================================
TGraphAsymmErrors* BuildDetectorBand(TH1D* h_nom61, const std::vector<TH1D*>& h_vars,
                                      TH1D* hMC)
{
    int Nvars = h_vars.size();
    TGraphAsymmErrors* gr = new TGraphAsymmErrors();
    int nB = hMC->GetNbinsX();
    for (int i = 1; i <= nB; i++) {
        double x   = hMC->GetBinCenter(i);
        double y65 = hMC->GetBinContent(i);
        double ex  = 0.5 * hMC->GetBinWidth(i);
        double y61 = h_nom61->GetBinContent(i);
        double sum = 0.0;
        for (int v = 0; v < Nvars; ++v) {
            double delta = h_vars[v]->GetBinContent(i) - y61;
            sum += delta * delta;
        }
        double sigma61  = (Nvars > 0) ? std::sqrt(sum / Nvars) : 0.0;
        double frac     = (y61 > 0.0) ? sigma61 / y61 : 0.0;
        gr->SetPoint(i-1, x, y65);
        gr->SetPointError(i-1, ex, ex, frac*y65, frac*y65);
    }
    gr->SetFillColorAlpha(kGreen+2, 0.35);
    gr->SetLineColor(kGreen+2);
    return gr;
}

// =============================================================
// Build Total band: GENIE ⊕ FSI ⊕ Detector in quadrature 
// =============================================================
TGraphAsymmErrors* BuildTotalBand(TGraphAsymmErrors* grGENIE,
                                   TGraphAsymmErrors* grFSI,
                                   TGraphAsymmErrors* grDet)
{
    TGraphAsymmErrors* gr = new TGraphAsymmErrors();
    int nB = std::min({grGENIE->GetN(), grFSI->GetN(), grDet->GetN()});
    for (int i = 0; i < nB; i++) {
    // int nB = grGENIE->GetN();
    // for (int i = 0; i < nB; i++) {
        double x, y;
        grGENIE->GetPoint(i, x, y);
        double sig_genie = grGENIE->GetErrorYhigh(i);
        double sig_fsi   = grFSI->GetErrorYhigh(i);
        double sig_det   = grDet->GetErrorYhigh(i);
        double sig_tot   = std::sqrt(sig_genie*sig_genie + sig_fsi*sig_fsi + sig_det*sig_det);
        double exl = grGENIE->GetErrorXlow(i);
        double exh = grGENIE->GetErrorXhigh(i);
        gr->SetPoint(i, x, y);
        gr->SetPointError(i, exl, exh, sig_tot, sig_tot);
    }
    gr->SetFillColorAlpha(kOrange+1, 0.50);
    gr->SetLineColor(kOrange+1);
    return gr;
}


// =============================================================
// Draw plot: Data, MC, Total band
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

    hData->SetLineColor(kBlack);  hData->SetLineWidth(3);
    hMC->SetLineColor(kRed);      hMC->SetLineWidth(3);

    hData->SetTitle(title);
    hData->GetXaxis()->SetTitle(xTitle);
    hData->GetYaxis()->SetTitle("Counts");

    double maxY = (ymax > 0) ? ymax : 1.4 * std::max(hData->GetMaximum(), hMC->GetMaximum());
    hData->GetYaxis()->SetRangeUser(0, maxY);

    hData->Draw("HIST E");
    grBand->Draw("2 SAME");
    hMC->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.45, 0.75, 0.82, 0.93);
    leg->SetTextFont(133); leg->SetTextSize(18);
    leg->SetFillStyle(0);  leg->SetBorderSize(0);
    leg->AddEntry(hData,   "Data Reco",               "l");
    leg->AddEntry(hMC,     "MC Reco (norm. to Data)", "l");
    // leg->AddEntry(grBand, "Total (GENIE #oplus FSI)",    "f");
    leg->AddEntry(grBand, "Total (GENIE #oplus FSI #oplus Det)",    "f");
    // leg->AddEntry(grBand,  bandLabel,                 "f");
    leg->Draw();

    TLatex tL; tL.SetNDC();
    tL.DrawLatex(0.20, 0.94, "#bf{DUNE:ND-LAr 2x2}");

    c->SaveAs(outPng.c_str());
    delete c;
}

// =============================================================
// Draw combined plot: Data, MC, GENIE band, FSI band, Total band
// =============================================================
void DrawCombinedPlot(TH1D* hData, TH1D* hMC,
                      TGraphAsymmErrors* grGENIE,
                      TGraphAsymmErrors* grFSI,
                      TGraphAsymmErrors* grTotal,
                      const char* xTitle,
                      const char* title,
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

    // Draw order: total (widest) first, then GENIE, then FSI, then MC, then data on top
    hData->Draw("HIST E");
    grTotal->Draw("2 SAME");
    grGENIE->Draw("2 SAME");
    grFSI->Draw("2 SAME");
    hMC->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.45, 0.68, 0.82, 0.93);
    leg->SetTextFont(133);
    leg->SetTextSize(18);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hData,   "Data Reco",                   "l");
    leg->AddEntry(hMC,     "MC Reco (norm. to Data)",     "l");
    leg->AddEntry(grGENIE, "GENIE Systematics",           "f");
    leg->AddEntry(grFSI,   "FSI Systematics",             "f");
    // leg->AddEntry(grTotal, "Total (GENIE #oplus FSI)",    "f");
    leg->AddEntry(grTotal, "Total (GENIE #oplus FSI #oplus Det)",    "f");
    leg->Draw();

    TLatex tL; tL.SetNDC();
    tL.DrawLatex(0.20, 0.94, "#bf{DUNE:ND-LAr 2x2}");

    c->SaveAs(outPng.c_str());
    delete c;
}

// =============================================================
// Main
// =============================================================
void plot_dataMC_proton_combined_band()
{
    gROOT->LoadMacro("protoDUNEStyle.C");
    gROOT->SetStyle("protoDUNEStyle");
    gROOT->ForceStyle();
    gStyle->SetTitleX(0.25);
    gSystem->mkdir(OUT_DIR.c_str(), true);

    // ── Load weight map ───────────────────────────────────────────────────
    std::map<WKey,W2> wmap_pat;
    LoadWeights(WFILE_PATCHED, wmap_pat);

    // ── Open files ────────────────────────────────────────────────────────
    TFile* fCAF  = TFile::Open(CAF_FILE.c_str());
    TFile* fData = TFile::Open(DATA_FILE.c_str());
    TFile* fSyst = TFile::Open(SYST_FILE.c_str());

    if (!fCAF || fCAF->IsZombie())  { std::cerr << "ERROR: CAF file\n";  return; }
    if (!fData || fData->IsZombie()) { std::cerr << "ERROR: Data file\n"; return; }
    if (!fSyst || fSyst->IsZombie()) { std::cerr << "ERROR: Syst file\n"; return; }

    TTree* tReco = (TTree*)fCAF->Get("syst_g4rw_proton_reco");
    if (!tReco) { std::cerr << "ERROR: Cannot find syst_g4rw_proton_reco\n"; return; }

    // ── Data histograms ───────────────────────────────────────────────────
    TH1D* hDataCos = RebinPreserveErrors((TH1D*)fData->Get("recoProtonCosLAll"),
                                          kNBinsC, kXMinC, kXMaxC, "hDataCos");
    TH1D* hDataLen = RebinPreserveErrors((TH1D*)fData->Get("recoProtonLengthAll"),
                                          kNBinsL, kXMinL, kXMaxL, "hDataLen");

    // ── MC histograms from CAF (for scaling and display) ──────────────────
    TH1D* hMC_cos_raw = RebinPreserveErrors((TH1D*)fCAF->Get("reco_matched_proton_cosL"),
                                             kNBinsC, kXMinC, kXMaxC, "hMC_cos_raw");
    TH1D* hMC_len_raw = RebinPreserveErrors((TH1D*)fCAF->Get("reco_matched_proton_length"),
                                             kNBinsL, kXMinL, kXMaxL, "hMC_len_raw");

    // ── Nominal + weighted histograms from tree (for FSI band) ────────────
    TH1D* h_nom_cos    = new TH1D("h_nom_cos",    "", kNBinsC, kXMinC, kXMaxC);
    TH1D* h_pat_cos_up = new TH1D("h_pat_cos_up", "", kNBinsC, kXMinC, kXMaxC);
    TH1D* h_pat_cos_dn = new TH1D("h_pat_cos_dn", "", kNBinsC, kXMinC, kXMaxC);
    TH1D* h_nom_len    = new TH1D("h_nom_len",    "", kNBinsL, kXMinL, kXMaxL);
    TH1D* h_pat_len_up = new TH1D("h_pat_len_up", "", kNBinsL, kXMinL, kXMaxL);
    TH1D* h_pat_len_dn = new TH1D("h_pat_len_dn", "", kNBinsL, kXMinL, kXMaxL);

    int cafFileIdx, G4ID;
    long long nu_id;
    double cosT, length;
    tReco->SetBranchAddress("cafFileIdx", &cafFileIdx);
    tReco->SetBranchAddress("nu_id",      &nu_id);
    tReco->SetBranchAddress("G4ID",       &G4ID);
    tReco->SetBranchAddress("cosTheta",   &cosT);
    tReco->SetBranchAddress("length",     &length);

    long N = tReco->GetEntries();
    std::cout << "Reco entries: " << N << "\n";
    for (long i = 0; i < N; i++) {
        tReco->GetEntry(i);
        WKey key = {cafFileIdx, (int)nu_id, G4ID};
        W2 wp = {1.0, 1.0}; if (wmap_pat.count(key)) wp = wmap_pat[key];

        h_nom_cos->Fill(cosT);
        h_pat_cos_up->Fill(cosT, wp[0]); h_pat_cos_dn->Fill(cosT, wp[1]);

        h_nom_len->Fill(length);
        h_pat_len_up->Fill(length, wp[0]); h_pat_len_dn->Fill(length, wp[1]);
    }

    // ── Scale MC to data ──────────────────────────────────────────────────
    double scaleCos = (hMC_cos_raw->Integral("width") > 0.0)
                    ? hDataCos->Integral("width") / hMC_cos_raw->Integral("width") : 1.0;
    TH1D* hMC_cos = (TH1D*)hMC_cos_raw->Clone("hMC_cos");
    hMC_cos->Scale(scaleCos);

    double scaleLen = (hMC_len_raw->Integral("width") > 0.0)
                    ? hDataLen->Integral("width") / hMC_len_raw->Integral("width") : 1.0;
    TH1D* hMC_len = (TH1D*)hMC_len_raw->Clone("hMC_len");
    hMC_len->Scale(scaleLen);

    // ── GENIE bands ───────────────────────────────────────────────────────
    TGraphAsymmErrors* grGENIE_cos_raw =
        (TGraphAsymmErrors*)fSyst->Get("gr_band_sumDials_pCosL");
    if (!grGENIE_cos_raw) { std::cerr << "ERROR: GENIE cosL band not found\n"; return; }
    grGENIE_cos_raw = (TGraphAsymmErrors*)grGENIE_cos_raw->Clone("grGENIE_cos_raw");

    // Note: update the length band object name below if it differs in your syst file
    TGraphAsymmErrors* grGENIE_len_raw =
        (TGraphAsymmErrors*)fSyst->Get("gr_band_sumDials_pLength");
    if (!grGENIE_len_raw) { std::cerr << "ERROR: GENIE length band not found — check object name in syst file\n"; return; }
    grGENIE_len_raw = (TGraphAsymmErrors*)grGENIE_len_raw->Clone("grGENIE_len_raw");

    auto* grGENIE_cos = BuildGENIEBand(grGENIE_cos_raw, hMC_cos);
    auto* grGENIE_len = BuildGENIEBand(grGENIE_len_raw, hMC_len);

    // ── FSI bands ─────────────────────────────────────────────────────────
    auto* grFSI_cos = BuildFSIBand(h_nom_cos, h_pat_cos_up, h_pat_cos_dn, hMC_cos);
    auto* grFSI_len = BuildFSIBand(h_nom_len, h_pat_len_up, h_pat_len_dn, hMC_len);

    // ── Detector bands ────────────────────────────────────────────────────
    TFile* fNom61 = TFile::Open(FILE_NOM61.c_str());
    if (!fNom61 || fNom61->IsZombie()) { std::cerr << "ERROR: MR6.1 file\n"; return; }

    TH1D* h_nom61_cos = RebinPreserveErrors((TH1D*)fNom61->Get("reco_matched_proton_cosL"),
                                             kNBinsC, kXMinC, kXMaxC, "h_nom61_cos"); h_nom61_cos->SetDirectory(0);
    TH1D* h_nom61_len = RebinPreserveErrors((TH1D*)fNom61->Get("reco_matched_proton_length"),
                                             kNBinsL, kXMinL, kXMaxL, "h_nom61_len"); h_nom61_len->SetDirectory(0);

    std::vector<TH1D*> h_vars_cos, h_vars_len;
    int Nvars = (int)VAR_FILES.size();
    for (int v = 0; v < Nvars; v++) {
        TFile* fVar = TFile::Open(VAR_FILES[v].c_str());
        if (!fVar || fVar->IsZombie()) { std::cerr << "ERROR: var file " << v << "\n"; return; }
        TH1D* hc = RebinPreserveErrors((TH1D*)fVar->Get("reco_matched_proton_cosL"),
                                        kNBinsC, kXMinC, kXMaxC, Form("h_var%d_cos", v));
        TH1D* hl = RebinPreserveErrors((TH1D*)fVar->Get("reco_matched_proton_length"),
                                        kNBinsL, kXMinL, kXMaxL, Form("h_var%d_len", v));
        hc->SetDirectory(0); hl->SetDirectory(0);
        h_vars_cos.push_back(hc);
        h_vars_len.push_back(hl);
        fVar->Close();
    }
    // Normalize variations to MR6.1 nominal (POT equalisation)
    double norm_cos = h_nom61_cos->Integral();
    double norm_len = h_nom61_len->Integral();
    for (int v = 0; v < Nvars; v++) {
        if (h_vars_cos[v]->Integral() > 0.0) h_vars_cos[v]->Scale(norm_cos / h_vars_cos[v]->Integral());
        if (h_vars_len[v]->Integral() > 0.0) h_vars_len[v]->Scale(norm_len / h_vars_len[v]->Integral());
    }

    auto* grDet_cos = BuildDetectorBand(h_nom61_cos, h_vars_cos, hMC_cos);
    auto* grDet_len = BuildDetectorBand(h_nom61_len, h_vars_len, hMC_len);

    // ── Total bands ───────────────────────────────────────────────────────
    auto* grTotal_cos = BuildTotalBand(grGENIE_cos, grFSI_cos, grDet_cos);
    auto* grTotal_len = BuildTotalBand(grGENIE_len, grFSI_len, grDet_len);

    // ── Draw ──────────────────────────────────────────────────────────────
    // DrawCombinedPlot(hDataCos, hMC_cos, grGENIE_cos, grFSI_cos, grTotal_cos,
    //                  "cos(#theta)", "Protons",
    //                  OUT_DIR + "/band_dataMC_cosTheta_combined.png", 50);

    // DrawCombinedPlot(hDataLen, hMC_len, grGENIE_len, grFSI_len, grTotal_len,
    //                  "Track Length [cm]", "Protons",
    //                  OUT_DIR + "/band_dataMC_length_combined.png", 70);

    DrawPlot(hDataCos, hMC_cos, grTotal_cos,
             "cos(#theta)", "Protons", "Total Systematics",
             OUT_DIR + "/band_dataMC_cosTheta_total.png", 50);

    DrawPlot(hDataLen, hMC_len, grTotal_len,
             "Track Length [cm]", "Protons", "Total Systematics",
             OUT_DIR + "/band_dataMC_length_total.png", 70);

    std::cout << "\nDone. Plots saved in:\n" << OUT_DIR << "\n";

    fCAF->Close();
    fNom61->Close();
    fData->Close();
    fSyst->Close();
}

int main()
{
    plot_dataMC_proton_combined_band();
    return 0;
}