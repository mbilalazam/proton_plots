///////////////////////////////////////////////////////////////////////////////
// plot_dataMC_proton_g4rw_band.C
//
// Data vs MC comparison with G4RW ±1σ (w=1.33 / w=0.67) unisim band
// for reco-matched protons (cosTheta and track length).
//
// Modelled on plot_dataMC_proton_cosTheta_GENIE.C — band is built as a
// fractional uncertainty from the G4RW weight files, then applied to the
// scaled MC histogram, identical to the GENIE band treatment.
//
// Four output plots:
//   band_reco_cosTheta_patched.png
//   band_reco_cosTheta_unpatched.png
//   band_reco_length_patched.png
//   band_reco_length_unpatched.png
//
// Compile and run:
//   g++ -std=c++17 -O2 -g plot_dataMC_proton_g4rw_band.C -o plot_dataMC_proton_g4rw_band `root-config --cflags --libs`
//   ./plot_dataMC_proton_g4rw_band
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

// =============================================================
// Paths
// =============================================================
const std::string CAF_FILE =
    // "/global/u1/m/mazam/MiniRuns/MR6.5/proton_final/mc_minirun6p5_minervaOff_geomContainRecoOnly_Lgt3cm_300cm_endZneg0to5cmVeto_cosZgtNeg0p9.root";
    "/global/u1/m/mazam/CAFs_Analysis/9_May_2026/slides/proton_systematics/output.root";

const std::string DATA_FILE =
    "/global/u1/m/mazam/sandbox/v11/sandbox_v11_minervaOff_geomContainRecoOnly_Lgt3cm_300cm_endZneg0to5cmVeto_cosZgtNeg0p9.root";

const std::string WFILE_PATCHED =
    "/global/cfs/cdirs/dune/users/mazam/G4RW/output_weights/patched_withRock_genieNu/weights_1000_all_LAr_P_cascade_0_4000_5weights.root";

const std::string WFILE_UNPATCHED =
    "/global/cfs/cdirs/dune/users/mazam/G4RW/output_weights/unpatched_withRock_genieNu/weights_1000_all_LAr_P_cascade_0_4000_5weights.root";

const std::string OUT_DIR =
    "/global/u1/m/mazam/CAFs_Analysis/9_May_2026/slides/proton_systematics/final_proton_codes";

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
// Rebin preserving errors (same as GENIE script)
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
// Build G4RW band as TGraphAsymmErrors with fractional uncertainty
// applied to scaled MC — identical treatment to GENIE band
//
// h_up and h_dn are filled with w=1.33 and w=0.67 weights (unscaled).
// We convert to fractional uncertainty bin-by-bin using h_nom_unscaled,
// then apply to h_mc_scaled.
// =============================================================
TGraphAsymmErrors* BuildG4RWBand(TH1D* h_nom_unscaled,
                                  TH1D* h_up_unscaled,
                                  TH1D* h_dn_unscaled,
                                  TH1D* h_mc_scaled)
{
    TGraphAsymmErrors* gr = new TGraphAsymmErrors();
    int nB = h_mc_scaled->GetNbinsX();
    for (int i = 1; i <= nB; i++) {
        double x    = h_mc_scaled->GetBinCenter(i);
        double yMC  = h_mc_scaled->GetBinContent(i);
        double ex   = 0.5 * h_mc_scaled->GetBinWidth(i);
        double yNom = h_nom_unscaled->GetBinContent(i);

        double eup = std::max(0.0, h_up_unscaled->GetBinContent(i) - yNom);
        double edn = std::max(0.0, yNom - h_dn_unscaled->GetBinContent(i));

        // convert to fractional, apply to scaled MC
        double frac_up = (yNom > 0.0) ? eup / yNom : 0.0;
        double frac_dn = (yNom > 0.0) ? edn / yNom : 0.0;

        gr->SetPoint(i-1, x, yMC);
        gr->SetPointError(i-1, ex, ex, frac_dn * yMC, frac_up * yMC);
    }
    gr->SetFillColorAlpha(kOrange+1, 0.50);
    gr->SetLineColor(kOrange+1);
    return gr;
}

// =============================================================
// Draw one data-MC plot with G4RW band
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
    leg->AddEntry(hData,   "Data Reco",                 "l");
    leg->AddEntry(hMC,     "MC Reco (norm. to Data)",   "l");
    leg->AddEntry(grBand,  bandLabel,                   "f");
    leg->Draw();

    TLatex tL; tL.SetNDC();
    tL.DrawLatex(0.20, 0.94, "#bf{DUNE:ND-LAr 2x2}");

    c->SaveAs(outPng.c_str());
    delete c;
}

// =============================================================
// Main
// =============================================================
void plot_dataMC_proton_g4rw_band()
{
    gROOT->LoadMacro("protoDUNEStyle.C");
    gROOT->SetStyle("protoDUNEStyle");
    gROOT->ForceStyle();
    gStyle->SetTitleX(0.25);
    gSystem->mkdir(OUT_DIR.c_str(), true);

    // ── Load weight maps ─────────────────────────────────────────────────
    std::map<WKey,W2> wmap_pat, wmap_unp;
    LoadWeights(WFILE_PATCHED,   wmap_pat);
    LoadWeights(WFILE_UNPATCHED, wmap_unp);

    // ── Open CAF file and get reco tree ──────────────────────────────────
    TFile* fCAF = TFile::Open(CAF_FILE.c_str());
    if (!fCAF || fCAF->IsZombie()) {
        std::cerr << "ERROR: Cannot open CAF file\n"; return;
    }
    TTree* tReco = (TTree*)fCAF->Get("syst_g4rw_proton_reco");
    if (!tReco) {
        std::cerr << "ERROR: Cannot find syst_g4rw_proton_reco tree\n"; return;
    }

    // ── Load data histograms ─────────────────────────────────────────────
    TFile* fData = TFile::Open(DATA_FILE.c_str());
    if (!fData || fData->IsZombie()) {
        std::cerr << "ERROR: Cannot open data file\n"; return;
    }
    TH1D* hDataCos_raw = (TH1D*)fData->Get("recoProtonCosLAll");
    TH1D* hDataLen_raw = (TH1D*)fData->Get("recoProtonLengthAll");
    if (!hDataCos_raw || !hDataLen_raw) {
        std::cerr << "ERROR: Missing data histograms\n"; return;
    }
    TH1D* hDataCos = RebinPreserveErrors(hDataCos_raw, kNBinsC, kXMinC, kXMaxC, "hDataCos");
    TH1D* hDataLen = RebinPreserveErrors(hDataLen_raw, kNBinsL, kXMinL, kXMaxL, "hDataLen");

    // ── Declare MC histograms (filled from reco tree) ────────────────────
    TH1D* h_nom_cos = new TH1D("h_nom_cos", "", kNBinsC, kXMinC, kXMaxC);
    TH1D* hMC_cos_raw = RebinPreserveErrors((TH1D*)fCAF->Get("reco_matched_proton_cosL"), kNBinsC, kXMinC, kXMaxC, "hMC_cos_raw");    
    TH1D* h_pat_cos_up   = new TH1D("h_pat_cos_up",  "", kNBinsC, kXMinC, kXMaxC);
    TH1D* h_pat_cos_dn   = new TH1D("h_pat_cos_dn",  "", kNBinsC, kXMinC, kXMaxC);
    TH1D* h_unp_cos_up   = new TH1D("h_unp_cos_up",  "", kNBinsC, kXMinC, kXMaxC);
    TH1D* h_unp_cos_dn   = new TH1D("h_unp_cos_dn",  "", kNBinsC, kXMinC, kXMaxC);

    TH1D* h_nom_len = new TH1D("h_nom_len", "", kNBinsL, kXMinL, kXMaxL);
    TH1D* hMC_len_raw = RebinPreserveErrors((TH1D*)fCAF->Get("reco_matched_proton_length"), kNBinsL, kXMinL, kXMaxL, "hMC_len_raw");
    TH1D* h_pat_len_up   = new TH1D("h_pat_len_up",  "", kNBinsL, kXMinL, kXMaxL);
    TH1D* h_pat_len_dn   = new TH1D("h_pat_len_dn",  "", kNBinsL, kXMinL, kXMaxL);
    TH1D* h_unp_len_up   = new TH1D("h_unp_len_up",  "", kNBinsL, kXMinL, kXMaxL);
    TH1D* h_unp_len_dn   = new TH1D("h_unp_len_dn",  "", kNBinsL, kXMinL, kXMaxL);

    // ── Fill from reco tree ───────────────────────────────────────────────
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
        // if (G4ID == -1) continue;
        WKey key = {cafFileIdx, (int)nu_id, G4ID};

        W2 wp = {1.0, 1.0}; if (wmap_pat.count(key)) wp = wmap_pat[key];
        W2 wu = {1.0, 1.0}; if (wmap_unp.count(key)) wu = wmap_unp[key];

        h_nom_cos->Fill(cosT);
        h_pat_cos_up->Fill(cosT, wp[0]); h_pat_cos_dn->Fill(cosT, wp[1]);
        h_unp_cos_up->Fill(cosT, wu[0]); h_unp_cos_dn->Fill(cosT, wu[1]);

        h_nom_len->Fill(length);
        h_pat_len_up->Fill(length, wp[0]); h_pat_len_dn->Fill(length, wp[1]);
        h_unp_len_up->Fill(length, wu[0]); h_unp_len_dn->Fill(length, wu[1]);
    }
        // std::cout << "G4RW nominal cos last bin: " << h_nom_cos->GetBinContent(kNBinsC) << "\n";


    // ── Scale MC to data ──────────────────────────────────────────────────
    // cosTheta
    double scaleCos = (hMC_cos_raw->Integral("width") > 0.0)
                    ? hDataCos->Integral("width") / hMC_cos_raw->Integral("width")
                    : 1.0;
    TH1D* hMC_cos = (TH1D*)hMC_cos_raw->Clone("hMC_cos");
    
    hMC_cos->Scale(scaleCos);

    // length
    double scaleLen = (hMC_len_raw->Integral("width") > 0.0)
                    ? hDataLen->Integral("width") / hMC_len_raw->Integral("width")
                    : 1.0;
    TH1D* hMC_len = (TH1D*)hMC_len_raw->Clone("hMC_len");
    
    hMC_len->Scale(scaleLen);

    // ── Build G4RW bands ──────────────────────────────────────────────────
    // Patched
    auto* gr_pat_cos = BuildG4RWBand(h_nom_cos, h_pat_cos_up, h_pat_cos_dn, hMC_cos);
    auto* gr_pat_len = BuildG4RWBand(h_nom_len, h_pat_len_up, h_pat_len_dn, hMC_len);
    // Unpatched
    auto* gr_unp_cos = BuildG4RWBand(h_nom_cos, h_unp_cos_up, h_unp_cos_dn, hMC_cos);
    auto* gr_unp_len = BuildG4RWBand(h_nom_len, h_unp_len_up, h_unp_len_dn, hMC_len);

    // ── Draw plots ────────────────────────────────────────────────────────
    DrawPlot(hDataCos, hMC_cos, gr_pat_cos,
             "cos(#theta)", "Protons",
             // "FSI (#pm1#sigma)",
            "FSI Systematics",
             OUT_DIR + "/band_dataMC_cosTheta_patched.png", 50);

    // DrawPlot(hDataCos, hMC_cos, gr_unp_cos,
    //          "cos(#theta)", "Reco Protons",
    //          "G4RW Unpatched (#pm1#sigma)",
    //          OUT_DIR + "/band_dataMC_cosTheta_unpatched.png", 50);

    DrawPlot(hDataLen, hMC_len, gr_pat_len,
             "Track Length [cm]", "Protons",
             "FSI Systematics",
        // "FSI Systematics (#pm1#sigma)",
             OUT_DIR + "/band_dataMC_length_patched.png", 70);

    // DrawPlot(hDataLen, hMC_len, gr_unp_len,
    //          "Track Length [cm]", "Reco Protons",
    //          "G4RW Unpatched (#pm1#sigma)",
    //          OUT_DIR + "/band_dataMC_length_unpatched.png", 70);

    std::cout << "\nDone. Plots saved in:\n" << OUT_DIR << "\n";
    fCAF->Close();
    fData->Close();
}

int main()
{
    plot_dataMC_proton_g4rw_band();
    return 0;
}