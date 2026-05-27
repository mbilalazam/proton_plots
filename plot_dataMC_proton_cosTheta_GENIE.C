#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGraphAsymmErrors.h"
#include <iostream>
#include <cmath>

// ------------------------------------------------------------
// Helper: Rebin preserving errors
// ------------------------------------------------------------
TH1D* RebinPreserveErrors(TH1D* h, int nbins, double xmin, double xmax)
{
	TH1D* hR = new TH1D("hR", "", nbins, xmin, xmax);
	hR->Sumw2();

	for (int i = 1; i <= h->GetNbinsX(); ++i)
	{
		double x = h->GetBinCenter(i);
		if (x < xmin || x >= xmax) continue;

		int b = hR->FindBin(x);

		double c_old = hR->GetBinContent(b);
		double e_old = hR->GetBinError(b);

		double c = h->GetBinContent(i);
		double e = h->GetBinError(i);

		hR->SetBinContent(b, c_old + c);
		hR->SetBinError(b, std::sqrt(e_old*e_old + e*e));
	}
	return hR;
}

// ------------------------------------------------------------
// MAIN
// ------------------------------------------------------------
void plot_dataMC_proton_cosTheta_GENIE()
{
	gROOT->LoadMacro("protoDUNEStyle.C");
	gROOT->SetStyle("protoDUNEStyle");
	gROOT->ForceStyle();
    gStyle->SetTitleX(0.25);


	const int    kNBins = 15;
	const double kXMin  = -1.0;
	const double kXMax  =  1.0;

	// --------------------------------------------------------
	// Open files
	// --------------------------------------------------------
	// TFile fMC("/global/u1/m/mazam/MiniRuns/MR6.5/proton_final/mc_minirun6p5_minervaOff_geomContainRecoOnly_Lgt3cm_300cm_endZneg0to5cmVeto_cosZgtNeg0p9.root");
	TFile fMC("/global/u1/m/mazam/CAFs_Analysis/9_May_2026/slides/proton_systematics/output.root");

	TFile fData("/global/u1/m/mazam/sandbox/v11/sandbox_v11_minervaOff_geomContainRecoOnly_Lgt3cm_300cm_endZneg0to5cmVeto_cosZgtNeg0p9.root");

	TFile fSyst("/global/u1/m/mazam/Geant4ReWeighting/2x2_G4RW/test/on_CAFS/step2_reco/nusyst_step2_reco_output.root");

	if (fSyst.IsZombie()) {
		std::cerr << "Cannot open nusyst file\n";
		return;
	}

	// --------------------------------------------------------
	// Load histograms
	// --------------------------------------------------------
	TH1D* hData_raw = (TH1D*)fData.Get("recoProtonCosLAll");
	TH1D* hMC_raw   = (TH1D*)fMC.Get("reco_matched_proton_cosL");

	if (!hData_raw || !hMC_raw) {
		std::cerr << "Missing Data or MC hist\n";
		return;
	}

	// Rebin
	TH1D* hData = RebinPreserveErrors(hData_raw, kNBins, kXMin, kXMax);
	TH1D* hMC   = RebinPreserveErrors(hMC_raw,   kNBins, kXMin, kXMax);

    // std::cout << "GENIE MC cos last bin: " << hMC->GetBinContent(kNBins) << "\n";

	// --------------------------------------------------------
	// Load GENIE sum-of-dials band
	// --------------------------------------------------------
	TGraphAsymmErrors* grBand =
		(TGraphAsymmErrors*)fSyst.Get("gr_band_sumDials_pCosL");

	if (!grBand) {
		std::cerr << "GENIE band not found\n";
		return;
	}

	// Clone locally
	grBand = (TGraphAsymmErrors*)grBand->Clone("grBand_local");

	// --------------------------------------------------------
	// Normalize MC to Data
	// --------------------------------------------------------
	double intData = hData->Integral("width");
	double intMC   = hMC->Integral("width");
	double scale   = (intMC > 0.0) ? intData / intMC : 1.0;

	hMC->Scale(scale);

	// --------------------------------------------------------
	// Convert absolute band → scaled band
	// --------------------------------------------------------
	TGraphAsymmErrors* grScaled = new TGraphAsymmErrors();

	for (int i = 0; i < grBand->GetN(); ++i)
	{
		double x, yNom;
		grBand->GetPoint(i, x, yNom);  // nominal absolute

		int b = hMC->FindBin(x);
		double yMC = hMC->GetBinContent(b);

		double sigma = grBand->GetErrorYhigh(i);

		// convert absolute → fractional
		double frac = (yNom > 0.0) ? sigma / yNom : 0.0;

		double exl = grBand->GetErrorXlow(i);
		double exh = grBand->GetErrorXhigh(i);

		grScaled->SetPoint(i, x, yMC);
		grScaled->SetPointError(i, exl, exh, frac*yMC, frac*yMC);
	}

    grScaled->SetFillColorAlpha(kOrange + 1, 0.50);
    grScaled->SetLineColor(kOrange + 1);

	// --------------------------------------------------------
	// Style
	// --------------------------------------------------------
	hData->SetLineColor(kBlack);
	hData->SetLineWidth(3);

	hMC->SetLineColor(kRed);
	hMC->SetLineWidth(3);

	hData->SetTitle("Protons");
	hData->GetXaxis()->SetTitle("cos(#theta)");
	hData->GetYaxis()->SetTitle("Counts");

	double maxY = std::max(hData->GetMaximum(), hMC->GetMaximum());
	// hData->GetYaxis()->SetRangeUser(0, 1.4*maxY);
    hData->GetYaxis()->SetRangeUser(0, 50);

	// --------------------------------------------------------
	// Draw
	// --------------------------------------------------------
	TCanvas* c = new TCanvas("c", "c", 800, 600);
	gPad->SetLeftMargin(0.18);
	gPad->SetRightMargin(0.15);

	hData->Draw("HIST E");
	grScaled->Draw("2 SAME");
	hMC->Draw("HIST SAME");

	TLegend* leg = new TLegend(0.50, 0.75, 0.80, 0.93);
	leg->SetTextFont(133);
	leg->SetTextSize(18);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);

	leg->AddEntry(hData, "Data Reco", "l");
	leg->AddEntry(hMC, "MC Reco (norm. to Data)", "l");
	// leg->AddEntry(grScaled, "GENIE syst. (sum-of-dials)", "f");
	leg->AddEntry(grScaled, "GENIE Systematics", "f");
    leg->Draw();

	TLatex tL;
	tL.SetNDC();
	tL.DrawLatex(0.20, 0.94, "#bf{DUNE:ND-LAr 2x2}");

	c->SaveAs("proton_cosTheta_data_reco_GENIE_band.png");
}