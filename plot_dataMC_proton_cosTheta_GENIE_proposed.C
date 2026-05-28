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
void plot_dataMC_proton_cosTheta_GENIE_proposed()
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

	hData->SetTitle("Proton Candidates");
	hData->GetXaxis()->SetTitle("Reconstructed cos(#theta)");
    hData->GetYaxis()->SetTitle("Counts");

	double maxY = std::max(hData->GetMaximum(), hMC->GetMaximum());
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

    TLegend* leg = new TLegend(0.49, 0.70, 0.77, 0.90);
	leg->SetTextFont(133);
	leg->SetTextSize(18);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);

    leg->SetHeader("NuMI RHC");
    leg->AddEntry(hData, "Validation Reco Data",  "l");
	leg->AddEntry(hMC, "MC Reco (norm. to Data)", "l");
	leg->AddEntry(grScaled, "Sample GENIE Systematics", "f");
    leg->Draw();

	TLatex tL;
	tL.SetNDC();
	tL.DrawLatex(0.20, 0.94, "#bf{DUNE:ND-LAr 2x2}");
    tL.DrawLatex(0.20, 0.84, "#bf{Work in Progress}");

	c->SaveAs("proton_cosTheta_data_reco_GENIE_band.png");


    // ── Shape normalized version ──────────────────────────────
    TH1D* hData_norm = (TH1D*)hData->Clone("hData_norm");
    TH1D* hMC_norm   = (TH1D*)hMC->Clone("hMC_norm");
    
    double normData = hData_norm->Integral();
    double normMC   = hMC_norm->Integral();
    
    if (normData > 0) hData_norm->Scale(1.0 / normData);
    if (normMC   > 0) hMC_norm->Scale(1.0 / normMC);
    
    // Scale GENIE band to unity normalization
    TGraphAsymmErrors* grScaled_norm = (TGraphAsymmErrors*)grScaled->Clone("grScaled_norm");
    for (int i = 0; i < grScaled_norm->GetN(); ++i)
    {
        double x, y;
        grScaled_norm->GetPoint(i, x, y);
        double eyl = grScaled_norm->GetErrorYlow(i);
        double eyh = grScaled_norm->GetErrorYhigh(i);
        grScaled_norm->SetPoint(i, x, y / normMC);
        grScaled_norm->SetPointError(i,
            grScaled_norm->GetErrorXlow(i),
            grScaled_norm->GetErrorXhigh(i),
            eyl / normMC,
            eyh / normMC);
    }
    
    hData_norm->GetYaxis()->SetTitle("Normalized Counts");
    // hData_norm->GetYaxis()->SetRangeUser(0, 1.4 * std::max(hData_norm->GetMaximum(), hMC_norm->GetMaximum()));
    hData_norm->GetYaxis()->SetRangeUser(0, 50.0 / normMC);
    
    TCanvas* c2 = new TCanvas("c2", "c2", 800, 600);
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(0.15);
    
    hData_norm->Draw("HIST E");
    grScaled_norm->Draw("2 SAME");
    hMC_norm->Draw("HIST SAME");
    
    TLegend* leg2 = new TLegend(0.49, 0.70, 0.77, 0.90);
    leg2->SetTextFont(133);
    leg2->SetTextSize(18);
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    // leg2->SetHeader("NuMI RHC");
    leg2->SetHeader("#splitline{NuMI RHC}{Data and MC normalized to unity}");
    leg2->AddEntry(hData_norm, "Validation Data Reco", "l");
    leg2->AddEntry(hMC_norm,     "MC Reco", "l");
    leg2->AddEntry(grScaled_norm, "Sample GENIE Systematics", "f");
    leg2->Draw();
    
    tL.DrawLatex(0.20, 0.94, "#bf{DUNE:ND-LAr 2x2}");
    tL.DrawLatex(0.20, 0.84, "#bf{Work in Progress}");
    
    c2->SaveAs("proton_cosTheta_data_reco_GENIE_band_norm.png");
}
