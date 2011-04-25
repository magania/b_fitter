/*
 *  Author: Ricardo Magana-Villalba
 *          magania@fnal.gov
 *
 *  September 2010
 */


#include "TH2.h"
#include "TF2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TVirtualFitter.h"
#include "TPaveLabel.h"
#include "TStyle.h"
#include "TMath.h"
#include "TROOT.h"
#include "TFrame.h"
#include "TFile.h"
#include "TTree.h"
#include "TEntryList.h"
#include "TCut.h"

#include "TRandom3.h"
#include "Math/SpecFuncMathMore.h"

//#include "Fit/FitConfig.h"


Double_t HarmonicSphericalY(int l, int m, Double_t _ctheta, Double_t _phi) {
        if (m==0)
                return TMath::Sqrt( ((2.0*l+1)*TMath::Gamma(l+1))/(4.0*TMath::Pi()*TMath::Gamma(l+1)) )*ROOT::Math::assoc_legendre(l,0,_ctheta);
        if (m>0)
                return TMath::Sqrt( ((2.0*l+1)*TMath::Gamma(l-m+1))/(2.0*TMath::Pi()*TMath::Gamma(l+m+1)) )*ROOT::Math::assoc_legendre(l,m,_ctheta)*TMath::Cos(m*_phi)*TMath::Power(-1.0,m);
        if (m<0)
                return TMath::Sqrt( ((2.0*l+1)*TMath::Gamma(l+m+1))/(2.0*TMath::Pi()*TMath::Gamma(l-m+1)) )*ROOT::Math::assoc_legendre(l,-m,_ctheta)*TMath::Sin(-m*_phi)*TMath::Power(-1.0,-m);
}

// Sum of background and peak function
Double_t fitFunction(Double_t *x, Double_t *par) {
	Double_t result=0;

	int n_par=0;
	for(int l=0; l<=2; l++)
		for ( int m=-l; m<=l; m++)
			result += par[n_par++]*HarmonicSphericalY(l,m,x[0],x[1]);

	return result;
}

int main() {
	gSystem->Load("libMathMore");
	gROOT->SetStyle("Plain");

	TFile acceptance_file("/work/elsanto-clued0/Z/Bs/tmva/bs_flat.root");
	TTree* acceptance_tree = (TTree*) acceptance_file.Get("tree");

	Double_t cpsi, ctheta, phi;
	acceptance_tree->SetBranchAddress("bs_angle_cpsi", &cpsi);
	acceptance_tree->SetBranchAddress("bs_angle_ctheta", &ctheta);
	acceptance_tree->SetBranchAddress("bs_angle_phi", &phi);

        Double_t muPlusPt, muMinusPt, jpsiPt, muPlusEta, muMinusEta;
	acceptance_tree->SetBranchAddress("mu_plus_pt", &muPlusPt);
	acceptance_tree->SetBranchAddress("mu_minus_pt", &muMinusPt);
	acceptance_tree->SetBranchAddress("mu_plus_eta", &muPlusEta);
	acceptance_tree->SetBranchAddress("mu_minus_eta", &muMinusEta);
	acceptance_tree->SetBranchAddress("jpsi_pt", &jpsiPt);

        TF1 *weightFcnCentral = new TF1("weightFcnCentral","landau(0)+pol1(3)",3.85,45);
        weightFcnCentral->SetParameters(9.76584e+00,1.13316e+01,3.69072e+00,-4.52973e-01,1.44476e-02);

        TF1 *weightFcnForward = new TF1("weightFcnCentral","landau(0)+pol1(3)",3.02,45);
        weightFcnForward->SetParameters(8.93988e+00,1.16932e+01,3.94829e+00,-3.07532e-01,2.27729e-02);

	TH2D* acceptance_histo = new TH2D("acceptance_histo", "acceptance", 20, -1.0, 1.0, 20, -TMath::Pi(), TMath::Pi());

        TCut QC("(k_minus_nSMT>1 && k_plus_nSMT>1 && mu_minus_nSMT>1 && mu_plus_nSMT>1 && ((mu_minus_nseg==1 && mu_plus_nseg==3) || (mu_minus_nseg==3 && mu_plus_nseg==1) || (mu_minus_nseg==3 && mu_plus_nseg==3)) && TMath::Abs(bs_vrt_z-bs_pv_vrt_z)<5.0)");
        TCut PRL("mu_plus_cpt>1.5 && mu_minus_cpt>1.5 && 2.9 < jpsi_mass  && jpsi_mass < 3.3  && k_plus_cpt>0.7 && k_minus_cpt>0.7 && phi_cpt>1.5 && 1.01 < phi_mass_corrected && phi_mass_corrected<1.03 && bs_pt>6 && bs_vrt_chi2<36 && TMath::Abs(bs_pv_vrt_z - bs_vrt_z) < 5 && k_plus_nSMT>1 && k_minus_nSMT>1 && mu_plus_nSMT>1 && mu_minus_nSMT>1  && k_plus_nSMT+k_plus_nCFT>7 && k_minus_nSMT+k_minus_nCFT>7  && mu_plus_nSMT+mu_plus_nCFT>7 && mu_minus_nSMT+mu_minus_nCFT>7  && (mu_plus_nseg==3||mu_plus_nseg==1) && (mu_minus_nseg==3||mu_minus_nseg==1) && ( jpsi_cpt>4 || (mu_plus_cpt>mu_minus_cpt && TMath::Abs(mu_plus_eta)>=1) || (mu_minus_cpt>mu_plus_cpt && TMath::Abs(mu_minus_eta)>=1) )");

	TCut cut0("mc_match == 1");

        TCut BDT10("inclusiveBDT>0.45 && promptBDT>0.42");
        TCut BDT12("inclusiveBDT>0.39 && promptBDT>0.35");
        TCut BDT14("inclusiveBDT>0.33 && promptBDT>0.28");
        TCut BDT16("inclusiveBDT>0.21 && promptBDT>0.29");
        TCut BDT18("inclusiveBDT>0.07 && promptBDT>0.29");
        TCut BDT20("inclusiveBDT>-0.05 && promptBDT>-0.01");

	TCut cut = BDT20 + QC + cut0;

	acceptance_tree->Draw(">>entry_list", cut, "entrylist");
	TEntryList *event_list = (TEntryList*)gDirectory->Get("entry_list");

	Long64_t n_entries = event_list->GetN();
	std::cout << "Selected: " << n_entries << " events." << std::endl;

	for (Long64_t i=0; i<n_entries; i++){
		acceptance_tree->GetEntry(event_list->GetEntry(i));
                 
		Double_t muLeadingPt = muPlusPt>muMinusPt?muPlusPt:muMinusPt;
		Double_t muLeadingEta = muPlusPt>muMinusPt?muPlusEta:muMinusEta;

                if (jpsiPt > 45)
		   jpsiPt = 45;
		Double_t weight = fabs(muLeadingEta)<1?weightFcnCentral->Eval(jpsiPt):weightFcnForward->Eval(jpsiPt);
                if (weight < 0)
		  weight = 0;

		acceptance_histo->Fill(ctheta,phi,weight);
//		std::cout << ctheta << ' ' << phi<< ' ' << weight << std::endl;
//		std::cout << muPlusPt << ' ' << muMinusPt << ' ' << muPlusEta << ' ' << muMinusEta << ' ' << weight << std::endl;
	}

//	acceptance_histo->Draw();

	TF2 *fitFcn = new TF2("fitFcn", fitFunction, -1.0, 1.0, -TMath::Pi(), TMath::Pi(), 9);
	fitFcn->SetParameters(0,0,0,0,0,0,0,0,0);
	//acceptance_histo->Fit(fitFcn,"EV");

	acceptance_histo->Fit(fitFcn,"EV0");

        gStyle->SetPalette(1);

	TCanvas canvas("canvas", "canvas", 600, 600);
	canvas.Divide(2,2);
	canvas.cd(1);
	acceptance_histo->Draw("COLZ");
//	acceptance_histo->ProjectionX()->Draw();
	acceptance_histo->Fit(fitFcn,"EV0");
//	fitFcn->Draw();
//
	TH2D* generated_histo = new TH2D("generated_histo", "acceptance", 20, -1.0, 1.0, 20, -TMath::Pi(), TMath::Pi());

        gRandom = new TRandom3();
	Double_t x, y;
	for(int i=0; i<n_entries*1000; i++){
		fitFcn->GetRandom2(x,y);
		generated_histo->Fill(x,y,0.001);
	}

	canvas.cd(2);
	generated_histo->Draw("COLZ");

	canvas.cd(3);
	acceptance_histo->ProjectionX()->Draw("E");
	generated_histo->ProjectionX()->Draw("same");

	canvas.cd(4);
	acceptance_histo->ProjectionY()->Draw("E");
	generated_histo->ProjectionY()->Draw("same");

	Double_t max_bin =  acceptance_histo->GetBinContent(acceptance_histo->GetMaximumBin());


	std::cout << acceptance_histo->Chi2Test(generated_histo) << std::endl;

	std::cout << "static const Double_t e_00p0 = " << fitFcn->GetParameter(0)/max_bin << ';' << std::endl;
	std::cout << "static const Double_t e_01n1 = " << fitFcn->GetParameter(1)/max_bin << ';' << std::endl;
	std::cout << "static const Double_t e_01p0 = " << fitFcn->GetParameter(2)/max_bin << ';' << std::endl;
	std::cout << "static const Double_t e_01p1 = " << fitFcn->GetParameter(3)/max_bin << ';' << std::endl;
	std::cout << "static const Double_t e_02n2 = " << fitFcn->GetParameter(4)/max_bin << ';' << std::endl;
	std::cout << "static const Double_t e_02n1 = " << fitFcn->GetParameter(5)/max_bin << ';' << std::endl;
	std::cout << "static const Double_t e_02p0 = " << fitFcn->GetParameter(6)/max_bin << ';' << std::endl;
	std::cout << "static const Double_t e_02p1 = " << fitFcn->GetParameter(7)/max_bin << ';' << std::endl;
	std::cout << "static const Double_t e_02p2 = " << fitFcn->GetParameter(8)/max_bin << ';' << std::endl;

	canvas.SaveAs("acceptance.eps");

	return 1;
}

