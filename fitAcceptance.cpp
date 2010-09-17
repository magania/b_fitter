#include "TH2.h"
#include "TF2.h"
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

#include "TRandom3.h"
#include "Math/SpecFuncMathMore.h"

//#include "Fit/FitConfig.h"

Double_t HarmonicSphericalY(int l, int m, Double_t ctheta, Double_t phi){
 	if (m==0)
		return TMath::Sqrt( ((2.0*l+1)*TMath::Gamma(l-m+1))/(4.0*TMath::Pi()*TMath::Gamma(l+m+1)) )*ROOT::Math::assoc_legendre(l,m,ctheta);
	if (m>0)
		return TMath::Sqrt( ((2.0*l+1)*TMath::Gamma(l-m+1))/(2.0*TMath::Pi()*TMath::Gamma(l+m+1)) )*ROOT::Math::assoc_legendre(l,m,ctheta)*TMath::Cos(m*phi);
	if (m<0)
	        return TMath::Sqrt( ((2.0*l+1)*TMath::Gamma(l-m+1))/(2.0*TMath::Pi()*TMath::Gamma(l+m+1)) )*ROOT::Math::assoc_legendre(l,-m,ctheta)*TMath::Sin(m*phi);
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

	TFile acceptance_file("bs.root");
	TTree* acceptance_tree = (TTree*) acceptance_file.Get("tree");

	Double_t cpsi, ctheta, phi;
	acceptance_tree->SetBranchAddress("bs_angle_cpsi", &cpsi);
	acceptance_tree->SetBranchAddress("bs_angle_ctheta", &ctheta);
	acceptance_tree->SetBranchAddress("bs_angle_phi", &phi);

	TH2D* acceptance_histo = new TH2D("acceptance_histo", "acceptance", 20, -1.0, 1.0, 20, -TMath::Pi(), TMath::Pi());

	const char* cut = "mc_match == 1 && mu_plus_nseg==3 && mu_minus_nseg==3";

	acceptance_tree->Draw(">>entry_list", cut, "entrylist");
	TEntryList *event_list = (TEntryList*)gDirectory->Get("entry_list");


	Long64_t n_entries = event_list->GetN();
	std::cout << "Selected: " << n_entries << " events." << std::endl;

	for (Long64_t i=0; i<n_entries; i++){
		acceptance_tree->GetEntry(event_list->GetEntry(i));
		acceptance_histo->Fill(ctheta,phi);
		//std::cout << ctheta << ' ' << phi << std::endl;
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

	std::cout << "static const Double_t e_00p0 = " << fitFcn->GetParameter(0)/max_bin << ';' << std::endl;
	std::cout << "static const Double_t e_01n1 = " << fitFcn->GetParameter(1)/max_bin << ';' << std::endl;
	std::cout << "static const Double_t e_01p0 = " << fitFcn->GetParameter(2)/max_bin << ';' << std::endl;
	std::cout << "static const Double_t e_01p1 = " << fitFcn->GetParameter(3)/max_bin << ';' << std::endl;
	std::cout << "static const Double_t e_02n2 = " << fitFcn->GetParameter(4)/max_bin << ';' << std::endl;
	std::cout << "static const Double_t e_02n1 = " << fitFcn->GetParameter(5)/max_bin << ';' << std::endl;
	std::cout << "static const Double_t e_02p0 = " << fitFcn->GetParameter(6)/max_bin << ';' << std::endl;
	std::cout << "static const Double_t e_02p1 = " << fitFcn->GetParameter(7)/max_bin << ';' << std::endl;
	std::cout << "static const Double_t e_02p2 = " << fitFcn->GetParameter(8)/max_bin << ';' << std::endl;

	canvas.SaveAs("acceptance.png");

	return 1;
}

