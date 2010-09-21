/*
 *  Author: Ricardo Magana-Villalba
 *          magania@fnal.gov
 *
 *  September 2010
 */


#include <iostream>

#include <TFile.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TEntryList.h>
#include <TTree.h>

#include <RooBsTimeAngle.hpp>

#include <RooWorkspace.h>
#include <RooDataSet.h>

#include "RooFit.h"
#include "RooGlobalFunc.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooPlot.h"
#include "RooCategory.h"
#include "RooFitResult.h"
#include "RooBinning.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"
#include "RooHistPdf.h"
#include "RooTruthModel.h"


int main() {
//	TFile fit_file("bs_fit.root");
//	RooDataSet* data = (RooDataSet*) fit_file.Get("data");
//	THDAcceptance* acceptance = (THDAcceptance*) fit_file.Get("acceptance");

	RooWorkspace ws;
//	ws.import(*data);

	ws.factory("{t[-3,15],et[0,1],cpsi[-1,1],ctheta[-1,1],phi[-3.15,3.15],D[-1,-1,1}");
	ws.factory("{A02[0.524],All2[0.231],Ap2[0.245],tau_L[1.4],tau_H[1.6],Dm[17.77],beta[0.25,0,1],delta_p[2.95],delta_l[0.2],delta_s[0.2],Fs[0.1]}");

        ws.factory("GaussModel::resol(t,meanRes[0.],sigmaRes[0.0001,0,1])");
//        ws.factory("GaussModel::resol(t,meanRes[0.],scaleRes[1.0,0.8,1.2],et)");
	RooTruthModel rtrue("rtrue","rtrue", *ws.var("t"));

	RooBsTimeAngle time_angle("time_angle", "time_angle",
		*ws.var("t"),
		*ws.var("cpsi"),
		*ws.var("ctheta"),
		*ws.var("phi"),
		*ws.var("A02"),
		*ws.var("All2"),
		*ws.var("Ap2"),
		*ws.var("tau_L"),
		*ws.var("tau_H"),
		*ws.var("Dm"),
		*ws.var("beta"),
		*ws.var("delta_p"),
		*ws.var("delta_l"),
		*ws.var("delta_s"),
		*ws.var("Fs"),
		*ws.var("D"),
		rtrue);
//		*((RooResolutionModel*)ws.allResolutionModels().find("resol")));

	RooDataSet *data = time_angle.generate(RooArgSet(*ws.var("t"), *ws.var("cpsi"), *ws.var("ctheta"), *ws.var("phi")),1000);


/*
        TFile fit_file("bs_fit.root");
        TTree* fit_tree = (TTree*) fit_file.Get("tree");

        Double_t d_pdl, d_epdl, d_cpsi, d_ctheta, d_phi;
        fit_tree->SetBranchAddress("bs_pdl", &d_pdl);
        fit_tree->SetBranchAddress("bs_epdl", &d_epdl);
        fit_tree->SetBranchAddress("bs_angle_cpsi", &d_cpsi);
        fit_tree->SetBranchAddress("bs_angle_ctheta", &d_ctheta);
        fit_tree->SetBranchAddress("bs_angle_phi", &d_phi);

        const char* cut = "mc_match == 1 && mu_plus_nseg==3 && mu_minus_nseg==3";

        fit_tree->Draw(">>entry_list", cut, "entrylist");
        TEntryList *event_list = (TEntryList*)gDirectory->Get("entry_list");

        Long64_t n_entries = event_list->GetN();
        std::cout << "Selected: " << n_entries << " events." << std::endl;

	RooRealVar *t = ws.var("t");
	RooRealVar *et = ws.var("et");
	RooRealVar *cpsi = ws.var("cpsi");
	RooRealVar *ctheta = ws.var("ctheta");
	RooRealVar *phi = ws.var("phi");
	RooRealVar *D = ws.var("D");

	RooDataSet *data = new RooDataSet("data", "data", RooArgSet(*t,*et,*cpsi,*ctheta,*phi));

        for (Long64_t i=0; i<n_entries; i++){
		fit_tree->GetEntry(event_list->GetEntry(i));
		*t = d_pdl/0.0299792458;
		*et = d_epdl/0.0299792458;
		*cpsi = d_cpsi;
		*ctheta = d_ctheta;
		*phi = d_phi;

		data->add(RooArgSet(*t,*et,*cpsi,*ctheta,*phi));
        }

*/
	ws.var("D")->setConstant(kTRUE); 
	time_angle.fitTo(*data);

	gROOT->SetStyle("Plain");

	TCanvas canvas("canvas", "canvas", 800,800);
	canvas.Divide(2,2);
	
	canvas.cd(1);
	RooPlot *t_frame = ws.var("t")->frame();
	data->plotOn(t_frame, RooFit::MarkerSize(0.2));
//	time_angle.plotOn(t_frame);
	gPad->SetLogy(1);
	t_frame->Draw();

	canvas.cd(2);
	gPad->SetLogy(0);
	RooPlot *cpsi_frame = ws.var("cpsi")->frame();
	data->plotOn(cpsi_frame,RooFit::MarkerSize(0.2));
	cpsi_frame->Draw();

	canvas.cd(3);
	RooPlot *ctheta_frame = ws.var("ctheta")->frame();
	data->plotOn(ctheta_frame,RooFit::MarkerSize(0.2));
	ctheta_frame->Draw();

	canvas.cd(4);
	RooPlot *phi_frame = ws.var("phi")->frame();
	data->plotOn(phi_frame,RooFit::MarkerSize(0.2));
	phi_frame->Draw();

	canvas.SaveAs("t.png");


}
