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

#include <RooBsTimeAngle.h>

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
#include "RooGaussModel.h"
#include "RooDataHist.h"

#include "RooMCStudy.h"


int main() {
//	TFile fit_file("bs_fit.root");
//	RooDataSet* data = (RooDataSet*) fit_file.Get("data");
//	THDAcceptance* acceptance = (THDAcceptance*) fit_file.Get("acceptance");

	RooWorkspace ws;
//	ws.import(*data);

	ws.factory("{t[-4,15],et[0.0,4.0],cpsi[-1,1],ctheta[-1,1],phi[-3.15,3.15],D[-1]");

//	ws.factory("{a0[0.712963,0,1],al[0.473377,0,1],ap[0.487511,0,1]}");

/*
	ws.factory("{a02[0.3,0,1],al2[0.3,0,1]}");
	ws.factory("expr:a0('TMath::Sqrt(@0)',a02)");
	ws.factory("expr:al('TMath::Sqrt(@0)',al2)");
	ws.factory("expr:ap('TMath::Sqrt(1-@0-@1)',a02,al2)");
*/

//	ws.factory("{Dm[17.77],beta[-1.58,1.58],delta_p[2.95,-3.15,3.15],delta_l[0.2,-3.15,3.15],delta_s[0.2,-3.15,3.15],Fs[0.0,0,1]}");
//	ws.factory("{DG[0.09],Dm[17.77],tau[1.5,1.0,2.0],beta[0.25],delta_l[0.2],delta_p[2.95],delta_s[0.2],fs[0.2],s[1.0]}");
	ws.factory("{DG[0.09,0,1],Dm[17.77],tau[1.5,1.0,2.0],beta[0.25,-1,1],delta_l[0.2,-3,3],delta_p[2.95,-3,3],delta_s[0.2,-3,3],fs[0.2,0,1],s[1.0,0.8,1.5]}");
//	ws.factory("{DG[0.09,0,1],Dm[17.77],tau[1.5,1.0,2.0],beta[0.25,-1,1],delta_l[0.2,-3,3],delta_p[2.95,-3,3],delta_s[0.5],fs[0.2],s[1.0,0.8,1.5]}");

	ws.factory("{A02[0.6,0,1],Al2[0.3,0,1]}");
	ws.factory("expr:Ap2('TMath::Sqrt(1-@0-@1)',A02,Al2)");

        ws.factory("expr:z('TMath::Cos(2.0*@0)*@1*@2/2.0',beta,DG,tau");
        ws.factory("expr:y('(1+@0)/(1-@0)',z");
        ws.factory("expr:ap('TMath::Sqrt(@0*@1/(@1+(1.0-@1)*@2))',Ap2,y,Ap2)");
        ws.factory("expr:al('TMath::Sqrt(@0*@1/(@1+(1.0-@1)*@2))',Al2,y,Ap2)");
        ws.factory("expr:a0('TMath::Sqrt(@0*@1/(@1+(1.0-@1)*@2))',A02,y,Ap2)");
	

	RooBsTimeAngle time_angle("time_angle", "time_angle",
		*ws.var("t"),
		*ws.var("et"),
		*ws.var("cpsi"),
		*ws.var("ctheta"),
		*ws.var("phi"),
		*ws.var("D"),
		*ws.function("a0"),
		*ws.function("al"),
		*ws.function("ap"),
		*ws.var("DG"),
		*ws.var("Dm"),
		*ws.var("tau"),
		*ws.var("beta"),
		*ws.var("delta_l"),
		*ws.var("delta_p"),
		*ws.var("delta_s"),
		*ws.var("fs"),
		*ws.var("s"));

//	RooMCStudy toy(time_angle,RooArgSet(*ws.var("t"), *ws.var("cpsi"), *ws.var("ctheta"), *ws.var("phi")),RooFit::FitOptions(RooFit::NumCPU(2)));
//	toy.generateAndFit(100,4000);

//	RooDataSet *data = time_angle.generate(RooArgSet(*ws.var("t"), *ws.var("cpsi"), *ws.var("ctheta"), *ws.var("phi")),5000,RooFit::Verbose(kTRUE));


        TFile fit_file("bs_fit.root");
        TTree* fit_tree = (TTree*) fit_file.Get("tree");

        Double_t d_pdl, d_epdl, d_cpsi, d_ctheta, d_phi;
        fit_tree->SetBranchAddress("bs_pdl", &d_pdl);
        fit_tree->SetBranchAddress("bs_epdl", &d_epdl);
        fit_tree->SetBranchAddress("bs_angle_cpsi", &d_cpsi);
        fit_tree->SetBranchAddress("bs_angle_ctheta", &d_ctheta);
        fit_tree->SetBranchAddress("bs_angle_phi", &d_phi);

        const char* cut = "mc_match == 1";

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
		*et = d_epdl/0.0299792458 ;
		*cpsi = d_cpsi;
		*ctheta = d_ctheta;
		*phi = d_phi;

/*
		if (d_epdl/0.0299792458 > 2) { 
			cout << "Event Dropped: et=" << d_epdl/0.0299792458 << endl;
			continue;
		}
*/


//		if ( (d_epdl/0.299792458) > 1.5 ) continue;

		data->add(RooArgSet(*t,*et,*cpsi,*ctheta,*phi));
        }


/*
	gROOT->SetStyle("Plain");

	TCanvas canvas("canvas", "canvas", 800,800);
	canvas.Divide(2,2);
	canvas.cd(1);
	gPad->SetMargin(1,1,1,1);
	toy.plotParam(*ws.var("beta"),RooFit::MarkerSize(0.2))->Draw();
	canvas.cd(2);
	gPad->SetMargin(1,1,1,1);
	toy.plotError(*ws.var("beta"),RooFit::MarkerSize(0.2))->Draw();
	canvas.cd(3);
	gPad->SetMargin(1,1,1,1);
	toy.plotPull(*ws.var("beta"),RooFit::MarkerSize(0.2),RooFit::FitGauss(kTRUE))->Draw();
	canvas.cd(4);
	gPad->SetMargin(1,1,1,1);
	toy.plotNLL(RooFit::MarkerSize(0.2))->Draw();
*/

	ws.factory("Landau::error_signal(et,mean[0.07,0,1],sigma[0.01,0,1])");
	ws.pdf("error_signal")->fitTo(*data);

//	time_angle.fitTo(*data,RooFit::Verbose(kTRUE), RooFit::ConditionalObservables(*ws.var("et")), RooFit::NumCPU(2));

	RooDataSet* et_data = ws.pdf("error_signal")->generate(*et,1000000);
//	*ws.var("fs") = 0.0;
	RooDataSet *datax = time_angle.generate(RooArgSet(*ws.var("t"), *ws.var("cpsi"), *ws.var("ctheta"), *ws.var("phi")), RooFit::ProtoData(*et_data));
	time_angle.fitTo(*datax,RooFit::Verbose(kTRUE), RooFit::NumCPU(2));

	return 0;
/*
	gROOT->SetStyle("Plain");

	TCanvas canvas("canvas", "canvas", 800,800);
	canvas.Divide(2);

	*ws.var("fs") = 0.0;
	RooDataSet et_data("et_data", "et_data", data, *et);
	RooDataSet *data2 = time_angle.generate(RooArgSet(*ws.var("t"), *ws.var("cpsi"), *ws.var("ctheta"), *ws.var("phi")), RooFit::ProtoData(et_data));

	canvas.cd(1);
	RooPlot *t_frame = ws.var("t")->frame();
	data->plotOn(t_frame, RooFit::MarkerSize(0.3), RooFit::Rescale(1));
	data2->plotOn(t_frame,
		RooFit::LineColor(kBlue), RooFit::DrawOption("L"));
	gPad->SetLogy(1);
	t_frame->Draw();

	canvas.cd(2);
	gPad->SetLogy(0);
	RooPlot *et_frame = ws.var("et")->frame();
	data->plotOn(et_frame,RooFit::MarkerSize(0.2));
	ws.pdf("error_signal")->plotOn(et_frame);
	gPad->SetLogy(1);
	et_frame->Draw();

	canvas.SaveAs("t.png");
*/	
/*
	canvas.cd(2);
	gPad->SetLogy(0);
	RooPlot *cpsi_frame = ws.var("cpsi")->frame();
	data->plotOn(cpsi_frame,RooFit::MarkerSize(0.2), RooFit::Rescale(1));
	data2->plotOn(cpsi_frame,
		RooFit::LineColor(kBlue), RooFit::DrawOption("L"));
	cpsi_frame->Draw();

	canvas.cd(3);
	RooPlot *ctheta_frame = ws.var("ctheta")->frame();
	data->plotOn(ctheta_frame,RooFit::MarkerSize(0.2), RooFit::Rescale(1));
	data2->plotOn(ctheta_frame,
		RooFit::LineColor(kBlue), RooFit::DrawOption("L"));
	ctheta_frame->Draw();

	canvas.cd(4);
	RooPlot *phi_frame = ws.var("phi")->frame();
	data->plotOn(phi_frame,RooFit::MarkerSize(0.2), RooFit::Rescale(1));
	data2->plotOn(phi_frame,
		RooFit::LineColor(kBlue), RooFit::DrawOption("L"));
	phi_frame->Draw();
*/
//	canvas.SaveAs("t.png");
}
