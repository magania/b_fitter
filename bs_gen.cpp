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

#include "TFile.h"

int main() {
	TFile root_file("toymc.root", "RECREATE");

	RooWorkspace ws;

	ws.factory("{t[-5.0,20.0],et[0,2],cpsi[-1,1],ctheta[-1,1],phi[-3.141592654,3.141592654],D[1]");
	
	ws.factory("{DG[0.09,0,1],Dm[17.77],tau[1.5,1.0,2.0],beta[0.25,-1.6,1.6],delta_l[0.2,-INF,INF],delta_p[2.95,-INF,INF],delta_s[0.2,-INF,INF],fs[0.2,0,1],s[1.0,0.8,1.5]}");
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

        ws.factory("Landau::error_signal(et,mean[0.0741194,0,1],sigma[0.0175981,0,1])");


	for (int i=0; i<10; i++){
		TString data_name = "data_";
		data_name+=i;
		TString fit_name = "fit_";
		fit_name+=i;
		RooDataSet* et_data = ws.pdf("error_signal")->generate(*ws.var("et"),5000);

		*ws.var("A02")= 0.6;
		*ws.var("Al2")= 0.3;
		*ws.var("DG")= 0.09;
		*ws.var("tau")= 1.5;
		*ws.var("beta")= 0.25;
		*ws.var("delta_l")= 0.2;
		*ws.var("delta_p")= 2.95;
		*ws.var("delta_s")= 0.2;
		*ws.var("fs")= 0.2;
		*ws.var("s")= 1.0;

		RooDataSet *data = time_angle.generate(RooArgSet(*ws.var("t"), *ws.var("cpsi"), *ws.var("ctheta"), *ws.var("phi")), RooFit::ProtoData(*et_data));
		data->Write(data_name);

		RooFitResult *fit_result = time_angle.fitTo(*data, RooFit::Save(kTRUE), RooFit::ConditionalObservables(*ws.var("et")), RooFit::NumCPU(2), RooFit::PrintLevel(-1));
		fit_result->Write(fit_name);
	}
/*
        gROOT->SetStyle("Plain");

        TCanvas canvas("canvas", "canvas", 800,400);
        canvas.Divide(2);

        canvas.cd(1);
        RooPlot *t_frame = ws.var("t")->frame();
        data->plotOn(t_frame, RooFit::MarkerSize(0.3));
        gPad->SetLogy(1);
        t_frame->Draw();

        canvas.cd(2);
        RooPlot *et_frame = ws.var("et")->frame();
        data->plotOn(et_frame,RooFit::MarkerSize(0.2));
        ws.pdf("error_signal")->plotOn(et_frame);
        gPad->SetLogy(1);
        et_frame->Draw();

        canvas.SaveAs("t.png"); 
*/

	root_file.Close();
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
//      canvas.SaveAs("t.png");


}
