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

	ws.factory("{t[-INF,INF],cpsi[-1,1],ctheta[-1,1],phi[-3.141592654,3.141592654],D[-1]");

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
		*ws.var("fs"));

//	RooMCStudy toy(time_angle,RooArgSet(*ws.var("t"), *ws.var("cpsi"), *ws.var("ctheta"), *ws.var("phi")),RooFit::FitOptions(RooFit::NumCPU(2)));
//	toy.generateAndFit(100,4000);

	RooDataSet *data = time_angle.generate(RooArgSet(*ws.var("t"), *ws.var("cpsi"), *ws.var("ctheta"), *ws.var("phi")),10000000,RooFit::Verbose(kTRUE));

}
