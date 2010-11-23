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
#include <RooAngle.h>
#include <RooErrPdf.h>

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
#include "RooHist.h"

#include "RooMCStudy.h"


int main() {
//	TFile fit_file("bs_fit.root");
//	RooDataSet* data = (RooDataSet*) fit_file.Get("data");
//	THDAcceptance* acceptance = (THDAcceptance*) fit_file.Get("acceptance");

	TFile data_file("dataset.root");
	RooWorkspace* ws = (RooWorkspace*) data_file.Get("rws");
	cout << ws->data("data") << endl;

//	ws.factory("{a0[0.712963,0,1],al[0.473377,0,1],ap[0.487511,0,1]}");

//	ws.factory("{Dm[17.77],beta[-1.58,1.58],delta_p[2.95,-3.15,3.15],delta_l[0.2,-3.15,3.15],delta_s[0.2,-3.15,3.15],Fs[0.0,0,1]}");
//	ws.factory("{DG[0.09],Dm[17.77],tau[1.5,1.0,2.0],beta[0.25],delta_l[0.2],delta_p[2.95],delta_s[0.2],fs[0.2],s[1.0]}");
	ws->factory("{DG[0.09,0,1],Dm[17.77],tau[1.5,1.0,2.0],beta[0.25,-1,1],delta_l[-2.93],delta_p[2.93],delta_s[0.2,-3,3],fs[0.2,0,1],scale[1.0,0.5,1.5]}");
//	ws.factory("{DG[0.09,0,1],Dm[17.77],tau[1.5,1.0,2.0],beta[0.25,-1,1],delta_l[0.2,-3,3],delta_p[2.95,-3,3],delta_s[0.5],fs[0.2],s[1.0,0.8,1.5]}");

	ws->factory("{A02[0.6,0,1],Al2[0.3,0,1]}");
	ws->factory("expr:Ap2('TMath::Sqrt(1-@0-@1)',A02,Al2)");

        ws->factory("expr:z('TMath::Cos(2.0*@0)*@1*@2/2.0',beta,DG,tau");
        ws->factory("expr:y('(1+@0)/(1-@0)',z");
        ws->factory("expr:ap('TMath::Sqrt(@0*@1/(@1+(1.0-@1)*@2))',Ap2,y,Ap2)");
        ws->factory("expr:al('TMath::Sqrt(@0*@1/(@1+(1.0-@1)*@2))',Al2,y,Ap2)");
        ws->factory("expr:a0('TMath::Sqrt(@0*@1/(@1+(1.0-@1)*@2))',A02,y,Ap2)");
	

	ws->factory("Gaussian::massSignal(m,mu[5.2,5.4],sigmaM[0.03,0.02,0.04]");

	RooBsTimeAngle timeAngle("timeAngle", "timeAngle",
		*ws->var("t"),
		*ws->var("et"),
		*ws->var("cpsi"),
		*ws->var("ctheta"),
		*ws->var("phi"),
		*ws->var("D"),
		*ws->function("a0"),
		*ws->function("al"),
		*ws->function("ap"),
		*ws->var("DG"),
		*ws->var("Dm"),
		*ws->var("tau"),
		*ws->var("beta"),
		*ws->var("delta_l"),
		*ws->var("delta_p"),
		*ws->var("delta_s"),
		*ws->var("fs"),
		*ws->var("scale"));

	ws->factory("Landau::errorSignal(et,meanErrSignal[0.0741194,0,1],sigmaErrSignal[0.0175981,0,1])");

	ws->import(timeAngle);

	ws->factory("PROD::signalTimeAngle(timeAngle|et,errorSignal");
	ws->factory("PROD::signal(massSignal,signalTimeAngle)");


	ws->factory("Landau::errLBkg(et,mean[0.072,0,1],sigma[0.011,0,1])");
	ws->factory("Landau::errRBkg(et,mean2[0.08,0,1],sigma2[0.011,0,1])");
//	ws->factory("Gaussian::errLBkg(et,mean[0.0741194,0,1],sigma[0.0175981,0,1])");
	ws->factory("{errA[3.8,0.5,10],errB[0.011,0.01,0.2],errS[0.035,0,0.05]}");
	ws->factory("{errA2[2.24,0.5,10],errB2[0.032,0.01,0.2],errS2[0.024,0,0.05]}");
	RooErrPdf errBkg("errBkg", "errBkg", *ws->var("et"), *ws->var("errA"), *ws->var("errB"),*ws->var("errS"));
	RooErrPdf err2Bkg("err2Bkg", "err2Bkg", *ws->var("et"), *ws->var("errA2"), *ws->var("errB2"),*ws->var("errS2"));
	ws->import(errBkg);
	ws->import(err2Bkg);

	ws->factory("SUM::errorBkg(xe[0.61,0,1]*errBkg,err2Bkg)");

        ws->factory("GaussModel::resolution(t,0,scale,et)");
//        ws->factory("GaussModel::resolution(t,0,sigma[0.02,0.01,0.1])");
/*        ws->factory("GaussModel::resW(t,0,sigma[0.02,0.01,0.5])");
        ws->factory("GaussModel::resG(t,0,sigma1[0.02,0.01,0.1])");
        ws->factory("AddModel::resolution({resW,resG},{f0[0,1])");
*/
        ws->factory("RSUM::tBkg(xr[0.73,0,1]*resolution, xn[0.15,0,1]*Decay(t,tauNeg[0.36,0.01,1],resolution,Flipped), xp[0.87,0,1]*Decay(t,tauPos[0.41,0.1,1],resolution,SingleSided), Decay(t,tauPosPos[1.9,0.5,3],resolution,SingleSided)");
//        ws->factory("RSUM::tBkg(xn[0.14,0,1]*Decay(t,tauNeg[0.13,0.01,1],resolution,Flipped), Decay(t,tauPosPos[0.78,0.5,2],resolution,SingleSided)");

        ws->factory("PROD::timeBkg(tBkg|et,errorBkg)");

//        ws->factory("{e01m1[-1,1],e01p0[-1,1],e01p1[-1,1],e02m2[-1,1],e02m1[-1,1],e02p0[-1,1],e02p1[-1,1],e02p2[-1,1]}");
//        ws->factory("{e20p0[-1,1],e21m1[-1,1],e21p0[-1,1],e21p1[-1,1],e22m2[-1,1],e22m1[-1,1],e22p0[-1,1],e22p1[-1,1],e22p2[-1,1]}");

        ws->factory("{e01m1[0],e01p0[0],e01p1[0],e02m2[0],e02m1[0],e02p0[-1,1],e02p1[0],e02p2[-1,1]}");
        ws->factory("{e20p0[-1,1],e21m1[0],e21p0[0],e21p1[0],e22m2[0],e22m1[0],e22p0[0],e22p1[0],e22p2[0]}");

	RooAngle angle("angle","angle", *ws->var("cpsi"), *ws->var("ctheta"), *ws->var("phi"),
	               *ws->var("e01m1"),
	               *ws->var("e01p0"),
	               *ws->var("e01p1"),
	               *ws->var("e02m2"),
	               *ws->var("e02m1"),
	               *ws->var("e02p0"),
	               *ws->var("e02p1"),
	               *ws->var("e02p2"),
	               *ws->var("e20p0"),
	               *ws->var("e21m1"),
	               *ws->var("e21p0"),
	               *ws->var("e21p1"),
	               *ws->var("e22m2"),
	               *ws->var("e22m1"),
	               *ws->var("e22p0"),
	               *ws->var("e22p1"),
	               *ws->var("e22p2")
		       );

        ws->import(angle);
	angle.Print();

	ws->pdf("errorBkg")->fitTo(*ws->data("dataBkg"), RooFit::NumCPU(2));
	ws->var("mean")->setConstant(kTRUE);
	ws->var("sigma")->setConstant(kTRUE);
	ws->var("errA")->setConstant(kTRUE);
	ws->var("errB")->setConstant(kTRUE);
	ws->var("errS")->setConstant(kTRUE);
	ws->var("xe")->setConstant(kTRUE);
	ws->var("errA2")->setConstant(kTRUE);
	ws->var("errB2")->setConstant(kTRUE);
	ws->var("errS2")->setConstant(kTRUE);

	ws->factory("Polynomial::massBkg(m,{slope[-0.125,-1,-0.1]})");

	ws->factory("PROD::background(massBkg,timeBkg,angle)");
//	ws->pdf("massBkg")->fitTo(*ws->data("data")/*, RooFit::ConditionalObservables(*ws->var("et"))*/, RooFit::NumCPU(2));
//	ws->pdf("angle")->fitTo(*ws->data("data")/*, RooFit::ConditionalObservables(*ws->var("et"))*/, RooFit::NumCPU(2));
//	ws->pdf("timeBkg")->fitTo(*ws->data("data")/*, RooFit::ConditionalObservables(*ws->var("et"))*/, RooFit::NumCPU(2), RooFit::Verbose(kTRUE));

	ws->factory("SUM::model(xs[0.05,0,1]*signal,background)");

	ws->pdf("background")->fitTo(*ws->data("dataBkg") ,RooFit::NumCPU(2));
	ws->var("tauNeg")->setConstant(kTRUE);
	ws->var("tauPos")->setConstant(kTRUE);
	ws->var("tauPosPos")->setConstant(kTRUE);
	ws->var("xn")->setConstant(kTRUE);
	ws->var("xp")->setConstant(kTRUE);
	ws->var("xr")->setConstant(kTRUE);
	ws->var("scale")->setConstant(kTRUE);
	ws->var("e02p0")->setConstant(kTRUE);
	ws->var("e02p2")->setConstant(kTRUE);
	ws->var("e20p0")->setConstant(kTRUE);
	ws->var("slope")->setConstant(kTRUE);

	ws->pdf("model")->fitTo(*ws->data("data"), RooFit::ConditionalObservables(*ws->var("D")) ,RooFit::NumCPU(2), RooFit::Verbose(kTRUE));

return 0;

        cout << "Printing ... " << endl;
	gROOT->SetStyle("Plain");

	TPad *plotPad, *resPad;
	TCanvas canvas("canvas", "canvas", 900,800);
	canvas.Divide(3,2);

//	TCanvas canvasMass("canvasMass", "canvas mass", 400,600);
	canvas.cd(1);
	plotPad = new TPad("plotPadMass", "Plot Pad", 0,0.2, 1,1);
	resPad = new TPad("resPadMass", "Res Pad", 0,0, 1,0.2);
	plotPad->Draw();
	resPad->Draw();

	plotPad->cd();
	RooPlot *m_frame = ws->var("m")->frame();
	ws->data("data")->plotOn(m_frame, RooFit::MarkerSize(0.3));
	ws->pdf("model")->plotOn(m_frame, RooFit::ProjWData(RooArgSet(*ws->var("D"), *ws->var("et")), *ws->data("data")), RooFit::NumCPU(2));
	m_frame->Draw();

	resPad->cd();
	RooHist *massResidual = m_frame->pullHist();
        RooPlot *massResidualPlot  = ws->var("m")->frame(RooFit::Title("Residual Distribution"));
	massResidualPlot->addPlotable(massResidual,"E");
	massResidualPlot->Draw();

//	TCanvas canvasTime("canvasTime", "canvas time", 400,600);
	canvas.cd(2);
	plotPad = new TPad("plotPadTime", "Plot Pad", 0,0.2, 1,1);
	resPad = new TPad("resPadTime", "Res Pad", 0,0, 1,0.2);
	plotPad->Draw();
	resPad->Draw();

	plotPad->cd();
	ws->var("t")->setMin(-1.5);
        ws->var("t")->setMax( 5.0);

	RooPlot *t_frame = ws->var("t")->frame();
	ws->data("data")->plotOn(t_frame, RooFit::MarkerSize(0.3));
//	ws->pdf("timeBkg")->plotOn(t_frame);
//	ws->pdf("tBkg")->plotOn(t_frame,RooFit::ProjWData(*ws->var("et"),*ws->data("data")), RooFit::NumCPU(2));
	ws->pdf("model")->plotOn(t_frame, RooFit::ProjWData(RooArgSet(*ws->var("D"), *ws->var("et")), *ws->data("data")), RooFit::NumCPU(2));
	gPad->SetLogy(1);
	t_frame->Draw();

	resPad->cd();
	RooHist *timeResidual = t_frame->pullHist();
        RooPlot *timeResidualPlot  = ws->var("t")->frame(RooFit::Title("Residual Distribution"));
	timeResidualPlot->addPlotable(timeResidual,"E");
	timeResidualPlot->Draw();

//	TCanvas canvasEt("canvasEt", "canvas et", 400,600);
	canvas.cd(3);
	plotPad = new TPad("plotPadEt", "Plot Pad", 0,0.2, 1,1);
	resPad = new TPad("resPadEt", "Res Pad", 0,0, 1,0.2);
	plotPad->Draw();
	resPad->Draw();

	plotPad->cd();
        ws->var("et")->setMax(0.6);
	gPad->SetLogy(0);
	RooPlot *et_frame = ws->var("et")->frame();
	ws->data("data")->plotOn(et_frame,RooFit::MarkerSize(0.2));
	ws->pdf("errorBkg")->plotOn(et_frame);
	gPad->SetLogy(1);
	et_frame->Draw();

	resPad->cd();
	RooHist *etResidual = et_frame->pullHist();
        RooPlot *etResidualPlot  = ws->var("et")->frame(RooFit::Title("Residual Distribution"));
	etResidualPlot->addPlotable(etResidual,"E");
	etResidualPlot->Draw();


//	TCanvas canvasCpsi("canvasCpsi", "canvas Cpsi", 400,600);
	canvas.cd(4);
	plotPad = new TPad("plotPadCpsi", "Plot Pad", 0,0.2, 1,1);
	resPad = new TPad("resPadCpsi", "Res Pad", 0,0, 1,0.2);
	plotPad->Draw();
	resPad->Draw();

	plotPad->cd();
	gPad->SetLogy(0);
	RooPlot *cpsi_frame = ws->var("cpsi")->frame();
	ws->data("data")->plotOn(cpsi_frame,RooFit::MarkerSize(0.2));
	ws->pdf("model")->plotOn(cpsi_frame, RooFit::ProjWData(RooArgSet(*ws->var("D"), *ws->var("et")), *ws->data("data")), RooFit::NumCPU(2));
	cpsi_frame->Draw();

	resPad->cd();
	RooHist *cpsiResidual = cpsi_frame->pullHist();
        RooPlot *cpsiResidualPlot  = ws->var("cpsi")->frame(RooFit::Title("Residual Distribution"));
	cpsiResidualPlot->addPlotable(cpsiResidual,"E");
	cpsiResidualPlot->Draw();


//	TCanvas canvasCtheta("canvasCtheta", "canvas Ctheta", 400,600);
	canvas.cd(5);
	plotPad = new TPad("plotPadCtheta", "Plot Pad", 0,0.2, 1,1);
	resPad = new TPad("resPadCtheta", "Res Pad", 0,0, 1,0.2);
	plotPad->Draw();
	resPad->Draw();

	plotPad->cd();
	RooPlot *ctheta_frame = ws->var("ctheta")->frame();
	ws->data("data")->plotOn(ctheta_frame,RooFit::MarkerSize(0.2));
	ws->pdf("model")->plotOn(ctheta_frame, RooFit::ProjWData(RooArgSet(*ws->var("D"), *ws->var("et")), *ws->data("data")), RooFit::NumCPU(2));
	ctheta_frame->Draw();

	resPad->cd();
	RooHist *cthetaResidual = ctheta_frame->pullHist();
        RooPlot *cthetaResidualPlot  = ws->var("ctheta")->frame(RooFit::Title("Residual Distribution"));
	cthetaResidualPlot->addPlotable(cthetaResidual,"E");
	cthetaResidualPlot->Draw();


//	TCanvas canvasPhi("canvasPhi", "canvas phi", 400,600);
	canvas.cd(6);
	plotPad = new TPad("plotPadPhi", "Plot Pad", 0,0.2, 1,1);
	resPad = new TPad("resPadPhi", "Res Pad", 0,0, 1,0.2);
	plotPad->Draw();
	resPad->Draw();

	plotPad->cd();
	RooPlot *phi_frame = ws->var("phi")->frame();
	ws->data("data")->plotOn(phi_frame,RooFit::MarkerSize(0.2));
	ws->pdf("model")->plotOn(phi_frame, RooFit::ProjWData(RooArgSet(*ws->var("D"), *ws->var("et")), *ws->data("data")), RooFit::NumCPU(2));
	phi_frame->Draw();

	resPad->cd();
	RooHist *phiResidual = phi_frame->pullHist();
        RooPlot *phiResidualPlot  = ws->var("phi")->frame(RooFit::Title("Residual Distribution"));
	phiResidualPlot->addPlotable(phiResidual,"E");
	phiResidualPlot->Draw();

	canvas.SaveAs("t.png");
}
