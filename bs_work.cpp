/*
 *  Author: Ricardo Magaña-Villalba
 *          magania@fnal.gov
 *
 *  September 2010
 */

#include <iostream>
#include <stdio.h>
#include <unistd.h>

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

int main (int argc, char **argv)
{
  const char* chOutFile = "ws.root";
  
  char option_char;
  while ( (option_char = getopt(argc,argv, "i:o:")) != EOF )
    switch (option_char)
      {  
         case 'o': chOutFile = optarg; break;
         case '?': fprintf (stderr, 
                            "usage: %s [o<output file>]\n", argv[0]);
      }

  cout << "Out File = " << chOutFile << endl;

  TFile outFile(chOutFile,"RECREATE");
  RooWorkspace* ws = new RooWorkspace();

  RooRealVar m("m", "Mass (#mu#mu KK)", 0, 5, 5.8);
  RooRealVar t("t", "t", 0,-5, 20);
  RooRealVar et("et", "#sigma(t))", 0.01, 0, 3);
  RooRealVar cpsi("cpsi", "cos(#psi)", -1, 1);
  RooRealVar ctheta("ctheta", "cos(#theta)", -1 , 1);
  RooRealVar phi("phi", "#phi", -TMath::Pi(), TMath::Pi());
  RooRealVar d("d", "d", -1, 1);
  ws->import(m);
  ws->import(t);
  ws->import(et);
  ws->import(cpsi);
  ws->import(ctheta);
  ws->import(phi);
  ws->import(d);

  ws->factory("{DG[0.09,0,1],Dm[17.77],tau[1.5,1.0,2.0],beta[0.25,-1,1],delta_l[-2.93],delta_p[2.93],delta_s[0.2,-3,3],fs[0.1,0,1],scale[1.0,0.5,1.5]}");

  ws->factory("{A02[0.6,0,1],A1[0.75,0,1]}");
//  ws->factory("{A02[0.6,0,1],Al2[0.3,0,1]}");
  ws->factory("expr:Al2('(1-@0)*@1',A02,A1)");
  ws->factory("expr:Ap2('TMath::Sqrt(1-@0-@1)',A02,Al2)");

  ws->factory("expr:z('TMath::Cos(2.0*@0)*@1*@2/2.0',beta,DG,tau");
  ws->factory("expr:y('(1+@0)/(1-@0)',z");
  ws->factory("expr:ap('TMath::Sqrt(@0*@1/(@1+(1.0-@1)*@2))',Ap2,y,Ap2)");
  ws->factory("expr:al('TMath::Sqrt(@0*@1/(@1+(1.0-@1)*@2))',Al2,y,Ap2)");
  ws->factory("expr:a0('TMath::Sqrt(@0*@1/(@1+(1.0-@1)*@2))',A02,y,Ap2)");

  ws->factory("expr:D('(@0/TMath::Abs(@0))*( 0.6625/(1 + TMath::Exp( (0.3104 - TMath::Abs(@0))/0.1184 ) ) - 0.6625/(1 + TMath::Exp(0.3104/0.1184) ) )',d)");
//  ws->factory("expr:D('0.0',d)");

  /* ************************************* Signal ********************************************** */
  /* Mass */
  ws->factory("Gaussian::massSignal(m,mu[5.2,5.4],sigmaM[0.03,0.02,0.04]");

  /* Time Angle */
  RooBsTimeAngle timeAngle("timeAngle", "timeAngle",
		*ws->var("t"),
		*ws->var("et"),
		*ws->var("cpsi"),
		*ws->var("ctheta"),
		*ws->var("phi"),
		*ws->function("D"),
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
  ws->import(timeAngle);

  /* Sigma(t) */
  ws->factory("Landau::errorSignal(et,meanErrSignal[0.0741194,0,1],sigmaErrSignal[0.0175981,0,1])");

  /* Full Signal PDF */
  ws->factory("PROD::signalTimeAngle(timeAngle|et,errorSignal");
  ws->factory("PROD::signal(massSignal,signalTimeAngle)");

  /* ************************************* Background ********************************************** */
  /* Mass */
  ws->factory("Polynomial::massBkg(m,{slope[-0.125,-1,-0.1]})");

  /* Time */
  ws->factory("GaussModel::resolution(t,0,scale,et)");
  ws->factory("Decay::negativeDecay(t,tauNeg[0.36,0.01,1],resolution,Flipped)");
  ws->factory("Decay::positiveDecay(t,tauPos[0.41,0.1,1],resolution,SingleSided)");
  ws->factory("Decay::positiveLongDecay(t,tauPosPos[1.9,0.5,3],resolution,SingleSided)");

  ws->factory("RSUM::tBkg(xr[0.73,0,1]*resolution,xn[0.15,0,1]*negativeDecay,xp[0.87,0,1]*positiveDecay,positiveLongDecay");

  /* Angle */
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

  /* Sigma(t) */
  //ws->factory("Landau::errLBkg(et,mean[0.072,0,1],sigma[0.011,0,1])");
  //ws->factory("Landau::errRBkg(et,mean2[0.08,0,1],sigma2[0.011,0,1])");
  ws->factory("{errA[3.8,0.5,10],errB[0.011,0.01,0.2],errS[0.035,0,0.05]}");
  ws->factory("{errA2[2.24,0.5,10],errB2[0.032,0.01,0.2],errS2[0.024,0,0.05]}");
  RooErrPdf errBkg("errBkg", "errBkg", *ws->var("et"), *ws->var("errA"), *ws->var("errB"),*ws->var("errS"));
  RooErrPdf err2Bkg("err2Bkg", "err2Bkg", *ws->var("et"), *ws->var("errA2"), *ws->var("errB2"),*ws->var("errS2"));
  ws->import(errBkg);
  ws->import(err2Bkg);

  ws->factory("SUM::errorBkg(xe[0.61,0,1]*errBkg,err2Bkg)");

  /* Full Background PDF */
  ws->factory("PROD::timeBkg(tBkg|et,errorBkg)");
  ws->factory("PROD::background(massBkg,timeBkg,angle)");

  /* Full PDF */
  ws->factory("SUM::model(xs[0.0436,0,1]*signal,background)");

  ws->Write("rws");
  outFile.Close();
}


