/*
 *  Author: Ricardo Magana-Villalba
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
#include "RooCategory.h"
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

  RooRealVar m("m", "Mass (#mu#mu KK)", 0, 5.17, 5.57);
  RooRealVar t("t", "t", 0,-5, 20);
  RooRealVar et("et", "#sigma(t))", 0.01, 0, 3);
  RooRealVar cpsi("cpsi", "cos(#psi)", -1, 1);
  RooRealVar ctheta("ctheta", "cos(#theta)", -1 , 1);
  RooRealVar phi("phi", "#phi", -TMath::Pi(), TMath::Pi());
  RooRealVar d("d", "d", 0, -1, 1);
  RooCategory dilution("dilution","dilution");
  dilution.defineType("NoTag",0);
  dilution.defineType("Tag",1);

  ws->import(m);
  ws->import(t);
  ws->import(et);
  ws->import(cpsi);
  ws->import(ctheta);
  ws->import(phi);
  ws->import(d);
  ws->import(dilution);

  ws->factory("{DZero[0],D10to20[0.102],D20to35[0.184],D35to45[0.325],D45to60[0.486],D60to100[0.543]}");
  ws->factory("{DG[0.09,0,1],Dm[17.77,17.4,18.2],tau[1.5,1.0,2.0],beta[0.25,-1,1],delta_l[-2.93],delta_p[2.93],delta_s[0.2,-3,3],fs[0.1,0,1],scale[1.0,0.5,1.5]}");
  ws->factory("Gaussian::DmConstraint(Dm,17.77,0.17");

  ws->factory("{A02[0.6,0,1],A1[0.75,0,1]}");
//  ws->factory("{A02[0.6,0,1],Al2[0.3,0,1]}");
  ws->factory("expr:Al2('(1-@0)*@1',A02,A1)");
  ws->factory("expr:Ap2('1-@0-@1',A02,Al2)");

  ws->factory("expr:z('TMath::Cos(2.0*@0)*@1*@2/2.0',beta,DG,tau");
  ws->factory("expr:y('(1+@0)/(1-@0)',z");
  ws->factory("expr:ap('TMath::Sqrt(@0*@1/(@1+(1.0-@1)*@2))',Ap2,y,Ap2)");
  ws->factory("expr:al('TMath::Sqrt(@0*@1/(@1+(1.0-@1)*@2))',Al2,y,Ap2)");
  ws->factory("expr:a0('TMath::Sqrt(@0*@1/(@1+(1.0-@1)*@2))',A02,y,Ap2)");

  ws->factory("expr:Dilution('(@0/TMath::Abs(@0))*( 0.612548/(1 + TMath::Exp( (0.3421 - TMath::Abs(@0))/0.119108 ) ) - 0.612548/(1 + TMath::Exp(0.3421/0.119108) ) )',d)");
////  ws->factory("expr:D('0.0',d)");

  /* ************************************* Signal ********************************************** */
  /* Mass */
  ws->factory("Gaussian::massSignal(m,mu[5.2,5.4],sigmaM[0.03,0.02,0.04]");

  /* Time Angle */
  RooBsTimeAngle timeAngleNoTag("timeAngleNoTag", "timeAngleNoTag",
		*ws->var("t"), *ws->var("et"), *ws->var("cpsi"), *ws->var("ctheta"), *ws->var("phi"),
		*ws->var("DZero"),
		*ws->function("a0"), *ws->function("al"), *ws->function("ap"),
		*ws->var("DG"), *ws->var("Dm"), *ws->var("tau"), *ws->var("beta"),
		*ws->var("delta_l"), *ws->var("delta_p"), *ws->var("delta_s"), *ws->var("fs"), *ws->var("scale"));

  RooBsTimeAngle timeAngleTag("timeAngleTag", "timeAngleTag",
		*ws->var("t"), *ws->var("et"), *ws->var("cpsi"), *ws->var("ctheta"), *ws->var("phi"),
		*ws->function("Dilution"),
		*ws->function("a0"), *ws->function("al"), *ws->function("ap"),
		*ws->var("DG"), *ws->var("Dm"), *ws->var("tau"), *ws->var("beta"),
		*ws->var("delta_l"), *ws->var("delta_p"), *ws->var("delta_s"), *ws->var("fs"), *ws->var("scale"));

  ws->import(timeAngleNoTag);
  ws->import(timeAngleTag);

  /* Sigma(t) */
  ws->factory("GaussModel::etGaussianS(et,meanGaussEtS[0.0614,0,0.2],sigmaGaussEtS[0.0116,0,0.2])");
  ws->factory("Decay::errorSignal(et,tauEtS[0.0481,0,0.2],etGaussianS,SingleSided]");

////  ws->factory("Landau::errorSignal(et,meanErrSignal[0.0741194,0,1],sigmaErrSignal[0.0175981,0,1])");

  /* Full Signal PDF */
  ws->factory("PROD::signalTimeAngleNoTag(timeAngleNoTag|et,errorSignal");
  ws->factory("PROD::signalTimeAngleTag(timeAngleTag|et,errorSignal");

  ws->factory("PROD::signalNoTag(massSignal,signalTimeAngleNoTag,DmConstraint)");
  ws->factory("PROD::signalTag(massSignal,signalTimeAngleTag,DmConstraint)");

  /* ************************************* Background ********************************************** */
  /* Mass */
  ws->factory("Polynomial::massBkgNP(m,{slopeNP[-0.1238,-0.1795,1]})");
  ws->factory("Polynomial::massBkgPR(m,{slopePR[-0.1283,-0.1795,1]})");

  /* Time */
  ws->factory("GaussModel::resolution(t,0,scale,et)");
////  ws->factory("GaussModel::resolutionWide(t,0,scaleWide[1,0,3],et)");
////  ws->factory("GaussModel::resolutionNarrow(t,0,scaleNarrow[1,0,1],et)");
////  ws->factory("AddModel::resolution({resolutionWide,resolutionNarrow},xResWide[0.9,0,1])");
  ws->factory("Decay::negativeDecay(t,tauNeg[0.387,0.01,1],resolution,Flipped)");
  ws->factory("Decay::positiveDecay(t,tauPos[0.151,0.01,1],resolution,SingleSided)");
  ws->factory("Decay::positiveLongDecay(t,tauLngPos[0.761,0.1,3],resolution,SingleSided)");
  ws->factory("Decay::positiveLongLongDecay(t,tauLngLngPos[0.761,0.1,5],resolution,SingleSided)");

  ws->factory("RSUM::tBkgNP(xn[0.0488,0,1]*negativeDecay,xp[0.767,0,1]*positiveDecay,xpp[0.1,0,1]*positiveLongDecay,positiveLongLongDecay");

  /* Angle */
  ws->factory("{pre01m1[0],pre01p0[0],pre01p1[0.006,-1,1],pre02m2[0.0003,-1,1],pre02m1[0],pre02p0[0.120,-1,1],pre02p1[0.003,-1,1],pre02p2[-0.200,-1,1]}");
  ws->factory("{pre20p0[0.057,-1,1],pre21m1[0],pre21p0[0],pre21p1[0],pre22m2[0],pre22m1[0],pre22p0[0],pre22p1[0],pre22p2[-0.013,-1,1]}");

  RooAngle anglePR("anglePR","anglePR", *ws->var("cpsi"), *ws->var("ctheta"), *ws->var("phi"),
	               *ws->var("pre01m1"),
	               *ws->var("pre01p0"),
	               *ws->var("pre01p1"),
	               *ws->var("pre02m2"),
	               *ws->var("pre02m1"),
	               *ws->var("pre02p0"),
	               *ws->var("pre02p1"),
	               *ws->var("pre02p2"),
	               *ws->var("pre20p0"),
	               *ws->var("pre21m1"),
	               *ws->var("pre21p0"),
	               *ws->var("pre21p1"),
	               *ws->var("pre22m2"),
	               *ws->var("pre22m1"),
	               *ws->var("pre22p0"),
	               *ws->var("pre22p1"),
	               *ws->var("pre22p2")
		       );

  ws->factory("{npe01m1[0],npe01p0[0],npe01p1[0],npe02m2[0],npe02m1[-0.006,-1,1],npe02p0[0.066,-1,1],npe02p1[0],npe02p2[-0.139,-1,1]}");
  ws->factory("{npe20p0[0],npe21m1[0],npe21p0[0],npe21p1[0],npe22m2[0],npe22m1[0],npe22p0[0],npe22p1[0.031,-1,1],npe22p2[0]}");

  RooAngle angleNP("angleNP","angleNP", *ws->var("cpsi"), *ws->var("ctheta"), *ws->var("phi"),
	               *ws->var("npe01m1"),
	               *ws->var("npe01p0"),
	               *ws->var("npe01p1"),
	               *ws->var("npe02m2"),
	               *ws->var("npe02m1"),
	               *ws->var("npe02p0"),
	               *ws->var("npe02p1"),
	               *ws->var("npe02p2"),
	               *ws->var("npe20p0"),
	               *ws->var("npe21m1"),
	               *ws->var("npe21p0"),
	               *ws->var("npe21p1"),
	               *ws->var("npe22m2"),
	               *ws->var("npe22m1"),
	               *ws->var("npe22p0"),
	               *ws->var("npe22p1"),
	               *ws->var("npe22p2")
		       );

  ws->import(anglePR);
  ws->import(angleNP);

  /* Sigma(t) */
//  ws->factory("Landau::errBkgPR(et,meanEtPR[0.07,0,1],sigmaEtPR[0.011,0,1])");
//  ws->factory("Landau::errBkgNP(et,meanEtNP[0.08,0,1],sigmaEtNP[0.011,0,1])");


/*
  ws->factory("{EtAPR[2.5,0.5,10],EtBPR[0.022,0.001,0.09],EtSPR[0.023,0,0.2]}");
  ws->factory("{EtANP[0.7,0.5,10],EtBNP[0.032,0.001,0.09],EtSNP[0.0189,0,0.2]}");
  RooErrPdf errBkgPR("errBkgPR", "errBkgPR", *ws->var("et"), *ws->var("EtAPR"), *ws->var("EtBPR"),*ws->var("EtSPR"));
  RooErrPdf errBkgNP("errBkgNP", "errBkgNP", *ws->var("et"), *ws->var("EtANP"), *ws->var("EtBNP"),*ws->var("EtSNP"));
  ws->import(errBkgPR);
  ws->import(errBkgNP);
*/

  ws->factory("GaussModel::etGaussianPR(et,meanGaussEtPR[0.0614,0,0.2],sigmaGaussEtPR[0.0116,0,0.2])");
  ws->factory("Decay::errBkgPR(et,tauEtPR[0.0481,0,0.2],etGaussianPR,SingleSided]");

  ws->factory("GaussModel::etGaussianNP(et,meanGaussEtNP[0.0614,0,0.2],sigmaGaussEtNP[0.0116,0,0.2])");
  ws->factory("Decay::errBkgNP(et,tauEtNP[0.0481,0,0.2],etGaussianNP,SingleSided]");

  //ws->factory("BifurGauss::errBkgPR(et,meanEtPR[0.07,0,1],sigmaLEtPR[0.011,0,1],sigmaREtPR[0.011,0,1])");

//  ws->factory("SUM::errorBkg(xe[0.61,0,1]*errBkg,err2Bkg)");

  //ws->factory("EXPR::errBkgPR('(@0<@3)*TMath::Exp(-0.5*(@0-@1)*(@0-@1)/(@2*@2)) + (@0>=@3)*( TMath::Exp(-(@0-@3)/@4)*TMath::Exp(-0.5*(@3-@1)*(@3-@1)/(@2*@2)) )',{et,EtMPR[0.07,0.05,0.15],EtSPR[0.013,0.005,0.020],EtM2PR[0.07,0.05,0.15],EtBPR[0.0475,0.01,0.1]})");
  //ws->factory("Landau::errBkgNP(et,meanEtNP[0.08,0,1],sigmaEtNP[0.011,0,1])");

/*               Promt and Non-Prompt                       */
   ws->factory("PROD::timeBkgNP(tBkgNP|et,errBkgNP)");
   ws->factory("PROD::timeBkgPR(resolution|et,errBkgPR)");

   ws->factory("PROD::Prompt(massBkgPR,timeBkgPR,anglePR)");
   ws->factory("PROD::NonPrompt(massBkgNP,timeBkgNP,angleNP)");

  ////ws->factory("PROD::Prompt(massBkgPR,resolution,anglePR)");
  ////ws->factory("PROD::NonPrompt(massBkgNP,tBkgNP,angleNP)");

  /* Full Background PDF */
  ws->factory("SUM::background(xprompt[0.853,0,1]*Prompt,NonPrompt)");

  /* Full PDF */
  ws->factory("SUM::modelNoTag(xsNoTag[0.0436,0,1]*signalNoTag,background)");
  ws->factory("SUM::modelTag(xsTag[0.1,0,1]*signalTag,background)");

  ws->factory("SIMUL:model(dilution,NoTag=modelNoTag,Tag=modelTag)");

  ws->Write("rws");
  outFile.Close();

  cout << endl << "Done." << endl;
}
