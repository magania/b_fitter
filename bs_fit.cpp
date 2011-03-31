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
  const char* chInFile = "ws_data.root";
  const char* chOutFile = "ws_data_fit.root";
  int nCPU = 4;
  bool justBkg = false;
  
  char option_char;
  while ( (option_char = getopt(argc,argv, "i:o:n:b")) != EOF )
    switch (option_char)
      {  
         case 'i': chInFile = optarg; break;
         case 'o': chOutFile = optarg; break;
         case 'n': nCPU = atoi(optarg); break;
         case 'b': justBkg = true; break;
         case '?': fprintf (stderr, 
                            "usage: %s [i<input file> o<output file>]\n", argv[0]);
      }

  cout << "DataSet = " << chInFile << endl;
  cout << "FitResult = " << chOutFile << endl;

  TFile inFile(chInFile);
  RooWorkspace* ws = (RooWorkspace*) inFile.Get("rws");
  cout << ws->data("data") << endl;
  TFile outFile(chOutFile,"RECREATE");


  RooAbsData *data = ws->data("dataGen");
  if (!data)
  {
     data = ws->data("data");
     cout << "Using data .. " << endl;
  }


  RooAbsData *dataBkg = ws->data("dataGenBkg");
  if (!dataBkg)
  {
     dataBkg = ws->data("dataBkg");
     cout << "Using dataBkg .. " << endl;
  }

  /* ************************************* Fit ********************************************** */

  cout << "Fitting Sigma(et) from sidebands" << endl;

/*  ws->var("meanEtNP")->setConstant(kTRUE);
  ws->var("sigmaEtNP")->setConstant(kTRUE);
  ws->var("xprompt")->setConstant(kTRUE);
  ws->var("xprompt")->setVal(1.0);
*/
  //ws->factory("SUM::tmpErrBkg(xprompt*errBkgPR,errBkgNP)");
  //ws->pdf("tmpErrBkg")->fitTo(*dataBkg, RooFit::NumCPU(nCPU));

  cout << "Fitting Background form sidebands" << endl;

  ws->factory("SUM::backgroundT(xprompt*timeBkgPR,timeBkgNP)");

  ws->pdf("backgroundT")->fitTo(*dataBkg,RooFit::NumCPU(nCPU),RooFit::Hesse(kFALSE),RooFit::Minos(kFALSE));

  ws->var("meanGaussEtNP")->setConstant(kTRUE);
  ws->var("meanGaussEtPR")->setConstant(kTRUE);
  ws->var("sigmaGaussEtNP")->setConstant(kTRUE);
  ws->var("sigmaGaussEtPR")->setConstant(kTRUE);
  ws->var("tauEtPR")->setConstant(kTRUE);
  ws->var("tauEtNP")->setConstant(kTRUE);

  ws->var("tauLngLngPos")->setConstant(kTRUE);
  ws->var("tauLngPos")->setConstant(kTRUE);
  ws->var("tauPos")->setConstant(kTRUE);
  ws->var("tauNeg")->setConstant(kTRUE);
  ws->var("xn")->setConstant(kTRUE);
  ws->var("xp")->setConstant(kTRUE);
  ws->var("xpp")->setConstant(kTRUE);
  ws->var("xprompt")->setConstant(kTRUE);
  ws->var("scale")->setConstant(kTRUE);

  ws->factory("Decay::signalTimeExp(t,tau,resolution,SingleSided)");
  ws->factory("PROD::signalTimeExpEt(signalTimeExp|et,errorSignal");
  ws->factory("PROD::signalMassTime(massSignal,signalTimeExpEt)");

  ws->factory("PROD::PromptTM(massBkgPR,timeBkgPR)");
  ws->factory("PROD::NonPromptTM(massBkgNP,timeBkgNP)");
  ws->factory("SUM::backgroundTM(xprompt*PromptTM,NonPromptTM)");

  ws->factory("SUM::backgroundMassTime(xs*signalMassTime,backgroundTM)");
  ws->factory("SUM::massTimeModel(xs*signalMassTime,backgroundMassTime)");

  ws->pdf("massTimeModel")->fitTo(*data,RooFit::NumCPU(nCPU),RooFit::Hesse(kFALSE),RooFit::Minos(kFALSE));
  ws->var("meanGaussEtS")->setConstant(kTRUE);
  ws->var("sigmaGaussEtS")->setConstant(kTRUE);
  ws->var("tauEtS")->setConstant(kTRUE);

  ws->var("mu")->setConstant(kTRUE);
  ws->var("sigmaM")->setConstant(kTRUE);
  ws->var("slopePR")->setConstant(kTRUE);
  ws->var("slopeNP")->setConstant(kTRUE);

  ////ws->var("xsNoTag")->setConstant(kTRUE);


  ws->factory("PROD::PromptNM(timeBkgPR,anglePR)");
  ws->factory("PROD::NonPromptNM(timeBkgNP,angleNP)");
  ws->factory("SUM::backgroundNM(xprompt*PromptNM,NonPromptNM)");

  ws->pdf("backgroundNM")->fitTo(*dataBkg,RooFit::NumCPU(nCPU),RooFit::Hesse(kFALSE),RooFit::Minos(kFALSE));
  ws->var("npe02m1")->setConstant(kTRUE);
  ws->var("npe02p0")->setConstant(kTRUE);
  ws->var("npe02p2")->setConstant(kTRUE);
  ws->var("npe22p1")->setConstant(kTRUE);

  ws->var("pre01p1")->setConstant(kTRUE);
  ws->var("pre02m2")->setConstant(kTRUE);
  ws->var("pre02p0")->setConstant(kTRUE);
  ws->var("pre02p1")->setConstant(kTRUE);
  ws->var("pre02p2")->setConstant(kTRUE);
  ws->var("pre20p0")->setConstant(kTRUE);
  ws->var("pre22p2")->setConstant(kTRUE);

  if (!justBkg)
  {
  cout << "Fitting signal only" << endl;
  ws->pdf("model")->fitTo(*data ,RooFit::NumCPU(nCPU),RooFit::Hesse(kFALSE),RooFit::Minos(kFALSE));

  ws->var("tauLngLngPos")->setConstant(kFALSE);
  ws->var("tauLngPos")->setConstant(kFALSE);
  ws->var("tauPos")->setConstant(kFALSE);
  ws->var("tauNeg")->setConstant(kFALSE);
  ws->var("xn")->setConstant(kFALSE);
  ws->var("xp")->setConstant(kFALSE);
  ws->var("xpp")->setConstant(kFALSE);
  ws->var("xprompt")->setConstant(kFALSE);
  ws->var("scale")->setConstant(kFALSE);

  ws->var("mu")->setConstant(kFALSE);
  ws->var("sigmaM")->setConstant(kFALSE);
  ws->var("slopePR")->setConstant(kFALSE);
  ws->var("slopeNP")->setConstant(kFALSE);

  ws->var("npe02m1")->setConstant(kFALSE);
  ws->var("npe02p0")->setConstant(kFALSE);
  ws->var("npe02p2")->setConstant(kFALSE);
  ws->var("npe22p1")->setConstant(kFALSE);

  ws->var("pre01p1")->setConstant(kFALSE);
  ws->var("pre02m2")->setConstant(kFALSE);
  ws->var("pre02p0")->setConstant(kFALSE);
  ws->var("pre02p1")->setConstant(kFALSE);
  ws->var("pre02p2")->setConstant(kFALSE);
  ws->var("pre20p0")->setConstant(kFALSE);
  ws->var("pre22p2")->setConstant(kFALSE);

  cout << "FullFit" << endl;
  RooFitResult *rfr = ws->pdf("model")->fitTo(*data, RooFit::NumCPU(nCPU), RooFit::Save(), RooFit::Verbose(kFALSE));
  rfr->Write("fitResult");
  }

  ws->Write("rws");
  outFile.Close();
  inFile.Close();

  cout << endl << "Done." << endl;
}
	
