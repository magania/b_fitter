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
  const char* chInFile = "ws_gen.root";
  const char* chOutFile = "ws_gen_fit.root";
  
  char option_char;
  while ( (option_char = getopt(argc,argv, "i:o:")) != EOF )
    switch (option_char)
      {  
         case 'i': chInFile = optarg; break;
         case 'o': chOutFile = optarg; break;
         case '?': fprintf (stderr, 
                            "usage: %s [i<input file> o<output file>]\n", argv[0]);
      }

  cout << "DataSet = " << chInFile << endl;
  cout << "FitResult = " << chOutFile << endl;

  TFile inFile(chInFile);
  RooWorkspace* ws = (RooWorkspace*) inFile.Get("rws");
  cout << ws->data("data") << endl;
  TFile outFile(chOutFile,"RECREATE");

  /* ************************************* Fit ********************************************** */

  cout << "Fitting Sigma(et) from sidebands" << endl;
  ws->pdf("errorBkg")->fitTo(*ws->data("dataBkg"), RooFit::NumCPU(2));
  ws->var("errA")->setConstant(kTRUE);
  ws->var("errB")->setConstant(kTRUE);
  ws->var("errS")->setConstant(kTRUE);
  ws->var("xe")->setConstant(kTRUE);
  ws->var("errA2")->setConstant(kTRUE);
  ws->var("errB2")->setConstant(kTRUE);
  ws->var("errS2")->setConstant(kTRUE);

  cout << "Fitting Background form sidebands" << endl;
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

  cout << "Fitting signal only" << endl;
  ws->pdf("model")->fitTo(*ws->data("data"), RooFit::ConditionalObservables(*ws->var("d")) ,RooFit::NumCPU(2));
  ws->var("tauNeg")->setConstant(kFALSE);
  ws->var("tauPos")->setConstant(kFALSE);
  ws->var("tauPosPos")->setConstant(kFALSE);
  ws->var("xn")->setConstant(kFALSE);
  ws->var("xp")->setConstant(kFALSE);
  ws->var("xr")->setConstant(kFALSE);
  ws->var("scale")->setConstant(kFALSE);
  ws->var("e02p0")->setConstant(kFALSE);
  ws->var("e02p2")->setConstant(kFALSE);
  ws->var("e20p0")->setConstant(kFALSE);
  ws->var("slope")->setConstant(kFALSE);

  cout << "FullFit" << endl;
  ws->pdf("model")->fitTo(*ws->data("data"), RooFit::ConditionalObservables(*ws->var("d")) ,RooFit::NumCPU(2));

  ws->Write("rws");
  outFile.Close();
  inFile.Close();		
}
	
