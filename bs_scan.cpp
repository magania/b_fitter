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
  int nCPU = 1;
  double DG = 0.0;
  double BETA = 0.0;
  
  char option_char;
  while ( (option_char = getopt(argc,argv, "i:o:n:b:d:")) != EOF )
    switch (option_char)
      {  
         case 'i': chInFile = optarg; break;
         case 'o': chOutFile = optarg; break;
         case 'n': nCPU = atoi(optarg); break;
         case 'd': DG = atof(optarg); break;
         case 'b': BETA = atof(optarg); break;
         case '?': fprintf (stderr, 
                            "usage: %s [i<input file> o<output file>]\n", argv[0]);
      }

  cout << "DataSet = " << chInFile << endl;
  cout << "FitResult = " << chOutFile << endl;
  cout << "Delta Gamma = " << DG << endl;
  cout << "Beta = " << BETA << endl;

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

  cout << "Fitting ..." << endl;

  ws->var("beta")->setVal(BETA);
  ws->var("DG")->setVal(DG);
  ws->var("beta")->setConstant(kTRUE);
  ws->var("DG")->setConstant(kTRUE);
 

  RooFitResult *rfr = ws->pdf("model")->fitTo(*data ,RooFit::ConditionalObservables(*ws->var("d")), RooFit::NumCPU(nCPU), RooFit::Save(), RooFit::Verbose(kFALSE)); 
  rfr->Write("fitResultScan");

  ws->Write("rws");
  outFile.Close();
  inFile.Close();

  cout << endl << "Done." << endl;
}
	
