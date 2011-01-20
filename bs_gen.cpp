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

int main (int argc, char **argv)
{
  const char* chInFile = "ws.root";
  const char* chOutFile = "ws_gen.root";
  int numSignal = 10000;
  int numBkg = 100000;

  char option_char;
  while ( (option_char = getopt(argc,argv, "i:o:s:b:")) != EOF )
    switch (option_char)
      {
         case 'i': chInFile = optarg; break;
         case 'o': chOutFile = optarg; break;
         case 's': numSignal = atoi(optarg); break;
         case 'b': numBkg = atoi(optarg); break;
         case '?': fprintf (stderr,
                            "usage: %s [i<input file> o<output file>]\n", argv[0]);
      }

  cout << "In File = " << chInFile << endl;
  cout << "Out File = " << chOutFile << endl;
  cout << "Signal Events = " << numSignal << endl;
  cout << "Bkg Events = " << numBkg << endl;

  TFile inFile(chInFile);
  RooWorkspace* ws = (RooWorkspace*) inFile.Get("rws");
  TFile outFile(chOutFile,"RECREATE");

  *ws->var("xs") = numSignal/(numSignal+numBkg);
//  int numSignal = numEvents * ws->var("xs")->getVal();
//  int numBkg = numEvents - numSignal;

  ws->factory("Gaussian::dilutionGauss(d,0,0.276)");
  ws->factory("SUM::dSignalPDF(xds[0.109]*dilutionGauss,TruthModel(d))");
  ws->factory("SUM::dBkgPDF(xdb[0.109]*dilutionGauss,TruthModel(d))");

  RooDataSet* dSignalData = ws->pdf("dSignalPDF")->generate(RooArgSet(*ws->var("d")),numSignal);
  RooDataSet *dataSignal = ws->pdf("signal")->generate(RooArgSet(*ws->var("m"),*ws->var("t"),*ws->var("et"),*ws->var("cpsi"),*ws->var("ctheta"),*ws->var("phi")), RooFit::ProtoData(*dSignalData));
 
  ws->factory("Decay::negativeDecay2(t,tauNeg,resolution,Flipped)");
  ws->factory("Decay::positiveDecay2(t,tauPos,resolution,SingleSided)");
  ws->factory("Decay::positiveLongDecay2(t,tauPosPos,resolution,SingleSided)");

  ws->factory("RSUM::tBkg2(xr*resolution,xn*negativeDecay2,xp*positiveDecay2,positiveLongDecay2");

  ws->factory("PROD::timeBkg2(tBkg2|et,errorBkg)");
  ws->factory("PROD::background2(massBkg,timeBkg2,angle,dBkgPDF)");

  RooDataSet* dBkgData = ws->pdf("dBkgPDF")->generate(RooArgSet(*ws->var("d")),numBkg);
  RooDataSet* dataBkg = ws->pdf("background2")->generate(RooArgSet(*ws->var("m"),*ws->var("t"),*ws->var("et"),*ws->var("cpsi"),*ws->var("ctheta"),*ws->var("phi")), numBkg);

  dataBkg->merge(dBkgData);
  dataBkg->SetName("dataBkg");
  ws->import(*dataBkg,RooFit::Rename("dataGenBkg"));

  dataSignal->append(*dataBkg);
  dataSignal->SetName("data");
  ws->import(*dataSignal,RooFit::Rename("dataGen"));

  //RooFitResult *fit_result = ws->pdf("model")->fitTo(*ws->data("data"), RooFit::Save(kTRUE), RooFit::ConditionalObservables(*ws->var("d")), RooFit::NumCPU(2), RooFit::PrintLevel(3));

/*
        gROOT->SetStyle("Plain");

        TCanvas canvas("canvas", "canvas", 800,400);
        canvas.Divide(2);

        canvas.cd(1);
        RooPlot *t_frame = ws->var("t")->frame();
        ws->data("data")->plotOn(t_frame, RooFit::MarkerSize(0.3));
        gPad->SetLogy(1);
        t_frame->Draw();

        canvas.cd(2);
        RooPlot *et_frame = ws->var("et")->frame();
        ws->data("data")->plotOn(et_frame,RooFit::MarkerSize(0.2));
        ws->pdf("errorSignal")->plotOn(et_frame);
        gPad->SetLogy(1);
        et_frame->Draw();

        canvas.SaveAs("t.png"); 


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

       canvas.SaveAs("t.png");

*/

  ws->Write("rws");
  outFile.Close();
  inFile.Close();
}
