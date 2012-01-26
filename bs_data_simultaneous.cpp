/*
 *  Author: Ricardo Magana-Villalba
 *          magania@fnal.gov
 *
 *  Dec 2010
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
#include <TCut.h>

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
  const char* chInFile = "ws.root";
  const char* chOutFile = "ws_data.root";
  const char* chRootFile = "BDT20.root";
  const char* chCut = "";

  char option_char;
  while ( (option_char = getopt(argc,argv, "i:o:r:c:")) != EOF )
    switch (option_char)
      {
         case 'i': chInFile = optarg; break;
         case 'o': chOutFile = optarg; break;
         case 'r': chRootFile = optarg; break;
         case 'c': chCut = optarg; break;
         case '?': fprintf (stderr,
                            "usage: %s [i<input file> o<output file>]\n", argv[0]);
      }

  cout << "In  Ws = " << chInFile << endl;
  cout << "Out Ws = " << chOutFile << endl;
  cout << "Data From = " << chRootFile << endl;
  cout << "Extra Cut = " << chCut << endl;


   TFile* in_file = new TFile(chInFile);
   RooWorkspace *rws = (RooWorkspace*) in_file->Get("rws");

   TFile* tree_file = new TFile(chRootFile);
   TTree* tree = (TTree*) tree_file->Get("tree");

   TFile* out_file = new TFile(chOutFile, "RECREATE");

   RooArgSet allVars(*rws->var("m"),*rws->var("t"),*rws->var("et"),*rws->var("cpsi"),*rws->var("ctheta"),*rws->var("phi"),*rws->var("d"),*rws->cat("dilution"));
   RooDataSet* data = new RooDataSet("data","data",allVars);
   RooDataSet* dataBkg = new RooDataSet("dataBkg","dataBkg",allVars);

   //TCut* cut = new TCut("5.17<bs_mass && bs_mass<5.57 && bs_pdl>-0.044 && bs_pdl<0.3 ");
   TCut* cut = new TCut("5.17<bs_mass && bs_mass<5.57 && bs_epdl<0.025 && bs_pdl<0.4 && bs_pdl>-0.44");
   *cut += chCut;
   tree->Draw(">>entry_list", *cut, "entrylist");
   TEntryList* event_list = (TEntryList*) out_file->Get("entry_list");

   Double_t dM, dT, dEt, dCpsi, dCtheta, dPhi, dd;
   Int_t ddDefined;
   tree->SetBranchAddress("bs_mass", &dM);
   tree->SetBranchAddress("bs_pdl", &dT);
   tree->SetBranchAddress("bs_epdl", &dEt);
   tree->SetBranchAddress("bs_angle_cpsi", &dCpsi);
   tree->SetBranchAddress("bs_angle_ctheta", &dCtheta);
   tree->SetBranchAddress("bs_angle_phi", &dPhi);
   tree->SetBranchAddress("newtag_ost", &dd);
   tree->SetBranchAddress("newtag_ost_defined", &ddDefined);


   for (Long_t i=0; i<event_list->GetN(); i++){
     tree->GetEntry(event_list->GetEntry(i));

       *rws->var("m")=dM;
       *rws->var("t")=dT/0.0299792458;
       *rws->var("et")=dEt/0.0299792458;
       *rws->var("cpsi")=dCpsi;
       *rws->var("ctheta")=dCtheta;
       *rws->var("phi")=dPhi;

       *rws->var("d")=0;
       rws->cat("dilution")->setIndex(0);
       if ( ddDefined==1 ){
    	   rws->cat("dilution")->setIndex(1);
    	   *rws->var("d")=dd;
       }

       data->add(allVars);
       if (dM<5.29 || dM>5.44)
           dataBkg->add(allVars);
   }

   rws->import(*data);
   rws->import(*dataBkg);
   rws->Write("rws");
   out_file->Close();
   in_file->Close();
   tree_file->Close();

   cout << endl << "Done." << endl;
}
