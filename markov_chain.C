#include "TMessage.h"
#include "RooFit.h"
#include "RooWorkspace.h"
#include "TTree.h"
#include "TMatrixDSym.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TSocket.h"
#include "TMath.h"
#include "TVectorD.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "TFile.h"
#include "RooFitResult.h"
#include "TServerSocket.h"
#include "TMatrixTSym.h"
#include "TDecompChol.h"
#include "RooDataSet.h"
#include "TMonitor.h"
#include "TRandom3.h"


//static int NSAMPLE = 10000;
//static const char* fitWorkspaceFile = "ws_BDT20_fit.root";
static const bool DEBUG = false;

class Cholesky : public TDecompChol{
        public:
        Cholesky(TMatrixDSym m): TDecompChol(m){
                Decompose();
        }

        double operator()(int i, int j){
                return fU(j,i);
        }
};

void markov_chain(const char* fitWorkspaceFile, int NSAMPLE, int numCPU, double beta_start, double dg_start) {
   gSystem->Load("lib/libBFitter.so");

   TFile *fitWorkspace = new TFile(fitWorkspaceFile);
   RooWorkspace *ws = (RooWorkspace *) fitWorkspace->Get("rws");
   RooDataSet *allData = (RooDataSet *) ws->data("data");
   RooFitResult *fitResult = (RooFitResult*) fitWorkspace->Get("fitResult");
   RooArgList fitParams = fitResult->floatParsFinal();
   Int_t nParams = fitParams.getSize();
   TMatrixDSym *sigma = (TMatrixDSym*) fitResult->covarianceMatrix().Clone();

   if (nParams > 50){
       cout << "EE: Increase x and mu arrays size." << endl;
       exit(1);
   }

   ws->var("beta")->setVal(beta_start);
   ws->var("DG")->setVal(dg_start);
   RooAbsReal *nll = ws->pdf("model")->createNLL(*allData, RooFit::Verbose(false), RooFit::NumCPU(numCPU), RooFit::ConditionalObservables(*ws->var("d")), RooFit::ExternalConstraints(*ws->pdf("DmConstraint")) );

/*
   double nll_mmm = nll->getVal();
   double dm=17.77-0.12;
   for (int i=0; i<100; i++){
        ws->var("Dm")->setVal(dm);
//        ws->var("Dm")->Print();
//        nll->Print();
	cout << nll->getVal() - nll_mmm << "  ";
        dm += 0.12/50;
   }
	
   return;
*/

   TObjArray names(nParams);
   for (int i=0; i<nParams; i++)
   	names.Add(new TObjString(fitParams[i].GetName()));

//   (*sigma) *= 1.0/nParams;
   (*sigma) *= 1.0/16.;
   Cholesky R(*sigma);

   if (DEBUG){
	   sigma->Print();
   }

   /* *** Create Output TTree ** */
   TVectorD mu(nParams);
   TVectorD mup(nParams);

   TFile *out_file = new TFile("mcmc.root","RECREATE");
   TTree *tree = new TTree("mcmc","mcmc");
   double *pmu = mu.GetMatrixArray();
   for (int i=0; i<nParams; i++){
//      RooRealVar *var = (RooRealVar*) &(fitParams[i]);
      RooRealVar *var = (RooRealVar*) ws->var(fitParams[i].GetName());
      pmu[i] =  var->getVal();
      TString data_type = var->GetName();
      data_type += "/D";
      tree->Branch(var->GetName(),&pmu[i],data_type.Data());
   }
   double old_nll = fitResult->minNll();
   tree->Branch("L",&old_nll,"L/D");
   int n_found=0;
   tree->Branch("n",&n_found,"n/I");
   tree->Fill();

   fitWorkspace->Close();

   TRandom3 *aleatorio = new TRandom3(0);
   double g[50], x[50];
   RooRealVar *vars[50];
   for (int i=0; i<nParams; i++)
           vars[i] = ws->var( ((TObjString*)names.At(i))->GetString().Data() );
   while (n_found<NSAMPLE) {
	   /* Generate sample of multivariate gaussian  g[i] = Gi(0,1)
	    * x = R*g    R = Low-Diagonal Cholesky Matrix SIGMA = R*R^T
	    */
           for (int i=0; i<nParams; i++){
                   g[i] = aleatorio->Gaus(0.0,1.0);
                   x[i] = 0.0;
           }

           for (int i=0; i<nParams; i++)
                   for (int j=0; j<=i; j++)
                           x[i] += R(i,j)*g[j];

          for (int i=0; i<nParams; i++)
                   mup(i) = mu(i) + x[i];
          
          bool outOfRange = false;
	  for (int i=0; i<nParams; i++)
		if ( mup[i] < vars[i]->getMin() || mup[i] > vars[i]->getMax() ) {
			outOfRange = true;
//			cout << "Out Of Range" << endl;
			break;
 	  	}
	  if (outOfRange)
		continue;
          
	  if (DEBUG){
		cout << "G: "<< endl;
		for (int i=0; i<nParams; i++)
			cout << g[i] << "\t";
		cout << endl << endl;
		cout << "X: "<< endl;
		for (int i=0; i<nParams; i++)
			cout << x[i] << "\t";
		cout << endl << endl;
		cout << "MU: "<< endl;
		for (int i=0; i<nParams; i++)
			cout << mu(i)<< "\t";
		cout << endl;
		cout << "MUP: "<< endl;
		for (int i=0; i<nParams; i++)
			cout << mup(i)<< "\t";
		cout << endl;
	  }

           /* Metropolis-Hasting */
           for (int i=0; i<nParams; i++)
	       ws->var( ((TObjString*)names.At(i))->GetString().Data() )->setVal(mup(i));
	   double new_nll = nll->getVal();
	
	   if (DEBUG) {
              for (int i=0; i<nParams; i++)
	          ws->var( ((TObjString*)names.At(i))->GetString().Data() )->Print();
              cout << " new_nll = " << new_nll << endl;
	   }

           double alpha = TMath::Exp(-new_nll + old_nll);
	   if (DEBUG)
                   cout << "alpha =" << alpha << endl;

           if (aleatorio->Uniform()<alpha){
                n_found++;
		old_nll = new_nll;
                for (int i=0; i<nParams; i++)
                    mu(i) = mup(i);
                tree->Fill();
                if (n_found%1000 == 1)
		    cout << "N FOUND = " << n_found << endl;
	   }
    }
		
   tree->Write();
   out_file->Close();

}
