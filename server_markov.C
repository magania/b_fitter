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


static const int CHILDREN = 64;
static int NSAMPLE = 100000;
static const char* fitWorkspaceFile = "ws_BDT20_fit.root";
static const bool DEBUG = false;

static bool EVOL = false;
static unsigned int kINIT  = 20000;
static unsigned int kPOINT = 20001;

class Cholesky : public TDecompChol{
        public:
        Cholesky(TMatrixDSym m): TDecompChol(m){
                Decompose();
        }

        double operator()(int i, int j){
                return fU(j,i);
        }
};

void server_markov() {
   TMessage::EnableSchemaEvolutionForAll(EVOL);

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


   RooAbsReal *nll;
   if (DEBUG)
       nll = ws->pdf("model")->createNLL(*allData, RooFit::Verbose(false), RooFit::NumCPU(4), RooFit::ConditionalObservables(*ws->var("d")));

   TObjArray names(nParams);
   for (int i=0; i<nParams; i++)
   	names.Add(new TObjString(fitParams[i].GetName()));

   /* create datasets */
   RooDataSet *data[CHILDREN];
   for (int i=0; i<CHILDREN; i++)
	data[i] = (RooDataSet*) allData->emptyClone();

   for (Int_t i=0; i<allData->numEntries(); i++)
        data[i%CHILDREN]->add(*allData->get(i));

   if (DEBUG){
	 cout << "DATASETS:" << endl;
         allData->Print();
	 for (int i=0; i<CHILDREN; i++)
		 data[i]->Print();
   }


   (*sigma) *= 1.0/nParams;
   Cholesky R(*sigma);

   if (DEBUG){
	   sigma->Print();
   }

   /* *** Start Server *** */
   cout << "Server Ready ..." << endl;
   TServerSocket *ss = new TServerSocket(9090, kTRUE);
   TMessage *message;

   // Accept a connection and return a full-duplex communication socket.
   TSocket *sock[CHILDREN];
   for(int i=0; i<CHILDREN; i++)
      sock[i] = ss->Accept();

   // tell the clients to start
   for(Int_t i=0; i<CHILDREN; i++){
      TString id = "j_";
      id+=i;
      TObjString oid;
      oid.SetString(id.Data());
      message = new TMessage(kINIT);
      message->WriteObject(ws);
      message->WriteObject(data[i]);
      message->WriteObject(&names);
      message->WriteObject(&oid);
      sock[i]->Send(*message);
      delete message;
   }

   // Close the server socket (unless we will use it later to wait for
   // another connection).
   ss->Close();

   // Check some options of socket.
   for(int i=0; i<CHILDREN; i++){
      printf("** Socket %d:\n",i);
      int val;
      sock[i]->GetOption(kSendBuffer, val);
      printf("sendbuffer size: %d\n", val);
      sock[i]->GetOption(kRecvBuffer, val);
      printf("recvbuffer size: %d\n", val);

      // Get the remote addresses (informational only).
      TInetAddress adr = sock[i]->GetInetAddress();
      adr.Print();
      printf("\n");
   }

   /* *** Create Output TTree ** */
   TVectorD mu(nParams);
   TVectorD mup(nParams);

   TFile *out_file = new TFile("mcmc.root","RECREATE");
   TTree *tree = new TTree("mcmc","mcmc");
   double *pmu = mu.GetMatrixArray();
   for (int i=0; i<nParams; i++){
      RooRealVar *var = (RooRealVar*) &(fitParams[i]);
      pmu[i] =  var->getVal();
      TString data_type = var->GetName();
      data_type += "/D";
      tree->Branch(var->GetName(),&pmu[i],data_type.Data());
   }
   tree->Fill();

   fitWorkspace->Close();

   TMonitor *monitor = new TMonitor;
   for(int i=0; i<CHILDREN; i++)
        monitor->Add(sock[i]);

   TRandom3 *aleatorio = new TRandom3(0);
   double L[CHILDREN];
   int n_found=0;
   double g[50], x[50];
   double old_nll = fitResult->minNll();
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

//	   cout << "sending ..." << endl;
	   /* Send point */
           message = new TMessage(kPOINT);
           message->WriteObject(&mup);
           for (int i=0; i<CHILDREN; i++){
               sock[i]->Send(*message);
	   }
           delete message;

           /* Receive L */
	   int n_complete = 0;
           while (n_complete<CHILDREN){
                 TSocket  *socket;
           
                 socket = monitor->Select();
                 socket->Recv(message);
 
                 TVectorD* likelihood = (TVectorD*) message->ReadObjectAny(TVectorD::Class());
                 L[n_complete] = (*likelihood)[0];
                 n_complete++;

	   	if (DEBUG)
			cout << "Received = " << n_complete << endl;
                 delete message;
           }
	   
	   if (DEBUG){
		for (int i=0; i<CHILDREN; i++)
			cout << "i = " << i << "  L = " << L[i] << endl;
	   }

           /* Metropolis-Hasting */
	   double new_nll = 0;
           for (int i=0; i<CHILDREN; i++)
                new_nll += L[i];
            
	   if (DEBUG){
   		for (int i=0; i<nParams; i++)
	           ws->var( ((TObjString*)names.At(i))->GetString().Data() )->setVal(mup(i));
		double test_nll = nll->getVal();
	
		cout << "TEST NLL = " << new_nll << " - " << test_nll << " = " << new_nll-test_nll << "  " << old_nll << endl;
//		if (  TMath::Abs(new_nll-test_nll) >  1e-5 )
//			exit(1);
	   }


	   if (DEBUG)
		   cout << " new_nll = " << new_nll << endl;

           double alpha = TMath::Exp(-new_nll + old_nll);
//           cout << "alpha =" << alpha << endl;

           if (aleatorio->Uniform()<alpha){
                n_found++;
		old_nll = new_nll;
                for (int i=0; i<nParams; i++)
                    mu(i) = mup(i);
                tree->Fill();
		cout << "N FOUND = " << n_found << endl;
	   }
    }

		
   for(int i=0; i<CHILDREN; i++)
      printf("Client %d: bytes recv = %d, bytes sent = %d\n",
		      i, sock[i]->GetBytesRecv(), sock[i]->GetBytesSent());

   // Close the socket.
   for(int i=0; i<CHILDREN; i++)
      sock[i]->Close();

   tree->Write();
   out_file->Close();

}
