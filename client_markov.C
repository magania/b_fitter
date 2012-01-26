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
#include "TDecompChol.h"
#include "TMatrixD.h"
#include "RooDataSet.h"

static int NUM_CPU = 1;
static bool VERBOSE = false;
//static const char* SERVER = "tlaloc-clued0.fnal.gov";

static bool EVOL = false;
static unsigned int kINIT  = 20000;
static unsigned int kPOINT = 20001;

int client_markov(const char* SERVER){
   TMessage::EnableSchemaEvolutionForAll(EVOL);

   gSystem->Load("lib/libBFitter.so");

   TMessage *message;
   RooWorkspace *ws;
   RooDataSet *data;
   TObjString *id;
   TObjArray *names;

   TSocket *socket = new TSocket(SERVER,9090);
   socket->Recv(message);

   if(message->What() == kINIT){
	   ws = (RooWorkspace*) message->ReadObjectAny(RooWorkspace::Class());
	   data = (RooDataSet*) message->ReadObjectAny(RooDataSet::Class());
	   names = (TObjArray*) message->ReadObjectAny(TObjArray::Class());
	   id = (TObjString*) message->ReadObjectAny(TObjString::Class());
	   cout << "Jod ID = " << id->GetString() << endl;
   } else {
	   cout << "EE: Failed to init." << endl;
	   exit(1);
   }

   data->Print();

   ws->Print();
   data->Print();
   names->Print();
//   id->Print();

   int N = names->GetSize();
   cout << "N = " << N << endl;

   if (N > 50){
  	cout << "EE: Increase x and mu arrays size." << endl;
	exit(1);
   }

   RooRealVar *params[50];
   for (int i=0; i<N; i++)
	params[i] = ws->var( ((TObjString*)names->At(i))->GetString().Data() );

   RooAbsReal *nll = ws->pdf("modelnc")->createNLL(*data, RooFit::Verbose(VERBOSE), RooFit::NumCPU(NUM_CPU), RooFit::ConditionalObservables(*ws->var("d")));

   TVectorD L(1);
   while(true){
	   /*Get point */
	   socket->Recv(message);
	   TVectorD *vmu = (TVectorD*) message->ReadObjectAny(TVectorD::Class());

/*
	   for (int i=0; i<N; i++)
		   cout << (*vmu)(i) << " ";
	   cout << endl;
*/

	   /* Evaluate Likelihood */
	   for (int i=0; i<N; i++)
		     params[i]->setVal((*vmu)[i]);
	   L[0] = nll->getVal();

//	   cout << "L = " << L[0] << endl;

	   delete message;
      
	   /* Send Partial Likelihood */
	   message = new TMessage(kPOINT);
	   message->WriteObject(&L);
	   socket->Send(*message);
	   delete message;
   }

   return 0;
}
