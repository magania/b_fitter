static int n_clients = 1;
static int NSAMPLE = 1;
const char* fitWorkspaceFile = "ws_CUT_fit.root";

static bool EVOL = false;
static int kINIT  = 20000;
static int kGET   = 20001;
static int kFOUND = 20002;

void markov_server() {
   TMessage::EnableSchemaEvolutionForAll(EVOL);

   gSystem->Load("libBFitter.so");

   TFile *fitWorkspace = new TFile(fitWorkspaceFile);
   RooWorkspace *ws = (RooWorkspace *) fitWorkspace->Get("rws");
   RooFitResult *fitResult = (RooFitResult*) fitWorkspace->Get("fitResult");
   RooArgList fitParams = fitResult->floatParsFinal();
   Int_t nParams = fitParams.getSize();
   TMatrixDSym *sigma = &((TMatrixDSym*) fitResult->covarianceMatrix());

   TObjArray names(nParams);
   for (int i=0; i<nParams; i++)
   	names.Add(new TObjString(fitParams[i].GetName()));

   /* *** Start Server *** */
   TServerSocket *ss = new TServerSocket(9090, kTRUE);

   // Accept a connection and return a full-duplex communication socket.
   TSocket *sock[200];
   for(int i=0; i<n_clients; i++)
      sock[i] = ss->Accept();

   // tell the clients to start
   for(Int_t i=0; i<n_clients; i++){
      TString id = "j_";
      id+=i;
      TObjString oid;
      oid.SetString(id.Data());
      message = new TMessage(kINIT);
      message->WriteObject(ws);
      message->WriteObject(sigma);
      message->WriteObject(&names);
      message->WriteObject(&oid);
      sock[i]->Send(*message);
      delete message;
   }

   // Close the server socket (unless we will use it later to wait for
   // another connection).
   ss->Close();

   // Check some options of socket.
   for(int i=0; i<n_clients; i++){
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
   for(int i=0; i<n_clients; i++)
	monitor->Add(sock[i]);

   int n_found=0;
   while (n_found<NSAMPLE) {
      TMessage *message;
      TSocket  *socket;

      socket = monitor->Select();
      socket->Recv(message);

      if (!message){
         printf("No mess :(\n");
         continue;
      }
     
      if (message->What() == kGET) {
	  cout << "GET" << endl;
	  TMessage response(kMESS_OBJECT);
	  response.WriteObject(&mu);
	  socket->Send(response);
      } else if (message->What() == kFOUND){
	  cout << "FOUND" << endl;
	  TObjArray *pair = (TObjArray*) message->ReadObjectAny(TObjArray::Class());
	  TVectorD *pold = pair->At(0);
	  TVectorD *pnew = pair->At(1);
	  bool from_current = true;
	  for (int i=0; i<nParams; i++)
		  if ( mu[i] != (*pold)[i]) 
			  from_current = false;
	  if (from_current);
	  for (int i=0; i<nParams; i++)
		  mu[i] = (*pnew)[i];
	  tree->Fill();
	  n_found++;
      } else if (message->What() == kMESS_STRING) {
         char str[64];
         message->ReadString(str, 64);
	 int i=0;
	 for(;i<n_clients;i++)
		 if (socket==sock[i]) break;
         printf("Client %d: %s\n", i, str);
      } else {
         printf("*** Unexpected message ***\n");
      }

      delete message;
   }

   for(int i=0; i<n_clients; i++)
      printf("Client %d: bytes recv = %d, bytes sent = %d\n",
		      i, sock[i]->GetBytesRecv(), sock[i]->GetBytesSent());

   // Close the socket.
   for(int i=0; i<n_clients; i++)
      sock[i]->Close();

   tree->Write();
   out_file->Close();

}
