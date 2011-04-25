static int NUM_CPU = 8;
static bool VERBOSE = false;
static char* SERVER = "localhost";

static bool EVOL = false;
static int kINIT  = 20000;
static int kGET   = 20001;
static int kFOUND = 20002;

void markov_client(){
   TMessage::EnableSchemaEvolutionForAll(EVOL);

   gSystem->Load("libBFitter.so");

   TMessage *message;
   RooWorkspace *ws;
   TMatrixDSym *sigma;
   TObjString *id;
   TObjArray *names;

   TSocket *socket = new TSocket(SERVER,9090);
   socket->Recv(message);

   if(message->What() == kINIT){
	   ws = (RooWorkspace*) message->ReadObjectAny(RooWorkspace::Class());
	   sigma = (TMatrixDSym*) message->ReadObjectAny(TMatrixDSym::Class());
	   names = (TObjArray*) message->ReadObjectAny(TObjArray::Class());
	   id = (TObjString*) message->ReadObjectAny(TObjString::Class());
	   cout << "Jod ID = " << id->GetString() << endl;
   } else {
	   cout << "EE: Failed to init." << endl;
	   exit(1);
   }
   /*
   ws->Print();
   sigma->Print();
   names->Print();
   id->Print();
   */
  
   Int_t N = sigma->GetNcols();
   if (N > 50){
  	cout << "EE: Increase x and mu arrays size." << endl;
	exit(1);
   }
   double x[50],ex[50],mu[50],mup[50];
   RooRealVar *params[50];
   for (int i=0; i<N; i++){
        ex[i] = TMath::Sqrt((*sigma)(i,i));
	params[i] = ws->var( ((TObjString*)names->At(i))->GetString().Data() );
   }
   sigma->Invert();
   RooAbsReal *nll = ws->pdf("model")->createNLL(*ws->data("data"), RooFit::Verbose(VERBOSE), RooFit::NumCPU(NUM_CPU));
   TRandom3 *aleatorio = new TRandom3(0);

   /* Generate sample of multivariate gaussian
      Now woks by accept-reject method with up boundary = 1
      TODO: check if its posible use as boundary G(x1)G(x2)..G(xn)
            where G is a gaussian distribution.
   */


   while(true){
	   /* Generate gausian point */
	   int tries_gauss = 0;
	   while(true){
		   tries_gauss++;
		   for (int i=0; i<N; i++)
			   x[i] = aleatorio->Uniform(-ex[i]/2.,ex[i]/2.);
			   
		   // Evaluate multivariate gaussian
		   double p=0;
		   for (int i=0; i<N; i++)
			   p += x[i]*(*sigma)(i,i)*x[i];
   		   for (int i=0; i<N; i++)
			   for (int j=i+1; j<N; j++)
				   p += 2.0*x[i]*(*sigma)(i,j)*x[j];
   
		   double P = TMath::Exp(-p/2.);       
		   if (aleatorio->Uniform()<P) break;
	   }
	   cout << "Multigaussian(" << tries_gauss << " tries.)"<< endl;

	   /*Get current point */
	   message = new TMessage(kGET);
	   socket->Send(*message);
	   delete message;
	   socket->Recv(message);
	   TVectorD *vmu = (TVectorD*) message->ReadObjectAny(TVectorD::Class());
	   for (int i=0; i<N; i++){
		   mu[i] = (*vmu)[i];
		   mup[i] = mu[i]+x[i];
	   }
	   delete message;

	   /*Evaluate Likelihood */
	   for (int i=0; i<N; i++)
		   mu[i] = params[i]->setVal(mu[i]);
	   double old_nll = nll->getVal();
      
	   for (int i=0; i<N; i++)
		   mu[i] = params[i]->setVal(mup[i]);
	   double new_nll = nll->getVal();

	   double alpha = TMath::Exp(-new_nll + old_nll);
	   cout << "alpha =" << alpha << endl;

	   /* If accepted Send */
	   if (aleatorio->Uniform()<alpha){
		   message = new TMessage(kFOUND);
		   TVectorD current(N,mu);
		   TVectorD found(N,mup);
		   TObjArray pair(2);
		   pair.AddAt(&current,0);
		   pair.AddAt(&found,1);
		   message->WriteObject(&pair);
		   socket->Send(*message);
		   delete message;
	   }
   }
}
