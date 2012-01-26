autocorrelation(){
	TFile file("mcmc.root");
	TTree *tree = (TTree*) file.Get("mcmc");

	double beta;
	tree->SetBranchAddress("beta",&beta);

	cout.precision(5);
	cout.setf(ios::fixed,ios::floatfield); 
	for ( int burnIn =0; burnIn<10000; burnIn += 1000){

	int n_entries = tree->GetEntries();
	double mean = 0;
	for(int i=burnIn; i<n_entries; i++){
		tree->GetEntry(i);
		mean += beta;
	}
	mean /= (n_entries - burnIn);
//	cout << "Mean = " << mean << endl;
	double sigma = 0;
	for(int i=burnIn; i<n_entries; i++){
		tree->GetEntry(i);
		sigma += (beta - mean)*(beta - mean);
	}
	sigma /= (n_entries-burnIn);
//	cout << "sigma^2 = " << sigma << endl;

	for (int k=1; k<=200; k+= 10){
		double r = 0;
		tree->GetEntry(0);
		double beta_old = beta;
		for(int i=k+burnIn; i<n_entries; i += k){
			tree->GetEntry(i);
			r += (beta - mean)*(beta_old - mean);
			beta_old = beta;
		}
		r /= (n_entries - k - burnIn)*sigma;
		cout <<  TMath::Abs(r) << "\t";
	}
	cout << endl;
   }
}
