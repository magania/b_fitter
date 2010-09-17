#include <iostream>

#include <THDAcceptance.h>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TMath.h>

#include <RooRealVar.h>
#include <RooDataSet.h>

using namespace std;

int main(int argc, char *argv[]){
    const char* flat_name = "bs_flat.root";
    const char* fit_name = "bs_fit.root";
    const char* cut = "";
    const char* dataset_name = "dataset.root";
 
    int c;
    while((c =  getopt(argc, argv, ":e:d:c:o:")) != EOF){
        switch (c){
             case 'e':
                 flat_name = optarg;
                 cout << "Efficiency file: "<< flat_name << endl;
                 break;
             case 'd':
                 fit_name = optarg;
                 cout << "Data file: " << fit_name << endl;
                 break;
             case 'c':
                 cut = optarg;
                 cout << "Cut: " << cut << endl;
                 break;
             case 'o':
                 dataset_name = optarg;
                 cout << "DataSet file: " << dataset_name << endl;
                 break;
             case ':':
                 cerr << "Missing option." << endl;
                 exit(1);
                 break;
        }
    }

   TFile flat_file(flat_name);

   TTree* bs = (TTree*) flat_file.Get("tree");
//   bs.Draw("cpsi");

   cout << "Flat: " << bs <<endl;

   Double_t cos0l_bins[] = {-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0};
   Double_t phi_bins[] = {-3.141592650,-2.513274119,-1.884955588,-1.256637057,-0.628318526,0.000000005,0.628318536,1.256637067,1.884955598,2.513274129,3.141592660};

   THDAcceptance acceptance("acceptance", "Acceptance histo",10,cos0l_bins,10,phi_bins); 

   Long64_t nentries = bs->GetEntries();

   Double_t _cpsi, _ctheta, _phi;
   bs->SetBranchAddress("bs_angle_cpsi", &_cpsi);
   bs->SetBranchAddress("bs_angle_ctheta", &_ctheta);
   bs->SetBranchAddress("bs_angle_phi", &_phi);

   Double_t cos0l_,cos0k_,phi_;
   for (Long64_t i=0; i<nentries; i++){
   	bs->GetEntry(i);
//	cout << "Trans= " << _cpsi << " " << _ctheta << " "  << _phi << endl;

	cos0l_ = -TMath::Sin(_phi)*_ctheta;
	cos0k_ = _cpsi;
	phi_ = TMath::ACos( TMath::Sin(_phi)*TMath::Sqrt(1-_ctheta*_ctheta)/TMath::Sqrt(1-TMath::Sin(_phi)*TMath::Sin(_phi)*_ctheta*_ctheta) );
	if (TMath::Cos(_phi)<0) 
		phi_ = -phi_;

//	cout << "Helic= " << cos0l_ << " " << cos0k_ << " "  << phi_ << endl;
	acceptance.Fill(cos0l_,phi_);
   }

   TCanvas canvas("canvas","canvas", 800, 600);
   acceptance.Draw();
   canvas.Print("acceptance.eps");

   for (int i=0; i<10; i++){
   	for(int j=0; j<10; j++){
	        cout.width(5);
		cout <<  acceptance.acceptance(cos0l_bins[j]+0.01,1,phi_bins[i]+0.01) << " ";
	}
	cout << endl;
   }
   //cout << "x=" << acceptance.acceptance(0.1,0.1,0.1);

   flat_file.Close();

   /*----------------------------------------------------------------------------------------------------------------------------------*/
   /*----------------------------------------------------------------------------------------------------------------------------------*/

   TFile fit_file(fit_name);

   bs = (TTree*) fit_file.Get("tree");
   cout << "Fit: " << bs <<endl;

   nentries = bs->GetEntries();

   Double_t _m, _pdl,_epdl,_d,D_;
   Int_t _defined;
   bs->SetBranchAddress("bs_mass", &_m);
   bs->SetBranchAddress("bs_pdl", &_pdl);
   bs->SetBranchAddress("bs_epdl", &_epdl);
   bs->SetBranchAddress("bs_angle_cpsi", &_cpsi);
   bs->SetBranchAddress("bs_angle_ctheta", &_ctheta);
   bs->SetBranchAddress("bs_angle_phi", &_phi);
   bs->SetBranchAddress("bs_angle_phi", &_phi);
   bs->SetBranchAddress("newtag_ost",&_d);
   bs->SetBranchAddress("newtag_ost_defined", &_defined);


   RooRealVar m("m","mass", 5.0, 5.8);
   RooRealVar t("t","time", -1.0, 3.0);
   RooRealVar et("et","#sigma(t)", 0.0, 1.0);
   RooRealVar cos0l("cos0l", "Cos(#Theta_{l})", -1, 1);
   RooRealVar cos0k("cos0k", "Cos(#Theta_{k})", -1, 1);
   RooRealVar phi("phi", "#phi", -TMath::Pi(), TMath::Pi());
   RooRealVar D("D", "dilution", -1.0, 1.0);

   RooDataSet data("data", "data", RooArgSet(m,t,et,cos0l,cos0k,phi,D));

   for (Long64_t i=0; i<nentries; i++){
        bs->GetEntry(i);
//      cout << "Trans= " << _cpsi << " " << _ctheta << " "  << _phi << endl;

        cos0l_ = -TMath::Sin(_phi)*_ctheta;
        cos0k_ = _cpsi;
        phi_ = TMath::ACos( TMath::Sin(_phi)*TMath::Sqrt(1-_ctheta*_ctheta)/TMath::Sqrt(1-TMath::Sin(_phi)*TMath::Sin(_phi)*_ctheta*_ctheta) );
        if (TMath::Cos(_phi)<0)
                phi_ = -phi_;

//      cout << "Helic= " << cos0l_ << " " << cos0k_ << " "  << phi_ << endl;
	m = _m;
	t = _pdl/0.0299792458;
	et = _epdl/0.0299792458;
	cos0l = cos0l_;
	cos0k = cos0k_;
	phi = phi_;

	D_=0;
        if ( _defined )
		D_ = (_d/TMath::Abs(_d))*( 0.6625/(1 + TMath::Exp( (0.3104 - TMath::Abs(_d))/0.1184 ) ) - 0.6625/(1 + TMath::Exp(0.3104/0.1184) ) );

	if (fabs(D_)>1){
		cout << "ERROR: Wrong Dilution" << D_ << endl;
		exit(1);
	}
	D = D_;

	data.add(RooArgSet(m,t,et,cos0l,cos0k,phi,D));
   }

   fit_file.Close();

   TFile dataset_file(dataset_name,"recreate");
   acceptance.Write();
   data.Write();
   dataset_file.Write();
   dataset_file.Close();
  
}
