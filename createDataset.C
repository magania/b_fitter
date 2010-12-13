createDataset(){
   TFile* tree_file = new TFile("PRL.root");
   
   TTree* tree = (TTree*) tree_file->Get("tree");
   TFile* out_file = new TFile("dataset.root", "RECREATE");

   RooWorkspace* rws = new RooWorkspace();
   RooRealVar m("m", "Mass (#mu#mu KK)", 0, 5, 5.8);
   RooRealVar t("t", "t", 0,-5, 20);
   RooRealVar et("et", "#sigma(t))", 0.01, 0, 3);
   RooRealVar cpsi("cpsi", "cos(#psi)", -1, 1);
   RooRealVar ctheta("ctheta", "cos(#theta)", -1 , 1);
   RooRealVar phi("phi", "#phi", -TMath::Pi(), TMath::Pi());
   RooRealVar d("d", "d", -1, 1);
  

   RooDataSet* data = new RooDataSet("data","data",RooArgSet(m,t,et,cpsi,ctheta,phi,d));
   RooDataSet* dataBkg = new RooDataSet("dataBkg","dataBkg",RooArgSet(m,t,et,cpsi,ctheta,phi,d));

   TCut* cut = new TCut("mu_plus_nseg==3 && mu_minus_nseg==3");
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

       m=dM;
       t=dT/0.0299792458;
       et=dEt/0.0299792458;
       cpsi=dCpsi;
       ctheta=dCtheta;
       phi=dPhi;
       d=ddDefined==1?dd:0.0;

       data->add(RooArgSet(m,t,et,cpsi,ctheta,phi,d));
       if (dM<5.2 || dM>5.6)
           dataBkg->add(RooArgSet(m,t,et,cpsi,ctheta,phi,d));
   }

   rws->import(*data);
   rws->import(*dataBkg);
   rws->Write("rws");
   out_file->Close();
}
