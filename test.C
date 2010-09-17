test(){
	TFile acceptance_file("bs.root");
        TTree* acceptance_tree = (TTree*) acceptance_file.Get("tree");

        Double_t cpsi, ctheta, phi;
        acceptance_tree->SetBranchAddress("bs_angle_cpsi", &cpsi);
        acceptance_tree->SetBranchAddress("bs_angle_ctheta", &ctheta);
        acceptance_tree->SetBranchAddress("bs_angle_phi", &phi);

        TH2D acceptance_histo("acceptance_histo", "acceptance", 10, -1.0, 1.0, 10, -TMath::Pi(), TMath::Pi());

        Long64_t n_entries = acceptance_tree->GetEntries();
        for (Long64_t i=0; i<n_entries; i++){
                acceptance_tree->GetEntry(i);
                acceptance_histo.Fill(ctheta,phi);
		//cout << cpsi << ' ' << ctheta << ' ' << phi << endl;
        }

	acceptance_histo.Print("all");
        acceptance_histo.Draw("acceptance.png");
}
