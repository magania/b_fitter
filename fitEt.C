/*
 *  Author: Ricardo Magana-Villalba
 *          magania@fnal.gov
 *
 *  September 2010
 */


#include "TH2.h"
#include "TF2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TVirtualFitter.h"
#include "TPaveLabel.h"
#include "TStyle.h"
#include "TMath.h"
#include "TROOT.h"
#include "TFrame.h"
#include "TFile.h"
#include "TTree.h"
#include "TEntryList.h"
#include "TCut.h"

#include "TRandom3.h"
#include "Math/SpecFuncMathMore.h"

//#include "Fit/FitConfig.h"


int fitEt() {
	gSystem->Load("libMathMore");
	gROOT->SetStyle("Plain");

        TRandom3 rnd(6438238139);

	TFile acceptance_file("/work/elsanto-clued0/Z/Bs/tmva/bs_flat.root");
	TTree* acceptance_tree = (TTree*) acceptance_file.Get("tree");

	Double_t epdl;
	acceptance_tree->SetBranchAddress("bs_epdl", &epdl);

        Double_t muPlusPt, muMinusPt, jpsiPt, muPlusEta, muMinusEta;
	acceptance_tree->SetBranchAddress("mu_plus_pt", &muPlusPt);
	acceptance_tree->SetBranchAddress("mu_minus_pt", &muMinusPt);
	acceptance_tree->SetBranchAddress("mu_plus_eta", &muPlusEta);
	acceptance_tree->SetBranchAddress("mu_minus_eta", &muMinusEta);
	acceptance_tree->SetBranchAddress("jpsi_pt", &jpsiPt);

        TF1 *weightFcnCentral = new TF1("weightFcnCentral","landau(0)+pol1(3)",3.85,45);
        weightFcnCentral->SetParameters(9.76584e+00,1.13316e+01,3.69072e+00,-4.52973e-01,1.44476e-02);

        TF1 *weightFcnForward = new TF1("weightFcnCentral","landau(0)+pol1(3)",3.02,45);
        weightFcnForward->SetParameters(8.93988e+00,1.16932e+01,3.94829e+00,-3.07532e-01,2.27729e-02);

	TH1D* et_histo = new TH1D("et_histo", "sigma(et)", 100, 0, 0.02 );
	RooRealVar et("et","et",0,0.4);
	RooDataSet dataSet("dataSet","dataSet",et);

        TCut QC("(k_minus_nSMT>1 && k_plus_nSMT>1 && mu_minus_nSMT>1 && mu_plus_nSMT>1 && ((mu_minus_nseg==1 && mu_plus_nseg==3) || (mu_minus_nseg==3 && mu_plus_nseg==1) || (mu_minus_nseg==3 && mu_plus_nseg==3)) && TMath::Abs(bs_vrt_z-bs_pv_vrt_z)<5.0)");
        TCut PRL("mu_plus_cpt>1.5 && mu_minus_cpt>1.5 && 2.9 < jpsi_mass  && jpsi_mass < 3.3  && k_plus_cpt>0.7 && k_minus_cpt>0.7 && phi_cpt>1.5 && 1.01 < phi_mass_corrected && phi_mass_corrected<1.03 && bs_pt>6 && bs_vrt_chi2<36 && TMath::Abs(bs_pv_vrt_z - bs_vrt_z) < 5 && k_plus_nSMT>1 && k_minus_nSMT>1 && mu_plus_nSMT>1 && mu_minus_nSMT>1  && k_plus_nSMT+k_plus_nCFT>7 && k_minus_nSMT+k_minus_nCFT>7  && mu_plus_nSMT+mu_plus_nCFT>7 && mu_minus_nSMT+mu_minus_nCFT>7  && (mu_plus_nseg==3||mu_plus_nseg==1) && (mu_minus_nseg==3||mu_minus_nseg==1) && ( jpsi_cpt>4 || (mu_plus_cpt>mu_minus_cpt && TMath::Abs(mu_plus_eta)>=1) || (mu_minus_cpt>mu_plus_cpt && TMath::Abs(mu_minus_eta)>=1) )");

	TCut cut0("mc_match == 1");

        TCut BDT10("inclusiveBDT>0.45 && promptBDT>0.42");
        TCut BDT12("inclusiveBDT>0.39 && promptBDT>0.35");
        TCut BDT14("inclusiveBDT>0.33 && promptBDT>0.28");
        TCut BDT16("inclusiveBDT>0.21 && promptBDT>0.29");
        TCut BDT18("inclusiveBDT>0.07 && promptBDT>0.29");
        TCut BDT20("inclusiveBDT>-0.05 && promptBDT>-0.01");

	TCut cut = BDT20 + QC + cut0;

	acceptance_tree->Draw(">>entry_list", cut, "entrylist");
	TEntryList *event_list = (TEntryList*)gDirectory->Get("entry_list");

	Long64_t n_entries = event_list->GetN();
	std::cout << "Selected: " << n_entries << " events." << std::endl;

	for (Long64_t i=0; i<n_entries; i++){
		acceptance_tree->GetEntry(event_list->GetEntry(i));
                 
		Double_t muLeadingPt = muPlusPt>muMinusPt?muPlusPt:muMinusPt;
		Double_t muLeadingEta = muPlusPt>muMinusPt?muPlusEta:muMinusEta;

                if (jpsiPt > 45)
		   jpsiPt = 45;
		Double_t weight = fabs(muLeadingEta)<1?weightFcnCentral->Eval(jpsiPt):weightFcnForward->Eval(jpsiPt);
                if (weight < 0)
		  weight = 0;

		et_histo->Fill(epdl,weight);
		et = epdl/0.0299792458;
		dataSet.add(et,weight);
//		std::cout << epdl <<  weight << std::endl;
//		std::cout << muPlusPt << ' ' << muMinusPt << ' ' << muPlusEta << ' ' << muMinusEta << ' ' << weight << std::endl;
	}

	cout << "All Entries = " << dataSet.sumEntries() << " " << et_histo->GetSumOfWeights() << endl;

	RooWorkspace rws;
	rws.import(et);
//	rws.import(dataSet);

//  rws.factory("GaussModel::etGaussianS(et,meanGaussEtS[0.0614,0,0.2],sigmaGaussEtS[0.0116,0,0.2])");
//  rws.factory("Decay::errorSignal(et,tauEtS[0.0481,0,0.2],etGaussianS,SingleSided]");

//	rws.factory("RSUM::errorSignal(x1[0,1]*Gaussian::g1(et,m1[0,1],s1[0,0.1]),x2[0,1]*Gaussian::g2(et,m2[0,1],s2[0,0.1]),x3[0,1]*Gaussian::g3(et,m3[0,1],s3[0,0.1]),x4[0,1]*Gaussian::g4(et,m4[0,1],s4[0,0.1]),Gaussian::g5(et,m5[0,1],s5[0,0.1]))");

//BDT
rws.factory("RSUM::errorSignal(x1[0.330954,0,1]*Gaussian::g1(et,m1[0.0712217,0,1],s1[0.0148496,0,1]),x2[0.550785,0,1]*Gaussian::g2(et,m2[0.101322,0,1],s2[0.0227363,0,1]),x3[0.0667897,0,1]*Gaussian::g3(et,m3[0.20006,0,1],s3[0.0767951,0,1]),x4[0.615575,0,1]*Gaussian::g4(et,m4[0.142648,0,1],s4[0.0366608,0,1]),Gaussian::g5(et,m5[0.0488805,0,1],s5[0.00966675,0,1]))");

rws.factory("RSUM::errorSignal1(y1[0.367333]*Gaussian::h1(et,n1[0.0714381],t1[0.0142117]),y2[0.545832]*Gaussian::h2(et,n2[0.0983138],t2[0.0221079]),y3[0.0819051]*Gaussian::h3(et,n3[0.186739],t3[0.0687823]),y4[0.545434]*Gaussian::h4(et,n4[0.141453],t4[0.0354564]),Gaussian::h5(et,n5[0.0485111],t5[0.00933244]))"); // variation 1

rws.factory("RSUM::errorSignal2(z1[0.347236]*Gaussian::i1(et,o1[0.0707638],u1[0.0153091]),z2[0.540398]*Gaussian::i2(et,o2[0.0992327],u2[0.0225856]),z3[0.0998333]*Gaussian::i3(et,o3[0.182831],u3[0.0718767]),z4[0.612224]*Gaussian::i4(et,o4[0.143124],u4[0.0347972]),Gaussian::i5(et,o5[0.0486886],u5[0.0092695]))"); // variation 2


//PRL
//rws.factory("RSUM::errorSignal(x1[0.0164413,0,1]*Gaussian::g1(et,m1[0.194133,0,1],s1[0.0754175,0,1]),x2[0.371492,0,1]*Gaussian::g2(et,m2[0.0968759,0,1],s2[0.0212262,0,1]),x3[0.185101,0,1]*Gaussian::g3(et,m3[0.0477364,0,1],s3[0.00920869,0,1]),x4[0.705251,0,1]*Gaussian::g4(et,m4[0.068875,0,1],s4[0.0141091,0,1]),Gaussian::g5(et,m5[0.135975,0,1],s5[0.033775,0,1]))");

//	rws.pdf("errorSignal")->fitTo(dataSet, RooFit::NumCPU(4));


        cout << "rws.factory(\"RSUM::errorSignal("
	     << "x1["<< rws.var("x1")->getVal() <<"]*Gaussian::g1(et,m1["<< rws.var("m1")->getVal() <<"],s1["<< rws.var("s1")->getVal() <<"]),"
             << "x2["<< rws.var("x2")->getVal() <<"]*Gaussian::g2(et,m2["<< rws.var("m2")->getVal() <<"],s2["<< rws.var("s2")->getVal() <<"]),"
             << "x3["<< rws.var("x3")->getVal() <<"]*Gaussian::g3(et,m3["<< rws.var("m3")->getVal() <<"],s3["<< rws.var("s3")->getVal() <<"]),"
             << "x4["<< rws.var("x4")->getVal() <<"]*Gaussian::g4(et,m4["<< rws.var("m4")->getVal() <<"],s4["<< rws.var("s4")->getVal() <<"]),"
             << "Gaussian::g5(et,m5["<< rws.var("m5")->getVal() <<"],s5["<< rws.var("s5")->getVal() <<"]))\");"
	     << endl<<endl;

        cout << "rws.factory(\"RSUM::errorSignal("
	     << "x1["<< rnd.Gaus(rws.var("x1")->getVal(),rws.var("x1")->getError()) <<"]*Gaussian::g1(et,m1["<< rnd.Gaus(rws.var("m1")->getVal(),rws.var("m1")->getError()) <<"],s1["<< rnd.Gaus(rws.var("s1")->getVal(),rws.var("s1")->getError()) <<"]),"
             << "x2["<< rnd.Gaus(rws.var("x2")->getVal(),rws.var("x2")->getError()) <<"]*Gaussian::g2(et,m2["<< rnd.Gaus(rws.var("m2")->getVal(),rws.var("m2")->getError()) <<"],s2["<< rnd.Gaus(rws.var("s2")->getVal(),rws.var("s2")->getError()) <<"]),"
             << "x3["<< rnd.Gaus(rws.var("x3")->getVal(),rws.var("x3")->getError()) <<"]*Gaussian::g3(et,m3["<< rnd.Gaus(rws.var("m3")->getVal(),rws.var("m3")->getError()) <<"],s3["<< rnd.Gaus(rws.var("s3")->getVal(),rws.var("s3")->getError()) <<"]),"
             << "x4["<< rnd.Gaus(rws.var("x4")->getVal(),rws.var("x4")->getError()) <<"]*Gaussian::g4(et,m4["<< rnd.Gaus(rws.var("m4")->getVal(),rws.var("m4")->getError()) <<"],s4["<< rnd.Gaus(rws.var("s4")->getVal(),rws.var("s4")->getError()) <<"]),"
             << "Gaussian::g5(et,m5["<< rnd.Gaus(rws.var("m5")->getVal(),rws.var("m5")->getError()) <<"],s5["<< rnd.Gaus(rws.var("s5")->getVal(),rws.var("s5")->getError()) <<"]))\");"
	     << endl<<endl;

        cout << "rws.factory(\"RSUM::errorSignal("
	     << "x1["<< rnd.Gaus(rws.var("x1")->getVal(),rws.var("x1")->getError()) <<"]*Gaussian::g1(et,m1["<< rnd.Gaus(rws.var("m1")->getVal(),rws.var("m1")->getError()) <<"],s1["<< rnd.Gaus(rws.var("s1")->getVal(),rws.var("s1")->getError()) <<"]),"
             << "x2["<< rnd.Gaus(rws.var("x2")->getVal(),rws.var("x2")->getError()) <<"]*Gaussian::g2(et,m2["<< rnd.Gaus(rws.var("m2")->getVal(),rws.var("m2")->getError()) <<"],s2["<< rnd.Gaus(rws.var("s2")->getVal(),rws.var("s2")->getError()) <<"]),"
             << "x3["<< rnd.Gaus(rws.var("x3")->getVal(),rws.var("x3")->getError()) <<"]*Gaussian::g3(et,m3["<< rnd.Gaus(rws.var("m3")->getVal(),rws.var("m3")->getError()) <<"],s3["<< rnd.Gaus(rws.var("s3")->getVal(),rws.var("s3")->getError()) <<"]),"
             << "x4["<< rnd.Gaus(rws.var("x4")->getVal(),rws.var("x4")->getError()) <<"]*Gaussian::g4(et,m4["<< rnd.Gaus(rws.var("m4")->getVal(),rws.var("m4")->getError()) <<"],s4["<< rnd.Gaus(rws.var("s4")->getVal(),rws.var("s4")->getError()) <<"]),"
             << "Gaussian::g5(et,m5["<< rnd.Gaus(rws.var("m5")->getVal(),rws.var("m5")->getError()) <<"],s5["<< rnd.Gaus(rws.var("s5")->getVal(),rws.var("s5")->getError()) <<"]))\");"
	     << endl<<endl;

	TCanvas canvas("canvas","Sigma (et)", 600,600);

	RooPlot *et_frame = et.frame();
	et_frame->SetTitle("");
	dataSet.plotOn(et_frame);
	rws.pdf("errorSignal1")->plotOn(et_frame,RooFit::LineColor(kRed));
	rws.pdf("errorSignal2")->plotOn(et_frame,RooFit::LineColor(kRed));
	rws.pdf("errorSignal")->plotOn(et_frame);
	et_frame->Draw();
	gPad->SetLogy(kTRUE);
	
/*
	TF2 *fitFcn = new TF2("fitFcn", fitFunction, -1.0, 1.0, -TMath::Pi(), TMath::Pi(), 9);
	fitFcn->SetParameters(0,0,0,0,0,0,0,0,0);
	//acceptance_histo->Fit(fitFcn,"EV");

	acceptance_histo->Fit(fitFcn,"EV0");

        gStyle->SetPalette(1);

	TCanvas canvas("canvas", "canvas", 600, 600);
	canvas.Divide(2,2);
	canvas.cd(1);
	acceptance_histo->Draw("COLZ");
//	acceptance_histo->ProjectionX()->Draw();
	acceptance_histo->Fit(fitFcn,"EV0");
//	fitFcn->Draw();
//
	TH2D* generated_histo = new TH2D("generated_histo", "acceptance", 20, -1.0, 1.0, 20, -TMath::Pi(), TMath::Pi());

        gRandom = new TRandom3();
	Double_t x, y;
	for(int i=0; i<n_entries*1000; i++){
		fitFcn->GetRandom2(x,y);
		generated_histo->Fill(x,y,0.001);
	}

	canvas.cd(2);
	generated_histo->Draw("COLZ");

	canvas.cd(3);
	acceptance_histo->ProjectionX()->Draw("E");
	generated_histo->ProjectionX()->Draw("same");

	canvas.cd(4);
	acceptance_histo->ProjectionY()->Draw("E");
	generated_histo->ProjectionY()->Draw("same");

	Double_t max_bin =  acceptance_histo->GetBinContent(acceptance_histo->GetMaximumBin());


	std::cout << acceptance_histo->Chi2Test(generated_histo) << std::endl;

	std::cout << "static const Double_t e_00p0 = " << fitFcn->GetParameter(0)/max_bin << ';' << std::endl;
	std::cout << "static const Double_t e_01n1 = " << fitFcn->GetParameter(1)/max_bin << ';' << std::endl;
	std::cout << "static const Double_t e_01p0 = " << fitFcn->GetParameter(2)/max_bin << ';' << std::endl;
	std::cout << "static const Double_t e_01p1 = " << fitFcn->GetParameter(3)/max_bin << ';' << std::endl;
	std::cout << "static const Double_t e_02n2 = " << fitFcn->GetParameter(4)/max_bin << ';' << std::endl;
	std::cout << "static const Double_t e_02n1 = " << fitFcn->GetParameter(5)/max_bin << ';' << std::endl;
	std::cout << "static const Double_t e_02p0 = " << fitFcn->GetParameter(6)/max_bin << ';' << std::endl;
	std::cout << "static const Double_t e_02p1 = " << fitFcn->GetParameter(7)/max_bin << ';' << std::endl;
	std::cout << "static const Double_t e_02p2 = " << fitFcn->GetParameter(8)/max_bin << ';' << std::endl;

*/
	canvas.SaveAs("et_signal.C");

	return 1;

}

