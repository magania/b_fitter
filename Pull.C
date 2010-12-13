Pull(){
	TFile *root_file = new TFile("toymc.root");

	RooWorkspace *ws2 = new RooWorkspace();

        ws2->factory("{DG[0.09,-0.1,0.3],tau[1.5,1.2,1.7],beta[0.25,0.0,0.5],delta_l[0.2,-2,2],delta_p[2.95,1,5.0],delta_s[0.2,-2,2],fs[0.2,0,0.4],s[1.0,0.8,1.2]}");
        ws2->factory("{A02[0.6,0.5,0.7],Al2[0.3,0.2,0.4]}");

	RooDataSet* data = new RooDataSet("data","data", ws2->argSet("DG,tau,beta,delta_l,delta_p,delta_s,fs,s,A02,Al2"));

	for(int i=0; i<10; i++){
		TString fit_name = "fit_";
		fit_name+=i;

		RooFitResult *fit_result = (RooFitResult*) root_file->Get(fit_name);
		if (fit_result->status() != 0){
			cout << fit_result->status() << endl;
			continue;
		}
		//cout << fit_result->minNll() << endl;
//		if ( fit_result->minNll() < 17500 )
			data->add(fit_result->floatParsFinal());
	}

        gROOT->SetStyle("Plain");

	ws2->factory("Gaussian::A02_gauss(A02,mean_A02[0.6,0.5,0.7],sigma_A02[0.15,0.0001,1.0])");
	ws2->pdf("A02_gauss")->fitTo(*data, RooFit::PrintLevel(-1));	

	TCanvas *canvas_A02 = new TCanvas("canvas_A02", "canvas_A02", 600,600);
	RooPlot *A02_plot = ws2->var("A02")->frame();
	data->plotOn(A02_plot,RooFit::MarkerSize(0.3));
//	ws2->pdf("A02_gauss")->plotOn(A02_plot);
//	ws2->pdf("A02_gauss")->paramOn(A02_plot);
	A02_plot->Draw();

        ws2->factory("Gaussian::Al2_gauss(Al2,mean_Al2[0.3,0.2,0.4],sigma_Al2[0.15,0.001,1.0])");
        ws2->pdf("Al2_gauss")->fitTo(*data, RooFit::PrintLevel(-1));

        TCanvas *canvas_Al2 = new TCanvas("canvas_Al2", "canvas_Al2", 600,600);
        RooPlot *Al2_plot = ws2->var("Al2")->frame();
        data->plotOn(Al2_plot,RooFit::MarkerSize(0.3));
//        ws2->pdf("Al2_gauss")->plotOn(Al2_plot);
//        ws2->pdf("Al2_gauss")->paramOn(Al2_plot);
        Al2_plot->Draw();


        ws2->factory("Gaussian::DG_gauss(DG,mean_DG[0.09,0.0,0.2],sigma_DG[0.15,0.01,1.0])");
        ws2->pdf("DG_gauss")->fitTo(*data, RooFit::PrintLevel(-1));

        TCanvas *canvas_DG = new TCanvas("canvas_DG", "canvas_DG", 600,600);
        RooPlot *DG_plot = ws2->var("DG")->frame();
        data->plotOn(DG_plot,RooFit::MarkerSize(0.3));
//        ws2->pdf("DG_gauss")->plotOn(DG_plot);
//        ws2->pdf("DG_gauss")->paramOn(DG_plot);
        DG_plot->Draw();

        ws2->factory("Gaussian::beta_gauss(beta,mean_beta[0.25,0.2,0.3],sigma_beta[0.15,0.001,1.0])");
        ws2->pdf("beta_gauss")->fitTo(*data, RooFit::PrintLevel(-1));

        TCanvas *canvas_beta = new TCanvas("canvas_beta", "canvas_beta", 600,600);
        RooPlot *beta_plot = ws2->var("beta")->frame();
        data->plotOn(beta_plot,RooFit::MarkerSize(0.3));
//        ws2->pdf("beta_gauss")->plotOn(beta_plot);
//        ws2->pdf("beta_gauss")->paramOn(beta_plot);
        beta_plot->Draw();


        ws2->factory("Gaussian::tau_gauss(tau,mean_tau[1.5,1,2],sigma_tau[0.1,0.01,1.0])");
        ws2->pdf("tau_gauss")->fitTo(*data, RooFit::PrintLevel(-1));

        TCanvas *canvas_tau = new TCanvas("canvas_tau", "canvas_tau", 600,600);
        RooPlot *tau_plot = ws2->var("tau")->frame();
        data->plotOn(tau_plot,RooFit::MarkerSize(0.3));
//        ws2->pdf("tau_gauss")->plotOn(tau_plot);
//        ws2->pdf("tau_gauss")->paramOn(tau_plot);
        tau_plot->Draw();

        ws2->factory("Gaussian::s_gauss(s,mean_s[1,0,2],sigma_s[0.15,0.01,1.0])");
        ws2->pdf("s_gauss")->fitTo(*data, RooFit::PrintLevel(-1));

        TCanvas *canvas_s = new TCanvas("canvas_s", "canvas_s", 600,600);
        RooPlot *s_plot = ws2->var("s")->frame();
        data->plotOn(s_plot,RooFit::MarkerSize(0.3));
//        ws2->pdf("s_gauss")->plotOn(s_plot);
//        ws2->pdf("s_gauss")->paramOn(s_plot);
       	s_plot->Draw();

        ws2->factory("Gaussian::fs_gauss(fs,mean_fs[0.2,0,2],sigma_fs[0.15,0.01,1.0])");
        ws2->pdf("fs_gauss")->fitTo(*data, RooFit::PrintLevel(-1));

        TCanvas *canvas_fs = new TCanvas("canvas_fs", "canvas_fs", 600,600);
        RooPlot *fs_plot = ws2->var("fs")->frame();
        data->plotOn(fs_plot,RooFit::MarkerSize(0.3));
//        ws2->pdf("fs_gauss")->plotOn(fs_plot);
//       ws2->pdf("fs_gauss")->paramOn(fs_plot);
        fs_plot->Draw();

//        ws2->factory("Gaussian::delta_s_gauss(delta_s,mean_delta_s[0.2,0,2],sigma_delta_s[0.15,0.1,1.0])");
//        ws2->pdf("delta_s_gauss")->fitTo(*data, RooFit::PrintLevel(-1));

        TCanvas *canvas_delta_s = new TCanvas("canvas_delta_s", "canvas_delta_s", 600,600);
        RooPlot *delta_s_plot = ws2->var("delta_s")->frame();
        data->plotOn(delta_s_plot,RooFit::MarkerSize(0.3));
//        ws2->pdf("delta_s_gauss")->plotOn(delta_s_plot);
//        ws2->pdf("delta_s_gauss")->paramOn(delta_s_plot);
        delta_s_plot->Draw();

//        ws2->factory("Gaussian::delta_p_gauss(delta_p,mean_delta_p[2.95,0,4],sigma_delta_p[0.15,0.1,1.0])");
//        ws2->pdf("delta_p_gauss")->fitTo(*data, RooFit::PrintLevel(-1));

        TCanvas *canvas_delta_p = new TCanvas("canvas_delta_p", "canvas_delta_p", 600,600);
        RooPlot *delta_p_plot = ws2->var("delta_p")->frame();
        data->plotOn(delta_p_plot,RooFit::MarkerSize(0.3));
//        ws2->pdf("delta_p_gauss")->plotOn(delta_p_plot);
//        ws2->pdf("delta_p_gauss")->paramOn(delta_p_plot);
        delta_p_plot->Draw();


//       ws2->factory("Gaussian::delta_l_gauss(delta_l,mean_delta_l[1,0,2],sigma_delta_l[0.25,0.2,2.0])");
//        ws2->pdf("delta_l_gauss")->fitTo(*data, RooFit::PrintLevel(-1));

        TCanvas *canvas_delta_l = new TCanvas("canvas_delta_l", "canvas_delta_l", 600,600);
        RooPlot *delta_l_plot = ws2->var("delta_l")->frame();
        data->plotOn(delta_l_plot,RooFit::MarkerSize(0.3));
//        ws2->pdf("delta_l_gauss")->plotOn(delta_l_plot);
//        ws2->pdf("delta_l_gauss")->paramOn(delta_l_plot);
        delta_l_plot->Draw();
}

