void plotFit(const char* file){
        cout << "Printing ... " << endl;
	gROOT->SetStyle("Plain");

	gSystem->Load("lib/libBFitter.so");

	Int_t nCPU = 4;

    TString rootFile(file);
	rootFile += ".root";
	TString pngFile(file);
	pngFile += ".eps";
	TFile* inFile = new TFile(rootFile);
	RooWorkspace* ws = (RooWorkspace*) inFile->Get("rws");
/*
	ws->factory("SUM::mm(xs*timeAngle,tBkg)");
        ws->factory("PROD::mmbkg(tBkg,angle)");
//	RooDataSet protoData("protoData","protoData", (RooDataSet*)ws->data("data"), RooArgSet(*ws->var("D"),*ws->var("et")) );

	RooDataSet* protoData = ws->pdf("errorBkg")->generate(*ws->var("et"), 1000);

        ws->factory("RSUM::tmptBkg(xr*resolution, xn*Decay(t,tauNeg,resolution,Flipped), xp*Decay(t,tauPos,resolution,SingleSided), Decay(t,tauPosPos,resolution,SingleSided)");
        ws->factory("PROD::tmptimeBkg(tBkg|et,errorBkg)");
	ws->factory("PROD::tmpbackground(massBkg,tmptimeBkg,angle)");
	ws->factory("SUM::tmpmodel(xs*signal,tmpbackground)");

//        ws->pdf("tmpbackground")->generate(RooArgSet(*ws->var("t"),*ws->var("cpsi"),*ws->var("ctheta"),*ws->var("phi")), RooFit::ProtoData(protoData) );
//        ws->pdf("angle")->generate(RooArgSet(*ws->var("cpsi"),*ws->var("ctheta"),*ws->var("phi")) );
//        ws->pdf("tmptimeBkg")->generate(RooArgSet(*ws->var("t")), RooFit::ProtoData(*protoData) );
*/

//	TPad *plotPad, *resPad;
//	TCanvas* canvas = new TCanvas("canvas", "canvas", 900,800);
//	canvas->Divide(3,2);


	TCanvas* canvasMass = new TCanvas("canvasMass", "canvas mass", 600,600);
//	canvas->cd(1);
//	plotPad = new TPad("plotPadMass", "Plot Pad", 0,0.2, 1,1);
//	resPad = new TPad("resPadMass", "Res Pad", 0,0, 1,0.2);
//	plotPad->Draw();
//	resPad->Draw();


	RooBinning binning(48,5.17,5.29);
//	binning.addUniform(50,5.44,5.57);

	ws->var("m")->setRange("left",5.17,5.29); 
	ws->var("m")->setRange("right",5.44,5.57); 

//	plotPad->cd();
	RooPlot *m_frame = ws->var("m")->frame();
	ws->data("dataBkg")->plotOn(m_frame, RooFit::MarkerSize(0.3), RooFit::Binning(binning));
	ws->pdf("background")->plotOn(m_frame, RooFit::NumCPU(nCPU), RooFit::Components("Prompt"), RooFit::LineColor(kBlack),RooFit::Range("left"),RooFit::NormRange("left") );
	ws->pdf("background")->plotOn(m_frame, RooFit::NumCPU(nCPU), RooFit::Components("NonPrompt"), RooFit::LineColor(kRed),RooFit::Range("left"), RooFit::NormRange("left"));
	ws->pdf("background")->plotOn(m_frame, RooFit::NumCPU(nCPU),RooFit::Range("left"), RooFit::NormRange("left"));
	RooHist *massResidualL = m_frame->pullHist();
	RooBinning binningR(52,5.44,5.57);
	ws->data("dataBkg")->plotOn(m_frame, RooFit::MarkerSize(0.3), RooFit::Binning(binningR));
	ws->pdf("background")->plotOn(m_frame, RooFit::NumCPU(nCPU), RooFit::Components("Prompt"), RooFit::LineColor(kBlack),RooFit::Range("right"),RooFit::NormRange("right") );
	ws->pdf("background")->plotOn(m_frame, RooFit::NumCPU(nCPU), RooFit::Components("NonPrompt"), RooFit::LineColor(kRed),RooFit::Range("right"), RooFit::NormRange("right"));
	ws->pdf("background")->plotOn(m_frame, RooFit::NumCPU(nCPU),RooFit::Range("right"), RooFit::NormRange("right"));
	RooHist *massResidualR = m_frame->pullHist();
	m_frame->SetTitle();
	m_frame->GetYaxis()->SetLabelSize(0.06);
	m_frame->GetYaxis()->SetTitleSize(0.06);
	m_frame->GetYaxis()->SetTitleOffset(1.3);
	m_frame->GetYaxis()->SetRangeUser(0,600);
	m_frame->GetXaxis()->SetNdivisions(505);
	m_frame->GetXaxis()->SetTitle("M(GeV)");
	m_frame->SetLabelSize(0.06);
        tex = new TLatex(5.37,500,"D0 Run II, 8fb^{-1}");
        tex->SetLineWidth(2);
        tex->Draw();
	m_frame->Draw();

/*
	resPad->cd();
//	RooHist *massResidual = m_frame->residHist();
        RooPlot *massResidualPlot  = ws->var("m")->frame();
	massResidualPlot->addPlotable(massResidualL,"E");
	massResidualPlot->addPlotable(massResidualR,"E");
	massResidualPlot->SetTitle();
	massResidualPlot->GetXaxis()->SetTitle();
	massResidualPlot->GetXaxis()->SetLabelSize(0);
	massResidualPlot->GetYaxis()->SetLabelSize(0.1);
	massResidualPlot->Draw();
*/


	TCanvas* canvasTime = new TCanvas("canvasTime", "canvas time", 600,600);
        canvasTime->Range(-2.726415,0.08873195,5.395178,4.216817);
        canvasTime->SetLeftMargin(0.1510067);
        canvasTime->SetRightMargin(0.04865772);
        canvasTime->SetTopMargin(0.06620209);
        canvasTime->SetBottomMargin(0.1341463);

//	canvas->cd(2);
/*	plotPad = new TPad("plotPadTime", "Plot Pad", 0,0.2, 1,1);
	resPad = new TPad("resPadTime", "Res Pad", 0,0, 1,0.2);
	plotPad->Draw();
	resPad->Draw();

	plotPad->cd();
*/
	ws->var("t")->setMin(-1.5);
        ws->var("t")->setMax( 5.0);

	RooPlot *t_frame = ws->var("t")->frame();
	ws->data("dataBkg")->plotOn(t_frame, RooFit::MarkerSize(0.3));
//	ws->pdf("timeBkg")->plotOn(t_frame);
	RooDataHist *projData = new RooDataHist("projData","projData",*ws->var("et"),*ws->data("dataBkg"));
	ws->pdf("background")->plotOn(t_frame, RooFit::ProjWData(*ws->var("et"), *projData), RooFit::Components("Prompt"), RooFit::LineColor(kBlack));
	ws->pdf("background")->plotOn(t_frame, RooFit::ProjWData(*ws->var("et"), *projData), RooFit::Components("NonPrompt"), RooFit::LineColor(kRed));
	ws->pdf("background")->plotOn(t_frame, RooFit::ProjWData(*ws->var("et"), *projData));
//	ws->pdf("tBkg")->plotOn(t_frame,  RooFit::NumCPU(nCPU));
	gPad->SetLogy(1);
	t_frame->SetTitle();
	t_frame->Draw();
	t_frame->GetYaxis()->SetLabelSize(0.06);
	t_frame->GetYaxis()->SetTitleSize(0.06);
	t_frame->GetXaxis()->SetTitle("t(ps)");
	t_frame->SetLabelSize(0.06);
	t_frame->GetYaxis()->SetTitleOffset(1.3);
	t_frame->GetXaxis()->SetNdivisions(505);
        TLatex *tex = new TLatex(2.138365,3836.544,"D0 Run II, 8fb^{-1}");
        tex->SetLineWidth(2);
        tex->Draw();

	canvasTime.SaveAs("time.C");	

/*
	resPad->cd();
	RooHist *timeResidual = t_frame->pullHist();
        RooPlot *timeResidualPlot  = ws->var("t")->frame(RooFit::Title("Residual Distribution"));
	timeResidualPlot->addPlotable(timeResidual,"E");
	timeResidualPlot->SetTitle();
	timeResidualPlot->GetXaxis()->SetTitle();
	timeResidualPlot->GetXaxis()->SetLabelSize(0);
	timeResidualPlot->GetYaxis()->SetLabelSize(0.1);
	timeResidualPlot->Draw();
*/


	TCanvas* canvasEt = new TCanvas("canvasEt", "canvas et", 600,600);
//	canvas->cd(3);
/*
	plotPad = new TPad("plotPadEt", "Plot Pad", 0,0.2, 1,1);
	resPad = new TPad("resPadEt", "Res Pad", 0,0, 1,0.2);
	plotPad->Draw();
	resPad->Draw();
*/

//	plotPad->cd();
        ws->var("et")->setMax(0.6);
	gPad->SetLogy(0);
	RooPlot *et_frame = ws->var("et")->frame();
	ws->data("dataBkg")->plotOn(et_frame,RooFit::MarkerSize(0.2));
	ws->pdf("background")->plotOn(et_frame, RooFit::NumCPU(nCPU), RooFit::Components("Prompt"), RooFit::LineColor(kBlack));
	ws->pdf("background")->plotOn(et_frame, RooFit::NumCPU(nCPU), RooFit::Components("NonPrompt"), RooFit::LineColor(kRed));
	ws->pdf("background")->plotOn(et_frame, RooFit::NumCPU(nCPU));
	gPad->SetLogy(1);
	et_frame->SetTitle();
	et_frame->GetYaxis()->SetLabelSize(0.06);
	et_frame->GetYaxis()->SetTitleSize(0.06);
	et_frame->GetXaxis()->SetTitle("#sigma(t) ps");
	et_frame->SetLabelSize(0.06);
	et_frame->GetYaxis()->SetTitleOffset(1.3);
//	et_frame->GetYaxis()->SetRangeUser(0,600);
	et_frame->GetXaxis()->SetNdivisions(505);
        tex = new TLatex(0,0,"D0 Run II, 8fb^{-1}");
        tex->SetLineWidth(2);
        tex->Draw();
	et_frame->Draw();

/*
	resPad->cd();
	RooHist *etResidual = et_frame->pullHist();
        RooPlot *etResidualPlot  = ws->var("et")->frame(RooFit::Title("Residual Distribution"));
	etResidualPlot->addPlotable(etResidual,"E");
	etResidualPlot->SetTitle();
	etResidualPlot->GetXaxis()->SetTitle();
	etResidualPlot->GetXaxis()->SetLabelSize(0);
	etResidualPlot->GetYaxis()->SetLabelSize(0.1);
	etResidualPlot->Draw();
*/

	TCanvas* canvasCpsi = new TCanvas("canvasCpsi", "canvas Cpsi", 600,600);
//	canvas->cd(4);
/*
	plotPad = new TPad("plotPadCpsi", "Plot Pad", 0,0.2, 1,1);
	resPad = new TPad("resPadCpsi", "Res Pad", 0,0, 1,0.2);
	plotPad->Draw();
	resPad->Draw();
*/


//	plotPad->cd();
	gPad->SetLogy(0);
	RooPlot *cpsi_frame = ws->var("cpsi")->frame();
	ws->data("dataBkg")->plotOn(cpsi_frame,RooFit::MarkerSize(0.2));
	ws->pdf("background")->plotOn(cpsi_frame, RooFit::NumCPU(nCPU), RooFit::Components("Prompt"), RooFit::LineColor(kBlack));
	ws->pdf("background")->plotOn(cpsi_frame, RooFit::NumCPU(nCPU), RooFit::Components("NonPrompt"), RooFit::LineColor(kRed));
	ws->pdf("background")->plotOn(cpsi_frame, RooFit::NumCPU(nCPU));
	cpsi_frame->SetTitle();
	cpsi_frame->GetYaxis()->SetLabelSize(0.06);
	cpsi_frame->GetYaxis()->SetTitleSize(0.06);
	cpsi_frame->GetXaxis()->SetTitle("cos(#psi)");
	cpsi_frame->SetLabelSize(0.06);
	cpsi_frame->GetYaxis()->SetTitleOffset(1.3);
	cpsi_frame->GetYaxis()->SetRangeUser(0,600);
	cpsi_frame->GetXaxis()->SetNdivisions(505);
        tex = new TLatex(0,500,"D0 Run II, 8fb^{-1}");
        tex->SetLineWidth(2);
        tex->Draw();
	cpsi_frame->Draw();

/*
	resPad->cd();
	RooHist *cpsiResidual = cpsi_frame->pullHist();
        RooPlot *cpsiResidualPlot  = ws->var("cpsi")->frame(RooFit::Title("Residual Distribution"));
	cpsiResidualPlot->addPlotable(cpsiResidual,"E");
	cpsiResidualPlot->SetTitle();
	cpsiResidualPlot->GetXaxis()->SetTitle();
	cpsiResidualPlot->GetXaxis()->SetLabelSize(0);
	cpsiResidualPlot->GetYaxis()->SetLabelSize(0.1);
	cpsiResidualPlot->Draw();
*/


	TCanvas* canvasCtheta = new TCanvas("canvasCtheta", "canvas Ctheta", 600,600);
/*
//	canvas->cd(5);
	plotPad = new TPad("plotPadCtheta", "Plot Pad", 0,0.2, 1,1);
	resPad = new TPad("resPadCtheta", "Res Pad", 0,0, 1,0.2);
	plotPad->Draw();
	resPad->Draw();
*/
//	plotPad->cd();
	RooPlot *ctheta_frame = ws->var("ctheta")->frame();
	ws->data("dataBkg")->plotOn(ctheta_frame,RooFit::MarkerSize(0.2));
	ws->pdf("background")->plotOn(ctheta_frame, RooFit::NumCPU(nCPU), RooFit::Components("Prompt"), RooFit::LineColor(kBlack));
	ws->pdf("background")->plotOn(ctheta_frame, RooFit::NumCPU(nCPU), RooFit::Components("NonPrompt"), RooFit::LineColor(kRed));
	ws->pdf("background")->plotOn(ctheta_frame, RooFit::NumCPU(nCPU));
	ctheta_frame->SetTitle();
	ctheta_frame->GetYaxis()->SetLabelSize(0.06);
	ctheta_frame->GetYaxis()->SetTitleSize(0.06);
	ctheta_frame->GetXaxis()->SetTitle("cos(#theta)");
	ctheta_frame->SetLabelSize(0.06);
	ctheta_frame->GetYaxis()->SetTitleOffset(1.3);
	ctheta_frame->GetYaxis()->SetRangeUser(0,600);
	ctheta_frame->GetXaxis()->SetNdivisions(505);
        tex = new TLatex(0,500,"D0 Run II, 8fb^{-1}");
        tex->SetLineWidth(2);
        tex->Draw();
	ctheta_frame->Draw();

/*
	resPad->cd();
	RooHist *cthetaResidual = ctheta_frame->pullHist();
        RooPlot *cthetaResidualPlot  = ws->var("ctheta")->frame(RooFit::Title("Residual Distribution"));
	cthetaResidualPlot->addPlotable(cthetaResidual,"E");
	cthetaResidualPlot->SetTitle();
	cthetaResidualPlot->GetXaxis()->SetTitle();
	cthetaResidualPlot->GetXaxis()->SetLabelSize(0);
	cthetaResidualPlot->GetYaxis()->SetLabelSize(0.1);
	cthetaResidualPlot->Draw();
*/


	TCanvas* canvasPhi = new TCanvas("canvasPhi", "canvas phi", 600,600);
/*
//	canvas->cd(6);
	plotPad = new TPad("plotPadPhi", "Plot Pad", 0,0.2, 1,1);
	resPad = new TPad("resPadPhi", "Res Pad", 0,0, 1,0.2);
	plotPad->Draw();
	resPad->Draw();
*/

//	plotPad->cd();
	RooPlot *phi_frame = ws->var("phi")->frame();
	ws->data("dataBkg")->plotOn(phi_frame,RooFit::MarkerSize(0.2));
	ws->pdf("background")->plotOn(phi_frame, RooFit::NumCPU(nCPU), RooFit::Components("Prompt"), RooFit::LineColor(kBlack));
	ws->pdf("background")->plotOn(phi_frame, RooFit::NumCPU(nCPU), RooFit::Components("NonPrompt"), RooFit::LineColor(kRed));
	ws->pdf("background")->plotOn(phi_frame, RooFit::NumCPU(nCPU));
	phi_frame->SetTitle();
	phi_frame->GetYaxis()->SetLabelSize(0.06);
	phi_frame->GetYaxis()->SetTitleSize(0.06);
	phi_frame->GetXaxis()->SetTitle("#varphi (rad)");
	phi_frame->SetLabelSize(0.06);
	phi_frame->GetYaxis()->SetTitleOffset(1.3);
	phi_frame->GetYaxis()->SetRangeUser(0,600);
	phi_frame->GetXaxis()->SetNdivisions(505);
        tex = new TLatex(0,500,"D0 Run II, 8fb^{-1}");
        tex->SetLineWidth(2);
        tex->Draw();
	phi_frame->Draw();

/*
	resPad->cd();
	RooHist *phiResidual = phi_frame->pullHist();
        RooPlot *phiResidualPlot  = ws->var("phi")->frame(RooFit::Title("Residual Distribution"));
	phiResidualPlot->addPlotable(phiResidual,"E");
	phiResidualPlot->SetTitle();
	phiResidualPlot->GetXaxis()->SetTitle();
	phiResidualPlot->GetXaxis()->SetLabelSize(0);
	phiResidualPlot->GetYaxis()->SetLabelSize(0.1);
	phiResidualPlot->Draw();
*/

//	canvas->SaveAs(pngFile);

}
