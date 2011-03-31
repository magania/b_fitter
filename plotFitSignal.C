void plotFitSignal(const char* file){
    cout << "Printing ... " << endl;
	gROOT->SetStyle("Plain");

	gSystem->Load("lib/libBFitter.so");

	Int_t nCPU = 4;

    TString rootFile(file);
	rootFile += ".root";
	TString pngFile(file);
	pngFile += ".png";
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

	TPad *plotPad, *resPad;
	TCanvas* canvas = new TCanvas("canvas", "canvas", 900,800);
	canvas->Divide(3,2);


//	TCanvas canvasMass("canvasMass", "canvas mass", 400,600);
	canvas->cd(1);
	plotPad = new TPad("plotPadMass", "Plot Pad", 0,0.2, 1,1);
	resPad = new TPad("resPadMass", "Res Pad", 0,0, 1,0.2);
	plotPad->Draw();
	resPad->Draw();

	plotPad->cd();
	RooPlot *m_frame = ws->var("m")->frame();
	ws->data("data")->plotOn(m_frame, RooFit::MarkerSize(0.3));
	ws->data("dataGen")->plotOn(m_frame, RooFit::Rescale(0.05), RooFit::LineColor(kBlue), RooFit::DrawOption("COLZX"));
	m_frame->Draw();

	resPad->cd();

//	TCanvas canvasTime("canvasTime", "canvas time", 400,600);
	canvas->cd(2);
	plotPad = new TPad("plotPadTime", "Plot Pad", 0,0.2, 1,1);
	resPad = new TPad("resPadTime", "Res Pad", 0,0, 1,0.2);
	plotPad->Draw();
	resPad->Draw();

	plotPad->cd();
	ws->var("t")->setMin(-1.5);
    ws->var("t")->setMax( 5.0);

	RooPlot *t_frame = ws->var("t")->frame();
	ws->data("data")->plotOn(t_frame, RooFit::MarkerSize(0.3));
	ws->data("dataGen")->plotOn(t_frame, RooFit::Rescale(0.05), RooFit::LineColor(kBlue), RooFit::DrawOption("COLZX"));
	gPad->SetLogy(1);
	t_frame->Draw();

	resPad->cd();


//	TCanvas canvasEt("canvasEt", "canvas et", 400,600);
	canvas->cd(3);
	plotPad = new TPad("plotPadEt", "Plot Pad", 0,0.2, 1,1);
	resPad = new TPad("resPadEt", "Res Pad", 0,0, 1,0.2);
	plotPad->Draw();
	resPad->Draw();

	plotPad->cd();
    ws->var("et")->setMax(0.6);
	gPad->SetLogy(0);
	RooPlot *et_frame = ws->var("et")->frame();
	ws->data("data")->plotOn(et_frame, RooFit::MarkerSize(0.3));
	ws->data("dataGen")->plotOn(et_frame, RooFit::Rescale(0.05), RooFit::LineColor(kBlue), RooFit::DrawOption("COLZX"));
	gPad->SetLogy(1);
	et_frame->Draw();

	resPad->cd();

//	TCanvas canvasCpsi("canvasCpsi", "canvas Cpsi", 400,600);
	canvas->cd(4);
	plotPad = new TPad("plotPadCpsi", "Plot Pad", 0,0.2, 1,1);
	resPad = new TPad("resPadCpsi", "Res Pad", 0,0, 1,0.2);
	plotPad->Draw();
	resPad->Draw();

	plotPad->cd();
	gPad->SetLogy(0);
	RooPlot *cpsi_frame = ws->var("cpsi")->frame();
	ws->data("data")->plotOn(cpsi_frame, RooFit::MarkerSize(0.3));
	ws->data("dataGen")->plotOn(cpsi_frame, RooFit::Rescale(0.05), RooFit::LineColor(kBlue), RooFit::DrawOption("COLZX"));
	cpsi_frame->Draw();

	resPad->cd();


//	TCanvas canvasCtheta("canvasCtheta", "canvas Ctheta", 400,600);
	canvas->cd(5);
	plotPad = new TPad("plotPadCtheta", "Plot Pad", 0,0.2, 1,1);
	resPad = new TPad("resPadCtheta", "Res Pad", 0,0, 1,0.2);
	plotPad->Draw();
	resPad->Draw();

	plotPad->cd();
	RooPlot *ctheta_frame = ws->var("ctheta")->frame();
	ws->data("data")->plotOn(ctheta_frame, RooFit::MarkerSize(0.3));
	ws->data("dataGen")->plotOn(ctheta_frame, RooFit::Rescale(0.05), RooFit::LineColor(kBlue), RooFit::DrawOption("COLZX"));
	ctheta_frame->Draw();

	resPad->cd();


//	TCanvas canvasPhi("canvasPhi", "canvas phi", 400,600);
	canvas->cd(6);
	plotPad = new TPad("plotPadPhi", "Plot Pad", 0,0.2, 1,1);
	resPad = new TPad("resPadPhi", "Res Pad", 0,0, 1,0.2);
	plotPad->Draw();
	resPad->Draw();

	plotPad->cd();
	RooPlot *phi_frame = ws->var("phi")->frame();
	ws->data("data")->plotOn(phi_frame, RooFit::MarkerSize(0.3));
	ws->data("dataGen")->plotOn(phi_frame, RooFit::Rescale(0.05), RooFit::LineColor(kBlue), RooFit::DrawOption("COLZX"));
	phi_frame->Draw();

	resPad->cd();

	canvas->SaveAs(pngFile);

}
