fitMass(){
	gROOT->SetStyle("Plain");
	const char* chInFile = "dataset.root";

	TFile* inFile = new TFile(chInFile);
	RooWorkspace* ws = (RooWorkspace*) inFile.Get("rws");

	ws->factory("Gaussian::massSignal(m,mu[5.2,5.4],sigmaM[0.03,0.02,0.04]");
	ws->factory("Polynomial::massBkg(m,{slope[-0.125,-1,-0.1]})");

	ws->factory("SUM::model(nSig[5000,500,15000]*massSignal, nBkg[200000,10000,1000000]*massBkg)");

	ws->pdf("model")->fitTo(*ws->data("data"),RooFit::Extended());

	TCanvas* canvasMass = new TCanvas("canvasMass", "canvas mass", 400,600);
	plotPad = new TPad("plotPadMass", "Plot Pad", 0,0.2, 1,1);
	resPad = new TPad("resPadMass", "Res Pad", 0,0, 1,0.2);
	plotPad->Draw();
	resPad->Draw();

	plotPad->cd();
	RooPlot *m_frame = ws->var("m")->frame();
	ws->data("data")->plotOn(m_frame, RooFit::MarkerSize(0.3));
	ws->pdf("model")->plotOn(m_frame, RooFit::NumCPU(2));
	ws->pdf("model")->paramOn(m_frame);
	m_frame->Draw();

	resPad->cd();
	RooHist *massResidual = m_frame->pullHist();
        RooPlot *massResidualPlot  = ws->var("m")->frame(RooFit::Title("Residual Distribution"));
	massResidualPlot->addPlotable(massResidual,"E");
	massResidualPlot->Draw();	

}
