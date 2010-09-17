#include "THDAcceptance.h"


THDAcceptance::THDAcceptance(const char* name, const char* title, Int_t nbins_cos0l, const Double_t* cos0l_bins, Int_t nbins_phi, const Double_t* phi_bins) :
  TH2D(name, title, nbins_cos0l, cos0l_bins, nbins_phi, phi_bins)
{
	
}


Double_t THDAcceptance::acceptance(Double_t cos0l, Double_t cos0k, Double_t phi) const {
	using namespace std;

	Int_t binx, biny, bin;
	binx = fXaxis.FindFixBin(cos0l);
	biny = fYaxis.FindFixBin(phi);
	if (binx <0 || biny <0) {
		cout << "ERROR: Bin not found. cos0l=" << cos0l << " phi=" << phi << endl;
		exit(1);
	};

	return GetBinContent(binx,biny);
}
