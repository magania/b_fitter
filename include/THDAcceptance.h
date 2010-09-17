#ifndef THDACCEPTANCE
#define THDACCEPTANCE

#include <stdlib.h>
#include <iostream>
#include <TH2D.h>

class THDAcceptance : public TH2D {
public:
	THDAcceptance(const char* name, const char* title, Int_t nbins_cos0l, const Double_t* cos0l_bins, Int_t nbins_phi, const Double_t* phi_bins);
	Double_t acceptance(Double_t cos0l, Double_t cos0k, Double_t phi) const;
};

#endif
