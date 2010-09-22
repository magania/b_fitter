/*
 *  Author: Ricardo Magana-Villalba
 *          magania@fnal.gov
 *
 *  September 2010
 */


#ifndef ROO_BS_TIME_ANGLE
#define ROO_BS_TIME_ANGLE

#include "THDAcceptance.h"

#include "TRandom3.h"
#include "TComplex.h"
#include "Math/SpecFuncMathMore.h"

#include "RooAbsAnaConvPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"

class RooBsTimeAngle : public RooAbsAnaConvPdf {
public:
    RooBsTimeAngle(const char *name, const char *title,
            RooRealVar& t,
            RooRealVar& cpsi,
            RooRealVar& ctheta,
            RooRealVar& phi,
            RooRealVar& A02,
            RooRealVar& All2,
            RooRealVar& Ap2,
            RooRealVar& tau_L,
            RooRealVar& tau_H,
            RooRealVar& Dm,
            RooRealVar& beta,
            RooRealVar& delta_p,
            RooRealVar& delta_l,
            RooRealVar& delta_s,
            RooRealVar& Fs,
            RooRealVar& D,
            const RooResolutionModel& model);

    RooBsTimeAngle(const RooBsTimeAngle& other, const char* name = 0);

    virtual TObject* clone(const char* newname) const ;
    virtual ~RooBsTimeAngle();

    virtual Double_t coefficient(Int_t basisIndex) const;

    virtual Int_t getCoefAnalyticalIntegral(Int_t coef, RooArgSet& allVars, RooArgSet& analVars, const char* rangeName = 0) const;
    virtual Double_t coefAnalyticalIntegral(Int_t coef, Int_t code, const char* rangeName = 0) const;

    Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK = kTRUE) const;
    void initGenerator(Int_t code);
    void generateEvent(Int_t code);

private:

Double_t HarmonicSphericalY(int l, int m, Double_t ctheta, Double_t phi) const;
Double_t acceptance() const;

Double_t rho_B (Int_t bi) const;    //8.19
Double_t rho_Bbar (Int_t bi) const;    //8.20
Double_t P_B (Int_t bi) const;    //5.5
Double_t P_Bbar (Int_t bi) const;    //5.6
Double_t Q_B (Int_t bi) const;    //8.15
Double_t Q_Bbar (Int_t bi) const;    //8.16

TComplex Aplus(int i) const;    // 5.2
TComplex Aminus(int i) const;    // 5.3
TComplex B(int i) const;    // 8.16
Double_t nhat(int i) const;    // 3.6

TComplex CrossDot(TComplex (RooBsTimeAngle::*A)(int)const ,TComplex (RooBsTimeAngle::*B)(int)const ) const;

Double_t f_sq(Double_t A, Double_t B, int bi) const;    //5.7  5.8
Double_t f_bar_plus_sq(int bi) const;    //5.7
Double_t f_bar_minus_sq(int bi) const;    //5.7
Double_t f_plus_sq(int bi) const;    //5.8
Double_t f_minus_sq(int bi) const;    //5.8
TComplex f_f_star(Double_t A, int bi) const;    //5.9 5.10
TComplex f_f_star_bar(int bi) const;    //5.9
TComplex f_f_star(int bi) const;    //5.10

TComplex I_mu() const;    //8.10


protected:
    RooRealProxy _t;
    RooRealProxy _cpsi;
    RooRealProxy _ctheta;
    RooRealProxy _phi;
    RooRealProxy _A02;
    RooRealProxy _All2;
    RooRealProxy _Ap2;
    RooRealProxy _tau_L;
    RooRealProxy _tau_H;
    RooRealProxy _Dm;
    RooRealProxy _beta;
    RooRealProxy _delta_p;
    RooRealProxy _delta_l;
    RooRealProxy _delta_s;
    RooRealProxy _Fs;
    RooRealProxy _D;

    RooFormulaVar _tau;

    TRandom3 _aleatorio;

    Int_t _basisExpGammaL;
    Int_t _basisExpGammaH;
    Int_t _basisExpCosDm;
    Int_t _basisExpSinDm;

	static const Double_t e_00p0 = 2.37858;
	static const Double_t e_01n1 = -0.00191028;
	static const Double_t e_01p0 = 0.00635797;
	static const Double_t e_01p1 = -0.0119643;
	static const Double_t e_02n2 = 2.71959e-06;
	static const Double_t e_02n1 = -0.00136709;
	static const Double_t e_02p0 = 0.215535;
	static const Double_t e_02p1 = -0.0313594;
	static const Double_t e_02p2 = -0.322357;

private:
	static const bool __debug = false;

//    ClassDef(RooBsTimeAngle, 1) // B0s Time and Angular decay PDF
};

#endif
