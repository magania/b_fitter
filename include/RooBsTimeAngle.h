/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROOBSTIMEANGLE
#define ROOBSTIMEANGLE

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

#include "TComplex.h"
#include "Math/SpecFuncMathMore.h"
#include "TRandom3.h"
 
class RooBsTimeAngle : public RooAbsPdf {
public:
  RooBsTimeAngle() {} ; 
  RooBsTimeAngle(const char *name, const char *title,
	      RooAbsReal& _t,
	      RooAbsReal& _et,
	      RooAbsReal& _cpsi,
	      RooAbsReal& _ctheta,
	      RooAbsReal& _phi,
	      RooAbsReal& _D,
	      RooAbsReal& _a0,
	      RooAbsReal& _al,
	      RooAbsReal& _ap,
	      RooAbsReal& _DG,
	      RooAbsReal& _Dm,
	      RooAbsReal& _tau,
	      RooAbsReal& _beta,
	      RooAbsReal& _delta_l,
	      RooAbsReal& _delta_p,
	      RooAbsReal& _delta_s,
	      RooAbsReal& _fs,
	      RooAbsReal& _s);

  RooBsTimeAngle(const char *name, const char *title,
              RooAbsReal& _t,
              RooAbsReal& _cpsi,
              RooAbsReal& _ctheta,
              RooAbsReal& _phi,
              RooAbsReal& _D,
              RooAbsReal& _a0,
              RooAbsReal& _al,
              RooAbsReal& _ap,
              RooAbsReal& _DG,
              RooAbsReal& _Dm,
              RooAbsReal& _tau,
              RooAbsReal& _beta,
              RooAbsReal& _delta_l,
              RooAbsReal& _delta_p,
              RooAbsReal& _delta_s,
              RooAbsReal& _fs);

  RooBsTimeAngle(const RooBsTimeAngle& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooBsTimeAngle(*this,newname); }
  inline virtual ~RooBsTimeAngle() { }

protected:

  RooRealProxy t ;
  RooRealProxy et ;
  RooRealProxy cpsi ;
  RooRealProxy ctheta ;
  RooRealProxy phi ;
  RooRealProxy D ;
  RooRealProxy a0 ;
  RooRealProxy al ;
  RooRealProxy ap ;
  RooRealProxy DG ;
  RooRealProxy Dm ;
  RooRealProxy tau ;
  RooRealProxy beta ;
  RooRealProxy delta_l ;
  RooRealProxy delta_p ;
  RooRealProxy delta_s ;
  RooRealProxy fs ;
  RooRealProxy s ;

 TRandom3 _aleatorio;
  
  Double_t evaluate() const ;
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /* rangeName*/) const ;
  Double_t analyticalIntegral(Int_t code, const char* range) const ;
  Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK) const ;
  void initGenerator(Int_t code);
  void generateEvent(Int_t code);

private:

Double_t HarmonicSphericalY(int l, int m, Double_t ctheta, Double_t phi) const;
void acceptance(Double_t& acceptance, Double_t ctheta, Double_t phi) const;
void rho_Bb(Double_t& rho_Bb, Double_t fs, Double_t& P_Bb, Double_t& Q_Bb, TComplex& CDAminusB, TComplex& CDAplusB, Double_t& fb_minus_sq, TComplex &fb_fb_star, TComplex& I_mu) const;    //8.19 8.20
void P_Bb(Double_t &P_Bb, Double_t &CDAplusAplus, Double_t &CDAminusAminus, TComplex &CDAplusAminus, Double_t &fb_plus_sq, Double_t &fb_minus_sq, TComplex &fb_fb_star) const ;    //5.5 5.6
void Q_Bb(Double_t &Q_Bb, Double_t &CDBB, Double_t &fb_minus_sq) const ;    //8.15 8.16
void AplusK(TComplex* AplusK, Double_t a0, Double_t al, Double_t delta_l) const ;    // 5.2
void AminusK(TComplex& AminusK, Double_t ap, Double_t delta_p) const;    // 5.3
void nhat(Double_t* nhat, Double_t ctheta, Double_t phi) const;    // 3.6
void CrossDot(TComplex &CD, TComplex* A, TComplex* B, Double_t* n) const;
void I_mu(TComplex &I_mu, Double_t fs, Double_t delta_s) const;    //8.10

TComplex Erf(TComplex z) const;
TComplex Erf(double x0, double y0) const;

static const Double_t e_00p0 = 2.78142;
static const Double_t e_01n1 = -0.00650078;
static const Double_t e_01p0 = 0.00290784;
static const Double_t e_01p1 = -0.000100291;
static const Double_t e_02n2 = 0.00157789;
static const Double_t e_02n1 = -0.00239414;
static const Double_t e_02p0 = 0.130073;
static const Double_t e_02p1 = 0.0139719;
static const Double_t e_02p2 = -0.211069;

bool _resolution;

static const bool __debug = false;
static const bool __debug2 = false;
static const bool __debug3 = false;
static const bool __debug4 = false;
static const bool __fast = true;

//  ClassDef(RooBsTimeAngle,1) // Your description goes here...
};
 
#endif