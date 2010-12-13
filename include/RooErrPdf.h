/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROOERRPDF
#define ROOERRPDF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class RooErrPdf : public RooAbsPdf {
public:
  RooErrPdf() {} ; 
  RooErrPdf(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _a,
	      RooAbsReal& _b,
	      RooAbsReal& _s);
  RooErrPdf(const RooErrPdf& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooErrPdf(*this,newname); }
  inline virtual ~RooErrPdf() { }
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName) const ;
  ClassDef(RooErrPdf,1) // PDL error PDF

protected:

  RooRealProxy x ;
  RooRealProxy a ;
  RooRealProxy b ;
  RooRealProxy s ;
  
  Double_t evaluate() const ;


};
 
#endif