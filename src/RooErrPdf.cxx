 /***************************************************************************** 
  * Project: RooFit                                                           * 
  *                                                                           * 
  * This code was autogenerated by RooClassFactory                            * 
  *****************************************************************************/ 

 // Your description goes here... 

 #include "Riostream.h" 

 #include "RooErrPdf.h" 
 #include "RooAbsReal.h" 
 #include "RooAbsCategory.h" 
 #include <math.h> 
 #include "TMath.h" 

 ClassImp(RooErrPdf) 

 RooErrPdf::RooErrPdf(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _a,
                        RooAbsReal& _b,
                        RooAbsReal& _s) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   a("a","a",this,_a),
   b("b","b",this,_b),
   s("s","s",this,_s)
 { 
//   std::cout << "ErrModel 1" << std::endl;
 } 


 RooErrPdf::RooErrPdf(const RooErrPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   a("a",this,other.a),
   b("b",this,other.b),
   s("s",this,other.s)
 { 
// std::cout << "ErrModel 2" << std::endl;
 } 

 Double_t RooErrPdf::evaluate() const 
 { 
//   std::cout << "ErrModel evaluate" << std::endl;
//   std::cout << x << " " << TMath::Power(x,a)*TMath::Exp(-x/b) << std::endl;
   if (x<s) return 0;
   return TMath::Power(x-s,a)*TMath::Exp(-(x-s)/b) ;
 } 

Int_t RooErrPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{
//  std::cout << "Analitical Integral ErrModel" << std::endl;
  if (matchArgs(allVars,analVars,x)) return 1 ;
  return 0 ;
}


//_____________________________________________________________________________
Double_t RooErrPdf::analyticalIntegral(Int_t code, const char* rangeName) const
{
//  std::cout << "ErrModel integral" << std::endl;
  const Double_t xmin = x.min();
  const Double_t xmax = x.max();

/*
  if(xmin != 0){
       std::cout << "Range Error xmin should be 0 " << xmin <<std::endl;
       exit(1);
  }
*/

//  std::cout << xmax << " " << TMath::Power(b,a+1)*TMath::Gamma(a+1,xmax/b)*TMath::Gamma(a+1) << std::endl;
//  return TMath::Power(b,a+1)*TMath::Gamma(a+1,xmax/b)*TMath::Gamma(a+1);
  return TMath::Power(b,a+1)*TMath::Gamma(a+1);
}

