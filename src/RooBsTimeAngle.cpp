/*
 *  Author: Ricardo Magana-Villalba
 *          magania@fnal.gov
 *
 *  September 2010
 */


#include <RooAbsRealLValue.h>

#include "RooFit.h"
#include "Riostream.h"
#include "TMath.h"
#include "RooRealVar.h"
#include "RooRandom.h"
#include "TRandom3.h"

#include "RooBsTimeAngle.hpp"

//ClassImp(RooBsTimeAngle)
//;

//_____________________________________________________________________________
RooBsTimeAngle::RooBsTimeAngle(const char *name, const char *title,
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
        const RooResolutionModel& model) :
  RooAbsAnaConvPdf(name, title, model, t),
  _t("_t", "time", this, t),
  _cpsi("_cpsi","Cos(#Psi)",this,cpsi),
  _ctheta("_ctheta","Cos(#Theta)",this,ctheta),
  _phi("_phi","#Phi",this,phi),
  _A02("_A02", "|A_0(0)|^2", this, A02),
  _All2("_All2", "|A_#parallell(0)|^2", this, All2),
  _Ap2("_All2", "|A_#perpendicular(0)|^2", this, Ap2),
  _tau_L("_tau_L", "#tau_{L}", this, tau_L),
  _tau_H("_tau_H", "#tau_{H}", this, tau_H),
  _Dm("_Dm", "#Delta m", this, Dm),
  _beta("_beta", "#beta", this, beta),
  _delta_p("_delta_p", "#delta_{#perpendicular}", this, delta_p),
  _delta_l("_delta_l", "#delta_{#parallell}", this, delta_l),
  _delta_s("_delta_s", "#delta_{s}", this, delta_s),
  _Fs("_Fs", "F_{s}", this, Fs),
  _D("_D", "Dilution", this, D),
  _aleatorio(0),
  _tau("tau","tau","2.0*@0*@1/(@0+@1)",RooArgList(tau_L,tau_H))
{
  // Constructor
  cout << "** Bs Time Anlgle **" << endl;
  _basisExpGammaL = declareBasis("exp(-@0/@1)", RooArgList(tau_L)) ;
  _basisExpGammaH = declareBasis("exp(-@0/@1)", RooArgList(tau_H)) ;
  _basisExpCosDm  = declareBasis("exp(-@0/@1)*cos(@0*@2)", RooArgList(_tau,Dm)) ;
  _basisExpSinDm  = declareBasis("exp(-@0/@1)*sin(@0*@2)", RooArgList(_tau,Dm)) ;
}

//_____________________________________________________________________________
RooBsTimeAngle::RooBsTimeAngle(const RooBsTimeAngle& other, const char* name) :
  RooAbsAnaConvPdf(other,name),
  _t("_t", this, other._t),
  _cpsi("_cpsi", this, other._cpsi),
  _ctheta("_ctheta", this, other._ctheta),
  _phi("_phi", this, other._phi),
  _A02("_A02", this, other._A02),
  _All2("_All2", this, other._All2),
  _Ap2("_Ap2", this, other._Ap2),
  _tau_L("_tau_L", this, other._tau_L),
  _tau_H("_tau_H", this, other._tau_H),
  _Dm("_Dm", this, other._Dm),
  _beta("_beta", this, other._beta),
  _delta_p("_delta_p", this, other._delta_p),
  _delta_l("_delta_l", this, other._delta_l),
  _delta_s("_delta_s", this, other._delta_s),
  _Fs("_Fs", this, other._Fs),
  _D("_D",this, other._D),
  _aleatorio(other._aleatorio),
  _tau(other._tau),
  _basisExpGammaL(other._basisExpGammaL),
  _basisExpGammaH(other._basisExpGammaH),
  _basisExpCosDm(other._basisExpCosDm),
  _basisExpSinDm(other._basisExpSinDm)
{
  cout << "** Bs Time Anlgle Clone **" << endl;
}

//_____________________________________________________________________________
RooBsTimeAngle::~RooBsTimeAngle()
{
  // Destructor
}

TObject* RooBsTimeAngle::clone(const char* newname) const {
        return new RooBsTimeAngle(*this, newname);
}

//_____________________________________________________________________________
Double_t RooBsTimeAngle::coefficient(Int_t basisIndex) const {
	if (__debug) {
		cout << "acceptance = " << acceptance() << endl;
		cout << "rho_B(" << basisIndex << ") = " << rho_B(basisIndex) << endl;
	}
	return  ( (1-_D/2.0)*rho_B(basisIndex) + (1+_D/2.0)*rho_Bbar(basisIndex) ) * acceptance();
}

Double_t RooBsTimeAngle::HarmonicSphericalY(int l, int m, Double_t ctheta, Double_t phi) const {
        if (m==0)
                return TMath::Sqrt( ((2.0*l+1)*TMath::Gamma(l+1))/(4.0*TMath::Pi()*TMath::Gamma(l+1)) )*ROOT::Math::assoc_legendre(l,0,ctheta);
        if (m>0)
                return TMath::Sqrt( ((2.0*l+1)*TMath::Gamma(l-m+1))/(2.0*TMath::Pi()*TMath::Gamma(l+m+1)) )*ROOT::Math::assoc_legendre(l,m,ctheta)*TMath::Cos(m*phi)*TMath::Power(-1.0,m);
        if (m<0)
                return TMath::Sqrt( ((2.0*l+1)*TMath::Gamma(l+m+1))/(2.0*TMath::Pi()*TMath::Gamma(l-m+1)) )*ROOT::Math::assoc_legendre(l,-m,ctheta)*TMath::Sin(m*phi);
}

Double_t RooBsTimeAngle::acceptance() const {
/*
	if (__debug){
		cout << "00 0 = " << e_00p0*HarmonicSphericalY(0,0 ,_ctheta,_phi) << endl;
		cout << "01-1 = " << e_01n1*HarmonicSphericalY(1,-1,_ctheta,_phi) << endl; 
		cout << "01 0 = " << e_01p0*HarmonicSphericalY(1,0 ,_ctheta,_phi) << endl;
		cout << "01 1 = " << e_01p1*HarmonicSphericalY(1,1 ,_ctheta,_phi) << endl;
		cout << "02-2 = " << e_02n2*HarmonicSphericalY(2,-2,_ctheta,_phi) << endl;
		cout << "02-1 = " << e_02n1*HarmonicSphericalY(2,-1,_ctheta,_phi) << endl;
		cout << "02 0 = " << e_02p0*HarmonicSphericalY(2,0 ,_ctheta,_phi) << endl;
		cout << "02 1 = " << e_02p1*HarmonicSphericalY(2,1 ,_ctheta,_phi) << endl;
		cout << "02 2 = " << e_02p2*HarmonicSphericalY(2,2 ,_ctheta,_phi) << endl;
	}
*/

        Double_t result=0;

	result += e_00p0*HarmonicSphericalY(0,0 ,_ctheta,_phi);
	
	result += e_01n1*HarmonicSphericalY(1,-1,_ctheta,_phi);
	result += e_01p0*HarmonicSphericalY(1,0 ,_ctheta,_phi);
	result += e_01p1*HarmonicSphericalY(1,1 ,_ctheta,_phi);

	result += e_02n2*HarmonicSphericalY(2,-2,_ctheta,_phi);
	result += e_02n1*HarmonicSphericalY(2,-1,_ctheta,_phi);
	result += e_02p0*HarmonicSphericalY(2,0 ,_ctheta,_phi);
	result += e_02p1*HarmonicSphericalY(2,1 ,_ctheta,_phi);
	result += e_02p2*HarmonicSphericalY(2,2 ,_ctheta,_phi);

/*
	if(__debug) {
		cout << "Acceptance = " << result << endl;
	}
*/

        return result;
}

//_____________________________________________________________________________
Double_t RooBsTimeAngle::rho_B (Int_t bi) const {    //8.19
	if (__debug) {
		cout << "P_B(" << bi << ") = " << P_B(bi) << endl; 
		cout << "Q_B(" << bi << ") = " << Q_B(bi) << endl; 
	}
	return	(1-_Fs)*P_B(bi) 
		+ _Fs*Q_B(bi) 
		+ 0.206748335783172033*(I_mu()*CrossDot(&RooBsTimeAngle::Aminus,&RooBsTimeAngle::B)).Re()*f_minus_sq(bi)
		+ 0.206748335783172033*(I_mu()*CrossDot(&RooBsTimeAngle::Aplus,&RooBsTimeAngle::B)*f_f_star(bi)).Re();
}

//_____________________________________________________________________________
Double_t RooBsTimeAngle::rho_Bbar (Int_t bi) const {    //8.20
	return	(1-_Fs)*P_Bbar(bi) 
		+ _Fs*Q_Bbar(bi) 
		+ 0.206748335783172033*(I_mu()*CrossDot(&RooBsTimeAngle::Aminus,&RooBsTimeAngle::B)).Re()*f_bar_minus_sq(bi)
		+ 0.206748335783172033*(I_mu()*CrossDot(&RooBsTimeAngle::Aplus,&RooBsTimeAngle::B)*f_f_star_bar(bi)).Re();
}

//_____________________________________________________________________________
Double_t RooBsTimeAngle::P_B (Int_t bi) const {    //5.5
/*
	if(__debug) {
		cout << "Aplus X Aplus = " << CrossDot(&RooBsTimeAngle::Aplus,&RooBsTimeAngle::Aplus) << endl;
		cout << "f_plus_sq(" << bi << ") = " << f_plus_sq(bi) << endl;
		cout << "Aminus X Aminus = " << CrossDot(&RooBsTimeAngle::Aminus,&RooBsTimeAngle::Aminus) << endl; 
		cout << "f_plus_sq(" << bi << ") = " << f_minus_sq(bi) << endl;
		cout << "Aplus X Aminus = " << CrossDot(&RooBsTimeAngle::Aplus,&RooBsTimeAngle::Aminus) << endl;
		cout << "f_f_star_sq(" << bi << ") = " << f_f_star(bi) << endl;
	}
*/

	return 0.179049310978382253*CrossDot(&RooBsTimeAngle::Aplus,&RooBsTimeAngle::Aplus).Re()*f_plus_sq(bi)
		+ 0.179049310978382253*CrossDot(&RooBsTimeAngle::Aminus,&RooBsTimeAngle::Aminus).Re()*f_minus_sq(bi)
		+ 0.358098621956764507*(CrossDot(&RooBsTimeAngle::Aplus,&RooBsTimeAngle::Aminus)*f_f_star(bi) ).Re();

}

//_____________________________________________________________________________
Double_t RooBsTimeAngle::P_Bbar (Int_t bi) const {    //5.6
	return 0.179049310978382253*CrossDot(&RooBsTimeAngle::Aplus,&RooBsTimeAngle::Aplus).Re()*f_bar_plus_sq(bi)
		+ 0.179049310978382253*CrossDot(&RooBsTimeAngle::Aminus,&RooBsTimeAngle::Aminus).Re()*f_bar_minus_sq(bi)
		+ 0.358098621956764507*(CrossDot(&RooBsTimeAngle::Aplus,&RooBsTimeAngle::Aminus)*f_f_star_bar(bi) ).Re();
}

//_____________________________________________________________________________
Double_t RooBsTimeAngle::Q_B (Int_t bi) const {    //8.15
	return 0.0596831036594607511*CrossDot(&RooBsTimeAngle::B,&RooBsTimeAngle::B).Re()*f_minus_sq(bi);
}

//_____________________________________________________________________________
Double_t RooBsTimeAngle::Q_Bbar (Int_t bi) const {    //8.16
	return 0.0596831036594607511*CrossDot(&RooBsTimeAngle::B,&RooBsTimeAngle::B).Re()*f_bar_minus_sq(bi);
}


//_____________________________________________________________________________
TComplex RooBsTimeAngle::Aplus(int i) const {    // 5.2
	Double_t DG = 1/_tau_L - 1/_tau_H;
	Double_t G = (1/_tau_L + 1/_tau_H)/2 ;
	Double_t z = TMath::Cos(2.0*_beta)*DG/(2.0*G);
	Double_t y = (1.0-z)/(1.0+z);

	Double_t xx = y/(y+(1.0-y)*_Ap2);
	Double_t a0 = TMath::Sqrt(_A02*xx);
	Double_t al = TMath::Sqrt(_All2*xx);

/*
	if (__debug) {
		cout << "DG = " << DG << endl;
		cout << "G = " << G << endl;
		cout << "z = " << z << endl;
		cout << "y = " << y << endl;
		cout << "xx = " << xx << endl;
		cout << "a0 = " << a0 << endl;
		cout << "al = " << al << endl;

		cout << "Aplus = " 
			<< TComplex(a0*_cpsi,0.0) << " , " 
			<< -0.707106781186547462*TMath::Sqrt(1-_cpsi*_cpsi)*al*(TMath::Cos(_delta_l) + TComplex::I()*TMath::Sin(_delta_l)) << " , "
			<< 0.0 << endl;
	}
*/

	if (i==1) return TComplex(a0*_cpsi,0.0);
	if (i==2) return -0.707106781186547462*TMath::Sqrt(1-_cpsi*_cpsi)*al*(TMath::Cos(_delta_l) + TComplex::I()*TMath::Sin(_delta_l));
	if (i==3) return 0.0;
	std::cout << "ERROR: Aplus("<< i <<")" << std::endl;
	exit(1);
}

//_____________________________________________________________________________
TComplex RooBsTimeAngle::Aminus(int i) const {    // 5.3
        Double_t DG = 1/_tau_L - 1/_tau_H;
        Double_t G = (1/_tau_L + 1/_tau_H)/2 ;
        Double_t z = TMath::Cos(2.0*_beta)*DG/(2.0*G);
        Double_t y = (1.0-z)/(1.0+z);

        Double_t xx = y/(y+(1.0-y)*_Ap2);
        Double_t ap = TMath::Sqrt(_Ap2*xx);

/*
	if(__debug){
		cout << "Aminus = 0 , 0 , " << 0.707106781186547462*TMath::Sqrt(1-_cpsi*_cpsi)*ap*(TComplex::I()*TMath::Cos(_delta_p) -TMath::Sin(_delta_p)) << endl;
	}
*/

	if (i==1) return 0.0;
	if (i==2) return 0.0; 
	if (i==3) return 0.707106781186547462*TMath::Sqrt(1-_cpsi*_cpsi)*ap*(TComplex::I()*TMath::Cos(_delta_p) -TMath::Sin(_delta_p));
	std::cout << "ERROR: Aminus("<< i <<")" << std::endl;
	exit(1);
}

//_____________________________________________________________________________
TComplex RooBsTimeAngle::B(int i) const {    // 8.16
	if (i==1) return TComplex(1.0,0.0);
	if (i==2) return 0.0; 
	if (i==3) return 0.0;
	std::cout << "ERROR: B("<< i <<")" << std::endl;
	exit(1);
}

//_____________________________________________________________________________
Double_t RooBsTimeAngle::nhat(int i) const {    // 3.6
/*	if(__debug) {
		cout << "n = " 
			<< TMath::Sqrt(1-_ctheta*_ctheta)*TMath::Cos(_phi) << ' ' 
			<< TMath::Sqrt(1-_ctheta*_ctheta)*TMath::Sin(_phi) << ' ' 
			<< _ctheta << endl;
	}*/
	if (i==1) return TMath::Sqrt(1-_ctheta*_ctheta)*TMath::Cos(_phi);
	if (i==2) return TMath::Sqrt(1-_ctheta*_ctheta)*TMath::Sin(_phi);
	if (i==3) return _ctheta; 
	std::cout << "ERROR: nhat("<< i <<")" << std::endl;
	exit(1);
}

//_____________________________________________________________________________
TComplex RooBsTimeAngle::CrossDot(TComplex (RooBsTimeAngle::*A)(int)const ,TComplex (RooBsTimeAngle::*B)(int)const ) const{
	return	-(this->*A)(2)*nhat(1)*TComplex::Conjugate(-(this->*B)(2)*nhat(1) + (this->*B)(1)*nhat(2)) + 
		 (this->*A)(1)*nhat(2)*TComplex::Conjugate(-(this->*B)(2)*nhat(1) + (this->*B)(1)*nhat(2)) + 
		 (this->*A)(3)*nhat(1)*TComplex::Conjugate( (this->*B)(3)*nhat(1) - (this->*B)(1)*nhat(3)) - 
		 (this->*A)(1)*nhat(3)*TComplex::Conjugate( (this->*B)(3)*nhat(1) - (this->*B)(1)*nhat(3)) - 
		 (this->*A)(3)*nhat(2)*TComplex::Conjugate(-(this->*B)(3)*nhat(2) + (this->*B)(2)*nhat(3)) + 
		 (this->*A)(2)*nhat(3)*TComplex::Conjugate(-(this->*B)(3)*nhat(2) + (this->*B)(2)*nhat(3)) ;
}

//_____________________________________________________________________________
Double_t RooBsTimeAngle::f_sq(Double_t A, Double_t B, int bi) const {    //5.7  5.8
	Double_t CB = TMath::Cos(2.0*_beta);
	Double_t L = 1+A*CB;
	Double_t H = 1-A*CB;
	Double_t C = 2.0*(_tau_L*L+_tau_H*H);
/*
	if(__debug) {
		cout << "CB = " << CB << endl;
		cout << "L = " << L << endl;
		cout << "H = " << H << endl;
		cout << "C = " << C << endl;

		cout << "f_sq = (" << A << "," << B << ") "  << L/C << " , " << H/C << " , " << A*B*2.0*TMath::Sin(2.0*_beta)/C << " , 0  = " << L/C + H/C + A*B*2.0*TMath::Sin(2.0*_beta)/C << endl;

	}
*/

	if (bi == _basisExpGammaL) return L/C; 
	if (bi == _basisExpGammaH) return H/C; 
	if (bi == _basisExpSinDm ) return A*B*2.0*TMath::Sin(2.0*_beta)/C; 
	if (bi == _basisExpCosDm ) return 0.0; 
	std::cout << "ERROR: f_plus_sq("<< bi<< ")" << std::endl;
	exit(1);
}

//_____________________________________________________________________________
Double_t RooBsTimeAngle::f_bar_plus_sq(int bi) const {    //5.7
	return f_sq(1.0,1.0,bi);
}
Double_t RooBsTimeAngle::f_bar_minus_sq(int bi) const {    //5.7
	return f_sq(-1.0,1.0,bi);
}
Double_t RooBsTimeAngle::f_plus_sq(int bi) const {    //5.8
	return f_sq(1.0,-1.0,bi);
}
Double_t RooBsTimeAngle::f_minus_sq(int bi) const {    //5.8
	return f_sq(-1.0,-1.0,bi);
}

//_____________________________________________________________________________
TComplex RooBsTimeAngle::f_f_star(Double_t A, int bi) const {    //5.9 5.10
	Double_t C = TMath::Sqrt(((_tau_L-_tau_H)*TMath::Sin(2.0*_beta))*((_tau_L-_tau_H)*TMath::Sin(2.0*_beta))+4.0*_tau_L*_tau_H);
/*
	if(__debug){
		cout << "C = " << C << endl;
	}
*/

	if (bi == _basisExpGammaL) return TComplex(0.0, TMath::Sin(2.0*_beta)/2.0/C);
	if (bi == _basisExpGammaH) return TComplex(0.0, -TMath::Sin(2.0*_beta)/2.0/C); 
	if (bi == _basisExpSinDm ) return TComplex(0.0, -A*TMath::Cos(2.0*_beta)/C); 
	if (bi == _basisExpCosDm ) return TComplex(-A/C, 0.0); 

	std::cout << "ERROR: f_f_star("<< bi <<")" << std::endl;
	exit(1);
}

//_____________________________________________________________________________
TComplex RooBsTimeAngle::f_f_star_bar(int bi) const {    //5.9
	return f_f_star(1.0,bi);
}
TComplex RooBsTimeAngle::f_f_star(int bi) const {    //5.10
	return f_f_star(-1.0,bi);
}

//_____________________________________________________________________________
TComplex RooBsTimeAngle::I_mu() const {    //8.10
	TComplex imu( TMath::Cos(_delta_s)*0.0032835 - TMath::Sin(_delta_s)*0.326332, 
	             -TMath::Cos(_delta_s)*0.326332 - TMath::Sin(_delta_s)*0.0032835 );
	
	if(__debug) {
		cout << "I_mu = " << imu << endl;
	}

	return TMath::Sqrt(_Fs*(1-_Fs))*imu; 
}


//_____________________________________________________________________________
Int_t RooBsTimeAngle::getCoefAnalyticalIntegral(Int_t /*code*/, RooArgSet& allVars, RooArgSet& analVars, const char* /* rangeName*/) const
{
//    cout << "Integrate .." << endl;
//    allVars.Print();
  if (matchArgs(allVars, analVars, _cpsi, _ctheta, _phi)) return 1;
  if (matchArgs(allVars, analVars,        _ctheta, _phi)) return 2;
  if (matchArgs(allVars, analVars, _cpsi,          _phi)) return 3;
  if (matchArgs(allVars, analVars, _cpsi, _ctheta      )) return 4;
  if (matchArgs(allVars, analVars, _cpsi               )) return 5;
  if (matchArgs(allVars, analVars,        _ctheta      )) return 6;
  if (matchArgs(allVars, analVars,                 _phi)) return 7;

 return 0;
}


//_____________________________________________________________________________
Double_t RooBsTimeAngle::coefAnalyticalIntegral(Int_t bi, Int_t code, const char* range) const
{
	if ( code == 1 ){
	        Double_t DG = 1/_tau_L - 1/_tau_H;
	        Double_t G = (1/_tau_L + 1/_tau_H)/2 ;
	        Double_t z = TMath::Cos(2.0*_beta)*DG/(2.0*G);
	        Double_t y = (1.0-z)/(1.0+z);
	
	        Double_t xx = y/(y+(1.0-y)*_Ap2);
	        Double_t a0 = TMath::Sqrt(_A02*xx);
	        Double_t al = TMath::Sqrt(_All2*xx);
	        Double_t ap = TMath::Sqrt(_Ap2*xx);

	
		Double_t IAplusXn2 = ((a0*a0+al*al)/3.54490770181103176) * e_00p0
					+ ((a0*a0+al*al)/15.8533091904240440) * e_02p0
					- ((a0*a0-al*al)*0.109254843059207907) * e_02p2 ;
		Double_t IAminusXn2 = (ap*ap/3.54490770181103176) * e_00p0
					+ (-ap*ap/7.92665459521202198) * e_02p0 ;
		TComplex IAplusXAminus = al*TComplex(TMath::Cos(_delta_l),TMath::Sin(_delta_l)) * ap*TComplex(TMath::Sin(_delta_p),TMath::Cos(_delta_p)) * 0.218509686118415813 * e_02n1 ;

		Double_t IBXn2 = e_00p0/3.54490770181103176
					+ e_02p0/15.8533091904240440
					- 0.109254843059207907 * e_02p2 ;
		TComplex IAminusXB = ap*TComplex(-TMath::Sin(_delta_p),TMath::Cos(_delta_p))* 1.01663295048409386 * e_02p1 ;
		TComplex IAplusXB  = al*TComplex(TMath::Cos(_delta_l),TMath::Sin(_delta_l)) * 1.01663295048409386 * e_02n2 ;


		Double_t IP_B = IAplusXn2*f_plus_sq(bi)
		                + IAminusXn2*f_minus_sq(bi)
		                + (IAplusXAminus*f_f_star(bi)).Re();

		Double_t IP_Bbar = IAplusXn2*f_bar_plus_sq(bi)
		                + IAminusXn2*f_bar_minus_sq(bi)
		                + (IAplusXAminus*f_f_star_bar(bi)).Re();

		Double_t IQ_B = IBXn2*f_minus_sq(bi);
		Double_t IQ_Bbar = IBXn2*f_bar_minus_sq(bi);


       		Double_t I_rho_B = (1-_Fs)*IP_B
			                + _Fs*IQ_B
			                + 0.206748335783172033*(I_mu()*IAminusXB).Re()*f_minus_sq(bi)
			                + 0.206748335783172033*(I_mu()*IAplusXB*f_f_star(bi)).Re();
      
		Double_t I_rho_Bbar = (1-_Fs)*IP_Bbar
			                + _Fs*IQ_Bbar
			                + 0.206748335783172033*(I_mu()*IAminusXB).Re()*f_bar_minus_sq(bi)
			                + 0.206748335783172033*(I_mu()*IAplusXB*f_f_star_bar(bi)).Re();
	
		if (__debug){
			cout << "a0 = "  << a0 << endl;
			cout << "al = "  << al*TComplex(TMath::Cos(_delta_l),TMath::Sin(_delta_l))  << " |al| = " << al << " delta_l = " << _delta_l << endl;
			cout << "ap = "  << ap*TComplex(TMath::Cos(_delta_p),TMath::Sin(_delta_p))  << " |ap| = " << ap << " delta_p = " << _delta_p << endl;

			cout << "IAplusXn2 = " << IAplusXn2 << endl;
			cout << "IAminusXn2 = " << IAminusXn2 << endl;
			cout << "IAplusXAminus = " << IAplusXAminus << endl;
			cout << "IBXn2 = " << IBXn2 << endl;
			cout << "IAplusXB = " << IAplusXB << endl;
			cout << "IAminusXB = " << IAminusXB << endl;
			
			cout << "INTEGRAL N(" << bi << ") = " << I_rho_B << endl;
		}

		return (1-_D/2.0)*I_rho_B + (1+_D/2.0)*I_rho_Bbar ;

		//Double_t I_rho_B = (1.0-_Fs)*( (a0*a0+al*al)*f_plus_sq(basisIndex) + ap*ap*f_minus_sq(basisIndex) )
		//			+ _Fs*f_minus_sq(basisIndex);
	}
	return 0;
}


//_____________________________________________________________________________
Int_t RooBsTimeAngle::getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK) const
{
  if (matchArgs(directVars, generateVars, _t, _cpsi, _ctheta, _phi)) return 1;
//  if (staticInitOK)
//      if( matchArgs(directVars, generateVars, _t) ) return 2;
  return 0;
}


//_____________________________________________________________________________
void RooBsTimeAngle::initGenerator(Int_t code)
{
}

//_____________________________________________________________________________
void RooBsTimeAngle::generateEvent(Int_t code) {
    Double_t max = 0.1;
    Double_t value = 0;

    while (1) {
//        _t = -_tau_H * log(RooRandom::uniform(&_aleatorio));
        _t = RooRandom::uniform(&_aleatorio)*15;
        switch (code) {
            case 1:
                _cpsi   = RooRandom::uniform(&_aleatorio) * (_cpsi.max() - _cpsi.min()) + _cpsi.min();
                _ctheta = RooRandom::uniform(&_aleatorio) * (_ctheta.max() - _ctheta.min()) + _ctheta.min();
                _phi    = RooRandom::uniform(&_aleatorio) * (_phi.max() - _phi.min()) + _phi.min();
		if (__debug) {
			_t = 1.28222;
			_cpsi = 0.918636;
			_ctheta = -0.4005;
			_phi = 0.485588;
		}

		if (__debug) cout << "t cpsi ctheta phi = " << _t << ' ' << _cpsi << ' ' << _ctheta << ' ' << _phi << endl;

                value = ( coefficient(_basisExpGammaL) * TMath::Exp(-_t/_tau_L)
                        + coefficient(_basisExpGammaH) * TMath::Exp(-_t/_tau_H)
                        + coefficient(_basisExpCosDm)  * TMath::Exp( -_t*(_tau_L+_tau_H)/(2.0*_tau_L*_tau_H) )*TMath::Cos(_Dm * _t)
                        + coefficient(_basisExpSinDm)  * TMath::Exp( -_t*(_tau_L+_tau_H)/(2.0*_tau_L*_tau_H) )*TMath::Sin(_Dm * _t));// * exp(-0.5 * (_tau_L-_tau_H) * _t);

		if (__debug) cout << "value = " << value << endl;
                
		break;
	    default:
		std::cout << "ERROR: generateEvent("<< code << ")" << std::endl;
		exit(1);

        }
	
	if(__debug) exit(0);

        if (value > max) {
        	cout << "ERROR: Value > max " << value << endl;
            max = value * 1.05;
        }


        Double_t rand = RooRandom::uniform(&_aleatorio) * max;
        //cout << rand << ' ' << value << endl;
        if (rand < value)
            break;
    }
}
