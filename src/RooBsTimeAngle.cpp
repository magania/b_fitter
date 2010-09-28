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
        RooAbsReal& A02,
        RooAbsReal& All2,
        RooAbsReal& Ap2,
        RooAbsReal& tau_L,
        RooAbsReal& tau_H,
        RooRealVar& Dm,
        RooRealVar& beta,
        RooRealVar& delta_p,
        RooRealVar& delta_l,
	RooRealVar& delta_s,
	RooRealVar& Fs,
	RooRealVar& D,
	RooAbsReal& tau,
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
//  _tau("tau","tau","2.0*@0*@1/(@0+@1)",RooArgList(tau_L,tau_H))
  _tau("_tau", "tau", this, tau)
{
  // Constructor
  cout << "** Bs Time Anlgle **" << endl;
  _basisExpGammaL = declareBasis("exp(-@0/@1)", RooArgList(tau_L)) ;
  _basisExpGammaH = declareBasis("exp(-@0/@1)", RooArgList(tau_H)) ;
  _basisExpCosDm  = declareBasis("exp(-@0/@1)*cos(@0*@2)", RooArgList(tau,Dm)) ;
  _basisExpSinDm  = declareBasis("exp(-@0/@1)*sin(@0*@2)", RooArgList(tau,Dm)) ;
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
		cout << endl << "coefficient " << endl;
		cout << "(cpsi, ctheta, phi) = (" << _cpsi << ", " << _ctheta << ", " << _phi << ")" << endl;
		cout << "rho_B = " << rho_B(_basisExpGammaL) << ", " << rho_B(_basisExpGammaH) << ", " << rho_B(_basisExpCosDm) << ", " << rho_B(_basisExpSinDm) << endl;
		cout << "rho_Bbar = " << rho_Bbar(_basisExpGammaL) << ", " << rho_Bbar(_basisExpGammaH) << ", " << rho_Bbar(_basisExpCosDm) << ", " << rho_Bbar(_basisExpSinDm) << endl;
	}

/*
	static Double_t __tau;
	if (__tau != _tau){
		__tau = _tau;
		cout << "(t, cpsi, ctheta, phi) = (" << _t << ", " << _cpsi << ", " << _ctheta << ", " << _phi << ")" << endl;
		cout << "rho_B = " << rho_B(_basisExpGammaL) << ", " << rho_B(_basisExpGammaH) << ", " << rho_B(_basisExpCosDm) << ", " << rho_B(_basisExpSinDm) << endl;
		cout << "rho_Bbar = " << rho_Bbar(_basisExpGammaL) << ", " << rho_Bbar(_basisExpGammaH) << ", " << rho_Bbar(_basisExpCosDm) << ", " << rho_Bbar(_basisExpSinDm) << endl;
	}
*/
	return ( (0.5-0.5*_D)*rho_B(basisIndex) + (0.5+0.5*_D)*rho_Bbar(basisIndex) ) * acceptance();
}

Double_t RooBsTimeAngle::HarmonicSphericalY(int l, int m) const {
        if (m==0)
                return TMath::Sqrt( ((2.0*l+1)*TMath::Gamma(l+1))/(4.0*TMath::Pi()*TMath::Gamma(l+1)) )*ROOT::Math::assoc_legendre(l,0,_ctheta);
	if (m>0)
                return TMath::Sqrt( ((2.0*l+1)*TMath::Gamma(l-m+1))/(2.0*TMath::Pi()*TMath::Gamma(l+m+1)) )*ROOT::Math::assoc_legendre(l,m,_ctheta)*TMath::Cos(m*_phi)*TMath::Power(-1.0,m);
	if (m<0)
		return TMath::Sqrt( ((2.0*l+1)*TMath::Gamma(l+m+1))/(2.0*TMath::Pi()*TMath::Gamma(l-m+1)) )*ROOT::Math::assoc_legendre(l,-m,_ctheta)*TMath::Sin(-m*_phi)*TMath::Power(-1.0,-m);
}

Double_t RooBsTimeAngle::acceptance() const {
	static Double_t __return, __cpsi, __ctheta;
	if ( !__optimize || __cpsi != _cpsi || __ctheta != _ctheta ){
		__cpsi = _cpsi; __ctheta = _ctheta;

		__return = 0;

		__return += e_00p0*HarmonicSphericalY(0, 0);
	
		__return += e_01n1*HarmonicSphericalY(1,-1);
		__return += e_01p0*HarmonicSphericalY(1, 0);
		__return += e_01p1*HarmonicSphericalY(1, 1);

		__return += e_02n2*HarmonicSphericalY(2,-2);
		__return += e_02n1*HarmonicSphericalY(2,-1);
		__return += e_02p0*HarmonicSphericalY(2, 0);
		__return += e_02p1*HarmonicSphericalY(2, 1);
		__return += e_02p2*HarmonicSphericalY(2, 2);

		if(__debug ) {
			cout << "00 0 = " << e_00p0*HarmonicSphericalY(0, 0) << endl;
			cout << "01-1 = " << e_01n1*HarmonicSphericalY(1,-1) << endl; 
			cout << "01 0 = " << e_01p0*HarmonicSphericalY(1, 0) << endl;
			cout << "01 1 = " << e_01p1*HarmonicSphericalY(1, 1) << endl;
			cout << "02-2 = " << e_02n2*HarmonicSphericalY(2,-2) << endl;
			cout << "02-1 = " << e_02n1*HarmonicSphericalY(2,-1) << endl;
			cout << "02 0 = " << e_02p0*HarmonicSphericalY(2, 0) << endl;
			cout << "02 1 = " << e_02p1*HarmonicSphericalY(2, 1) << endl;
			cout << "02 2 = " << e_02p2*HarmonicSphericalY(2, 2) << endl;

			cout << "acceptance = " << __return << endl;
		}
	}
        return __return;
}

//_____________________________________________________________________________
Double_t RooBsTimeAngle::rho_B (Int_t bi) const {    //8.19
	if (__debug) {
		cout << "Fs = " << _Fs << endl;
		cout << "P_B = " << P_B(_basisExpGammaL) << ", " << P_B(_basisExpGammaH) << ", " << P_B(_basisExpCosDm) << ", " << P_B(_basisExpSinDm)  << endl; 
		cout << "Q_B = " << Q_B(_basisExpGammaL) << ", " << Q_B(_basisExpGammaH) << ", " << Q_B(_basisExpCosDm) << ", " << Q_B(_basisExpSinDm)  << endl; 
	}

	return	(1-_Fs)*P_B(bi) 
		+ _Fs*Q_B(bi) 
		+ 0.206748335783172033*(I_mu()*CrossDot(&RooBsTimeAngle::Aminus,&RooBsTimeAngle::B)).Re()*f_minus_sq(bi)
		+ 0.206748335783172033*(I_mu()*CrossDot(&RooBsTimeAngle::Aplus,&RooBsTimeAngle::B)*f_f_star(bi)).Re();
}

//_____________________________________________________________________________
Double_t RooBsTimeAngle::rho_Bbar (Int_t bi) const {    //8.20
	if (__debug) {
		cout << "P_Bbar(" << bi << ") = " << P_B(bi) << endl; 
		cout << "Q_Bbar(" << bi << ") = " << Q_B(bi) << endl; 
	}

	return	(1-_Fs)*P_Bbar(bi) 
		+ _Fs*Q_Bbar(bi) 
		+ 0.206748335783172033*(I_mu()*CrossDot(&RooBsTimeAngle::Aminus,&RooBsTimeAngle::B)).Re()*f_bar_minus_sq(bi)
		+ 0.206748335783172033*(I_mu()*CrossDot(&RooBsTimeAngle::Aplus,&RooBsTimeAngle::B)*f_f_star_bar(bi)).Re();
}

//_____________________________________________________________________________
Double_t RooBsTimeAngle::P_B (Int_t bi) const {    //5.5
	if(__debug) {
		cout << "Aplus X Aplus = " << CrossDot(&RooBsTimeAngle::Aplus,&RooBsTimeAngle::Aplus) << endl;
		cout << "Aminus X Aminus = " << CrossDot(&RooBsTimeAngle::Aminus,&RooBsTimeAngle::Aminus) << endl; 
		cout << "Aplus X Aminus = " << CrossDot(&RooBsTimeAngle::Aplus,&RooBsTimeAngle::Aminus) << endl;
	}


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
	static TComplex __return1, __return2;
	static Double_t __tau_L, __tau_H, __beta, __A02, __Ap2, __All2,  __delta_l;

	if( !__optimize || __tau_L != _tau_L || __tau_H != _tau_H || __beta != _beta || __A02 != _A02 ||__Ap2 != _Ap2 || __All2 != _All2 || __delta_l != _delta_l){
		__tau_L = _tau_L;  __tau_H = _tau_H; __beta = _beta; __A02 = _A02; __Ap2 = _Ap2; __All2 = _All2; __delta_l = _delta_l;
	
		Double_t DG = 1/_tau_L - 1/_tau_H;
		Double_t G = (1/_tau_L + 1/_tau_H)/2 ;
		Double_t z = TMath::Cos(2.0*_beta)*DG/(2.0*G);
		Double_t y = (1.0-z)/(1.0+z);

		Double_t xx = y/(y+(1.0-y)*_Ap2);
		Double_t a0 = TMath::Sqrt(_A02*xx);
		Double_t al = TMath::Sqrt(_All2*xx);

		__return1 = TComplex(a0, 0.0);
		__return2 = -0.707106781186547462*al*TComplex(TMath::Cos(_delta_l), TMath::Sin(_delta_l));

		if (__debug) {
			cout << "DG = " << DG << endl;
			cout << "G = " << G << endl;
			cout << "z = " << z << endl;
			cout << "y = " << y << endl;
			cout << "xx = " << xx << endl;
			cout << "a0 = " << a0 << endl;
			cout << "al = " << al << endl;

			cout << "Aplus = " << __return1*_cpsi << " , " << __return2*TMath::Sqrt(1.0-_cpsi*_cpsi) << " , 0.0"  << endl;
		}
	}

	if (i==1) return __return1*_cpsi;
	if (i==2) return __return2*TMath::Sqrt(1.0-_cpsi*_cpsi);
	if (i==3) return TComplex(0.0, 0.0);

	std::cout << "ERROR: Aplus("<< i <<")" << std::endl;
	exit(1);
}

//_____________________________________________________________________________
TComplex RooBsTimeAngle::Aminus(int i) const {    // 5.3
	static TComplex __return;
	static Double_t __tau_L, __tau_H, __beta, __Ap2, __delta_p;

	if( !__optimize || __tau_L != _tau_L || __tau_H != _tau_H || __beta != _beta || __Ap2 != _Ap2 || __delta_p != _delta_p){
		 __tau_L = _tau_L;  __tau_H = _tau_H; __beta = _beta; __Ap2 = _Ap2; __delta_p = _delta_p;
	
	        Double_t DG = 1/_tau_L - 1/_tau_H;
	        Double_t G = (1/_tau_L + 1/_tau_H)/2 ;
	        Double_t z = TMath::Cos(2.0*_beta)*DG/(2.0*G);
	        Double_t y = (1.0-z)/(1.0+z);

	        Double_t xx = y/(y+(1.0-y)*_Ap2);
	        Double_t ap = TMath::Sqrt(_Ap2*xx);

		__return = 0.707106781186547462*ap*TComplex(-TMath::Sin(_delta_p), TMath::Cos(_delta_p));

		
		if(__debug){
			cout << "Aminus = 0, 0, " << __return*TMath::Sqrt(1-_cpsi*_cpsi) << endl;
		}
	}

	if (i==1) return TComplex(0.0, 0.0);
	if (i==2) return TComplex(0.0, 0.0); 
	if (i==3) return __return*TMath::Sqrt(1-_cpsi*_cpsi);

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
	static Double_t __return1, __return2, __return3, __ctheta, __phi;
	if ( !__optimize || __ctheta != _ctheta || __phi != _phi ){
		__ctheta = _ctheta; __phi = _phi;

		__return1 = TMath::Sqrt(1-_ctheta*_ctheta)*TMath::Cos(_phi);
		__return2 = TMath::Sqrt(1-_ctheta*_ctheta)*TMath::Sin(_phi);
		__return3 = _ctheta;
	
		if(__debug) 
			cout << "n = " << __return1 << ' ' << __return2 << ' ' << __return3 << endl;
	}
	
	if (i==1) return __return1;
	if (i==2) return __return2;
	if (i==3) return __return3; 
	
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
	static Double_t __return1p, __return2p, __return3p;
	static Double_t __return1m, __return2m, __return3m;
	static Double_t __tau_L, __tau_H, __beta;
	if ( !__optimize || __tau_L != _tau_L || __tau_H != _tau_H || __beta != _beta ){
		__tau_L = _tau_L; __tau_H = _tau_H, __beta = _beta;

		Double_t CB = TMath::Cos(2.0*_beta);
		Double_t L = 1+CB;
		Double_t H = 1-CB;
		Double_t Cp = 2.0*(_tau_L*L+_tau_H*H);
		Double_t Cm = 2.0*(_tau_L*H+_tau_H*L);


		__return1p = L/Cp;
		__return2p = H/Cp;
		__return3p = 2.0*TMath::Sin(2.0*_beta)/Cp;

		__return1m = H/Cm;
		__return2m = L/Cm;
		__return3m = -2.0*TMath::Sin(2.0*_beta)/Cm;

		if(__debug) {
			cout << "CB = " << CB << endl;
			cout << "L = " << L << endl;
			cout << "H = " << H << endl;
			cout << "Cp = " << Cp << endl;
			cout << "Cm = " << Cm << endl;

			cout << "f_sq + =  "  << __return1p << " , " << __return2p << " , " << __return3p << " , " << 0.0  << endl;
			cout << "f_sq - =  "  << __return1m << " , " << __return2m << " , " << __return3m << " , " << 0.0  << endl;
		}

	}


	if (bi == _basisExpCosDm) return 0.0; 
	if (bi == _basisExpSinDm)
		if (A>0)
			return B*__return3p;
		else
			return B*__return3m;

	if (bi == _basisExpGammaL) 
		if ( A>0 )
			return __return1p;
		else
			return __return1m;

	if (bi == _basisExpGammaH) 
		if ( A>0 )
			return __return2p;
		else
			return __return2m;

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
	static TComplex __return1, __return2, __return3, __return4;
	static Double_t __tau_L, __tau_H, __beta;

	if ( !__optimize || __tau_L != _tau_L || __tau_H != _tau_H || __beta != _beta ){
		__tau_L = _tau_L; __tau_H = _tau_H, __beta = _beta;

		Double_t C = TMath::Sqrt(((_tau_L-_tau_H)*TMath::Sin(2.0*_beta))*((_tau_L-_tau_H)*TMath::Sin(2.0*_beta))+4.0*_tau_L*_tau_H);

		__return1 = TComplex(0.0, TMath::Sin(2.0*_beta)/2.0/C);
		__return2 = TComplex(0.0, -TMath::Sin(2.0*_beta)/2.0/C);
		__return3 = TComplex(0.0, -TMath::Cos(2.0*_beta)/C);
		__return4 = TComplex(-1.0/C, 0.0);
	
		if(__debug){
			cout << "f_f_star = " << __return1 << " " << __return2 << " " << __return3 << " " << __return4 << endl;
		}
	}

	if (bi == _basisExpGammaL) return __return1;
	if (bi == _basisExpGammaH) return __return2;
	if (bi == _basisExpSinDm ) return A*__return3;
	if (bi == _basisExpCosDm ) return A*__return4;

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
	static TComplex __return;
	static Double_t __delta_s;

	if ( !__optimize || __delta_s != _delta_s ){
		__delta_s = _delta_s;

		__return = TMath::Sqrt(_Fs*(1-_Fs))*TComplex(  TMath::Cos(_delta_s)*0.0032835 - TMath::Sin(_delta_s)*0.326332, 
		                                              -TMath::Cos(_delta_s)*0.326332 - TMath::Sin(_delta_s)*0.0032835  );

		if (__debug){
			cout << "I_mu = " << __return << endl;
		}
	}

	return __return;
}


//_____________________________________________________________________________
Int_t RooBsTimeAngle::getCoefAnalyticalIntegral(Int_t /*code*/, RooArgSet& allVars, RooArgSet& analVars, const char* /* rangeName*/) const
{
    cout << "Integrate .." << endl;
    allVars.Print();
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
	
		if (__debug2){
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

		return (0.5-0.5*_D)*I_rho_B + (0.5+0.5*_D)*I_rho_Bbar ;

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
    Double_t max = 0.08;
    Double_t value = 0;
	
    Double_t G;
    while (1) {
        _t = -_tau_H * log(RooRandom::uniform(&_aleatorio));
//        _t = RooRandom::uniform(&_aleatorio)*15;
        switch (code) {
            case 1:
                _cpsi   = RooRandom::uniform(&_aleatorio) * (_cpsi.max() - _cpsi.min()) + _cpsi.min();
                _ctheta = RooRandom::uniform(&_aleatorio) * (_ctheta.max() - _ctheta.min()) + _ctheta.min();
                _phi    = RooRandom::uniform(&_aleatorio) * (_phi.max() - _phi.min()) + _phi.min();
		if (__debug) {
			_t = 0;
			_cpsi = 0.918636;
			_ctheta = -0.4005;
			_phi = 0.485588;
		}


//		if (__debug) cout << "t cpsi ctheta phi = " << _t << ' ' << _cpsi << ' ' << _ctheta << ' ' << _phi << endl;

		G = (1/_tau_L + 1/_tau_H)/2 ;
                value = ( coefficient(_basisExpGammaL) * TMath::Exp(-_t/_tau_L)
                        + coefficient(_basisExpGammaH) * TMath::Exp(-_t/_tau_H)
                        + coefficient(_basisExpCosDm)  * TMath::Exp( -_t*G )*TMath::Cos(_Dm * _t)
                        + coefficient(_basisExpSinDm)  * TMath::Exp( -_t*G )*TMath::Sin(_Dm * _t))* acceptance() * exp(_t/ _tau_H);

		if (__debug) cout << "coefficients = " << coefficient(_basisExpGammaL) << ", " << coefficient(_basisExpGammaH) << ", " << coefficient(_basisExpCosDm) << ", " << coefficient(_basisExpSinDm) << endl;
		if (__debug) cout << "value = " << value << endl;
                
		if(__debug) exit(0);

		break;
	    default:
		std::cout << "ERROR: generateEvent("<< code << ")" << std::endl;
		exit(1);
        }
	
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
