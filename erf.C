#include "TMath.h"

erf(){
	double x0 = 4.048113054268301;
	double y0 = 12.565287501684949;
//	double x0 = 2.0;
	//double y0 = 2.0;

//	Erf(x0,y0);

	double v0 = x0*y0;
	double x02 = x0*x0;
	long double gamma_v02n = TMath::Erfc(x0);
	double Emx02 = TMath::Exp(-x0*x0);

	long double H1 = gamma_v02n;
	long double nfact = 1;
	double k = 1.12837916709551256;
	long double a = Emx02*(1.0/nfact)*(1.0/x0);
	long double b = 1.77245385090551588*gamma_v02n;
	gamma_v02n = k*(a-b);
	
	long double h2 = gamma_v02n;

	long double v0n = v0;
	long double x02np1 = x0*x02;

	double _H1, _h2;
	double n =0;
	while (1){
		n++;
		cout << "n = " << n << endl;
		cout << "nfact = " << nfact << endl;
		cout << "gamma_v02n = " << gamma_v02n << endl;
		cout << "H1 = " << H1 << endl;
		cout << "h2 = " << h2 << endl;
		cout << "x02np1 = " << x02np1 << endl;
		cout << "v0n = " << v0n << endl;
		cout << endl;

		gamma_v02n *= v0*v0;
		_H1=H1;
		H1 += gamma_v02n;
		
		nfact = (n+1)*nfact;
		k = 2.0/(1.77245385090551588*(2.0*n + 1));
		a = Emx02*(v0n/nfact)*(v0n/x02np1);
		b = 1.77245385090551588*gamma_v02n/(n+1);
		gamma_v02n = k*(a-b);

		_h2 = h2;
		h2 += (n+1)*gamma_v02n;

		if (_H1 == H1 && _h2 == h2)
			break;

		if ( nfact > 1e4900){
			cout << "Erf: ERROR z = (" << x0 << ", " << y0 << ")" << endl;
			exit(1);
		}
	
		x02np1 *= x02;
		v0n *= v0;
	}

	double H2 = x0*h2;

	double erfc_re =  H1*TMath::Cos(2*v0) - y0*H2*TMath::Sin(2*v0);
	double erfc_im = -H1*TMath::Sin(2*v0) - y0*H2*TMath::Cos(2*v0);

	cout << "erfc = " << erfc_re << " + i " << erfc_im << endl;
}

TComplex Erf(double x0, double y0) const{
		cout << "Erf: " << x0 << " " << y0 << "I" << endl;

	double gamma_n = TMath::Erfc(x0);
	double nfact = 1;
	double x02np1 = x0;

	double v0 = x0*y0;
	double v02n = 1;

	double H1=0;
	double h2=0;

	double _H1, _h2;

	double a, b;
	for (int n=0; n<5; n++){

			cout << "Erf: n = " << n << endl;
			cout << "Erf: nfact = " << nfact << endl;
			cout << "Erf: gamma_n = " << gamma_n << endl;
			cout << "Erf: H1 = " << H1 << endl;
			cout << "Erf: h2 = " << h2 << endl;
			cout << "Erf: a = " << a << endl;
			cout << "Erf: b = " << b  << endl;

		_H1 = H1;
		H1 += gamma_n*v02n;

		double k = 2.0/((2.0*n + 1)*TMath::Sqrt(TMath::Pi()));
		nfact = (n+1)*nfact;
		a = TMath::Exp(-x0*x0)/(nfact*x02np1);
		b = TMath::Sqrt(TMath::Pi())*gamma_n/(n+1);
		gamma_n = k*(a - b);

		_h2 = h2;
		h2 += (n+1)*gamma_n*v02n;

		if (_H1 == H1 && _h2 == h2)
			break;

		x02np1 *= x0*x0;
		v02n *= v0*v0;
	}

	double H2 = x0*h2;

	double erfc_re =  H1*TMath::Cos(2*v0) - y0*H2*TMath::Sin(2*v0);
	double erfc_im = -H1*TMath::Sin(2*v0) - y0*H2*TMath::Cos(2*v0);

		cout << "erfc = " << erfc_re << " + i " << erfc_im << endl;
		cout << "erf = " << 1.0-erfc_re << " + i " << - erfc_im << endl;

	return TComplex(1.0-erfc_re, -erfc_im);
}

