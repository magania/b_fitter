#include <iostream>
#include <TComplex.h>

int main(void){
	using namespace std;

	TComplex z = TComplex(2.0,3.0);

	TComplex log_z = TComplex::Log(z);
	TComplex log_erf;

	double n = 0;
	double nfact = 1;
	double k;
	double m1n = 1;
	while (1) {
		n++;
		nfact = n*nfact;
		k = nfact*(2.0*n+1);
		m1n = -1*m1n;

		log_erf += m1n*n*log_z/k;

		cout << "m1n: " << m1n << " " << log_erf << " k: " << k << endl;

                if (n==10) 
                        break;
	}
	TComplex erf = 1.1283791670955126*TComplex::Exp(log_erf);

	cout << "erf = " << erf << endl;
}
