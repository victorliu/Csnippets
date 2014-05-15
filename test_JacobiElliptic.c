#include <stdio.h>
extern void JacobiElliptic(double x, double k, double *cn, double *sn, double *dn);
extern void JacobiElliptic_complex(
	double zr, double zi, double k,
	double *cnr, double *cni, double *snr, double *sni, double *dnr, double *dni
);
int main(){
	double k = 0.5;
	
	/*
	double x;
	for(x = 0; x < 10; x += 0.1){
		double cn, sn, dn;
		JacobiElliptic(x, k, &cn, &sn, &dn);
		printf("%g\t%g\t%g\t%g\t%g\n", x, k, cn, sn, dn);
	}
	*/
	
	double zr = 0.3;
	double zi = 0.4;
	double cnr, cni, snr, sni, dnr, dni;
	JacobiElliptic_complex(zr, zi, k, &cnr, &cni, &snr, &sni, &dnr, &dni);
	printf("%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", zr, zi, k, cnr, cni, snr, sni, dnr, dni);
}