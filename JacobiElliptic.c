#include <math.h>
#include <float.h>
#include <stdio.h>

// Returns the values of the Jacobi elliptic functions Cn, Sn and Dn,
// evaluated for argument x and parameter 0 <= k <= 1.
// Uses the method of the arithmetic-geometric mean described in [1].
// 
//    [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical
//        Functions" Dover Publications", 1965, Ch. 16-17.6.
void JacobiElliptic(double x, double k, double *cn, double *sn, double *dn){
	const double tol = DBL_EPSILON;
	double a[16], b[16], c[16], phi[16];
	int i;
	if(1 == k){
		*cn = 1./cosh(x);
		*sn = tanh(x);
		*dn = *cn;
		return;
	}else if(0 == k){
		*cn = cos(x);
		*sn = sin(x);
		*dn = 1;
		return;
	}
	a[0] = 1;
	c[0] = sqrt(k);
	b[0] = sqrt(1-k);
	i = 0;
	while(fabs(c[i]) > tol){
		a[i+1] = 0.5 * a[i] + 0.5 * b[i];
		b[i+1] = sqrt(a[i] * b[i]);
		c[i+1] = 0.5 * a[i] - 0.5 * b[i];
		++i;
	}
	phi[i] = pow(2.,i) * a[i] * x;
	while(i > 0){
		phi[i-1] = 0.5 * (phi[i] + asin(c[i]/a[i] * sin(phi[i])));
		--i;
	}
	*cn = cos(phi[0]);
	*sn = sin(phi[0]);
	*dn = sqrt(1. - k * (*sn)*(*sn));
}

void JacobiElliptic_complex(
	double zr, double zi, double k,
	double *cnr, double *cni, double *snr, double *sni, double *dnr, double *dni
){
	double kp = 1.-k;
	double cr, sr, dr;
	double ci, si, di;
	double idenom;
	JacobiElliptic(zr, k , &cr, &sr, &dr);
	JacobiElliptic(zi, kp, &ci, &si, &di);
	idenom = 1./ (ci*ci + k*sr*sr*si*si);
	*snr = idenom * (sr*di);
	*sni = idenom * (cr*dr*si*ci);
	*cnr = idenom * (cr*ci);
	*cni = idenom * (-sr*dr*si*di);
	*dnr = idenom * (dr*ci*di);
	*dni = idenom * (-k*sr*cr*si);
}

double EllipticK(double x){
	if(0 == x){ return 0.5*M_PI; }
	double a = 1;
	double b = cos(x);
	double c = sin(x);
	while(fabs(c) > DBL_EPSILON){
		double an = 0.5 * a + 0.5 * b;
		double bn = sqrt(a * b);
		c = 0.5 * a - 0.5 * b;
	}
	return 0.5*M_PI / a;
}
