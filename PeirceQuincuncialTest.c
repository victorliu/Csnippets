#include <stdio.h>
#include <math.h>

void JacobiElliptic_complex(
	double zr, double zi, double k,
	double *cnr, double *cni, double *snr, double *sni, double *dnr, double *dni
);

#define KHALF 1.8540746773013719184

// maps x, y in [-1,1] to phi in [0,pi/2] and theta in [-pi,pi]
void PeirceQuincuncial(double x, double y, double *phi, double *theta){
	double cnr, cni, dum;
	double z[2] = { 0.5*(x-y) + 1., 0.5*(x+y) };
	z[0] *= KHALF;
	z[1] *= KHALF;
	JacobiElliptic_complex(z[0], z[1], 0.5, &cnr, &cni, &dum, &dum, &dum, &dum);
	*theta = atan2(cni,cnr);
	*phi = 2 * atan(hypot(cnr, cni));
}

int main(){
	const int nx = 10;
	const int ny = 10;
	int ix, iy;
	
	printf("PVF(1)\n");
	
	for(ix = 0; ix < nx; ++ix){
		const double tx = ((double)ix+0.5) / (double)nx;
		const double x = 2*tx - 1.;
		for(iy = 0; iy < ny; ++iy){
			const double ty = ((double)iy+0.5) / (double)ny;
			const double y = 2*ty - 1.;
			
			double phi; // in 0 to pi/2
			double theta; // in -pi to pi
			PeirceQuincuncial(x, y, &phi, &theta);
			
			printf("p{%g, %g, %g}\n", sin(phi)*cos(theta), sin(phi)*sin(theta), cos(phi));
		}
	}
}
