#include "strdata.h"

void m1_prim2conc(double * S_prim,double * S_cons) {

	double rho;

	// rho
	S_cons[0] = S_prim[0];
	rho       = S_prim[0];
	
	// ru rv rw
	S_cons[1] = rho*S_prim[1];
	S_cons[2] = rho*S_prim[2];
	S_cons[3] = rho*S_prim[3];
		
}
