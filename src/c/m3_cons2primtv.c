#include "strdata.h"

void m3_cons2primtv(double * S_cons,double * S_prim) {

	int i,i_prim,i_cons;
	double avfsum,rho;
	double avf[10];

	// avf
	avfsum = 0.0;
	i_cons = eqtypi[3];
	for(i=0;i<eqtypn[3];++i){
		S_prim[i] = S_cons[i_cons +i];
		avf[i] 	  = S_cons[i_cons +i];
		avfsum   += avf[i]; 
	}
	avf[eqtypn[3]] = 1.0 - avfsum;

	// rho_1 ...
	rho			=	0.0;
	i_prim 	= eqtypn[3];
	for(i=0;i<eqtypn[0];++i){
		S_prim[i_prim + i] = S_cons[i];
		rho								+= S_cons[i];
	}

	i_prim = eqtypn[3] + eqtypn[0];
	i_cons = eqtypi[1];

	S_prim[i_prim + 0] = S_cons[i_cons + 0]/rho; 
	S_prim[i_prim + 1] = S_cons[i_cons + 1]/rho; 
	S_prim[i_prim + 2] = S_cons[i_cons + 2]/rho; 

	/* Needs fixing
	i_prim 	= eqtypn[3] + eqtypn[0] + eqtypn[1];
	i_cons	= eqtypi[4];
	for(i=0;i<eqtypn[4];++i) {
		S_prim[i_prim + i] = ((MATERGAMA[i]-1.0)*(S_cons[i_cons + i])/avf[i]) - MATERGAMA[i]*MATERPINF[i];
	}
	*/

}