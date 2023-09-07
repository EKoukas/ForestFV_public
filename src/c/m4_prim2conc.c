#include "strdata.h"

void m4_prim2conc(double * S_prim,double * S_cons) {

	int i,j;
	int i_prim;
	int i_cons;
	double avfsum,rho,sumare,avf_dif,avf_per;
	double avf[10];
	
	// avf
	//i_prim = 0;
  avfsum=0.0;
	for(i=0;i<eqtypn[3];++i){
		avf[i] = S_prim[i]; // i_prim + i
	
		if (avf[i]<= AMIN/10.0)              { avf[i] =           AMIN/10.0;   } 
		if (avf[i]>=(1.0-(2.0*(AMIN/10.0)))) { avf[i] = 1.0-(2.0*(AMIN/10.0)); }    
		avfsum  += avf[i];
	}

	if (avfsum>=(1.0-(AMIN/100.0))) { 
    avf_dif = avfsum - (1.0-(AMIN/100.0));
    for(j=0;j<eqtypn[3];++j){
      avf_per = avf[j]/avfsum;
      avf[j]  = avf[j] - avf_per*avf_dif;
    }
    avfsum = 1.0-(AMIN/100.0);
	}
	avf[eqtypn[3]] = 1.0 - avfsum;


	// ar
	rho = 0.0;
	i_prim = eqtypn[3];
	//i_cons = eqtypi[0]; =0
	for(i=0;i<eqtypn[0];++i){ 
		S_cons[i]  = S_prim[i_prim + i];
		rho       += S_prim[i_prim + i];
	}


	// ru rv rw
	i_prim = eqtypn[3] + eqtypn[0];
	i_cons = eqtypi[1];
	for(i=0;i<eqtypn[1];++i){
		S_cons[i_cons + i] = rho*S_prim[i_prim + i];
	}


	//  vf
	// i_prim = 0
	i_cons = eqtypi[3];
	for(i=0;i<eqtypn[3];++i){
		S_cons[i_cons + i] = avf[i];
	}


	// are
	sumare = 0.0;
	i_cons = eqtypi[4];
	i_prim = eqtypn[3] + eqtypn[0] + eqtypn[1];
	for(i=0;i<eqtypn[4];++i){
		S_cons[i_cons + i] = (avf[i]*(S_prim[i_prim + i] + MATERGAMA[i]*MATERPINF[i]) / (MATERGAMA[i]-1.0));	
		sumare 						+= S_cons[i_cons + i];
	}


	// E
	i_prim = eqtypn[3] + eqtypn[0];
	i_cons = eqtypi[2];
	S_cons[i_cons] = sumare + 0.5*rho*(pow(S_prim[i_prim+0],2.0) + pow(S_prim[i_prim+1],2.0) + pow(S_prim[i_prim+2],2.0));
			
}