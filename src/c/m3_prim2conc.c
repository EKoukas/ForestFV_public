#include "strdata.h"

void m3_prim2conc(double * S_prim,double * S_cons) {

	int iv,j,ivv,v;
	double avfsum,rho,sumare,avf_dif,avf_per;
	double avf[10];
	
	// avf
  avfsum=0.0;
	for(iv=0;iv<eqtypn[3];++iv){
		avf[iv] = S_prim[iv];
	
		if (avf[iv]<= AMIN/10.0)              { avf[iv] =           AMIN/10.0;   } 
		if (avf[iv]>=(1.0-(2.0*(AMIN/10.0)))) { avf[iv] = 1.0-(2.0*(AMIN/10.0)); }    
		avfsum  += avf[iv];
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
	ivv = eqtypn[3];
	v   = eqtypi[0];
	for(iv=0;iv<eqtypn[0];++iv){ 
		S_cons[iv + v] = S_prim[ivv + iv];
		rho           += S_prim[ivv + iv];
	}


	// ru rv rw
	ivv =eqtypn[3] + eqtypn[0];
	v = eqtypi[1];
	for(iv=0;iv<eqtypn[1];++iv){
		S_cons[iv + v] = rho*S_prim[iv + ivv];
	}


	//  vf
	ivv=eqtypi[3];
	for(iv=0;iv<eqtypn[3];++iv){
		S_cons[iv + ivv] = S_prim[iv];
	}


	// E
	v = eqtypn[3] + eqtypn[0];
	S_cons[eqtypi[2]] = sumare + 0.5*rho*(pow(S_prim[v+0],2.0) + pow(S_prim[v+1],2.0) + pow(S_prim[v+2],2.0));
			
}