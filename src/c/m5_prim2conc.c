#include "strdata.h"

void m5_prim2conc(double * S_prim,double * S_cons) {

	int i_prim,i_cons,i,j;
	double avfsum,rho,sumare,are_elast,avf_dif,avf_per,Gmag;
	double avf[10];
	double Amat[9];
	double Gmat[9];

	//===============================================
	// vf -> vf
	// i_prim = 0;
	i_cons = eqtypi[3];
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
	
	for(i=0;i<eqtypn[3];++i){
		S_cons[i_cons + i] = avf[i];
	}
	//===============================================


	//===============================================
	// rho_1 ...
	rho = 0.0;
	i_prim = eqtypn[3];
	// i_cons = 0
	for(i=0;i<eqtypn[0];++i){ 
		S_cons[i]  = S_prim[i_prim + i];
		rho       += S_prim[i_prim + i];
	}
	//===============================================

	//===============================================
	// u,v,w -> ru rv rw
	i_prim = eqtypn[3] + eqtypn[0];
	i_cons = eqtypi[1];
	S_cons[i_cons + 0] = rho*S_prim[i_prim + 0];
	S_cons[i_cons + 1] = rho*S_prim[i_prim + 1];
	S_cons[i_cons + 2] = rho*S_prim[i_prim + 2];
	//===============================================

	//===============================================
	// Amat -> Amat
	i_cons = eqtypi[5];
	i_prim = eqtypn[3] + eqtypn[0] + eqtypn[1] + eqtypn[4];
	for(i=0;i<eqtypn[5];++i){
		S_cons[i_cons + i] = S_prim[i_prim + i];
		Amat[i] 					 = S_prim[i_prim + i];
	}
	//===============================================

	//===============================================
	// p->are
	Gmat[0] = Amat[0]*Amat[0] + Amat[1]*Amat[1] + Amat[2]*Amat[2];
  Gmat[1] = Amat[0]*Amat[3] + Amat[1]*Amat[4] + Amat[2]*Amat[5];
  Gmat[2] = Amat[0]*Amat[6] + Amat[1]*Amat[7] + Amat[2]*Amat[8];

  Gmat[3] = Amat[3]*Amat[0] + Amat[4]*Amat[1] + Amat[5]*Amat[2];
  Gmat[4] = Amat[3]*Amat[3] + Amat[4]*Amat[4] + Amat[5]*Amat[5];
  Gmat[5] = Amat[3]*Amat[6] + Amat[4]*Amat[7] + Amat[5]*Amat[8];

  Gmat[6] = Amat[6]*Amat[0] + Amat[7]*Amat[1] + Amat[8]*Amat[2];
  Gmat[7] = Amat[6]*Amat[3] + Amat[7]*Amat[4] + Amat[8]*Amat[5];
  Gmat[8] = Amat[6]*Amat[6] + Amat[7]*Amat[7] + Amat[8]*Amat[8];

	Gmag = Gmat[0]*(Gmat[4]*Gmat[8] - Gmat[5]*Gmat[7]) - Gmat[1]*(Gmat[3]*Gmat[8] - Gmat[5]*Gmat[6]) + Gmat[2]*(Gmat[3]*Gmat[7] - Gmat[4]*Gmat[6]);
  
	sumare = 0.0;
	i_cons = eqtypi[4];
	i_prim = eqtypn[3] + eqtypn[0] + eqtypn[1];
	for(i=0;i<eqtypn[4];++i){
		
		if (MATERMUSH[i]!=0.0) { are_elast=m5_elastic_energy_solid_V2(Gmat,Gmag,MATERMUSH[i]); } 
		else {                   are_elast=0.0;    																	      		 }
	
		S_cons[i_cons + i] = (avf[i]*(S_prim[i_prim + i] + MATERGAMA[i]*MATERPINF[i]) / (MATERGAMA[i]-1.0)) + are_elast;		
		sumare 						+= S_cons[i_cons + i];
	}
	//===============================================

	//===============================================
	// E
	i_prim = eqtypn[3] + eqtypn[0];
	i_cons = eqtypi[2];
	S_cons[i_cons] = sumare + 0.5*rho*(pow(S_prim[i_prim+0],2.0) + pow(S_prim[i_prim+1],2.0) + pow(S_prim[i_prim+2],2.0));
	//===============================================

}