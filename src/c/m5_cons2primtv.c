#include "strdata.h"

void m5_cons2primtv(double * S_cons,double * S_prim) {

	int i_prim,i_cons,i,i_blakeia;
	double rho,avfsum,are_elast,Gmag;
	double avf[10];
	double Amat[9];
	double Gmat[9];

	//===============================================
	// vf -> vf
	// i_prim = 0;
	avfsum = 0.0;
	i_cons = eqtypi[3];
	for(i=0;i<eqtypn[3];++i){
		S_prim[i] = S_cons[i_cons +i];
		avf[i] 	  = S_cons[i_cons +i];
		avfsum   += avf[i]; 
	}
	avf[eqtypn[3]] = 1.0 - avfsum;
	//===============================================

	//===============================================
	// rho_1 ...
	rho		 =	0.0;
	i_prim = eqtypn[3];
	// i_cons = 0
	for(i=0;i<eqtypn[0];++i){
		S_prim[i_prim + i] = S_cons[i];
		rho								+= S_cons[i];
	}
	//===============================================

	//===============================================
	// ru rv rw -> u,v,w
	i_prim = eqtypn[3] + eqtypn[0];
	i_cons = eqtypi[1];
	S_prim[i_prim + 0] = S_cons[i_cons + 0]/rho; 
	S_prim[i_prim + 1] = S_cons[i_cons + 1]/rho; 
	S_prim[i_prim + 2] = S_cons[i_cons + 2]/rho; 
	//===============================================

	//===============================================
	// Amat -> Amat
	i_cons = eqtypi[5];
	i_prim = eqtypn[3] + eqtypn[0] + eqtypn[1] + eqtypn[4];
	for(i=0;i<eqtypn[5];++i) {
		S_prim[i_prim + i] = S_cons[i_cons + i];
		Amat[i] 					 = S_cons[i_cons + i];
	}
	//===============================================
		
	//===============================================
	// are -> p
	Gmat[0] = Amat[0]*Amat[0] + Amat[1]*Amat[1] + Amat[2]*Amat[2];
  Gmat[1] = Amat[0]*Amat[3] + Amat[1]*Amat[4] + Amat[2]*Amat[5];
  Gmat[2] = Amat[0]*Amat[6] + Amat[1]*Amat[7] + Amat[2]*Amat[8];

  Gmat[3] = Amat[3]*Amat[0] + Amat[4]*Amat[1] + Amat[5]*Amat[2];
  Gmat[4] = Amat[3]*Amat[3] + Amat[4]*Amat[4] + Amat[5]*Amat[5];
  Gmat[5] = Amat[3]*Amat[6] + Amat[4]*Amat[7] + Amat[5]*Amat[8];

  Gmat[6] = Amat[6]*Amat[0] + Amat[7]*Amat[1] + Amat[8]*Amat[2];
  Gmat[7] = Amat[6]*Amat[3] + Amat[7]*Amat[4] + Amat[8]*Amat[5];
  Gmat[8] = Amat[6]*Amat[6] + Amat[7]*Amat[7] + Amat[8]*Amat[8];

	Gmag = Gmat[0]*(Gmat[4]*Gmat[8] - Gmat[5]*Gmat[7]) - 
				 Gmat[1]*(Gmat[3]*Gmat[8] - Gmat[5]*Gmat[6]) + 
				 Gmat[2]*(Gmat[3]*Gmat[7] - Gmat[4]*Gmat[6]);
  
	if (isnan(Gmag)==1) { 
		printf("nan %f \n",Gmag); 
		for(i=0;i<eqtypn[5];++i) {
			printf("Gmat,Amat %d %f %f \n",i,Gmat[i],Amat[i]); 
		}
		exit(0); 
	}

	i_prim 	= eqtypn[3] + eqtypn[0] + eqtypn[1];
	i_cons	= eqtypi[4];
	for(i=0;i<eqtypn[4];++i) {


		if (MATERMUSH[i]!=0.0) { are_elast=m5_elastic_energy_solid_V2(Gmat,Gmag,MATERMUSH[i]); } 
		else {                   are_elast=0.0;    																             }

		S_prim[i_prim + i] = ((MATERGAMA[i]-1.0)*(S_cons[i_cons + i] - are_elast)/avf[i]) - MATERGAMA[i]*MATERPINF[i];
	}
	//===============================================

}