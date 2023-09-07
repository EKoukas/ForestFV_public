#include "strdata.h"

void m5_RK(int irk,struct RUN *run) {

	int e,f,v,flag_rk;
	double buf,as_new;
  double phi_old[4];
	double phi_new[4];
	struct BRANCH * crnt;
	
	crnt=run->topo->locl;  	// RUN->FOREST->TREE (points to current tree)
	while (crnt!=NULL){ 	
		
		for(f=0; f<(crnt->nlfc); f++) {      
			crnt->el->flux_flag[f]=0;
		}
		
		for(v=0;v<eqtypn[3];++v) {   
    	phi_old[v] = crnt->el->S[eqtypi[3]+v] - AVF_LIM;  
    }

		switch(irk){

			case 0:
			for(e=0; e<NEQ; e++){ 	
			  crnt->el->SN[e] = crnt->el->S[e]; // SN:Solution vector 2, S:Solution vector
			}
			break;

			case 1:

			flag_rk=0;
			for(v=0;v<eqtypn[3];++v) {
				if (MATERMUSH[v]!=0.0) {
					e = eqtypi[3] + v;
					buf = (1.0/crnt->cl->Vol) * crnt->el->RHS[e];
					as_new = crnt->el->SN[e] + DT*buf;
					phi_new[v] = as_new - AVF_LIM;
					if (phi_old[v]*phi_new[v]<0.0) {
						flag_rk=1;
					}
				}
			}

			if (flag_rk==1) { 					
				for(e=0; e<NEQ; e++) {
					if (e>=eqtypi[2]) { 
						buf = (1.0/crnt->cl->Vol) * crnt->el->RHS[e];
						crnt->el->S[e] = crnt->el->SN[e] + DT * buf;
					}
				}
				m5_gfm(run,crnt);
					
			} else {

				for(e=0; e<NEQ; e++) {
					buf = (1.0/crnt->cl->Vol) * crnt->el->RHS[e];
					crnt->el->S[e] = crnt->el->SN[e] + DT * buf;							
				}
						
			}
			break;

			case 2:
			for(e=0;e<NEQ;e++){
	      crnt->el->S[e]=(3.0/4.0)*crnt->el->SN[e]+(1.0/4.0)*crnt->el->S[e]+(1.0/4.0)*DT*crnt->el->RHS[e]/crnt->cl->Vol;
			}
			break;

			case 3:
			for(e=0;e<NEQ;e++){
			  crnt->el->S[e]=(1.0/3.0)*crnt->el->SN[e]+(2.0/3.0)*crnt->el->S[e]+(2.0/3.0)*DT*crnt->el->RHS[e]/crnt->cl->Vol;
			}
			break;

			case 4:
			for(e=0;e<NEQ;e++){
			  crnt->el->S[e]=(1.0/2.0)*crnt->el->SN[e]+(1.0/2.0)*crnt->el->S[e]+(1.0/2.0)*DT*crnt->el->RHS[e]/crnt->cl->Vol;
			}
			break;

			default:
			break;
		}

		for(v=0;v<eqtypn[0];v++) {
			e = eqtypi[0] + v;
			if (crnt->el->S[e]<0.0) {
				crnt->el->S[e]=crnt->el->SN[e];
			}
		}

		if (ZEROU0==1) {
			e = 0 + eqtypi[1];
			crnt->el->S[e]=crnt->el->SN[e];

			e = eqtypi[5];
			crnt->el->S[e+1]=crnt->el->SN[e+1]; // A01, A2
			crnt->el->S[e+2]=crnt->el->SN[e+2]; // A02, A3
			crnt->el->S[e+3]=crnt->el->SN[e+3]; // A10, B1
			crnt->el->S[e+6]=crnt->el->SN[e+6]; // A20, C1
		}
		if (ZEROU1==1) {
			e = 1 + eqtypi[1];
			crnt->el->S[e]=crnt->el->SN[e];

			e = eqtypi[5];
			crnt->el->S[e+1]=crnt->el->SN[e+1]; // A01, A2
			crnt->el->S[e+3]=crnt->el->SN[e+3]; // A10, B1 
			crnt->el->S[e+5]=crnt->el->SN[e+5]; // A12, B3 
			crnt->el->S[e+7]=crnt->el->SN[e+7]; // A21, C2
		}
		if (ZEROU2==1) {
			e = 2 + eqtypi[1];
			crnt->el->S[e]=crnt->el->SN[e];

			e = eqtypi[5];
			crnt->el->S[e+2]=crnt->el->SN[e+2]; // A02, A3
			crnt->el->S[e+5]=crnt->el->SN[e+5]; // A12, B3
			crnt->el->S[e+6]=crnt->el->SN[e+6]; // A20, C1
			crnt->el->S[e+7]=crnt->el->SN[e+7]; // A21, C2
		}

  	crnt=crnt->lnxt; 
	}

}