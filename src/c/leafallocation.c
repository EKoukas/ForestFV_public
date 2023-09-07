#include "strdata.h"

void leafallocation(struct RUN * run,struct BRANCH *crnt){

	int m,ord,typ,iv,ifc,ing,idr;

  crnt->el=malloc(sizeof(struct LEAF));
	crnt->el->crith=0;
  crnt->el->inner=-1;
	crnt->el->SN=malloc(NEQ*sizeof(double));
	crnt->el->S =malloc(NEQ*sizeof(double));
  crnt->el->SI=malloc(2*sizeof(int)); 	// +1 to track changes
	for(iv=0;iv<NEQ;iv++){ 
	  crnt->el->S[iv]=0.0; 
	}
  for(iv=0;iv<NEQ;iv++){ 
	  crnt->el->SN[iv]=0.0; 
	}
  for(iv=0;iv<2;iv++) { crnt->el->SI[iv] = 0.0; }
     

  
  crnt->el->WG = malloc(NEQ*sizeof(double *));
  for(iv=0;iv<NEQ;++iv){
    crnt->el->WG[iv] = malloc(3*sizeof(double));
    for(idr=0;idr<6;++idr){
      crnt->el->WG[iv][idr]=0.0;
    }
  }

  crnt->el->SF        = malloc(6*sizeof(double *));
  crnt->el->flux_face = malloc(6*sizeof(double *));
  for(ifc=0;ifc<6;++ifc){
    crnt->el->SF[ifc]        = malloc((NEQ)*sizeof(double ));
    crnt->el->flux_face[ifc] = malloc((NEQ)*sizeof(double ));
    for(iv=0;iv<NEQ;++iv){
      crnt->el->SF[ifc][iv]        = 0.0;
      crnt->el->flux_face[ifc][iv] = 0.0;  
    }
  }
   
	crnt->el->flux_flag = malloc(4*sizeof(int));
	for(ifc=0;ifc<6;++ifc){
	  for(m=0;m<3;++m){   
      crnt->el->flux_flag[ifc] = 0;   
		}
	}
  

  crnt->el->RHS=malloc(NEQ*sizeof(double ));
  for(iv=0;iv<NEQ;++iv){
    crnt->el->RHS[iv] = 0.0;
  }

  crnt->el->gfm_face = malloc(crnt->nlfc*sizeof(double));   
}