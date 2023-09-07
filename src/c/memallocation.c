#include "strdata.h"

void memallocation(RUN *run){

	struct BRANCH * crnt;;

	solallocation(run->sol_L,run);
	solallocation(run->sol_R,run);

	crnt=run->topo->locl; 
	while (crnt!=NULL){
		leafallocation(run,crnt);
		crnt=crnt->lnxt;
	} 

}

void solallocation(SOL * sol, RUN *run){
    
	int v,i,j;

	sol->vec     = malloc(NMATERIALS*sizeof(double ));
	sol->wec     = malloc(NPRIMITIV*sizeof(double ));
	sol->rho     = malloc(NMATERIALS*sizeof(double ));
	sol->ra      = malloc(NMATERIALS*sizeof(double ));
	sol->Y       = malloc(NMATERIALS*sizeof(double ));
	sol->avf     = malloc(NMATERIALS*sizeof(double ));
	sol->are     = malloc(NMATERIALS*sizeof(double ));
	sol->p_hydro = malloc(NMATERIALS*sizeof(double ));
	sol->c2Y     = malloc(NMATERIALS*sizeof(double ));
	sol->c2      = malloc(NMATERIALS*sizeof(double ));
	sol->phi     = malloc(NMATERIALS*sizeof(double));

	sol->Amat = malloc(9*sizeof(double));
	sol->st   = malloc(3*sizeof(double *));
	
	sol->avf_stress=malloc(NMATERIALS*sizeof(double **));

	for (i=0;i<3;i++){
		sol->st[i]   = malloc(3*sizeof(double));
	}

	
	for (v=0;v<NMATERIALS;v++){
		sol->avf_stress[v]=malloc(3*sizeof(double *));
		for (i=0;i<3;i++){
			sol->avf_stress[v][i]=malloc(3*sizeof(double));
		}
	}

}