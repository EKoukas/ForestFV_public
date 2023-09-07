
#include "strdata.h"

void createutility(struct RUN *run){

  int i,j,v,f,ieq;
    
  g_flux_adv       = malloc(NEQ*sizeof(double *));
  g_flux_viscous   = malloc(NEQ*sizeof(double *));
  g_flux_HLLC      = malloc(NEQ*sizeof(double *));
	for (v=0;v<NEQ;++v){
    g_flux_adv[v]     = malloc(3*sizeof(double));
    g_flux_viscous[v] = malloc(3*sizeof(double));
  }

  g_S_cons_L     = malloc(NEQ       * sizeof(double));
  g_S_prim_L     = malloc(NPRIMITIV * sizeof(double));
  g_S_cons_R     = malloc(NEQ       * sizeof(double));
  g_S_prim_R     = malloc(NPRIMITIV * sizeof(double));
  g_delta_back   = malloc(NEQ_TEMP  * sizeof(double));
  g_delta_front  = malloc(NEQ_TEMP  * sizeof(double));

  g_grad_limited = malloc(NEQ_TEMP  * sizeof(double));
  g_S_prim_rec_L = malloc(NPRIMITIV * sizeof(double));
  g_S_cons_rec_L = malloc(NEQ       * sizeof(double));

  g_flux_sym_adv = malloc(NEQ       * sizeof(double *));
  for (v=0;v<NEQ;++v){
    g_flux_sym_adv[v] = malloc(3    * sizeof(double));
  }

  for (v=0;v<NEQ;++v){
    g_S_cons_L[v]     = 0.0;
    g_S_cons_R[v]     = 0.0;
    g_S_cons_rec_L[v] = 0.0;
  }

  for (v=0;v<NPRIMITIV;++v){
    g_S_prim_L[v]     = 0.0;
    g_S_prim_R[v]     = 0.0;
    g_S_prim_rec_L[v] = 0.0;
  }

  for (v=0;v<NEQ_TEMP;++v){
    g_delta_back[v]     = 0.0;
    g_delta_front[v]    = 0.0;
    g_grad_limited[v]   = 0.0;
  }


  // Auxilliary variables for residual
	run->Amat 				= malloc(3*sizeof(double *));
	run->Gmat 				= malloc(3*sizeof(double *));
	run->stress_solid	= malloc(3*sizeof(double *));
	run->sstensor 		= malloc(3*sizeof(double *));
	for(i=0;i<3;++i){
		run->Gmat[i]				 = malloc(3*sizeof(double));
		run->Amat[i] 				 = malloc(3*sizeof(double));
		run->stress_solid[i] = malloc(3*sizeof(double));
		run->sstensor[i] 		 = malloc(3*sizeof(double));
	}

  run->ra    = malloc(NMATERIALS*sizeof(double));
	run->avf   = malloc(NMATERIALS*sizeof(double));
	run->avf_c = malloc(NMATERIALS*sizeof(double));
	run->are   = malloc(NMATERIALS*sizeof(double));

	
  run->source_avf=malloc((NMATERIALS-1)*sizeof(double )); 
  for(v=0;v<(NMATERIALS-1);v++){
    run->source_avf[v]=0.0;
  }

  run->source_st=malloc(NMATERIALS*sizeof(double **)); 
  for(v=0;v<NMATERIALS;v++){
    run->source_st[v]=malloc(3*sizeof(double *));
    for(i=0;i<3;i++){
      run->source_st[v][i]=malloc(3*sizeof(double));
      for(j=0;j<3;j++){
        run->source_st[v][i][j]=0.0;
      }
    }
  }

    

  run->vel_temp    = malloc(3*sizeof(double)); // static
  run->tensor_temp = malloc(9*sizeof(double));
 
  // Auxilliary variables for hllc 
  run->Us_temp     = malloc(NEQ*sizeof(double));
  run->ar_S        = malloc(NMATERIALS*sizeof(double));
  run->rm_S        = malloc(NMATERIALS*sizeof(double));
  run->rm          = malloc(NMATERIALS*sizeof(double));
  run->e_S      = malloc(NMATERIALS*sizeof(double));
  run->avf_S = malloc(NMATERIALS*sizeof(double));
	run->are_S = malloc(NMATERIALS*sizeof(double));

  run->Amat_S = malloc(3*sizeof(double *));
  for (i=0;i<3;i++){
    run->Amat_S[i] = malloc(3*sizeof(double));
  }

  for(v=0;v<NCONSEQ;++v){    
		run->Us_temp[v] = 0.0;        
	}	

	for(v=0;v<NMATERIALS;++v){
		run->ar_S[v]  = 0.0;
		run->rm_S[v]  = 0.0;
		run->rm[v]    = 0.0;
		run->e_S[v]   = 0.0;
		run->avf_S[v] = 0.0;
		run->are_S[v] = 0.0;
	}
  
  // Auxilliary variables for relaxation
	run->Vec_temp_R = malloc(NEQ*sizeof(double));
	run->rhoKs  = malloc(NMATERIALS*sizeof(double));  //nmaterial
	run->pkinit = malloc(NMATERIALS*sizeof(double));  //nmaterial  // ok
	run->ar_R  	= malloc(NMATERIALS*sizeof(double));  //nmaterial  // ok
	run->avf_R 	= malloc(NMATERIALS*sizeof(double));  //nmaterial  // ok
	run->are_R 	= malloc(NMATERIALS*sizeof(double));  //nmaterial  // ok
	
	run->t0_m = malloc(NMATERIALS*sizeof(double));  //nmaterial // ok
	run->t1_m = malloc(NMATERIALS*sizeof(double));  //nmaterial // ok
	run->p0_m = malloc(NMATERIALS*sizeof(double));  //nmaterial // ok
	run->r_m  = malloc(NMATERIALS*sizeof(double));  //nmaterial // ok

  run->delta_temp=malloc(6 * sizeof(double *));
  for (f=0;f<6;f++) {
    run->delta_temp[f]=malloc(NEQ_TEMP * sizeof(double));
  }



  run->amat = malloc(3*sizeof(double *));
  run->bmat = malloc(3*sizeof(double *));
  run->cmat = malloc(3*sizeof(double *));
  run->xmat = malloc(3*sizeof(double *));
  
  for (i=0;i<3;i++){
    run->amat[i] = malloc(3*sizeof(double));
    run->bmat[i] = malloc(3*sizeof(double));
    run->cmat[i] = malloc(NEQ*sizeof(double));
    run->xmat[i] = malloc(NEQ*sizeof(double));
  }

}