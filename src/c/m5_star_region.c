#include "strdata.h"

void m5_star_region(struct RUN *run, struct BRANCH * crnt, struct SOL * sol,
                            int ifc,int flag_ifs,double u_S,double SL,double SR,
                            double uL,double vL,double wL,double uR,double vR,double wR,
                            double s11L,double s12L,double s13L,double s11R,double s12R,double s13R) {

	int i,j,iv,v;
	double v_S,w_S;
	double r_S,re,Em_S;
	double s11_S,s12_S,s13_S;
	double p_S,p_k;  
  double uK,vK,wK;
  double SK,s11K,s12K,s13K;
  double are_elast_S_K,e_elast_S_K;
  double * tensor_temp;
  struct SOL * t_sol_L; 
  struct SOL * t_sol_R;
  
  t_sol_L = run->sol_L;
  t_sol_R = run->sol_R;

  tensor_temp = malloc(9*sizeof(double));

  if (flag_ifs==1) {
    SK = SL; uK = uL; vK = vL; wK = wL;
    s11K = s11L; s12K = s12L; s13K = s13L;
  } else if (flag_ifs==2) {
    SK = SR; uK = uR; vK = vR; wK = wR;
    s11K = s11R; s12K = s12R; s13K = s13R;
  }

  r_S = 0.0;
  for(v=0;v<eqtypn[0];++v){
    run->ar_S[v]  = sol->ra[v] * (SK-uK)/(SK-u_S); 
    r_S          += run->ar_S[v];   
  }

  s11_S = ( (uR-SR)*t_sol_R->r*(s11L)  -  (uL-SL)*t_sol_L->r*(s11R)  +  (uL-SL)*t_sol_L->r * (uR-SR)*t_sol_R->r * (uR-uL) ) / 
          ( (uR-SR)*t_sol_R->r         -  (uL-SL)*t_sol_L->r );

  if ((t_sol_L->phi[0]>0.0)&&(t_sol_R->phi[0]>0.0)) {
    s12_S = ( (uR-SR)*t_sol_R->r*(s12L)  -  (uL-SL)*t_sol_L->r*(s12R)  +  (uL-SL)*t_sol_L->r * (uR-SR)*t_sol_R->r * (vR-vL) ) /
            ( (uR-SR)*t_sol_R->r         -  (uL-SL)*t_sol_L->r );

    s13_S = ( (uR-SR)*t_sol_R->r*(s13L)  -  (uL-SL)*t_sol_L->r*(s13R)  +  (uL-SL)*t_sol_L->r * (uR-SR)*t_sol_R->r * (wR-wL) ) /
            ( (uR-SR)*t_sol_R->r         -  (uL-SL)*t_sol_L->r );
  } else {
    s12_S = 0.0;
    s13_S = 0.0;
  }

  v_S = vK + (s12_S-s12K) / ((uK-SK) * sol->r);
  w_S = wK + (s13_S-s13K) / ((uK-SK) * sol->r);

  rotate_vector_ltg_v2(run,crnt,&ifc,&(u_S),&(v_S),&(w_S));
    
  re = sol->r * sol->e;
  Em_S = ( re*(uK - SK) + s11K*uK - s12K*vK - s13K*wK + s11_S*u_S + s12_S*v_S + s13_S*w_S ) / (r_S*(u_S + SK));
  
  // Amat, global to local
  v=0;
  for(i=0;i<3;++i){ 
    for(j=0;j<3;++j){ 
      run->Amat[i][j] = sol->Amat[v];
      ++v;
    }
  }

  rotate_tensor(tensor_temp,run->Amat,crnt->cl->nx[ifc]  ,crnt->cl->ny[ifc]  ,crnt->cl->nz[ifc],
                                      crnt->cl->nxt1[ifc],crnt->cl->nyt1[ifc],crnt->cl->nzt1[ifc],
                                      crnt->cl->nxt2[ifc],crnt->cl->nyt2[ifc],crnt->cl->nzt2[ifc],
                                      VERBOSE,run->con->rank);
  
  v=0;
  for(i=0;i<3;++i){ 
    for(j=0;j<3;++j){ 
      run->Amat[i][j] = tensor_temp[v];
      ++v;
    }
  }
  
  run->Amat_S[0][0] = (run->Amat[0][0]*(uK-SK) - run->Amat[1][0]*(vK-v_S) + run->Amat[2][0]*(wK-w_S) ) / (u_S-SK);
  run->Amat_S[0][1] = (run->Amat[0][1]*(uK-SK) - run->Amat[1][1]*(vK-v_S) + run->Amat[2][1]*(wK-w_S) ) / (u_S-SK);
  run->Amat_S[0][2] = (run->Amat[0][2]*(uK-SK) - run->Amat[1][2]*(vK-v_S) + run->Amat[2][2]*(wK-w_S) ) / (u_S-SK);
  
  run->Amat_S[1][0] = run->Amat[1][0];
  run->Amat_S[1][1] = run->Amat[1][1];
  run->Amat_S[1][2] = run->Amat[1][2];

  run->Amat_S[2][0] = run->Amat[2][0];
  run->Amat_S[2][1] = run->Amat[2][1];
  run->Amat_S[2][2] = run->Amat[2][2];

  // Amat, local to global
  rotate_tensor(tensor_temp,run->Amat_S,crnt->cl->nx[ifc], crnt->cl->nxt1[ifc], crnt->cl->nxt2[ifc],
          		                          crnt->cl->ny[ifc], crnt->cl->nyt1[ifc], crnt->cl->nyt2[ifc],
          				                      crnt->cl->nz[ifc], crnt->cl->nzt1[ifc], crnt->cl->nzt2[ifc],
          				                      VERBOSE,run->con->rank);
  
  v=0;
  for(i=0;i<3;++i){  
    for(j=0;j<3;++j){ 
      run->Amat_S[i][j] = tensor_temp[v];
      ++v;
    }
  }

  for (i=0;i<3;++i){
    for (j=0;j<3;++j){
      run->Gmat[i][j] = run->Amat_S[i][0]*run->Amat_S[j][0] + 
                        run->Amat_S[i][1]*run->Amat_S[j][1] + 
                        run->Amat_S[i][2]*run->Amat_S[j][2];
    }
  }
  
  // Internal Energies
  for(v=0;v<eqtypn[4];++v){
    run->rm[v]    = (sol->ra[v]/sol->avf[v]);
    run->rm_S[v]  = (sol->ra[v]/sol->avf[v]) * (uK-SK)/(u_S-SK);

    p_k = -sol->p_hydro[v];
    p_S = (p_k+MATERPINF[v])*( ( (MATERGAMA[v]-1.0)*run->rm[v]   + (MATERGAMA[v]+1.0)*run->rm_S[v] ) / 
                               ( (MATERGAMA[v]-1.0)*run->rm_S[v] + (MATERGAMA[v]+1.0)*run->rm[v]   )  ) - MATERPINF[v];
    
    are_elast_S_K = m5_elastic_energy_solid(run->Gmat,MATERMUSH[v]); 
      e_elast_S_K = are_elast_S_K/(sol->avf[v]*run->rm_S[v]);
    
    run->e_S[v]=(p_S + MATERGAMA[v]*MATERPINF[v]) / (run->rm_S[v]*(MATERGAMA[v]+1.0)) + e_elast_S_K;

  }

  // Compute * vector
  for(v=0;v<eqtypn[0];++v){ run->Us_temp[eqtypi[0]+v] = run->ar_S[v];           }
  for(v=0;v<eqtypn[1];++v){ run->Us_temp[eqtypi[1]+v] = r_S * run->vel_temp[v]; }
  for(v=0;v<eqtypn[2];++v){ run->Us_temp[eqtypi[2]+v] = r_S * Em_S;             }
 
  //___________________________________________________________
  run->u_nc = run->vel_temp[0];  // For residual and fluxes
  run->v_nc = run->vel_temp[1];
  run->w_nc = run->vel_temp[2];
  
  for(v=0;v<eqtypn[3];++v){ run->avf_S[v]  = sol->avf[v];                }
  for(v=0;v<eqtypn[4];++v){ run->are_S[v]  = run->ar_S[v] * run->e_S[v]; }
  // Amat already has the value from tensor_temp ltg rotation, see line 92
  //___________________________________________________________

	free(tensor_temp);
  
  for(iv=0;iv<NCONSEQ;++iv){
    if ((isnan(run->Us_temp[iv])==1)){ 
      printf(" \n");
      printf("Compute Star Region nan | %d | %f %f %f \n",iv,run->Us_temp[iv],run->ar_S[iv],run->e_S[iv]);
      printf(" \n");
      exit(0);
    }
  }

}