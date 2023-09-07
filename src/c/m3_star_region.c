#include "strdata.h"

void m3_star_region(struct RUN *run, struct BRANCH * crnt, struct SOL * sol,
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
 
  s11_S = ( (uR-SR)*run->sol_R->r*(s11L)  -  (uL-SL)*run->sol_L->r*(s11R)  +  (uL-SL)*run->sol_L->r * (uR-SR)*run->sol_R->r * (uR-uL) ) / 
          ( (uR-SR)*run->sol_R->r         -  (uL-SL)*run->sol_L->r );

  s12_S = 0.0;
  s13_S = 0.0;

  v_S = vK + (s12_S-(s12K)) / ((uK-SK) * sol->r);
  w_S = wK + (s13_S-(s13K)) / ((uK-SK) * sol->r);
  
  rotate_vector_ltg_v2(run,crnt,&ifc,&(u_S),&(v_S),&(w_S));

  re = sol->r * sol->e;
  Em_S = ( re*(uK - SK) - s11K*uK - s12K*vK - s13K*wK + s11_S*u_S + s12_S*v_S + s13_S*w_S ) / (r_S*(u_S - SK));
    
  // Compute * vector
  for(v=0;v<eqtypn[0];++v){ run->Us_temp[eqtypi[0]+v] = run->ar_S[v];           }
  for(v=0;v<eqtypn[1];++v){ run->Us_temp[eqtypi[1]+v] = r_S * run->vel_temp[v]; }
  for(v=0;v<eqtypn[2];++v){ run->Us_temp[eqtypi[2]+v] = r_S * Em_S;             }

  //___________________________________________________________
  run->u_nc = run->vel_temp[0];  // For residual and fluxes
  run->v_nc = run->vel_temp[1];
  run->w_nc = run->vel_temp[2];

  for(v=0;v<eqtypn[3];++v){ run->avf_S[v] = sol->avf[v]; }
  // Amat already has the value from tensor_temp ltg rotation, see line 149
  //___________________________________________________________

	free(tensor_temp);

  for(iv=0;iv<NCONSEQ;++iv){
    if ((isnan(run->Us_temp[iv])==1)){ 
      printf(" \n");
      printf("Compute Star Region nan | %d | %f %f \n",iv,run->Us_temp[iv],run->ar_S[iv]);
      printf(" \n");
      exit(0);
    }
  }

}