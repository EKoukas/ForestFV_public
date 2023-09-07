#include "strdata.h"

void m2_star_region(struct RUN *run, struct BRANCH * crnt, struct SOL * sol,
                    int ifc,int flag_ifs,double S_S,double SL,double SR,
                    double uL,double vL,double wL,double uR,double vR,double wR) {

	double u_S,v_S,w_S,coNSt_S;
  double uK,vK,wK;
  double u,v,w;
  double SK;

  if (flag_ifs==1) {
    SK = SL; uK = uL; vK = vL; wK = wL;  
  } else if (flag_ifs==2) {
    SK = SR; uK = uR; vK = vR; wK = wR; 
  }

	coNSt_S = sol->r * (SK-uK)/(SK-S_S); 

  u= S_S;
  v= vK;
  w= wK;

  // Compute star vector
  // Toro-page 325
  run->Us_temp[0] = coNSt_S;           
  run->Us_temp[1] = coNSt_S *( u*crnt->cl->nx[ifc] + v*crnt->cl->nxt1[ifc] + w*crnt->cl->nxt2[ifc]);
	run->Us_temp[2] = coNSt_S *( u*crnt->cl->ny[ifc] + v*crnt->cl->nyt1[ifc] + w*crnt->cl->nyt2[ifc]);
	run->Us_temp[3] = coNSt_S *( u*crnt->cl->nz[ifc] + v*crnt->cl->nzt1[ifc] + w*crnt->cl->nzt2[ifc]);
  run->Us_temp[4] = coNSt_S *sol->ymass;
  	            
}