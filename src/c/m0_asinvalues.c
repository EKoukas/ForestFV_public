#include "strdata.h"

void m0_asinvalues(struct RUN *run,struct BRANCH * crnt,struct SOL * r_sol,double rho,double ru,double rv, double rw, double e) {

  int iv,i;

  r_sol->r = rho;
  r_sol->u = ru / rho;
  r_sol->v = rv / rho;
  r_sol->w = rw / rho;
  
  r_sol->vec[0] = rho;
  r_sol->vec[1] = ru;
  r_sol->vec[2] = rv;
  r_sol->vec[3] = rw;
  
  r_sol->e          = e;
  r_sol->vec[4]     = e;
  r_sol->e_internal = (r_sol->e/r_sol->r) - 0.5*(pow(r_sol->u,2.0) + pow(r_sol->v,2.0) + pow(r_sol->w,2.0));
  r_sol->p_hydro[0] = r_sol->r*(MATERGAMA[0]-1.0)*r_sol->e_internal - MATERGAMA[0]*MATERPINF[0];
  r_sol->c          = sqrt(MATERGAMA[0]*(r_sol->p_hydro[0] + MATERPINF[0])/r_sol->r);
  
  r_sol->st[0][0] = r_sol->p_hydro[0]; 
  r_sol->st[0][1] = 0.0;
  r_sol->st[0][2] = 0.0;
           
  r_sol->st[1][0] = 0.0;
  r_sol->st[1][1] = r_sol->p_hydro[0]; 
  r_sol->st[1][2] = 0.0;
  
  r_sol->st[2][0] = 0.0;
  r_sol->st[2][1] = 0.0;
  r_sol->st[2][2] = r_sol->p_hydro[0];

}