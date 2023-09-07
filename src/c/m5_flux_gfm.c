#include "strdata.h"

void m5_flux_gfm(struct RUN *run,struct SOL* sol,int flag_ifs) {

  int i,j,v,iv;
  
  //-------------------------------------------------------
  // Vector F1 of the coNServative fluxes
  //-------------------------------------------------------

  // Continuity         
  for(v=0; v <eqtypn[0]; ++v){
    iv = v + eqtypi[0]; 
    g_flux_adv[iv][0] = sol->u_temp * sol->ra[v]; // (a*rho)*u
    g_flux_adv[iv][1] = sol->v_temp * sol->ra[v]; // (a*rho)*v
    g_flux_adv[iv][2] = sol->w_temp * sol->ra[v]; // (a*rho)*w
  }
	
  // Momentum
  iv = 0 + eqtypi[1];
  g_flux_adv[iv][0] = sol->u_temp * sol->u_temp * sol->r - sol->st[0][0]; // rho*u^2-sigma_00
  g_flux_adv[iv][1] = sol->v_temp * sol->u_temp * sol->r - sol->st[0][1]; // rho*v*u-sigma_01
  g_flux_adv[iv][2] = sol->w_temp * sol->u_temp * sol->r - sol->st[0][2]; // rho*w*u-sigma_02
  
  iv = 1 + eqtypi[1]; 
  g_flux_adv[iv][0] = sol->u_temp * sol->v_temp * sol->r - sol->st[1][0]; // rho*u*v-sigma_10
  g_flux_adv[iv][1] = sol->v_temp * sol->v_temp * sol->r - sol->st[1][1]; // rho*v^2-sigma_11
  g_flux_adv[iv][2] = sol->w_temp * sol->v_temp * sol->r - sol->st[1][2]; // rho*w*v-sigma_12

  iv = 2 + eqtypi[1]; 
  g_flux_adv[iv][0] = sol->u_temp * sol->w_temp * sol->r - sol->st[2][0]; // rho*u*w-sigma_20
  g_flux_adv[iv][1] = sol->v_temp * sol->w_temp * sol->r - sol->st[2][1]; // rho*v*w-sigma_21
  g_flux_adv[iv][2] = sol->w_temp * sol->w_temp * sol->r - sol->st[2][2]; // rho*w^2-sigma_22
	
  // Total Energy  
  iv=eqtypi[2]; 
  g_flux_adv[iv][0] = sol->r*sol->e * sol->u_temp - sol->st[0][0]*sol->u_temp - sol->st[0][1]*sol->v_temp - sol->st[0][2]*sol->w_temp;  // r*e*u - (u,v,w)(sigma_00,sigma_01,sigma_02)
  g_flux_adv[iv][1] = sol->r*sol->e * sol->v_temp - sol->st[1][0]*sol->u_temp - sol->st[1][1]*sol->v_temp - sol->st[1][2]*sol->w_temp;  // r*e*v - (u,v,w)(sigma_10,sigma_11,sigma_12)
  g_flux_adv[iv][2] = sol->r*sol->e * sol->w_temp - sol->st[2][0]*sol->u_temp - sol->st[2][1]*sol->v_temp - sol->st[2][2]*sol->w_temp;  // r*e*w - (u,v,w)(sigma_20,sigma_21,sigma_22)

  //-------------------------------------------------------
  // Vector F2 of the non-conservative fluxes
  //-------------------------------------------------------

  // Only for right and left state, differently it  
  // gets the values from m5_star_region.c
  if ((flag_ifs==0) || (flag_ifs==3))  {
    run->u_nc = sol->u_temp; 
    run->v_nc = sol->v_temp;  
    run->w_nc = sol->w_temp;
    
    for(v=0;v<eqtypn[3];++v){ run->avf_S[v]  = sol->avf[v]; }
    for(v=0;v<eqtypn[4];++v){ run->are_S[v]  = sol->are[v]; }
    for (i=0;i<3;++i){ 
      for (j=0;j<3;++j){ 
        run->Amat_S[i][j] = sol->Amat[i*3+j]; 
      } 
    }
  }
  
  // Volume fractions;
  for(v=0; v<eqtypn[3]; ++v){
    iv = v + eqtypi[3];             
    g_flux_adv[iv][0] = run->u_nc * run->avf_S[v];  
    g_flux_adv[iv][1] = run->v_nc * run->avf_S[v];  
    g_flux_adv[iv][2] = run->w_nc * run->avf_S[v];  
  }

  // Specific energy;
  for(v=0; v<eqtypn[4]; ++v){  
    iv = v + eqtypi[4];           
    g_flux_adv[iv][0] = run->u_nc * run->are_S[v];
    g_flux_adv[iv][1] = run->v_nc * run->are_S[v];
    g_flux_adv[iv][2] = run->w_nc * run->are_S[v];
  }
  
  // Deformation matrix;
  iv = eqtypi[5];
  for (i=0;i<3;++i){ 
    for (j=0;j<3;++j){
      g_flux_adv[iv][0] = run->u_nc * run->Amat_S[i][j]; 
      g_flux_adv[iv][1] = run->v_nc * run->Amat_S[i][j]; 
      g_flux_adv[iv][2] = run->w_nc * run->Amat_S[i][j];
      iv++;
    }
  }

}