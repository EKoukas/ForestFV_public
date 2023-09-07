#include "strdata.h"

void m3_flux(struct RUN *run,struct SOL* sol,int flag_ifs) {

  int iv,v;
  
  //-------------------------------------------------------
  // Vector F1 of the coNServative fluxes
  //-------------------------------------------------------

  // Continuity          
  for(v=0; v <eqtypn[0]; v++){
    iv = v + eqtypi[0]; 
    g_flux_adv[iv][0] = sol->u * sol->ra[v]; // (a*rho)*u
    g_flux_adv[iv][1] = sol->v * sol->ra[v]; // (a*rho)*v
    g_flux_adv[iv][2] = sol->w * sol->ra[v]; // (a*rho)*w
  }
	
  // Momentum
  iv = 0 + eqtypi[1]; // 0+2
  g_flux_adv[iv][0] = sol->u * sol->u * sol->r - sol->st[0][0]; // rho*u^2-sigma_00
  g_flux_adv[iv][1] = sol->v * sol->u * sol->r - sol->st[0][1]; // rho*v*u-sigma_01
  g_flux_adv[iv][2] = sol->w * sol->u * sol->r - sol->st[0][2]; // rho*w*u-sigma_02
 
  iv = 1 + eqtypi[1]; // 1+2
  g_flux_adv[iv][0] = sol->u * sol->v * sol->r - sol->st[1][0]; // rho*u*v-sigma_10
  g_flux_adv[iv][1] = sol->v * sol->v * sol->r - sol->st[1][1]; // rho*v^2-sigma_11
  g_flux_adv[iv][2] = sol->w * sol->v * sol->r - sol->st[1][2]; // rho*w*v-sigma_12

  iv = 2 + eqtypi[1]; // 2+2
  g_flux_adv[iv][0] = sol->u * sol->w * sol->r - sol->st[2][0]; // rho*u*w-sigma_20
  g_flux_adv[iv][1] = sol->v * sol->w * sol->r - sol->st[2][1]; // rho*v*w-sigma_21
  g_flux_adv[iv][2] = sol->w * sol->w * sol->r - sol->st[2][2]; // rho*w^2-sigma_22
	
  // Total Energy
  iv=eqtypi[2]; // 0+5 
  g_flux_adv[iv][0] = sol->r*sol->e * sol->u - sol->st[0][0] * sol->u - sol->st[0][1] * sol->v - sol->st[0][2] * sol->w;  // r*e*u - (u,v,w)(sigma_00,sigma_01,sigma_02)
  g_flux_adv[iv][1] = sol->r*sol->e * sol->v - sol->st[1][0] * sol->u - sol->st[1][1] * sol->v - sol->st[1][2] * sol->w;  // r*e*v - (u,v,w)(sigma_10,sigma_11,sigma_12)
  g_flux_adv[iv][2] = sol->r*sol->e * sol->w - sol->st[2][0] * sol->u - sol->st[2][1] * sol->v - sol->st[2][2] * sol->w;  // r*e*w - (u,v,w)(sigma_20,sigma_21,sigma_22)

  //-------------------------------------------------------
  // Vector F2 of the non-coNServative fluxes
  //-------------------------------------------------------

  // Only for right and left state, differently it  
  // gets the values from m4_star_region.c
  if ((flag_ifs==0) || (flag_ifs==3))  {
    run->u_nc = sol->u; 
    run->v_nc = sol->v;  
    run->w_nc = sol->w;
    for(v=0;v<eqtypn[3];++v){ run->avf_S[v] = sol->avf[v]; }
  } else if (flag_ifs==4) {
    for(v=0;v<eqtypn[3];++v){ run->avf_S[v] = sol->avf[v]; }
  }

  // Volume fractioNS;
  for(v=0; v<eqtypn[3]; ++v){
    iv = v + eqtypi[3];            
    g_flux_adv[iv][0] = run->u_nc * run->avf_S[v];  
    g_flux_adv[iv][1] = run->v_nc * run->avf_S[v];  
    g_flux_adv[iv][2] = run->w_nc * run->avf_S[v];  
  }

}