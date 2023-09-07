#include "strdata.h"

void m1_flux(struct SOL* sol) {
  
  // Continuity 
  g_flux_adv[0][0] = sol->u * sol->r; // (a*rho)*u   F
  g_flux_adv[0][1] = sol->v * sol->r; // (a*rho)*v   G
  g_flux_adv[0][2] = sol->w * sol->r; // (a*rho)*w   H
	
  // Momentum
  g_flux_adv[1][0] = sol->r*sol->u*sol->u + sol->st[0][0]; // rho*u^2 + p  F
  g_flux_adv[1][1] = sol->r*sol->u*sol->v + sol->st[0][1]; // rho*v*u      G
  g_flux_adv[1][2] = sol->r*sol->u*sol->w + sol->st[0][2]; // rho*w*u      H
  
  g_flux_adv[2][0] = sol->r*sol->v*sol->u + sol->st[1][0]; // rho*u*v
  g_flux_adv[2][1] = sol->r*sol->v*sol->v + sol->st[1][1]; // rho*v^2 + p
  g_flux_adv[2][2] = sol->r*sol->v*sol->w + sol->st[1][2]; // rho*w*v
  
  g_flux_adv[3][0] = sol->r*sol->w*sol->u + sol->st[2][0]; // rho*u*w
  g_flux_adv[3][1] = sol->r*sol->w*sol->v + sol->st[2][1]; // rho*v*w
  g_flux_adv[3][2] = sol->r*sol->w*sol->w + sol->st[2][2]; // rho*w^2 + p
    
}