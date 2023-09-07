//-------------------------------------------------------
// Viscous Fluxes, page 16 Toro's book
//-------------------------------------------------------

#include "strdata.h"

void flux_viscous(struct RUN *run,struct BRANCH * crnt,struct SOL* sol_L,struct SOL* sol_R) {

  double temp,visc;

  double ux,uy,uz;
  double vx,vy,vz;
  double wy,wx,wz;

	double txx,txy,txz;
	double tyx,tyy,tyz;
	double tzx,tzy,tzz;

  double uf,vf,wf;
  
  // velocites on the face
  uf = 0.5*(sol_L->u + sol_R->u);
  vf = 0.5*(sol_L->v + sol_R->v);
  wf = 0.5*(sol_L->w + sol_R->w);
  
  // velocity gradients on the face
  ux = 0.5*(sol_L->ux + sol_R->ux);
  uy = 0.5*(sol_L->uy + sol_R->uy);
  uz = 0.5*(sol_L->uz + sol_R->uz);

  vx = 0.5*(sol_L->vx + sol_R->vx);
  vy = 0.5*(sol_L->vy + sol_R->vy);
  vz = 0.5*(sol_L->vz + sol_R->vz);

  wx = 0.5*(sol_L->wx + sol_R->wx);
  wy = 0.5*(sol_L->wy + sol_R->wy);
  wz = 0.5*(sol_L->wz + sol_R->wz);
  
  visc = MATERVISC[0];
  
	txx = visc*( 2.0*ux - (2.0/3.0)* (ux+vy+wz) );
	tyy = visc*( 2.0*vy - (2.0/3.0)* (ux+vy+wz) );
	tzz = visc*( 2.0*wz - (2.0/3.0)* (ux+vy+wz) );

	txy = visc * (uy+vx);
	txz = visc * (uz+wx);
	tyz = visc * (vz+wy);

	tyx = txy;
	tzy = tyz;
	tzx = txz;

  // Continuity 
  g_flux_viscous[0][0] = 0.0; 
  g_flux_viscous[0][1] = 0.0;
  g_flux_viscous[0][2] = 0.0;
  
  // Momentum
  g_flux_viscous[1][0] = txx; // F
  g_flux_viscous[1][1] = tyx; // G
  g_flux_viscous[1][2] = tzx; // H
  
  g_flux_viscous[2][0] = txy; // 
  g_flux_viscous[2][1] = tyy; //
  g_flux_viscous[2][2] = tzy; // 

  g_flux_viscous[3][0] = txz; // 
  g_flux_viscous[3][1] = tyz; // 
  g_flux_viscous[3][2] = tzz; // 

  if(MODEL==1){
    // Total Energy
    g_flux_viscous[4][0] = uf*txx + vf*txy + wf*txz; 
    g_flux_viscous[4][1] = uf*tyx + vf*tyy + wf*tyz;  
    g_flux_viscous[4][2] = uf*tzx + vf*tzy + wf*tzz;  
  }

}