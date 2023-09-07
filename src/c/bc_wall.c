#include "strdata.h"

void bc_wall(struct BRANCH * crnt,int f,struct RUN *run) {

  double u,v,w,u_sol_L_lc,v_sol_L_lc,w_sol_L_lc;
  double u_sol_R_lc,v_sol_R_lc,w_sol_R_lc;

  // Rotate global velocities to local coordinate system
  u = run->sol_L->u;
  v = run->sol_L->v;
  w = run->sol_L->w;

  u_sol_L_lc = u*crnt->cl->nx[f]   + v*crnt->cl->ny[f]   + w*crnt->cl->nz[f];
  v_sol_L_lc = u*crnt->cl->nxt1[f] + v*crnt->cl->nyt1[f] + w*crnt->cl->nzt1[f];
  w_sol_L_lc = u*crnt->cl->nxt2[f] + v*crnt->cl->nyt2[f] + w*crnt->cl->nzt2[f];

  // Opposite velocity 
  u_sol_R_lc = -u_sol_L_lc;
  v_sol_R_lc = -v_sol_L_lc;
  w_sol_R_lc = -w_sol_L_lc;

  // Rotate back to global
  run->sol_R->u = u_sol_R_lc*crnt->cl->nx[f] + v_sol_R_lc*crnt->cl->nxt1[f] + w_sol_R_lc*crnt->cl->nxt2[f];
  run->sol_R->v = u_sol_R_lc*crnt->cl->ny[f] + v_sol_R_lc*crnt->cl->nyt1[f] + w_sol_R_lc*crnt->cl->nyt2[f];
  run->sol_R->w = u_sol_R_lc*crnt->cl->nz[f] + v_sol_R_lc*crnt->cl->nzt1[f] + w_sol_R_lc*crnt->cl->nzt2[f];
  
	return;
}
