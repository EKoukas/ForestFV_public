#include "strdata.h"

void bc_reflective(struct BRANCH * crnt,int f,struct RUN *run) {

  double u,v,w,u_sol_L_lc ,v_sol_L_lc ,w_sol_L_lc;
  double u_sol_R_lc,v_sol_R_lc,w_sol_R_lc;
  struct CELL * crnt_cl;

  crnt_cl = crnt->cl;

  // Rotate global velocities to local coordinate system
  u = run->sol_L->u;
  v = run->sol_L->v;
  w = run->sol_L->w;

  u_sol_L_lc = u*crnt_cl->nx[f]   + v*crnt_cl->ny[f]   + w*crnt_cl->nz[f];
  v_sol_L_lc = u*crnt_cl->nxt1[f] + v*crnt_cl->nyt1[f] + w*crnt_cl->nzt1[f];
  w_sol_L_lc = u*crnt_cl->nxt2[f] + v*crnt_cl->nyt2[f] + w*crnt_cl->nzt2[f];

  // Opposite normal velocity 
  u_sol_R_lc = -u_sol_L_lc;
  v_sol_R_lc =  v_sol_L_lc;
  w_sol_R_lc =  w_sol_L_lc;

  // Rotate back to global
  run->sol_R->u = u_sol_R_lc*crnt_cl->nx[f] + v_sol_R_lc*crnt_cl->nxt1[f] + w_sol_R_lc*crnt_cl->nxt2[f];
  run->sol_R->v = u_sol_R_lc*crnt_cl->ny[f] + v_sol_R_lc*crnt_cl->nyt1[f] + w_sol_R_lc*crnt_cl->nyt2[f];
  run->sol_R->w = u_sol_R_lc*crnt_cl->nz[f] + v_sol_R_lc*crnt_cl->nzt1[f] + w_sol_R_lc*crnt_cl->nzt2[f];

}