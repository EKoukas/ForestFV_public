#include "strdata.h"

void m1_facevalues(struct RUN *run,struct BRANCH * crnt,int in,int ifc,int bc_temp,int imode) {

  int iv,i;
  struct SOL * r_sol;
  struct LEAF * temp_el;
  
  //==============================================================================================
  if       (imode==0)               { temp_el = crnt->el;                  r_sol = run->sol_L; } 
  else if ((imode==1)&&(bc_temp==0)){ temp_el = crnt->neigtr[ifc][in]->el; r_sol = run->sol_R; }  
  else if ((imode==1)&&(bc_temp!=0)){ temp_el = crnt->el;                  r_sol = run->sol_R; }
  //==============================================================================================

  r_sol->r = temp_el->S[0];
  r_sol->u = temp_el->S[1] / r_sol->r;
  r_sol->v = temp_el->S[2] / r_sol->r;
  r_sol->w = temp_el->S[3] / r_sol->r;
  
  r_sol->vec[0] = r_sol->r;
  r_sol->vec[1] = temp_el->S[1];
  r_sol->vec[2] = temp_el->S[2];
  r_sol->vec[3] = temp_el->S[3];
  
  m1_barotropic(r_sol); // p_hydro & c

  r_sol->st[0][0] = r_sol->p_hydro[0]; 
  r_sol->st[0][1] = 0.0;
  r_sol->st[0][2] = 0.0;
           
  r_sol->st[1][0] = 0.0;
  r_sol->st[1][1] = r_sol->p_hydro[0]; 
  r_sol->st[1][2] = 0.0;
  
  r_sol->st[2][0] = 0.0;
  r_sol->st[2][1] = 0.0;
  r_sol->st[2][2] = r_sol->p_hydro[0];

  if (NS==1) {  // Navier-stokes
    r_sol->ux = temp_el->WG[1][0];
    r_sol->uy = temp_el->WG[1][1];
    r_sol->uz = temp_el->WG[1][2];

    r_sol->vx = temp_el->WG[2][0];
    r_sol->vy = temp_el->WG[2][1];
    r_sol->vz = temp_el->WG[2][2];
    
    r_sol->wx = temp_el->WG[3][0];
    r_sol->wy = temp_el->WG[3][1];
    r_sol->wz = temp_el->WG[3][2];
  }
  
  if ((imode==1)&&(bc_temp!=0)){
    m1_bc(crnt,bc_temp,ifc,run);
  }
  
} 