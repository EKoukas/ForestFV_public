#include "strdata.h"

void m0_facercnstr(struct RUN *run,struct BRANCH * crnt,int in,int ifc,int bc_temp,int imode) {

  int iv,ifc_opp,ifc_temp;
  struct SOL * r_sol;
  struct SOL * r_sol1;
  struct LEAF * temp_el;

  //==============================================================================================================================
  if       (imode==0)               { temp_el = crnt->el;                  r_sol = run->sol_L; ifc_temp = ifc;               } 
  else if ((imode==1)&&(bc_temp==0)){ temp_el = crnt->neigtr[ifc][in]->el; r_sol = run->sol_R; ifc_temp = crnt->neigfc[ifc]; }  
  else if ((imode==1)&&(bc_temp!=0)){ temp_el = crnt->el;                  r_sol = run->sol_R; ifc_temp = ifc;               }
  //==============================================================================================================================

  r_sol->r = temp_el->SF[ifc_temp][0];
  r_sol->u = temp_el->SF[ifc_temp][1] / r_sol->r;
  r_sol->v = temp_el->SF[ifc_temp][2] / r_sol->r;
  r_sol->w = temp_el->SF[ifc_temp][3] / r_sol->r;
  
  r_sol->vec[0] = r_sol->r;
  r_sol->vec[1] = temp_el->SF[ifc_temp][1];
  r_sol->vec[2] = temp_el->SF[ifc_temp][2];
  r_sol->vec[3] = temp_el->SF[ifc_temp][3];
  
  r_sol->e          = temp_el->SF[ifc_temp][4];
  r_sol->vec[4]     = temp_el->SF[ifc_temp][4];
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
    m0_bc(crnt,bc_temp,ifc,run);
  }

}