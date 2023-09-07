#include "strdata.h"

void m3_facercnstr(struct RUN *run,struct BRANCH * crnt,int in,int ifc,int bc_temp,int imode) {

  int iv,v,i,j,ifc_temp;
  int ifc_opp;
  double avfsum,avf_per,avf_dif,press;
  double avf_temp[10];
  double avf_c_temp[10];
  struct SOL * r_sol;
  struct LEAF * temp_el;

  //==============================================================================================================================
  if       (imode==0)               { temp_el = crnt->el;                  r_sol = run->sol_L; ifc_temp = ifc;               } 
  else if ((imode==1)&&(bc_temp==0)){ temp_el = crnt->neigtr[ifc][in]->el; r_sol = run->sol_R; ifc_temp = crnt->neigfc[ifc]; }  
  else if ((imode==1)&&(bc_temp!=0)){ temp_el = crnt->el;                  r_sol = run->sol_R; ifc_temp = ifc;               }
  //==============================================================================================================================

  // density
  r_sol->r=0.0;
	for(v=0;v<eqtypn[0];++v){
	  r_sol->r     += temp_el->SF[ifc_temp][v];
	  r_sol->ra[v]  = temp_el->SF[ifc_temp][v];
	  r_sol->vec[v] = temp_el->SF[ifc_temp][v];
	}

  // mass fraction
	for(v=0;v<eqtypn[0];++v){
	  r_sol->Y[v] = r_sol->ra[v] / r_sol->r;
	}

  // velocities + total energy
  v = eqtypi[1];
	r_sol->u = temp_el->SF[ifc_temp][v+0] / r_sol->r;
	r_sol->v = temp_el->SF[ifc_temp][v+1] / r_sol->r;
	r_sol->w = temp_el->SF[ifc_temp][v+2] / r_sol->r;

  r_sol->vec[v+0] = temp_el->SF[ifc_temp][v+0];
  r_sol->vec[v+1] = temp_el->SF[ifc_temp][v+1];
  r_sol->vec[v+2] = temp_el->SF[ifc_temp][v+2];

  r_sol->e = temp_el->SF[ifc_temp][eqtypi[2]] / r_sol->r;
  r_sol->vec[eqtypi[2]] = temp_el->SF[ifc_temp][eqtypi[2]];
  
  // AMIN check
  for(iv=0;iv<eqtypn[3];++iv){ 
    avf_temp[iv]   = temp_el->SF[ifc_temp][iv+eqtypi[3]]; 
    avf_c_temp[iv] = temp_el->S[           iv+eqtypi[3]]; 
  }

  avfsum=0.0; 
  for(v=0;v<eqtypn[3];++v){
    if (avf_temp[v]  <=      AMIN/100.0)   { avf_temp[v]   =      AMIN/100.0;  } 
    if (avf_c_temp[v]<=      AMIN/100.0)   { avf_c_temp[v] =      AMIN/100.0;  } 
    if (avf_temp[v]  >=(1.0-(AMIN/100.0))) { avf_temp[v]   = 1.0-(AMIN/100.0); }   
    if (avf_c_temp[v]>=(1.0-(AMIN/100.0))) { avf_c_temp[v] = 1.0-(AMIN/100.0); }
    avfsum   += avf_temp[v]; 
  }  
  
  if (avfsum>(1.0-(AMIN/100.0))) { 
    avf_dif = avfsum - (1.0-(AMIN/100.0));  
    for(j=0;j<eqtypn[3];++j){
      avf_per = avf_temp[j]/avfsum;
      avf_temp[j]  = avf_temp[j] - avf_per*avf_dif;
    }
    avfsum = 1.0-(AMIN/100.0); 
  }
  avf_temp[eqtypn[3]] = 1.0 - avfsum;
  
  avfsum=0.0;
  iv=eqtypi[3];
	for(v=0;v<eqtypn[3];++v){
    r_sol->avf[v]    = avf_temp[v];
    avfsum          += avf_temp[v];
    r_sol->vec[iv+v] = avf_temp[v];
	}
  r_sol->avf[eqtypn[3]] = 1.0 - avfsum;
    
  for(v=0;v<eqtypn[3];++v){
    r_sol->phi[v] = avf_c_temp[v]-AVF_LIM;
  } 

		
  r_sol->st[0][0] = 0.0;
  r_sol->st[0][1] = 0.0;
  r_sol->st[0][2] = 0.0;

  r_sol->st[1][0] = 0.0;
  r_sol->st[1][1] = 0.0;
  r_sol->st[1][2] = 0.0;

  r_sol->st[2][0] = 0.0;
  r_sol->st[2][1] = 0.0;
  r_sol->st[2][2] = 0.0;

  for (v=0;v<eqtypn[5];++v) {   
    
    r_sol->c2Y[v] = (r_sol->avf[v]/r_sol->r)*MATERGAMA[v]*(press+MATERPINF[v]);

    r_sol->st[0][0] += -r_sol->avf[v]*press;
    r_sol->st[1][1] += -r_sol->avf[v]*press;
    r_sol->st[2][2] += -r_sol->avf[v]*press;
    
    r_sol->p_hydro[v] = press;
  } 

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
    m3_bc(crnt,bc_temp,ifc,run);
  }
  
}