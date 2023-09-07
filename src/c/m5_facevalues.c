#include "strdata.h"

void m5_facevalues(struct RUN *run,struct BRANCH * crnt,int in,int ifc,int bc_temp,int imode){

  int i,j,k,i_cons;
  double press,avfsum,avf_dif,avf_per,Gmag,are_elast;
  double Amat[9];
  double Gmat[9];
  double avf_temp[10];
  double avf_c_temp[10];
  struct SOL * r_sol;
  struct LEAF * temp_el;
  
  //==============================================================================================
  if       (imode==0)               { temp_el = crnt->el;                  r_sol = run->sol_L; } 
  else if ((imode==1)&&(bc_temp==0)){ temp_el = crnt->neigtr[ifc][in]->el; r_sol = run->sol_R; }  
  else if ((imode==1)&&(bc_temp!=0)){ temp_el = crnt->el;                  r_sol = run->sol_R; }
  //==============================================================================================

  // density
  r_sol->r=0.0;
	for(i=0;i<eqtypn[0];++i){
	  r_sol->r     += temp_el->S[i];
	  r_sol->ra[i]  = temp_el->S[i];
	  r_sol->vec[i] = temp_el->S[i];
	}

  // mass fraction
	for(i=0;i<eqtypn[0];++i){
	  r_sol->Y[i] = r_sol->ra[i] / r_sol->r;
	}
  
  // velocities + total energy
  i_cons = eqtypi[1];
	r_sol->u = temp_el->S[i_cons + 0] / r_sol->r;
	r_sol->v = temp_el->S[i_cons + 1] / r_sol->r;
	r_sol->w = temp_el->S[i_cons + 2] / r_sol->r;

  r_sol->vec[i_cons + 0] = temp_el->S[i_cons + 0];
  r_sol->vec[i_cons + 1] = temp_el->S[i_cons + 1];
  r_sol->vec[i_cons + 2] = temp_el->S[i_cons + 2];

  r_sol->e = temp_el->S[eqtypi[2]] / r_sol->r;
  r_sol->vec[eqtypi[2]] = temp_el->S[eqtypi[2]];

  // AMIN check
  i_cons = eqtypi[3];
  for(i=0;i<eqtypn[3];++i){ 
    avf_temp[i]   = temp_el->S[i_cons + i]; 
    avf_c_temp[i] = temp_el->S[i_cons + i]; 
  }

  avfsum=0.0; 
  for(i=0;i<eqtypn[3];++i){
    if      (avf_temp  [i]<=      AMIN/100.0)   { avf_temp  [i] =      AMIN/100.0;  } 
    else if (avf_temp  [i]>=(1.0-(AMIN/100.0))) { avf_temp  [i] = 1.0-(AMIN/100.0); }

    if      (avf_c_temp[i]<=      AMIN/100.0)   { avf_c_temp[i] =      AMIN/100.0;  }    
    else if (avf_c_temp[i]>=(1.0-(AMIN/100.0))) { avf_c_temp[i] = 1.0-(AMIN/100.0); }
    avfsum += avf_temp[i]; 
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
  i_cons=eqtypi[3];
	for(i=0;i<eqtypn[3];++i){
    r_sol->avf[i]          = avf_temp[i];
    r_sol->vec[i_cons + i] = avf_temp[i];
    avfsum                += avf_temp[i];
	}
  r_sol->avf[eqtypn[3]] = 1.0 - avfsum;
  
  // phi
  for(i=0;i<eqtypn[3];++i){
    r_sol->phi[i] = avf_c_temp[i]-AVF_LIM;
  } 

  if (n_solids==1) {
    r_sol->phi[0] = avf_c_temp[0]-AVF_LIM;
  } else if (n_solids==2) {
    r_sol->phi[0] = max(avf_c_temp[0],avf_c_temp[1])-AVF_LIM;
  } else if (n_solids==3) {
    r_sol->phi[0] = max(max(avf_c_temp[0],avf_c_temp[1]),avf_c_temp[2])-AVF_LIM;
  }

  for(i=0;i<eqtypn[3];++i){
    if (MATERMUSH[i]==0.0) {
      r_sol->phi[i] = avf_c_temp[i]-AVF_LIM;
    }
  }

  // are
  i_cons=eqtypi[4];
	for(i=0;i<eqtypn[4];++i){
	  r_sol->are[i]          = temp_el->S[i_cons + i];
    r_sol->vec[i_cons + i] = temp_el->S[i_cons + i];
	}

  // Amat
	k=0;
  i_cons = eqtypi[5];
  for (i=0;i<3;++i){
    for (j=0;j<3;++j){
      Amat[i*3+j]            = temp_el->S[i_cons + k];
      r_sol->Amat[i*3+j]     = temp_el->S[i_cons + k]; 
      r_sol->vec[i_cons + k] = temp_el->S[i_cons + k]; 
      k++;
    }
  }

  Gmat[0] = Amat[0]*Amat[0] + Amat[1]*Amat[1] + Amat[2]*Amat[2];
  Gmat[1] = Amat[0]*Amat[3] + Amat[1]*Amat[4] + Amat[2]*Amat[5];
  Gmat[2] = Amat[0]*Amat[6] + Amat[1]*Amat[7] + Amat[2]*Amat[8];

  Gmat[3] = Amat[3]*Amat[0] + Amat[4]*Amat[1] + Amat[5]*Amat[2];
  Gmat[4] = Amat[3]*Amat[3] + Amat[4]*Amat[4] + Amat[5]*Amat[5];
  Gmat[5] = Amat[3]*Amat[6] + Amat[4]*Amat[7] + Amat[5]*Amat[8];

  Gmat[6] = Amat[6]*Amat[0] + Amat[7]*Amat[1] + Amat[8]*Amat[2];
  Gmat[7] = Amat[6]*Amat[3] + Amat[7]*Amat[4] + Amat[8]*Amat[5];
  Gmat[8] = Amat[6]*Amat[6] + Amat[7]*Amat[7] + Amat[8]*Amat[8];

	Gmag = Gmat[0]*(Gmat[4]*Gmat[8] - Gmat[5]*Gmat[7]) - Gmat[1]*(Gmat[3]*Gmat[8] - Gmat[5]*Gmat[6]) + Gmat[2]*(Gmat[3]*Gmat[7] - Gmat[4]*Gmat[6]);
  
  r_sol->st[0][0] = 0.0;
  r_sol->st[0][1] = 0.0;
  r_sol->st[0][2] = 0.0;

  r_sol->st[1][0] = 0.0;
  r_sol->st[1][1] = 0.0;
  r_sol->st[1][2] = 0.0;

  r_sol->st[2][0] = 0.0;
  r_sol->st[2][1] = 0.0;
  r_sol->st[2][2] = 0.0;

  for (i=0;i<eqtypn[4];++i) {   
    
    if (MATERMUSH[i]!=0.0) { are_elast=m5_elastic_energy_solid_V2(Gmat,Gmag,MATERMUSH[i]); } 
    else {                   are_elast=0.0;                                                }

    press = ((MATERGAMA[i]-1.0)*(r_sol->are[i] - are_elast)/avf_temp[i]) - MATERGAMA[i]*MATERPINF[i];  // Stiffened-gas / ideal gas

    //r_sol->c2Y[i]     = (avf_temp[i]/r_sol->r)*(MATERGAMA[i]*(press+MATERPINF[i]) + (4.0/3.0)*MATERMUSH[i]);
    r_sol->c2Y[i]     = (avf_temp[i]/r_sol->r)*MATERGAMA[i]*(press+MATERPINF[i]);
    r_sol->p_hydro[i] = press;
    
    m5_eos_stress_V2(r_sol->avf_stress[i],r_sol->avf[i],r_sol->p_hydro[i],Gmat,Gmag,MATERMUSH[i]);
    
    
    r_sol->st[0][0] += r_sol->avf_stress[i][0][0];
    r_sol->st[0][1] += r_sol->avf_stress[i][0][1];
    r_sol->st[0][2] += r_sol->avf_stress[i][0][2];

    r_sol->st[1][0] += r_sol->avf_stress[i][1][0];
    r_sol->st[1][1] += r_sol->avf_stress[i][1][1];
    r_sol->st[1][2] += r_sol->avf_stress[i][1][2];

    r_sol->st[2][0] += r_sol->avf_stress[i][2][0];
    r_sol->st[2][1] += r_sol->avf_stress[i][2][1];
    r_sol->st[2][2] += r_sol->avf_stress[i][2][2];
    
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
    m5_bc(crnt,bc_temp,ifc,run);
  }

}