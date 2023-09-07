#include "strdata.h"

void m4_bc(struct BRANCH * crnt, int p, int f, struct RUN *run){
  
  switch(p){

    case 0: // no bc  
    break;
    
    case 1: // Symmetry  
      bc_symmetry(crnt,f,run);
    break;
    
    case 2: // outlet 
      m4_bc_outlet(crnt,f,run);
    break;

    case 3: // inlet 
      m4_bc_inlet(crnt,f,run);
    break;

    case 4: // Reflective (slip wall)
      bc_reflective(crnt,f,run);
    break;

    case 5: // No slip wall  
      bc_wall(crnt,f,run);
    break;

    default:
    break;
  }

}

void m4_bc_inlet(struct BRANCH * crnt,int f,struct RUN *run) {

  int iv;
  double rho,u,v,w,ru,rv,rw,re;
  struct SOL * r_sol_R;

  r_sol_R = run->sol_R;

  if ((CASE==3400)|| (CASE==3410)) {
    
    double press_inlet = 150.0; //MATERPINI[0]; 

    /*
    // ar0/ar1, rho
    r_sol_R->ra[0]  = MATERRINI[0]*(1.0-AMIN);
    r_sol_R->vec[0] = MATERRINI[0]*(1.0-AMIN);
    
    r_sol_R->ra[1]  = MATERRINI[1]*(AMIN);
    r_sol_R->vec[1] = MATERRINI[1]*(AMIN);
    
    r_sol_R->r = r_sol_R->ra[0] + r_sol_R->ra[1];
      
    // u,v,w
    r_sol_R->u = 0.0;
    r_sol_R->v = 0.0;
    r_sol_R->w = 0.0;

    iv = eqtypi[1];
    r_sol_R->vec[iv+0] = r_sol_R->r*r_sol_R->u;
    r_sol_R->vec[iv+1] = r_sol_R->r*r_sol_R->v;
    r_sol_R->vec[iv+2] = r_sol_R->r*r_sol_R->w;

    // avf
    r_sol_R->avf[0]  = (1.0-AMIN);
    r_sol_R->avf[1]  = 			AMIN;

    iv=eqtypi[3];  
    r_sol_R->vec[iv] = r_sol_R->avf[0];
    */

    // are
    r_sol_R->p_hydro[0] = press_inlet*pow(10.0,6.0); // 100 MPa
    r_sol_R->p_hydro[1] = press_inlet*pow(10.0,6.0); // 100 MPa

    r_sol_R->st[0][0] = - press_inlet*pow(10.0,6.0); // 100 MPa
    r_sol_R->st[0][1] = 0.0;
    r_sol_R->st[0][2] = 0.0;

    r_sol_R->st[1][0] = 0.0;
    r_sol_R->st[1][1] = - press_inlet*pow(10.0,6.0); // 100 MPa
    r_sol_R->st[1][2] = 0.0;

    r_sol_R->st[2][0] = 0.0;
    r_sol_R->st[2][1] = 0.0;
    r_sol_R->st[2][2] = - press_inlet*pow(10.0,6.0); // 100 MPa

    r_sol_R->are[0]  = r_sol_R->avf[0] * (r_sol_R->p_hydro[0] + MATERGAMA[0]*MATERPINF[0])/(MATERGAMA[0]-1.0);   
    r_sol_R->are[1]  = r_sol_R->avf[1] * (r_sol_R->p_hydro[1] + MATERGAMA[1]*MATERPINF[1])/(MATERGAMA[1]-1.0);  

    r_sol_R->vec[eqtypi[4]+0] = r_sol_R->are[0];
    r_sol_R->vec[eqtypi[4]+1] = r_sol_R->are[1];

    r_sol_R->c2Y[0]=(r_sol_R->avf[0]/r_sol_R->r)*MATERGAMA[0]*(r_sol_R->p_hydro[0]+MATERPINF[0]);
    r_sol_R->c2Y[1]=(r_sol_R->avf[1]/r_sol_R->r)*MATERGAMA[1]*(r_sol_R->p_hydro[1]+MATERPINF[1]);

    // r£
    re = r_sol_R->are[0] + r_sol_R->are[1] + 0.5*r_sol_R->r*(pow(r_sol_R->u,2.0) + 
                                                             pow(r_sol_R->v,2.0) +
                                                             pow(r_sol_R->w,2.0));  // rE

    r_sol_R->e = re / r_sol_R->r;
    r_sol_R->vec[eqtypi[2]+0] = re;

  } 

}

void m4_bc_outlet(struct BRANCH * crnt,int f,struct RUN *run) {

  int iv;
  double rho,u,v,w,ru,rv,rw,re;
  struct SOL * r_sol_R;

  r_sol_R = run->sol_R;

  if ((CASE==3400)|| (CASE==3410)) {

    double press_outlet = 1.0; //  MPa; 
    /*
    // ar0/ar1, rho
    r_sol_R->ra[0]  = MATERRINI[0]*(1.0-AMIN);
    r_sol_R->vec[0] = MATERRINI[0]*(1.0-AMIN);
    
    r_sol_R->ra[1]  = MATERRINI[1]*(AMIN);
    r_sol_R->vec[1] = MATERRINI[1]*(AMIN);
    
    r_sol_R->r = r_sol_R->ra[0] + r_sol_R->ra[1];
      
    // u,v,w
    r_sol_R->u = 0.0;
    r_sol_R->v = 0.0;
    r_sol_R->w = 0.0;

    iv = eqtypi[1];
    r_sol_R->vec[iv+0] = r_sol_R->r*r_sol_R->u;
    r_sol_R->vec[iv+1] = r_sol_R->r*r_sol_R->v;
    r_sol_R->vec[iv+2] = r_sol_R->r*r_sol_R->w;

    // avf
    r_sol_R->avf[0]  = (1.0-AMIN);
    r_sol_R->avf[1]  = 			AMIN;

    iv=eqtypi[3];  
    r_sol_R->vec[iv] = r_sol_R->avf[0];
    */
    // are
    r_sol_R->p_hydro[0] = press_outlet*pow(10.0,6.0); // 100 MPa
    r_sol_R->p_hydro[1] = press_outlet*pow(10.0,6.0); // 100 MPa

    r_sol_R->st[0][0] = - press_outlet*pow(10.0,6.0); // 100 MPa
    r_sol_R->st[0][1] = 0.0;
    r_sol_R->st[0][2] = 0.0;

    r_sol_R->st[1][0] = 0.0;
    r_sol_R->st[1][1] = - press_outlet*pow(10.0,6.0); // 100 MPa
    r_sol_R->st[1][2] = 0.0;

    r_sol_R->st[2][0] = 0.0;
    r_sol_R->st[2][1] = 0.0;
    r_sol_R->st[2][2] = - press_outlet*pow(10.0,6.0); // 100 MPa

    r_sol_R->are[0]  = r_sol_R->avf[0] * (r_sol_R->p_hydro[0] + MATERGAMA[0]*MATERPINF[0])/(MATERGAMA[0]-1.0);   
    r_sol_R->are[1]  = r_sol_R->avf[1] * (r_sol_R->p_hydro[1] + MATERGAMA[1]*MATERPINF[1])/(MATERGAMA[1]-1.0);  

    r_sol_R->vec[eqtypi[4]+0] = r_sol_R->are[0];
    r_sol_R->vec[eqtypi[4]+1] = r_sol_R->are[1];

    r_sol_R->c2Y[0]=(r_sol_R->avf[0]/r_sol_R->r)*MATERGAMA[0]*(r_sol_R->p_hydro[0]+MATERPINF[0]);
    r_sol_R->c2Y[1]=(r_sol_R->avf[1]/r_sol_R->r)*MATERGAMA[1]*(r_sol_R->p_hydro[1]+MATERPINF[1]);

    // r£
    re = r_sol_R->are[0] + r_sol_R->are[1] + 0.5*r_sol_R->r*(pow(r_sol_R->u,2.0) + 
                                                             pow(r_sol_R->v,2.0) +
                                                             pow(r_sol_R->w,2.0));  // rE

    r_sol_R->e = re / r_sol_R->r;
    r_sol_R->vec[eqtypi[2]+0] = re;
  }

}