#include "strdata.h"

void m0_bc(struct BRANCH * crnt, int p, int f, struct RUN *run){
  
  switch(p){

    case 0: // no bc  
    break;
    
    case 1: // Symmetry  
      bc_reflective(crnt,f,run);
    break;
    
    case 2: // outlet 
    break;

    case 3: // inlet 
      m0_bc_inlet(crnt,f,run);
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

void m0_bc_inlet(struct BRANCH * crnt,int f,struct RUN *run) {

  double rho,u,v,w,ru,rv,rw,p,e,e_internal;

	if (CASE==200) { // Forward faceing step, https://amroc.sourceforge.net/examples/euler/2d/html/ffstep_n.htm
  
		rho = 1.4;
		u   = 3.0;
		v   = 0.0;
		w   = 0.0;
		ru  = rho*u;
		rv  = rho*v;
		rw  = rho*w;
		p   = 1.0;

    e_internal = (p + MATERGAMA[0]*MATERPINF[0])/(rho*(MATERGAMA[0]-1.0)); 
		e = rho*( e_internal + 0.5*(pow(u,2.0) + pow(v,2.0) + pow(w,2.0)));
		
		m0_asinvalues(run,crnt,run->sol_R,rho,ru,rv,rw,e);

  } 
  else if (CASE==203) { // boundary layer
  
		rho = 1.0;
		u   = 40.0;
		v   = 0.0;
		w   = 0.0;
		ru  = rho*u;
		rv  = rho*v;
		rw  = rho*w;
		p   = 101325.0;

    e_internal = (p + MATERGAMA[0]*MATERPINF[0])/(rho*(MATERGAMA[0]-1.0)); 
		e = rho*( e_internal + 0.5*(pow(u,2.0) + pow(v,2.0) + pow(w,2.0)));
		
		m0_asinvalues(run,crnt,run->sol_R,rho,ru,rv,rw,e);

  }

}