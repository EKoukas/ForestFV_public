#include "strdata.h"

void m1_bc(struct BRANCH * crnt, int p, int f, struct RUN *run){
  
  switch(p){

    case 0: // no bc  
    break;
    
    case 1: // Symmetry  
      bc_reflective(crnt,f,run);
    break;
    
    case 2: // outlet 
    break;

    case 3: // inlet 
      m1_bc_inlet(crnt,f,run);
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

void m1_bc_inlet(struct BRANCH * crnt,int f,struct RUN *run) {

  double rho,u,v,w,ru,rv,rw;

	if (CASE==200) { // Forward faceing step, https://amroc.sourceforge.net/examples/euler/2d/html/ffstep_n.htm
  

  } 
  else if (CASE==203) { // boundary layer
  

  }

}