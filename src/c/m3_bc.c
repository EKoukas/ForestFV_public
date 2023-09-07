#include "strdata.h"

void m3_bc(struct BRANCH * crnt, int p, int f, struct RUN *run){
  
  switch(p){

    case 0: // no bc  
    break;
    
    case 1: // Symmetry  
      bc_symmetry(crnt,f,run);
    break;
    
    case 2: // outlet 
    break;

    case 3: // inlet 
      m3_bc_inlet(crnt,f,run);
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

void m3_bc_inlet(struct BRANCH * crnt,int f,struct RUN *run) {

  double rho,u,v,w,ru,rv,rw;

	if (CASE==200) { 
  

  } 
  else if (CASE==203) { 
  

  }

}