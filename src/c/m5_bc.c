#include "strdata.h"

void m5_bc(struct BRANCH * crnt, int p, int f, struct RUN *run){
  
  switch(p){

    case 0: // no bc  
    break;
    
    case 1: // Symmetry  
      m5_bc_symmetry(crnt,f,run);
    break;
    
    case 2: // outlet 
    break;

    case 3: // inlet 
      m5_bc_inlet(crnt,f,run);
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

void m5_bc_inlet(struct BRANCH * crnt,int f,struct RUN *run) {

  double rho,u,v,w,ru,rv,rw;

	if (CASE==200) { 
  

  } 
  else if (CASE==203) { 
  

  }

}


void m5_bc_symmetry(struct BRANCH * crnt,int f,struct RUN *run) {

 int i,j,v;
  double u_sol_lc ,v_sol_lc ,w_sol_lc;
  double u_sol1_lc,v_sol1_lc,w_sol1_lc;
  
  double * tensor_temp; 
  tensor_temp = malloc(9*sizeof(double));

  // Rotate global velocities to local coordinate system
  rotate_vector_gtl_v2(run,crnt,&f,&(run->sol_L->u),&(run->sol_L->v),&(run->sol_L->w));

  u_sol_lc = run->vel_temp[0];
  v_sol_lc = run->vel_temp[1]; 
  w_sol_lc = run->vel_temp[2];

  // Opposite normal velocity 
  u_sol1_lc = -u_sol_lc;
  v_sol1_lc =  v_sol_lc;
  w_sol1_lc =  w_sol_lc;

  // Rotate back to global
  rotate_vector_ltg_v2(run,crnt,&f,&(u_sol1_lc),&(v_sol1_lc),&(w_sol1_lc));

  run->sol_R->u = run->vel_temp[0];
	run->sol_R->v = run->vel_temp[1];
	run->sol_R->w = run->vel_temp[2];

  
  
  free(tensor_temp);
	return;

}