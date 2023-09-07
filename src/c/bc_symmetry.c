#include "strdata.h"

void bc_symmetry(struct BRANCH * crnt,int f,struct RUN *run) {

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

  /*
  // ==========================================================================================
  rotate_tensor(tensor_temp,run->sol_L->Amat,crnt->cl->nx[f]  ,crnt->cl->ny[f]  ,crnt->cl->nz[f],
                                            crnt->cl->nxt1[f],crnt->cl->nyt1[f],crnt->cl->nzt1[f],
                                            crnt->cl->nxt2[f],crnt->cl->nyt2[f],crnt->cl->nzt2[f],
                                            VERBOSE,run->con->rank);

   v=0;
  for(i=0;i<3;++i){ 
    for(j=0;j<3;++j){ 
      run->Amat[i*3+j] = tensor_temp[v];
      ++v;
    }
  }
  
  for (i=0;i<3;++i){ 
    for (j=0;j<3;++j){ 
      if ((j==0) && (i!=j)){  
        run->Amat[i*3+j] = -run->Amat[i*3+j]; // A01 A02 local     
      } else {
        run->Amat[i*3+j] = run->Amat[i*3+j];
      }
      
    } 
  }
  
  rotate_tensor(tensor_temp,run->Amat,crnt->cl->nx[f], crnt->cl->nxt1[f], crnt->cl->nxt2[f],
          		                        crnt->cl->ny[f], crnt->cl->nyt1[f], crnt->cl->nyt2[f],
          				                    crnt->cl->nz[f], crnt->cl->nzt1[f], crnt->cl->nzt2[f],
          				                    VERBOSE,run->con->rank);
  
  v=0;
  for(i=0;i<3;++i){  
    for(j=0;j<3;++j){ 
      run->sol_R->Amat[i*3+j] = tensor_temp[v];
      ++v;
    }
  }
  // ==========================================================================================
  
  // ==========================================================================================
  rotate_tensor(tensor_temp,run->sol_L->st,crnt->cl->nx[f]  ,crnt->cl->ny[f]  ,crnt->cl->nz[f],
                                         crnt->cl->nxt1[f],crnt->cl->nyt1[f],crnt->cl->nzt1[f],
                                         crnt->cl->nxt2[f],crnt->cl->nyt2[f],crnt->cl->nzt2[f],
                                         VERBOSE,run->con->rank);

   v=0;
  for(i=0;i<3;++i){ 
    for(j=0;j<3;++j){ 
      run->Amat[i*3+j] = tensor_temp[v]; //st
      ++v;
    }
  }
  
  for (i=0;i<3;++i){ 
    for (j=0;j<3;++j){ 
      if ((j==0) && (i!=j)){  
        run->Amat[i*3+j] = -run->Amat[i*3+j]; // A01 A02 local     
      } else {
        run->Amat[i*3+j] = run->Amat[i*3+j];
      }
      
    } 
  }
  
  rotate_tensor(tensor_temp,run->Amat,crnt->cl->nx[f], crnt->cl->nxt1[f], crnt->cl->nxt2[f],
          		                        crnt->cl->ny[f], crnt->cl->nyt1[f], crnt->cl->nyt2[f],
          				                    crnt->cl->nz[f], crnt->cl->nzt1[f], crnt->cl->nzt2[f],
          				                    VERBOSE,run->con->rank);
  
  v=0;
  for(i=0;i<3;++i){  
    for(j=0;j<3;++j){ 
      run->sol_R->st[i][j] = tensor_temp[v];
      ++v;
    }
  }
  // ==========================================================================================
  */

  free(tensor_temp);
	return;
}
