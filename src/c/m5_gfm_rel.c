#include "strdata.h"

void m5_gfm_rel(struct BRANCH * crnt,struct RUN *run) {

  int f,iv,v,f_inter,ineig,f_temp,bc_temp;
  double rho_old,rho_new,phi,flag_sols,rho;
  double u_flow;
  double e_internal,Y_temp;

  double u_temp_gb,v_temp_gb,w_temp_gb;
  double u_temp_lc,v_temp_lc,w_temp_lc;

  double u_new_gb,v_new_gb,w_new_gb;
  double u_new_lc,v_new_lc,w_new_lc;

  double u_sol1_inter,v_sol1_inter,w_sol1_inter;
  
  double norm,total_area;
  double nx_inter ,ny_inter ,nz_inter;
  double nx1_inter,ny1_inter,nz1_inter;
  double nx2_inter,ny2_inter,nz2_inter;

  double u_old_gb,v_old_gb,w_old_gb;

  double * u_sol1_inter_gb;
  double * v_sol1_inter_gb;
  double * w_sol1_inter_gb;

  double * u_sol1_inter_lc;
  double * v_sol1_inter_lc;
  double * w_sol1_inter_lc;

  double * vel_temp2;
  
  vel_temp2 = malloc(3*sizeof(double));
  
  u_sol1_inter_lc = malloc(crnt->nlfc*sizeof(double));
  v_sol1_inter_lc = malloc(crnt->nlfc*sizeof(double));
  w_sol1_inter_lc = malloc(crnt->nlfc*sizeof(double));

  u_sol1_inter_gb = malloc(crnt->nlfc*sizeof(double));
  v_sol1_inter_gb = malloc(crnt->nlfc*sizeof(double));
  w_sol1_inter_gb = malloc(crnt->nlfc*sizeof(double));
    
  // ---------------------------------------------------------------------
  // Normal speed from step 2
    rho = 0.0;
    for(v=0;v<eqtypn[0];++v) {
      rho += crnt->el->S[v];
    }

    iv = eqtypi[1];
    u_temp_gb = crnt->el->S[iv+0]/rho;
    v_temp_gb = crnt->el->S[iv+1]/rho;
    w_temp_gb = crnt->el->S[iv+2]/rho;
  // ---------------------------------------------------------------------

  nx_inter = 0.0; ny_inter = 0.0; nz_inter = 0.0;
  f_inter  = 0; total_area = 0.0;

  for(f=0; f<(crnt->nlfc); f++){
    
    ineig=0;
    // ---------------------------------------------------------
    // sol, sol1 
      bc_temp=crnt->cl->fc[f].bc;
      m5_facevalues(run,crnt,ineig,f,bc_temp,0);
      m5_facevalues(run,crnt,ineig,f,bc_temp,1);

      /*
      if(crnt->cl->fc[f].bc==0){
        m5_facevalues(run,crnt,ineig,f,0,0);
        m5_facevalues(run,crnt,ineig,f,1,1);
      }
      
      if(crnt->cl->fc[f].bc!=0){
        m5_facevalues(run,crnt,ineig,f,0,0);
        m5_facevalues(run,crnt,ineig,f,0,1);
        m5_bc(crnt,crnt->cl->fc[f].bc,f,run);
      }
      */
    //---------------------------------------------------------

    if (crnt->el->gfm_face[f]==1.0) { // Interface surface

      phi = (run->sol_L->phi[0] + run->sol_R->phi[0])/2.0; 
      
      nx_inter += phi*crnt->cl->Area[f][ineig]*crnt->cl->nx[f];
      ny_inter += phi*crnt->cl->Area[f][ineig]*crnt->cl->ny[f];
      nz_inter += phi*crnt->cl->Area[f][ineig]*crnt->cl->nz[f];

      u_sol1_inter_gb[f_inter] = run->sol_R->u;
      v_sol1_inter_gb[f_inter] = run->sol_R->v;
      w_sol1_inter_gb[f_inter] = run->sol_R->w;

      total_area += crnt->cl->Area[f][ineig];

      /*
      if (run->debug_gfm==1) {
        printf("Interface face %d | %f %f \n",f,run->sol_L->phi[0],run->sol_R->phi[0]);
        printf("Normal vector  %d | %f %f %f \n",f,crnt->cl->nx[f],crnt->cl->ny[f],crnt->cl->nz[f]);
        printf("phi %f \n",phi);
        printf("Interface sol1 vel gl %d | %f %f %f | %e %e \n",f_inter,u_sol1_inter_gb[f_inter], v_sol1_inter_gb[f_inter], w_sol1_inter_gb[f_inter],
        total_area,crnt->cl->Area[f][ineig]);
      }
      */

      f_inter++;
    }
  }

  if (f_inter!=0) {

    norm = sqrt(pow(nx_inter,2.0) + pow(ny_inter,2.0) + pow(nz_inter,2.0));
    
    nx_inter = -nx_inter/norm;
    ny_inter = -ny_inter/norm;
    nz_inter = -nz_inter/norm;

    if (fabs(nz_inter)<0.99) {   // (nx_inter,ny_inter,nz_inter) x (0,0,-1)
      nx1_inter =  ny_inter;
      ny1_inter = -nx_inter;
      nz1_inter =  0.0;
    } else {                    // (nx_inter,ny_inter,nz_inter) x (0,-1,0)
      nx1_inter =  nz_inter;
      ny1_inter =  0.0;
      nz1_inter = -nx_inter;
    }

    // (nx_inter,ny_inter,nz_inter)x(nx1_inter,ny1_inter,nz1_inter)
    nx2_inter = ny_inter*nz1_inter - nz_inter*ny1_inter;
    ny2_inter = nz_inter*nx1_inter - nx_inter*nz1_inter; 
    nz2_inter = nx_inter*ny1_inter - ny_inter*nx1_inter;  
    
    /*
    if (run->debug_gfm==1) {
      printf("Normal vector %f %f %f | %f \n",nx_inter,ny_inter,nz_inter,norm);
    }
    */

    // ---------------------------------------------------------------------------------
    // Rotate updated velocities on interface surface
    rotate_vector(vel_temp2,u_temp_gb,v_temp_gb,w_temp_gb,nx_inter , ny_inter , nz_inter,
                                                          nx1_inter, ny1_inter, nz1_inter,
                                                          nx2_inter, ny2_inter, nz2_inter,
                                                          VERBOSE,run->con->rank);
    u_temp_lc = vel_temp2[0];
    v_temp_lc = vel_temp2[1];
    w_temp_lc = vel_temp2[2];
    
    /*
    if (run->debug_gfm==1) {
      printf("vel temp gb %f %f %f \n",u_temp_gb,v_temp_gb,w_temp_gb);
      printf("vel temp lc %f %f %f \n",u_temp_lc,v_temp_lc,w_temp_lc);
    }
    */
    // ---------------------------------------------------------------------------------

    // Rotate velocities on interface surface
    for(f=0; f<f_inter; f++){ 
      rotate_vector(vel_temp2,u_sol1_inter_gb[f],v_sol1_inter_gb[f],w_sol1_inter_gb[f],
                                                      nx_inter , ny_inter , nz_inter,
                                                      nx1_inter, ny1_inter, nz1_inter,
                                                      nx2_inter, ny2_inter, nz2_inter,
                                                      VERBOSE,run->con->rank);
      u_sol1_inter_lc[f] = vel_temp2[0];
      v_sol1_inter_lc[f] = vel_temp2[1];
      w_sol1_inter_lc[f] = vel_temp2[2];

      /*
      if (run->debug_gfm==1) {
        printf("Interface f_inter lc %d | %f %f %f | \n",f,u_sol1_inter_lc[f], v_sol1_inter_lc[f], w_sol1_inter_lc[f]);
      }
      */

    }
    
    // ----------------------------------------------------------------
    u_new_lc = u_temp_lc;
    v_new_lc = 0.0;
    w_new_lc = 0.0;
    
    f_temp=0;
    for(f=0; f<(crnt->nlfc); f++){
      if (crnt->el->gfm_face[f]==1.0) { // Interface surface
        //u_new_lc = u_sol1_lc[f];
        v_new_lc = v_sol1_inter_lc[f_temp]*(crnt->cl->Area[f][ineig]/total_area);
        w_new_lc = w_sol1_inter_lc[f_temp]*(crnt->cl->Area[f][ineig]/total_area);

        f_temp++;

        /*
        if (run->debug_gfm==1) {
          printf("vel new summation: %d %d %d| %f %f %f | \n",f,f_inter,f_temp,u_new_lc,v_new_lc,w_new_lc);
        }
        */
      }

    }
    
    rotate_vector(vel_temp2,u_new_lc,v_new_lc,w_new_lc,nx_inter,nx1_inter,nx2_inter,
                                                       ny_inter,ny1_inter,ny2_inter,
                                                       nz_inter,nz1_inter,nz2_inter,
                                                       VERBOSE,run->con->rank);                                             
    u_new_gb = vel_temp2[0];
    v_new_gb = vel_temp2[1];
    w_new_gb = vel_temp2[2];  

    /*
    if (run->debug_gfm==1) {
      printf("vel new gb %f %f %f \n",u_new_gb,v_new_gb,w_new_gb);
    } 
    */

    // ----------------------------------------------------------------
    iv = eqtypi[1];
    crnt->el->S[iv+0] = rho*u_new_gb;
    crnt->el->S[iv+1] = rho*v_new_gb;
    crnt->el->S[iv+2] = rho*w_new_gb;
    
    /*
    if (run->debug_gfm==1) {
      printf("s new %f | %f %f %f \n",rho,crnt->el->S[iv+0]/rho,crnt->el->S[iv+1]/rho,crnt->el->S[iv+2]/rho);
      printf("\n");
    }
    */

    iv = eqtypi[2];
    crnt->el->S[iv + 0] -= rho*0.5*(pow(u_temp_gb,2.0) + pow(v_temp_gb,2.0) + pow(w_temp_gb,2.0));
    crnt->el->S[iv + 0] += rho*0.5*(pow(u_new_gb, 2.0) + pow(v_new_gb, 2.0) + pow(w_new_gb, 2.0));   

  } else {
    printf("GFM_rel with no interfaces \n");
    printf("\n");
  }
  
  free(vel_temp2);

  free(u_sol1_inter_gb);
  free(v_sol1_inter_gb);
  free(w_sol1_inter_gb);

  free(u_sol1_inter_lc);
  free(v_sol1_inter_lc);
  free(w_sol1_inter_lc);

}