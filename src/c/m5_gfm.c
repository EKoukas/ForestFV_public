#include "strdata.h"

void m5_gfm(struct RUN *run,struct BRANCH * crnt) {

  int f,iv,v,f_inter,ing,f_temp,ind,bc_temp;
  double rho_old,rho_new,phi,flag_sols;
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

  double rx,ry,rz;

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

  /*
  if (run->debug_gfm==1) {
    printf("\n");
    printf("GFM NOW %d %d \n",crnt->root,g_istep_tot);
  }
  */

  // ---------------------------------------------------------------------
  // Normal speed from step 2
    rho_old = 0.0;
    rho_new = 0.0;
    iv = eqtypi[0];
    for(v=0;v<eqtypn[0];++v) {
      rho_old += crnt->el->S[iv+v];
      rho_new += crnt->el->S[iv+v] + (DT/crnt->cl->Vol)*crnt->el->RHS[iv+v];
    }

    iv = eqtypi[1];
    u_temp_gb = crnt->el->S[iv+0] + (DT/crnt->cl->Vol)*crnt->el->RHS[iv+0];
    v_temp_gb = crnt->el->S[iv+1] + (DT/crnt->cl->Vol)*crnt->el->RHS[iv+1];
    w_temp_gb = crnt->el->S[iv+2] + (DT/crnt->cl->Vol)*crnt->el->RHS[iv+2];
    
    u_temp_gb = u_temp_gb/rho_new;
    v_temp_gb = v_temp_gb/rho_new;
    w_temp_gb = w_temp_gb/rho_new;
  // ---------------------------------------------------------------------

  nx_inter = 0.0; ny_inter = 0.0; nz_inter = 0.0;
  f_inter  = 0; total_area = 0.0;
  
  for(f=0; f<(crnt->nlfc); f++){
    
    ing=0;
    // ---------------------------------------------------------
    // sol, sol1 
    m5_facevalues(run,crnt,ing,f,0,0);          // crnt->el->S  to sol_L

    bc_temp = crnt->cl->fc[f].bc;
    m5_facevalues(run,crnt,ing,f,bc_temp,1);  // sol_R
    //---------------------------------------------------------

    if (crnt->el->gfm_face[f]==1.0) { // Interface surface

      // face value
      phi = (run->sol_L->phi[0] + run->sol_R->phi[0])/2.0; 
     
      nx_inter = phi*crnt->cl->Area[f][ing]*crnt->cl->nx[f];
      ny_inter = phi*crnt->cl->Area[f][ing]*crnt->cl->ny[f];
      nz_inter = phi*crnt->cl->Area[f][ing]*crnt->cl->nz[f];

      u_sol1_inter_gb[f_inter] = run->sol_R->u;
      v_sol1_inter_gb[f_inter] = run->sol_R->v;
      w_sol1_inter_gb[f_inter] = run->sol_R->w;

      total_area += crnt->cl->Area[f][ing];

      /*
      if (run->debug_gfm==1) {
        printf("Interface face %d | %f %f \n",f,run->sol_L->phi[0],run->sol_R->phi[0]);
        printf("Normal vector  %d | %f %f %f \n",f,crnt->cl->nx[f],crnt->cl->ny[f],crnt->cl->nz[f]);
        printf("phi %f \n",phi);
        printf("Interface sol1 vel gl %d | %f %f %f | %e %e \n",f_inter,u_sol1_inter_gb[f_inter], v_sol1_inter_gb[f_inter], w_sol1_inter_gb[f_inter],
        total_area,crnt->cl->Area[f][ing]);
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
    

    // ---------------------------------------------------------------------------------
    // Rotate updated velocities on interface surface
    rotate_vector(vel_temp2,u_temp_gb,v_temp_gb,w_temp_gb,nx_inter , ny_inter , nz_inter,
                                                          nx1_inter, ny1_inter, nz1_inter,
                                                          nx2_inter, ny2_inter, nz2_inter,
                                                          VERBOSE,run->con->rank);
    u_temp_lc = vel_temp2[0];
    v_temp_lc = vel_temp2[1];
    w_temp_lc = vel_temp2[2];
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

    }
    
    // ----------------------------------------------------------------
    u_new_lc = u_temp_lc;
    v_new_lc = 0.0;
    w_new_lc = 0.0;
    
    f_temp=0;
    for(f=0; f<(crnt->nlfc); f++){
      if (crnt->el->gfm_face[f]==1.0) { // Interface surface
        //u_new_lc = u_sol1_lc[f];
        v_new_lc = v_sol1_inter_lc[f_temp]*(crnt->cl->Area[f][ing]/total_area);
        w_new_lc = w_sol1_inter_lc[f_temp]*(crnt->cl->Area[f][ing]/total_area);

        f_temp++;
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
    iv = eqtypi[0];
    for(v=0;v<eqtypn[0];++v) {
      crnt->el->S[iv+v] = crnt->el->S[iv+v] + (DT/crnt->cl->Vol)*crnt->el->RHS[iv+v]; //avf_rho
    }

    iv = eqtypi[1];
    crnt->el->S[iv+0] = rho_new*u_new_gb;
    crnt->el->S[iv+1] = rho_new*v_new_gb;
    crnt->el->S[iv+2] = rho_new*w_new_gb;
    
    /*
    if (run->debug_gfm==1) {
      printf("s new %f | %f %f %f \n",rho_new,crnt->el->S[iv+0]/rho_new,crnt->el->S[iv+1]/rho_new,crnt->el->S[iv+2]/rho_new);
      printf("\n");
    }
    */

    iv = eqtypi[2];
    crnt->el->S[iv + 0] -= rho_new*0.5*(pow(u_temp_gb,2.0) + pow(v_temp_gb,2.0) + pow(w_temp_gb,2.0));
    crnt->el->S[iv + 0] += rho_new*0.5*(pow(u_new_gb,2.0)  + pow(v_new_gb,2.0)  + pow(w_new_gb,2.0));   

  } else {
    
    iv = eqtypi[0];
    for(v=0;v<eqtypn[0];++v) {
      crnt->el->S[iv+v] = crnt->el->S[iv+v] + (DT/crnt->cl->Vol)*crnt->el->RHS[iv+v];
    }

    iv = eqtypi[1];
    crnt->el->S[iv+0] = crnt->el->S[iv+0] + (DT/crnt->cl->Vol)*crnt->el->RHS[iv+0];
    crnt->el->S[iv+1] = crnt->el->S[iv+1] + (DT/crnt->cl->Vol)*crnt->el->RHS[iv+1];
    crnt->el->S[iv+2] = crnt->el->S[iv+2] + (DT/crnt->cl->Vol)*crnt->el->RHS[iv+2];
    
    rx=0.0;ry=0.0;rz=0.0;
    for (ind=0;ind<crnt->nlnd;ind++){
      rx+=crnt->cl->nd[ind].x/((double)crnt->nlnd);
      ry+=crnt->cl->nd[ind].y/((double)crnt->nlnd);
      rz+=crnt->cl->nd[ind].z/((double)crnt->nlnd);
    }

    printf("GFM with no interfaces \n");
    printf("coordinates: %e %e %e \n",rx,ry,rz);
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
