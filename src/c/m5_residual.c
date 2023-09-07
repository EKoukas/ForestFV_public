//-------------------------------------------------------
//      Computes Right Hand Side of the equatioNS
//-------------------------------------------------------

#include "strdata.h"

void m5_residual(struct RUN * run, struct BRANCH * crnt) {

  int iv,v,i,j,f,bc_temp;
  int ing,f_oppo,ifc_opp,ineig,Nneig;
  int flag_con;
  double A1,A2,A3,B1,B2,B3,C1,C2,C3;
  double u_nc,v_nc,w_nc;
  double nx,ny,nz,Area_f,as;
  double ** Area;
  int flag_compute;
  struct BRANCH * oppo;
  
  Area=malloc(crnt->nlfc*sizeof(double *));
  for(f=0; f<(crnt->nlfc); f++) {
    Area[f]=malloc(4*sizeof(double));
  }
  
  for(iv=0;iv<NEQ;iv++){ 
    crnt->el->RHS[iv] = 0.0;
    g_flux_viscous[iv][0] = 0.0;
    g_flux_viscous[iv][1] = 0.0;
    g_flux_viscous[iv][2] = 0.0;
  }

  m5_facevalues(run,crnt,0,f,0,0);

  for(iv=0;iv<eqtypn[3];++iv){  
    run->source_avf[iv] = run->sol_L->avf[iv];  
  }

  as=0.0;
  if (n_solids==1) {
    as = run->sol_L->avf[0];
  } else if (n_solids==2) {
    as = max(run->sol_L->avf[0],run->sol_L->avf[1]);
  } else if (n_solids==3) {
    as = max(max(run->sol_L->avf[0],run->sol_L->avf[1]),run->sol_L->avf[2]);
  }
  
  A1 = run->sol_L->Amat[0]; B1 = run->sol_L->Amat[3]; C1 = run->sol_L->Amat[6];
  A2 = run->sol_L->Amat[1]; B2 = run->sol_L->Amat[4]; C2 = run->sol_L->Amat[7];
  A3 = run->sol_L->Amat[2]; B3 = run->sol_L->Amat[5]; C3 = run->sol_L->Amat[8];

 
  for(iv=0;iv<eqtypn[4];++iv){         
    for(i=0;i<3;++i){
      for(j=0;j<3;++j){
        run->source_st[iv][i][j] = run->sol_L->avf_stress[iv][i][j];
      }
    }
  }
  
  for(f=0; f<(crnt->nlfc); f++) {  // Loop for faces  

    for(iv=0;iv<NEQ;++iv){  
      crnt->el->flux_face[f][iv] = 0.0;
    } 

    crnt->el->gfm_face[f] = 0.0;

    crnt->el->u_face[f][0] = 0.0;
    crnt->el->u_face[f][1] = 0.0;
    crnt->el->u_face[f][2] = 0.0;
    
    if (crnt->cl->fc[f].bc==0) {        // Face is not boundary

      if (crnt->lfc_neigs[f]==1) {  // 1 neighbor 

        ing=0;
        Nneig=1;
        oppo=crnt->neigtr[f][ing];
        ifc_opp=crnt->neigfc[f];

        if (oppo->lfc_neigs[ifc_opp]==1) { 	// 1:1  connectivity
          flag_con=1;
          if ((oppo->part==crnt->part) && (crnt->el->SI[0]==0)) {  // Opposite part in current processor and its not interface cell
            
            if (crnt->el->flux_flag[f]==1) {  // Flux has been computed from the opposite face                                                  
              f_oppo=crnt->neigfc[f];
   
              crnt->el->u_face[f][0] = -oppo->el->u_face[f_oppo][0];
              crnt->el->u_face[f][1] = -oppo->el->u_face[f_oppo][1];
              crnt->el->u_face[f][2] = -oppo->el->u_face[f_oppo][2];               
              
              for(iv=0;iv<NEQ;++iv){    
                crnt->el->flux_face[f][iv] = -oppo->el->flux_face[f_oppo][iv];
              }
              flag_compute=0;                // Do not compute again                
            } else {
              f_oppo=crnt->neigfc[f];
              oppo->el->flux_flag[f_oppo]=1;    // Flag for flux calculation 
              flag_compute=1;  
            }
            
          } else {
            flag_compute=1; 
          }
        
        } else { // 2:1 connectivity
          flag_con=2;
          flag_compute=1;                             
        }
         
      }

      if (crnt->lfc_neigs[f]!=1) {  // 1:2 connectivity
        flag_con=3;
        Nneig=crnt->lfc_neigs[f];
        flag_compute=1;
      }
        
    } else { // Face is boundary
      Nneig=1;
      flag_compute=1;
    }
    
    Area_f = 0.0;
    for (ineig=0;ineig<Nneig;++ineig) {
      Area[f][ineig] = crnt->cl->Area[f][ineig];
      Area_f += crnt->cl->Area[f][ineig];
    }

    if (flag_compute==1) {
      
      ineig=0;
      bc_temp=crnt->cl->fc[f].bc;
      if      (ORDER==1){ m5_facevalues(run,crnt,ineig,f,bc_temp,0); }  // sol_L
      else if (ORDER==2){ m5_facercnstr(run,crnt,ineig,f,bc_temp,0); }  // sol_L

      for (ineig=0;ineig<Nneig;++ineig) {
        
        if      (ORDER==1){ m5_facevalues(run,crnt,ineig,f,bc_temp,1); }  // sol_R
        else if (ORDER==2){ m5_facercnstr(run,crnt,ineig,f,bc_temp,1); }  // sol_R
        
        if (NS==1) {flux_viscous(run,crnt,run->sol_L,run->sol_R);}
        
        if (run->sol_L->phi[0]*run->sol_R->phi[0]<0.0) { // Solid-fluid interface found
          crnt->el->gfm_face[f] = 1.0;

          if (flag_con!=1) { 
            printf("Problem m5_hllc_gfm %d %d | %d | %f %f\n",flag_con,g_istep,crnt->root,run->sol_L->phi[0],run->sol_R->phi[0]); 
          }
          m5_hllc_gfm(run,crnt,ineig,f);
        } else {
          m5_hllc(run,crnt,ineig,f);
        }     
       

        for(iv=0;iv<NEQ;++iv){  
          crnt->el->flux_face[f][iv] +=  (Area[f][ineig])* 
                                         (((g_flux_adv[iv][0] - g_flux_viscous[iv][0]) * crnt->cl->nx[f]  +
                                           (g_flux_adv[iv][1] - g_flux_viscous[iv][1]) * crnt->cl->ny[f]  +
                                           (g_flux_adv[iv][2] - g_flux_viscous[iv][2]) * crnt->cl->nz[f]) + 
                                            g_flux_HLLC[iv]);
        }
        
        crnt->el->u_face[f][0] = run->u_nc/((double) Nneig);
        crnt->el->u_face[f][1] = run->v_nc/((double) Nneig);
        crnt->el->u_face[f][2] = run->w_nc/((double) Nneig);
        
      } // if compute==1
    
    } // neig loop

    for(iv=0;iv<NEQ;++iv){
      crnt->el->RHS[iv] += crnt->el->flux_face[f][iv]; 
    }

    // --------------------------------------------------------------
    //                     Non-CoNSevative part
    // --------------------------------------------------------------   
    u_nc = crnt->el->u_face[f][0]; nx = crnt->cl->nx[f];
    v_nc = crnt->el->u_face[f][1]; ny = crnt->cl->ny[f];
    w_nc = crnt->el->u_face[f][2]; nz = crnt->cl->nz[f];

    v=eqtypi[3];
    for(iv=0;iv<eqtypn[3];++iv){
      crnt->el->RHS[v + iv] += Area_f*( run->source_avf[iv] * u_nc*nx +
                                        run->source_avf[iv] * v_nc*ny +
                                        run->source_avf[iv] * w_nc*nz);   
    }

    
    v=eqtypi[4];
    for(iv=0; iv<eqtypn[4]; ++iv){
      crnt->el->RHS[v + iv] += Area_f*(  run->source_st[iv][0][0]*u_nc*nx
                                        +run->source_st[iv][0][1]*v_nc*nx
                                        +run->source_st[iv][0][2]*w_nc*nx
                                        +run->source_st[iv][1][0]*u_nc*ny
                                        +run->source_st[iv][1][1]*v_nc*ny
                                        +run->source_st[iv][1][2]*w_nc*ny
                                        +run->source_st[iv][2][0]*u_nc*nz
                                        +run->source_st[iv][2][1]*v_nc*nz
                                        +run->source_st[iv][2][2]*w_nc*nz);
    } 
    
    v=eqtypi[5];
    crnt->el->RHS[v + 0] += Area_f*((A1*v_nc*ny + A1*w_nc*nz) + B1*v_nc*nx - C1*w_nc*nx); // A1 00 1
    crnt->el->RHS[v + 1] += Area_f*((A2*v_nc*ny + A2*w_nc*nz) + B2*v_nc*nx - C2*w_nc*nx); // A2 01 0 
    crnt->el->RHS[v + 2] += Area_f*((A3*v_nc*ny + A3*w_nc*nz) + B3*v_nc*nx - C3*w_nc*nx); // A3 02 0

    crnt->el->RHS[v + 3] += Area_f*((B1*u_nc*nx + B1*w_nc*nz) + A1*u_nc*ny - C1*w_nc*ny); // B1 10 0
    crnt->el->RHS[v + 4] += Area_f*((B2*u_nc*nx + B2*w_nc*nz) + A2*u_nc*ny - C2*w_nc*ny); // B2 11 1
    crnt->el->RHS[v + 5] += Area_f*((B3*u_nc*nx + B3*w_nc*nz) + A3*u_nc*ny - C3*w_nc*ny); // B3 12 0
    
    crnt->el->RHS[v + 6] += Area_f*((C1*u_nc*nx + C1*v_nc*ny) + A1*u_nc*nz - B1*v_nc*nz); // C1 20 0
    crnt->el->RHS[v + 7] += Area_f*((C2*u_nc*nx + C2*v_nc*ny) + A2*u_nc*nz - B2*v_nc*nz); // C2 21 0
    crnt->el->RHS[v + 8] += Area_f*((C3*u_nc*nx + C3*v_nc*ny) + A3*u_nc*nz - B3*v_nc*nz); // C3 22 1  
            
  }

 
  // ======================
  //      Nan - Check
  // ======================
  

  for(iv=0;iv<NEQ;iv++){
    if ((isnan(crnt->el->RHS[iv])==1)){ 
      printf(" \n");
      printf("RHS NaN, eq: %d, element: %d, time-step: %d \n",iv,crnt->root,g_istep_tot);
      printf(" \n");
      exit(0);
    }
  }
  
  for(f=0; f<(crnt->nlfc); f++) {
    free(Area[f]);
  }
  free(Area);

}
