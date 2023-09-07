//-------------------------------------------------------
//      Computes Right Hand Side of the equations
//-------------------------------------------------------

#include "strdata.h"

void m4_residual(struct RUN * run, struct BRANCH * crnt) {

  int iv,v,i,j,f;
  int ing,f_oppo,ifc_opp,ineig,Nneig,bc_temp;
  int flag_con;
  double u_nc,v_nc,w_nc;
  double nx,ny,nz,Area_f;
  double ** Area;
  int flag_compute;
  struct BRANCH * oppo;
  
  Area=malloc(crnt->nlfc*sizeof(double *));
  for(f=0; f<(crnt->nlfc); f++) {
    Area[f]=malloc(4*sizeof(double));
  }
  
  for(iv=0;iv<NEQ;++iv){ 
    crnt->el->RHS[iv] = 0.0;
    g_flux_viscous[iv][0] = 0.0;
    g_flux_viscous[iv][1] = 0.0;
    g_flux_viscous[iv][2] = 0.0;
  }

  m4_facevalues(run,crnt,0,0,0,0);

  for(iv=0;iv<eqtypn[3];++iv){  
    run->source_avf[iv] = run->sol_L->avf[iv];  
  }

  for(iv=0;iv<eqtypn[4];++iv){         
    for(i=0;i<3;++i){  
      for(j=0;j<3;++j){
        run->source_st[iv][i][j] = 0.0;
      }
      run->source_st[iv][i][i] = -run->sol_L->avf[iv]*run->sol_L->p_hydro[iv];
    }
  }
 
  for(f=0; f<(crnt->nlfc); f++) {  // Loop for faces  

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
      
      for(iv=0;iv<NEQ;++iv){  
        crnt->el->flux_face[f][iv] = 0.0;
      } 

      ineig=0;
      bc_temp=crnt->cl->fc[f].bc;
      if      (ORDER==1){ m4_facevalues(run,crnt,ineig,f,bc_temp,0); }  // sol_L
      else if (ORDER==2){ m4_facercnstr(run,crnt,ineig,f,bc_temp,0); }  // sol_L

      for (ineig=0;ineig<Nneig;++ineig) {
         
        if      (ORDER==1){ m4_facevalues(run,crnt,ineig,f,bc_temp,1); }  // sol_R
        else if (ORDER==2){ m4_facercnstr(run,crnt,ineig,f,bc_temp,1); }  // sol_R
                
        if (NS==1) {flux_viscous(run,crnt,run->sol_L,run->sol_R);}

        m4_hllc(run,crnt,ineig,f); 
        
        for(iv=0;iv<NEQ;++iv){  
          crnt->el->flux_face[f][iv] += (Area[f][ineig])* 
                                         (((g_flux_adv[iv][0] - g_flux_viscous[iv][0]) * crnt->cl->nx[f]  +
                                           (g_flux_adv[iv][1] - g_flux_viscous[iv][1]) * crnt->cl->ny[f]  +
                                           (g_flux_adv[iv][2] - g_flux_viscous[iv][2]) * crnt->cl->nz[f]) + 
                                            g_flux_HLLC[iv]);
        }
        
        crnt->el->u_face[f][0] += run->u_nc/((double) Nneig);
        crnt->el->u_face[f][1] += run->v_nc/((double) Nneig);
        crnt->el->u_face[f][2] += run->w_nc/((double) Nneig);
        
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
  }

  // ======================
  //      Nan - Check
  // ======================
  for(v=0;v<NEQ;v++){
        
    if ((isnan(crnt->el->RHS[v])==1)){ 
      printf(" \n");
      printf("RHS NaN, eq: %d, element: %d, time-step: %d \n",v,crnt->root,g_istep_tot);
      printf(" \n");
      exit(0);
    }
  }
  
  for(f=0; f<(crnt->nlfc); f++) {
    free(Area[f]);
  }
  free(Area);

}
