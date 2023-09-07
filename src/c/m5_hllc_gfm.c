#include "strdata.h"

void m5_hllc_gfm(struct RUN *run,struct BRANCH * crnt,int in,int ifc) {

  int i,iv,flag_ifs,flag_end; 
  double u_S,SL,SR,uL,vL,wL,uR,vR,wR;
  double s11L,s11R,s12L,s12R,s13L,s13R;

  // sol to sol1: reverse

  // Left, Right local velocities 
  rotate_vector_gtl_v2(run,crnt,&ifc,&(run->sol_L->u),&(run->sol_L->v),&(run->sol_L->w));
  uL = run->vel_temp[0]; 
  vL = run->vel_temp[1]; 
  wL = run->vel_temp[2]; 

  rotate_vector_gtl_v2(run,crnt,&ifc,&(run->sol_R->u),&(run->sol_R->v),&(run->sol_R->w));
  uR = run->vel_temp[0]; 
  vR = run->vel_temp[1]; 
  wR = run->vel_temp[2];

  //----------------------------------------------
  // changes to the tangestial velocities
  // Local

  vR = vL; // sol to sol1
  wR = wL; // sol to sol1
  
  //----------------------------------------------
  // changes to global
  run->sol_L->u_temp = run->sol_L->u; // for flux_v2 
  run->sol_L->v_temp = run->sol_L->v; 
  run->sol_L->w_temp = run->sol_L->w;

  rotate_vector_ltg_v2(run,crnt,&ifc,&(uR),&(vR),&(wR));
  run->sol_R->u_temp = run->vel_temp[0]; // for flux_v2 
  run->sol_R->v_temp = run->vel_temp[1]; 
  run->sol_R->w_temp = run->vel_temp[2];

  run->sol_R->vec[eqtypi[1]+0] = run->sol_R->u_temp*run->sol_R->r; // for HLLC: csr_v3
  run->sol_R->vec[eqtypi[1]+1] = run->sol_R->v_temp*run->sol_R->r;
  run->sol_R->vec[eqtypi[1]+2] = run->sol_R->w_temp*run->sol_R->r;

  run->sol_R->vec[eqtypi[2]] = run->sol_R->r*run->sol_R->e;
  //----------------------------------------------

  m5_wave_speeds_gfm(run,crnt,ifc,&u_S,&SL,&SR,&uL,&vL,&wL,&s11L,&s12L,&s13L,
                                                   &uR,&vR,&wR,&s11R,&s12R,&s13R);
  
  run->temp_u_S= u_S;
  run->temp_SL = SL;
  run->temp_SR = SR;
  
  for(iv=0;iv<NEQ;++iv){ g_flux_HLLC[iv] = 0.0; }

  flag_ifs = 0; flag_end = 0;

  if((0.0<=SL) && (flag_end==0)){ // Left state
    flag_ifs = 0; flag_end = 1;        
    m5_flux_gfm(run,run->sol_L,flag_ifs);                            
  }

  if(((SL<0.0)&&(0.0<=u_S)) && (flag_end==0)) { // Left star state  
    flag_ifs = 1; flag_end = 1;

    m5_star_region(run,crnt,run->sol_L,ifc,flag_ifs,u_S,SL,SR,uL,vL,wL,uR,vR,wR,s11L,s12L,s13L,s11R,s12R,s13R);
    m5_flux_gfm(run,run->sol_L,flag_ifs);  

    for(iv=0;iv<NCONSEQ;++iv){
      g_flux_HLLC[iv] = SL*(run->Us_temp[iv]-run->sol_L->vec[iv]); 
    }

  }

  if(((u_S<0.0)&&(0.0<SR)) && (flag_end==0)) { // Right star state
    flag_ifs = 2; flag_end = 1;
    m5_star_region(run,crnt,run->sol_R,ifc,flag_ifs,u_S,SL,SR,uL,vL,wL,uR,vR,wR,s11L,s12L,s13L,s11R,s12R,s13R);    
    m5_flux_gfm(run,run->sol_R,flag_ifs); 

    for(iv=0;iv<NCONSEQ;++iv){
      g_flux_HLLC[iv] = SR*(run->Us_temp[iv]-run->sol_R->vec[iv]); 
    }
    
  }

  if((SR<=0.0) && (flag_end==0)){  // Right state
    flag_ifs = 3; flag_end = 1;                 
    m5_flux_gfm(run,run->sol_R,flag_ifs);              
  }
  
  // ======================================================================================
  

  for(iv=0;iv<NEQ;++iv){
    if ((isnan(g_flux_HLLC[iv])==1)){ 
      printf("HLLC nan | %d %d | %d | %f | %f \n",iv,crnt->root,flag_ifs,g_flux_HLLC[iv],run->Us_temp[iv]);
      exit(0);
    }
  }

}