#include "strdata.h"

void m3_hllc(struct RUN *run,struct BRANCH * crnt,int in,int ifc) {

  int i,iv,ieq,flag_ifs,flag_end; 
  double u_S,SL,SR,uL,vL,wL,uR,vR,wR;
  double s11L,s11R,s12L,s12R,s13L,s13R;

  flag_ifs = 0; flag_end = 0; 
 
  m3_wave_speeds(run,crnt,ifc,&u_S,&SL,&SR,&uL,&vL,&wL,&s11L,&s12L,&s13L,
                                           &uR,&vR,&wR,&s11R,&s12R,&s13R);

  if (crnt->cl->fc[ifc].bc==5){u_S=0.0;}

  run->temp_u_S = u_S;
  run->temp_SL  = SL;
  run->temp_SR  = SR;
  
  if((0.0<=SL) && (flag_end==0)){ // Left state
    flag_ifs = 0; flag_end = 1;
    m3_flux(run,run->sol_L,flag_ifs);                            
  }

  if(((SL<0.0)&&(0.0<=u_S)) && (flag_end==0)) { // Left star state  
    flag_ifs = 1; flag_end = 1;
    m3_star_region(run,crnt,run->sol_L,ifc,flag_ifs,
                   u_S,SL,SR,uL,vL,wL,uR,vR,wR,s11L,s12L,s13L,s11R,s12R,s13R);
    m3_flux(run,run->sol_L,flag_ifs); 

    for(iv=0;iv<NCONSEQ;++iv){
      g_flux_HLLC[iv] = SL*(run->Us_temp[iv]-run->sol_L->vec[iv]); 
    }
  }

  if(((u_S<0.0)&&(0.0<SR)) && (flag_end==0)) { // Right star state
    flag_ifs = 2; flag_end = 1;
    m3_star_region(run,crnt,run->sol_R,ifc,flag_ifs,
                           u_S,SL,SR,uL,vL,wL,uR,vR,wR,s11L,s12L,s13L,s11R,s12R,s13R);
    m3_flux(run,run->sol_R,flag_ifs); 

    for(iv=0;iv<NCONSEQ;++iv){
      g_flux_HLLC[iv] = SR*(run->Us_temp[iv]-run->sol_R->vec[iv]);  
    }   
  }

  if((SR<=0.0) && (flag_end==0)){  // Right state
    flag_ifs = 3; flag_end = 1;
    m3_flux(run,run->sol_R,flag_ifs);              
  }
  
 
  // ======================================================================================
  if (crnt->cl->fc[ifc].bc==1){ // Symmetry

    flag_ifs = 4;
    for(iv=0;iv<NEQ;++iv){         
      g_flux_HLLC[iv] =  0.0; 
      for(i=0;i<3;i++){ g_flux_sym_adv[iv][i] = 0.0;}
    }

    m3_facevalues(run,crnt,in,ifc,0,0);
    m3_facevalues(run,crnt,in,ifc,0,1);
    m3_bc(crnt,crnt->cl->fc[ifc].bc,ifc,run);

    run->u_nc = 0.5*(run->sol_L->u+run->sol_R->u);  
    run->v_nc = 0.5*(run->sol_L->v+run->sol_R->v);
    run->w_nc = 0.5*(run->sol_L->w+run->sol_R->w);  
    
    m3_flux(run,run->sol_L,flag_ifs);
    for(iv=0;iv<NEQ;iv++){       // 0->17
      for(i=0;i<3;i++){ g_flux_sym_adv[iv][i] += 0.5*g_flux_adv[iv][i]; } 
    }

    m3_flux(run,run->sol_R,flag_ifs);
    for(iv=0;iv<NEQ;iv++){       // 0->17
      for(i=0;i<3;i++){ g_flux_sym_adv[iv][i] += 0.5*g_flux_adv[iv][i]; } 
    }

    for(iv=0;iv<NEQ;iv++){       // 0->17
      for(i=0;i<3;i++){ g_flux_adv[iv][i] = g_flux_sym_adv[iv][i]; } 
    }
  }
  
  // ======================================================================================

  for(iv=0;iv<NEQ;++iv){
    if ((isnan(g_flux_HLLC[iv])==1)){ 
      printf("HLLC nan %d \n",iv);
      exit(0);
    }
  }

}