#include "strdata.h"

void m4_hllc (struct RUN *run,struct BRANCH * crnt,int in,int ifc) {

  int i,iv,ieq,flag_ifs,flag_end; 
  double u_S,SL,SR,uL,vL,wL,uR,vR,wR;
  double s11L,s11R,s12L,s12R,s13L,s13R;

  flag_ifs = 0; flag_end = 0; 
 
  m4_wave_speeds(run,crnt,ifc,&u_S,&SL,&SR,&uL,&vL,&wL,&s11L,&s12L,&s13L,
                                           &uR,&vR,&wR,&s11R,&s12R,&s13R);
  
  if (crnt->cl->fc[ifc].bc==4){u_S=0.0;}

  run->temp_u_S= u_S;
  run->temp_SL = SL;
  run->temp_SR = SR;
  
  for(iv=0;iv<NEQ;++iv){ g_flux_HLLC[iv] =  0.0; }

  if((0.0<=SL) && (flag_end==0)){ // Left state
    flag_ifs = 0; flag_end = 1;
    m4_flux(run,run->sol_L,flag_ifs);                                
  }

  if(((SL<=0.0)&&(0.0<=u_S)) && (flag_end==0)) { // Left star state  
    flag_ifs = 1; flag_end = 1;
    m4_star_region(run,crnt,run->sol_L,ifc,flag_ifs,
                   u_S,SL,SR,uL,vL,wL,uR,vR,wR,s11L,s12L,s13L,s11R,s12R,s13R);
    m4_flux(run,run->sol_L,flag_ifs); 

    for(iv=0;iv<NCONSEQ;++iv){
      g_flux_HLLC[iv] = SL*(run->Us_temp[iv]-run->sol_L->vec[iv]); 
    }
  }

  if(((u_S<=0.0)&&(0.0<=SR)) && (flag_end==0)) { // Right star state
    flag_ifs = 2; flag_end = 1;
    m4_star_region(run,crnt,run->sol_R,ifc,flag_ifs,
                   u_S,SL,SR,uL,vL,wL,uR,vR,wR,s11L,s12L,s13L,s11R,s12R,s13R);
    m4_flux(run,run->sol_R,flag_ifs); 

    for(iv=0;iv<NCONSEQ;++iv){
      g_flux_HLLC[iv] = SR*(run->Us_temp[iv]-run->sol_R->vec[iv]);  
    }
  }

  if((SR<=0.0) && (flag_end==0)){  // Right state
    flag_ifs = 3; flag_end = 1;
    m4_flux(run,run->sol_R,flag_ifs);              
  }

  for(iv=0;iv<NEQ;++iv){
    if ((isnan(g_flux_HLLC[iv])==1)){ 
      printf("m4_hllc nan %d \n",iv);
      exit(0);
    }
  }

}