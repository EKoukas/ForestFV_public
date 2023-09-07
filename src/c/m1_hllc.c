#include "strdata.h"

void m1_hllc(struct RUN *run,struct BRANCH * crnt,int in,int ifc) {

  int i,iv,flags_if;
  double S_S,SL,SR,uL,vL,wL,uR,vR,wR;
  double s11L,s11R;

  m1_wave_speeds(run,crnt,ifc,&S_S,&SL,&SR,&uL,&vL,&wL,&s11L,
                                           &uR,&vR,&wR,&s11R);
    
  flags_if=0;

  for(iv=0;iv<NEQ;++iv){ g_flux_HLLC[iv] = 0.0; }

  if (crnt->cl->fc[ifc].bc==5){S_S=0.0;}

  // ==================================================================================
  if ((0.0<=SL)&&(flags_if==0)) { // Left state
    flags_if=1;
    m1_flux(run->sol_L);                                    
  }  
  
  if (( (SL<=0.0)&&(0.0<=S_S) ) && (flags_if==0)) { // Left star state  
    flags_if=2;
    m1_star_region(run,crnt,run->sol_L,ifc,1,S_S,SL,SR,uL,vL,wL,uR,vR,wR);
    m1_flux(run->sol_L);          
    for(iv=0;iv<NEQ;++iv){
      g_flux_HLLC[iv] = SL*(run->Us_temp[iv]-run->sol_L->vec[iv]); 
    }
  } 
  
  if (( (S_S<=0.0)&&(0.0<=SR) ) && (flags_if==0)) { // Right star state
    flags_if=3;
    m1_star_region(run,crnt,run->sol_R,ifc,2,S_S,SL,SR,uL,vL,wL,uR,vR,wR);
    m1_flux(run->sol_R);          
    for(iv=0;iv<NEQ;++iv){
      g_flux_HLLC[iv] = SR*(run->Us_temp[iv]-run->sol_R->vec[iv]); 
    }
  } 
  
  if ((SR<=0.0)&&(flags_if==0)) {  // Right state
    flags_if=4;
    m1_flux(run->sol_R);                       
  }
  // ==================================================================================

  if (flags_if==0){
    printf("m1_hllc: stop if %e %e %e \n",S_S,SL,SR);
    exit(0);
  }
  
  for(iv=0;iv<NEQ;++iv){
    if ((isnan(g_flux_HLLC[iv])==1)){ 
      printf("HLLC nan %d \n",iv);
      exit(0);
    }
  }

}