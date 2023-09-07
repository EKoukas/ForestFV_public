#include "strdata.h"

void m2_deltavalues_V0(struct RUN *run,struct BRANCH *crnt,double * delta_temp,int f_temp) {

  int iv,bc_temp;
  
  bc_temp = crnt->cl->fc[f_temp].bc;
  m2_facevalues_MUSCL(run,crnt,f_temp,0,bc_temp,g_S_cons_L);  // gives value to g_S_cons_L: matrix with conservative values from the center of the cell 
  m2_facevalues_MUSCL(run,crnt,f_temp,1,bc_temp,g_S_cons_R);  // gives value to g_S_cons_R: matrix with conservative values from the neigboring cell or BC
  
  if (PRIMTV==1){
    m2_cons2primtv(g_S_cons_L,g_S_prim_L);  // g_S_cons_L  to g_S_prim_L
    m2_cons2primtv(g_S_cons_R,g_S_prim_R);  // g_S_cons_R  to g_S_prim_R
    for(iv=0;iv<NPRIMITIV;++iv) { delta_temp[iv]=g_S_prim_R[iv]-g_S_prim_L[iv]; }
  } else {
    for(iv=0;iv<NEQ;++iv)       { delta_temp[iv]=g_S_cons_R[iv]-g_S_cons_L[iv]; }
  }
    
  return;

}

void m2_deltavalues_V1(struct RUN *run,struct BRANCH *crnt,int f,int fo) {

  int iv,bc_temp,f_temp;
  
  //===================================================================================================
  f_temp = f;

  bc_temp = crnt->cl->fc[f_temp].bc;
  m2_facevalues_MUSCL(run,crnt,f_temp,1,bc_temp,g_S_cons_R);  // gives value to g_S_cons_R: matrix with conservative values from the neigboring cell or BC
  
  if (PRIMTV==1){
    m2_cons2primtv(g_S_cons_R,g_S_prim_R);  // g_S_cons_R  to g_S_prim_R
    for(iv=0;iv<NPRIMITIV;++iv) { g_delta_front[iv]=g_S_prim_R[iv]-g_S_prim_L[iv]; }
  } else {
    for(iv=0;iv<NEQ;++iv)       { g_delta_front[iv]=g_S_cons_R[iv]-g_S_cons_L[iv]; }
  }
  //===================================================================================================

  //===================================================================================================
  f_temp = fo;

  bc_temp = crnt->cl->fc[f_temp].bc;
  m2_facevalues_MUSCL(run,crnt,f_temp,1,bc_temp,g_S_cons_R);  // gives value to g_S_cons_R: matrix with conservative values from the neigboring cell or BC
  
  if (PRIMTV==1){
    m2_cons2primtv(g_S_cons_R,g_S_prim_R);  // g_S_cons_R  to g_S_prim_R
    for(iv=0;iv<NPRIMITIV;++iv) { g_delta_back[iv]=g_S_prim_R[iv]-g_S_prim_L[iv]; }
  } else {
    for(iv=0;iv<NEQ;++iv)       { g_delta_back[iv]=g_S_cons_R[iv]-g_S_cons_L[iv]; }
  }
  //===================================================================================================

  return;

}