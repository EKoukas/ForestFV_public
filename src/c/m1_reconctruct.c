#include "strdata.h"

void m1_reconctruct_V0(struct RUN *run,struct BRANCH *crnt,int f,double dist_p1_cf,double dist_p1_cc,
                                                                 double dist_m1_cf,double dist_m1_cc) {
  
  int iv;
  double product;

  // Variable extrapolation
  
  for(iv=0;iv<NEQ_TEMP;++iv){    //variables
    g_grad_limited[iv]=0.0;
    if (crnt->cl->fc[f].bc==0){  
      if (fabs(g_delta_front[iv])>pow(10.0,-10.0)) {   

				product = (-g_delta_back[iv]/dist_m1_cc)*(g_delta_front[iv]/dist_p1_cc);

        if (product>0.0) {
					if ( fabs(-g_delta_back[iv]/dist_m1_cc) > fabs(g_delta_front[iv]/dist_p1_cc)) {
						g_grad_limited[iv] =  g_delta_front[iv]/dist_p1_cc;
					} else {
						g_grad_limited[iv] = -g_delta_back[iv]/dist_m1_cc;
					}
					
        } else {
          g_grad_limited[iv] = 0.0;
        }
      }
    }

		if (PRIMTV==1){
			g_S_prim_rec_L[iv]   = g_S_prim_L[iv] + dist_p1_cf*g_grad_limited[iv];
		} else {
      crnt->el->SF[f][iv]  = g_S_cons_L[iv] + dist_p1_cf*g_grad_limited[iv];
		}
        
  }

  if (PRIMTV==1){
    m1_prim2conc(g_S_prim_rec_L,g_S_cons_rec_L);  // g_S_prim_rec_L to g_S_cons_rec_L
    for(iv=0;iv<NEQ;++iv){        
      crnt->el->SF[f][iv] = g_S_cons_rec_L[iv];
    }
  }
 
  
  return;

}


void m1_reconctruct_V1(struct RUN *run,struct BRANCH *crnt,int f, int fo) {
  
  int iv;
  double dist_p1_cc,dist_p1_cf,dist_m1_cc;
  double product;

  // Variable extrapolation
  
  for(iv=0;iv<NEQ_TEMP;++iv) {    //variables
    g_grad_limited[iv]=0.0;
    if (crnt->cl->fc[f].bc==0) {  
      if (fabs(g_delta_front[iv])>pow(10.0,-10.0)) {   
        
        dist_p1_cc = crnt->cl->dist_cc[f];
        dist_m1_cc = crnt->cl->dist_cc[fo];

        dist_p1_cf = crnt->cl->dist_cf[f]; 
        
				product = (-g_delta_back[iv]/dist_m1_cc)*(g_delta_front[iv]/dist_p1_cc);  // ((i) - (i-1))/d * ((i+1) - (i))/d

        if (product>0.0) {
					if ( fabs(g_delta_back[iv]/dist_m1_cc) > fabs(g_delta_front[iv]/dist_p1_cc)) {
						g_grad_limited[iv] =  g_limitedvars[iv]*g_delta_front[iv]/dist_p1_cc;
					} else {
						g_grad_limited[iv] = -g_limitedvars[iv]*g_delta_back[iv]/dist_m1_cc;
					}
        } else {
          g_grad_limited[iv] = 0.0;
        }
      }
    }

		if (PRIMTV==1){
			g_S_prim_rec_L[iv]   = g_S_prim_L[iv] + dist_p1_cf*g_grad_limited[iv];
		} else {
      crnt->el->SF[f][iv]  = g_S_cons_L[iv] + dist_p1_cf*g_grad_limited[iv];
		}
        
  }

  if (PRIMTV==1) {
    m1_prim2conc(g_S_prim_rec_L,g_S_cons_rec_L);  // g_S_prim_rec_L to g_S_cons_rec_L
    for(iv=0;iv<NEQ;++iv){        
      crnt->el->SF[f][iv] = g_S_cons_rec_L[iv];
    }
  }
 
  
  return;

}