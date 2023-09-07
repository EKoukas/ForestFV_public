#include "strdata.h"

void m4_reconctruct_V0(struct RUN *run,struct BRANCH *crnt,int f,double dist_p1_cf,double dist_p1_cc,
                                                                 double dist_m1_cf,double dist_m1_cc) {
  
  
  int j,iv;
  double product,avf_dif,avfsum,avf_per;
  double avf[10];

  // Variable extrapolation
  
  for(iv=0;iv<NEQ_TEMP;++iv){    //variables
    g_grad_limited[iv]=0.0;
		
    if (crnt->cl->fc[f].bc==0){  
      if (fabs(g_delta_front[iv])>pow(10.0,-10.0)) {   

				product = (-g_delta_back[iv]/dist_m1_cc)*(g_delta_front[iv]/dist_p1_cc);

        if (product>0.0) {
					if ( fabs(-g_delta_back[iv]/dist_m1_cc) > fabs(g_delta_front[iv]/dist_p1_cc)) {
						g_grad_limited[iv] =  g_limitedvars[iv]* (g_delta_front[iv]/dist_p1_cc);
					} else {
						g_grad_limited[iv] =  g_limitedvars[iv]* (-g_delta_back[iv]/dist_m1_cc);
					}
					
        } 

      }
    }

		if (PRIMTV==1){
			g_S_prim_rec_L[iv]  = g_S_prim_L[iv] + dist_p1_cf*g_grad_limited[iv];
		} else {
      crnt->el->SF[f][iv] = g_S_cons_L[iv] + dist_p1_cf*g_grad_limited[iv];
		}
        
  }


  // ================================================================================
	// VF limit
  avfsum=0.0;
  for(iv=0;iv<eqtypn[3];++iv){
    
    if      (PRIMTV==0){ avf[iv] = crnt->el->SF[f][iv+eqtypi[3]]; } 
    else if (PRIMTV==1){ avf[iv] = g_S_prim_rec_L[iv];            }
 
    if      (avf[iv]<= AMIN/100.0)              {  avf[iv] =           AMIN/100.0;   } 
    else if (avf[iv]>=(1.0-(2.0*(AMIN/100.0)))) {  avf[iv] = 1.0-(2.0*(AMIN/100.0)); }    
    avfsum  += avf[iv];

  }

  if (avfsum>=(1.0-(AMIN/100.0))) { 
    avf_dif = avfsum - (1.0-(AMIN/100.0));
    for(j=0;j<eqtypn[3];++j){
      avf_per       = avf[j]/avfsum;
      avf[j] = avf[j] - avf_per*avf_dif;
    }
    avfsum  = 1.0-(AMIN/100.0);
  }
  avf[eqtypn[3]] = 1.0 - avfsum;
  
  if (PRIMTV==1) {
    for(iv=0;iv<eqtypn[3];++iv){
      g_S_prim_rec_L[iv] = avf[iv];
    }
  } else {
    for(iv=0;iv<eqtypn[3];++iv){
      crnt->el->SF[f][iv+ eqtypi[3]] = avf[iv];
    }
  }
   // ================================================================================

  if (PRIMTV==1){
    m4_prim2conc(g_S_prim_rec_L,g_S_cons_rec_L);  // g_S_prim_rec_L to g_S_cons_rec_L
    for(iv=0;iv<NEQ;++iv){ run->Vec_temp_R[iv] = g_S_cons_rec_L[iv];  }
    m4_relaxationvec(run,crnt);
    for(iv=0;iv<NEQ;++iv){ crnt->el->SF[f][iv] = run->Vec_temp_R[iv]; }
  } else {
    for(iv=0;iv<NEQ;++iv){ run->Vec_temp_R[iv] = crnt->el->SF[f][iv];  }
    m4_relaxationvec(run,crnt);
    for(iv=0;iv<NEQ;++iv){ crnt->el->SF[f][iv] = run->Vec_temp_R[iv]; }
  }
  
  return;

}


void m4_reconctruct_V1(struct RUN *run,struct BRANCH *crnt,int f, int fo) {
  
  int i,j,iv;
  double dist_p1_cc,dist_p1_cf,dist_m1_cc;
  double product,avf_dif,avfsum,avf_per;
  double avf[10];

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

  // ================================================================================
	// VF limit
  avfsum=0.0;
  for(iv=0;iv<eqtypn[3];++iv){
    
    if      (PRIMTV==0){ avf[iv] = crnt->el->SF[f][iv+eqtypi[3]]; } 
    else if (PRIMTV==1){ avf[iv] = g_S_prim_rec_L[iv];            }
 
    if      (avf[iv]<= AMIN/100.0)              {  avf[iv] =           AMIN/100.0;   } 
    else if (avf[iv]>=(1.0-(2.0*(AMIN/100.0)))) {  avf[iv] = 1.0-(2.0*(AMIN/100.0)); }    
    avfsum  += avf[iv];

  }

  if (avfsum>=(1.0-(AMIN/100.0))) { 
    avf_dif = avfsum - (1.0-(AMIN/100.0));
    for(j=0;j<eqtypn[3];++j){
      avf_per       = avf[j]/avfsum;
      avf[j] = avf[j] - avf_per*avf_dif;
    }
    avfsum  = 1.0-(AMIN/100.0);
  }
  avf[eqtypn[3]] = 1.0 - avfsum;
  
  if (PRIMTV==1) {
    for(iv=0;iv<eqtypn[3];++iv){
      g_S_prim_rec_L[iv] = avf[iv];
    }
  } else {
    for(iv=0;iv<eqtypn[3];++iv){
      crnt->el->SF[f][iv+ eqtypi[3]] = avf[iv];
    }
  }
  // ================================================================================

  if (PRIMTV==1){
    m4_prim2conc(g_S_prim_rec_L,g_S_cons_rec_L);  // g_S_prim_rec_L to g_S_cons_rec_L
    for(iv=0;iv<NEQ;++iv){ run->Vec_temp_R[iv] = g_S_cons_rec_L[iv];  }
    m4_relaxationvec(run,crnt);
    for(iv=0;iv<NEQ;++iv){ crnt->el->SF[f][iv] = run->Vec_temp_R[iv]; }
  } else {
    for(iv=0;iv<NEQ;++iv){ run->Vec_temp_R[iv] = crnt->el->SF[f][iv];  }
    m4_relaxationvec(run,crnt);
    for(iv=0;iv<NEQ;++iv){ crnt->el->SF[f][iv] = run->Vec_temp_R[iv]; }
  }
 
  
  return;

}