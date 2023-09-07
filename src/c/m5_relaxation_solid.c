#include "strdata.h"

void m5_relaxation_solid(struct RUN * run) {

  int iv,e;
  struct BRANCH * crnt;
  double phi_old[4],phi_new[4];

  run->Trelax=timecpu(run->Trelax,0);

  crnt=run->topo->locl;  
  while (crnt!=NULL){    

    //if (PLASTIC==1) {plastic_rlx(crnt,run);} 

    for(iv=0;iv<eqtypn[3];++iv) {   
  	  phi_old[iv] = crnt->el->S[eqtypi[3]+iv]-AVF_LIM;  // not relaxaxed volume fraction
    }

    for(iv=0;iv<NEQ;++iv){ run->Vec_temp_R[iv] = crnt->el->S[iv]; }
    m5_relaxationvec(run,crnt); 
    for(iv=0;iv<NEQ;++iv){ crnt->el->S[iv] = run->Vec_temp_R[iv]; }

    
    for(iv=0;iv<eqtypn[3];++iv) {
      if (MATERMUSH[iv]!=0.0) {
        e = eqtypi[3] + iv;
        phi_new[iv] = crnt->el->S[e] - AVF_LIM; // relaxaxed volume fraction
      }
    }

    if (phi_old[0]*phi_new[0]<0.0) {
      m5_gfm_rel(crnt,run);
    }
    

    crnt=crnt->lnxt;            
  }

  run->Trelax=timecpu(run->Trelax,1);
}