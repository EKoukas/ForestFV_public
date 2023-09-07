#include "strdata.h"

void m4_relaxation_fluid(struct RUN * run) {

  int iv;
  struct BRANCH * crnt;

  run->Trelax=timecpu(run->Trelax,0);

  crnt=run->topo->locl;  
  while (crnt!=NULL){    

    for(iv=0;iv<NEQ;++iv){ run->Vec_temp_R[iv] = crnt->el->S[iv]; }
    m4_relaxationvec(run,crnt); 
    for(iv=0;iv<NEQ;++iv){ crnt->el->S[iv] = run->Vec_temp_R[iv]; }

    crnt=crnt->lnxt;            
  }

   run->Trelax=timecpu(run->Trelax,1);
}