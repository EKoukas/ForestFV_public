#include "strdata.h"

void m3_init_domain(struct RUN *run){

  int i,iv,ieq,j;
  double xw,yw,zw;
  struct BRANCH * crnt;
  
  crnt=run->topo->locl;
  while (crnt!=NULL){
    
    xw=crnt->cl->xc;
    yw=crnt->cl->yc;
    zw=crnt->cl->zc;

    for(iv=0;iv<NEQ;iv++){
      crnt->el->S[iv] = 0.0;
    }

    if (CASE==0) { // Interface advection [Liquid Gas] 

        if((xw)>0.5) {
          m3_init_S(run,crnt,1); // air
        } else {
          m3_init_S(run,crnt,0); // water
        }        

    }
    else if (CASE==1) {   // Shock tube case [Liquid Gas] 
      
      if((xw)>0.5) {
        m3_init_S(run,crnt,1); // air
      } else {
        m3_init_S(run,crnt,0);  // water
      } 

    } 
    else {
      printf("No case initialization found: Exiting");
      exit(0);
    }


    crnt=crnt->lnxt;
  }
}