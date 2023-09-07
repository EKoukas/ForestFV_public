#include "strdata.h"

void m1_init_domain(struct RUN *run){

  int i,iv,ieq,j;
  double xw,yw,zw;
  double rho,u,v,w,p,dist_wall,R0;
  struct BRANCH * crnt;
  
  crnt=run->topo->locl;
  while (crnt!=NULL){
    
    xw=crnt->cl->xc;
    yw=crnt->cl->yc;
    zw=crnt->cl->zc;

    for(iv=0;iv<NEQ;iv++){
      crnt->el->S[iv] = 0.0;
    }
    if (CASE==1) {          // Model 1, shocktube case barotropic
     
     if(xw<0.0) {
        rho = 1000.0;
        u   = 0.0;
        v   = 0.0;
        w   = 0.0;
        p   = 3810000.0;
        m1_init_S(run,crnt,rho,u,v,w);  
      } else {
        rho =  0.998;
        u   =  435;
        v   =  0.0;
        w   =  0.0;
        p   =  887.8;
        m1_init_S(run,crnt,rho,u,v,w);
      }

    }
    else if (CASE==1) {          // Model 1, shocktube case barotropic
     
     if(xw<0.0) {
        rho = 1000.0;
        u   = 0.0;
        v   = 0.0;
        w   = 0.0;
        p   = 3810000.0;
        m1_init_S(run,crnt,rho,u,v,w);  
      } else {
        rho =  0.998;
        u   =  435;
        v   =  0.0;
        w   =  0.0;
        p   =  887.8;
        m1_init_S(run,crnt,rho,u,v,w);
      }

    }
    else if (CASE==2500) {  // Model 1, bubble-wall collapse 
      R0 = 0.4/1000.0;
      dist_wall = 0.016/1000.0;
      if ( (pow((xw),2.0) + pow((yw-R0-dist_wall),2.0)) <= pow(R0,2.0)) {     // atm pressure gas
        p=2311.5;   // low pressure
        rho=49.91;    // low deNSity
        u   = 0.0;
        v   = 0.0;
        w   = 0.0;
        m1_init_S(run,crnt,rho,u,v,w);
      }
      else{
        p=pow(10.0,7); // high pressure
        rho=1002.89;               // high deNSity
        u   = 0.0;
        v   = 0.0;
        w   = 0.0;
        m1_init_S(run,crnt,rho,u,v,w);
      }
    } 
    else { 
      printf("No case initialization found: Exiting");
      exit(0);
    }

    crnt=crnt->lnxt;
  }

}