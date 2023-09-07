#include "strdata.h"

void m2_init_domain(struct RUN *run){

  int i,iv,ieq,j;
  double xw,yw,zw;
  double rho,u,v,w,p;
  double R0,R1,dist_wall,ymass;
  struct BRANCH * crnt;

  crnt=run->topo->locl;
  while (crnt!=NULL){
    
    xw=crnt->cl->xc;
    yw=crnt->cl->yc;
    zw=crnt->cl->zc;

    for(iv=0;iv<NEQ;iv++){
      crnt->el->S[iv] = 0.0;
    }

    if (CASE==1) {          // Model 2, Case A (page 16) https://dx.doi.org/10.6084/m9.figshare.c.4388069.
      
      if(crnt->root<500) {  // left

        rho =  0.0170;
        //rho =  1.055;
        u   =  0.0; //100.0;
        v   =  0.0;
        w   =  0.0;
        //p   =  140.0;
        ymass = 1-YMIN;
        
       // ymass = 1.0;
        m2_init_S(run,crnt,rho,u,v,w,ymass);

        
      } else {  // Liquid
       rho = 998.200;
       // rho = 1500.202;
        u   = 0.0; //100.0;
        v   = 0.0;
        w   = 0.0;

        ymass = YMIN;
       // ymass = 0.0;
        m2_init_S(run,crnt,rho,u,v,w,ymass);  
      }

    }
    else if (CASE==2) {     // Model 3, Case C (page 16) https://dx.doi.org/10.6084/m9.figshare.c.4388069.
     
      if(crnt->root<500) {  // left
        rho = 1000.0;
        u   = 0.0;
        v   = 0.0;
        w   = 0.0;
        p   = 3810000.0;
        ymass=0.0;
        m2_init_S(run,crnt,rho,u,v,w,ymass);   
      } else {
        rho =  0.998;
        u   =  435;
        v   =  0.0;
        w   =  0.0;
        p   =  887.8;
        ymass=1.0;
        m2_init_S(run,crnt,rho,u,v,w,ymass); 
      }
    }
    else if (CASE==2100) {  // Model 3, needlesfree injection
      R0 = 0.24/1000.0;
      R1 = 0.08/1000.0;
      dist_wall = 1.64/1000.0;
      if ( (pow((xw-dist_wall),2.0) + pow((yw),2.0)) <= pow(R0,2.0) || xw >=dist_wall ){ // gas
        p=101325.0;
        rho=1.225;
        u  = 0.0;
        v  = 0.0;
        w  = 0.0;
        ymass = 1.0;
        m2_init_S(run,crnt,rho,u,v,w,ymass);
      } else if ( (pow((xw),2.0) + pow((yw),2.0)) <= pow(R1,2.0) && (xw<=(0.08/1000))){ // 
        p=5.0*pow(10.0,7);
        rho=303.96; // 300 degrees
        u = 0.0;
        v = 0.0;
        w = 0.0;
        ymass = 1.0;
        m2_init_S(run,crnt,rho,u,v,w,ymass);
      } else {
        p=101325.0;
        rho=998.17;
        u = 0.0;
        v = 0.0;
        w = 0.0;
        ymass = 0.0;
        m2_init_S(run,crnt,rho,u,v,w,ymass);
      }
    }
    else if (CASE==2900) {  // injector

      if (xw <= 0.056) {    
        rho   =998.20;    // hight deNSity
        u     = 0.0;
        v     = 0.0;
        w     = 0.0;
        ymass = 0.0;
        m2_init_S(run,crnt,rho,u,v,w,ymass);
      } else {
        //p=pow(10.0,7); // high pressure
        rho   = 1.225;               // low deNSity
        u     = 0.0;
        v     = 0.0;
        w     = 0.0;
        ymass = 1.0;
        m2_init_S(run,crnt,rho,u,v,w,ymass);
      }
    } else {
      printf("No case initialization found: Exiting");
      exit(0);
    }
    crnt=crnt->lnxt;
      
  }

}