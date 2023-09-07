#include "strdata.h"

void m0_init_domain(struct RUN *run){

  int i,iv,ieq,j;
  double xw,yw,zw;
  double rho,u,v,w,p;
  struct BRANCH * crnt;
  
  crnt=run->topo->locl;
  while (crnt!=NULL){
    
    xw=crnt->cl->xc;
    yw=crnt->cl->yc;
    zw=crnt->cl->zc;

    for(iv=0;iv<NEQ;iv++){
      crnt->el->S[iv] = 0.0;
    }

    if (CASE==1) {        // Model 0, Test 1, Air, gamma=1.4, p_inf=0, Toro pg 225, table 6.2 

      if(crnt->root<700) {  // left
        rho = 1.0;
        u   = 0.75;
        v   = 0.0;
        w   = 0.0;
        p   = 1.0;
        m0_init_S(run,crnt,rho,u,v,w,p);
      } else {  // right
        rho = 0.125;
        u   = 0.0;
        v   = 0.0;
        w   = 0.0;
        p   = 0.1;
        m0_init_S(run,crnt,rho,u,v,w,p);
      }        

    }
    else if (CASE==2) {   // Model 0, Test 2, Air, gamma=1.4, p_inf=0, Toro pg 225, table 6.2

      if(crnt->root>500) {  // left
        
        rho = 1.0;
        u   = -2.0;
        v   = 0.0;
        w   = 0.0;
        p   = 0.4;
        m0_init_S(run,crnt,rho,u,v,w,p);

      } else {  // right

        rho = 1.0;
        u   = 2.0;
        v   = 0.0;
        w   = 0.0;
        p   = 0.4;
        m0_init_S(run,crnt,rho,u,v,w,p);

      }        

    }
    else if (CASE==3) {   // Model 0, Test 3, Air, gamma=1.4, p_inf=0, Toro pg 225, table 6.2

      if(crnt->root<500) {  // left
        rho = 1.0;
        u   = 0.0;
        v   = 0.0;
        w   = 0.0;
        p   = 0.01;
        m0_init_S(run,crnt,rho,u,v,w,p);
      } else {  // right
        rho = 1.0;
        u   = 0.0;
        v   = 0.0;
        w   = 0.0;
        p   = 1000.0;
        m0_init_S(run,crnt,rho,u,v,w,p);
      }  
    }      
    else if (CASE==4) {   // Model 0, Test 4, Air, gamma=1.4, p_inf=0, Toro pg 225, table 6.2

      if(crnt->root<600) {  // left
         rho = 5.99924;
        u   = -6.19633;
        v   = 0.0;
        w   = 0.0;
        p   = 46.095;
        m0_init_S(run,crnt,rho,u,v,w,p);
      } else {  // right
        rho = 5.99924;
        u   = 19.5975;
        v   = 0.0;
        w   = 0.0;
        p   = 460.894;
        m0_init_S(run,crnt,rho,u,v,w,p);
      }        

    }
    else if (CASE==200) { // Model 0, Forward faceing step, https://amroc.sourceforge.net/examples/euler/2d/html/ffstep_n.htm
  
      rho = 1.4;
      u   = 3.0;
      v   = 0.0;
      w   = 0.0;
      p   = 1.0;
      m0_init_S(run,crnt,rho,u,v,w,p);

    }
    else if (CASE==201) { // Model 0, Backward facing step, https://amroc.sourceforge.net/examples/euler/2d/html/bfstep_n.htm
  
      rho = 1.0;
      u   = 0.0;
      v   = 0.0;
      w   = 0.0;
      p   = 1.0;
      m0_init_S(run,crnt,rho,u,v,w,p);

    }
    else if (CASE==202) { // Model 0, Ramp, https://amroc.sourceforge.net/examples/euler/2d/html/ramp_n.htm

      if(xw+yw>0.0) { 
        rho = 1.4;
        u   = 0.0;
        v   = 0.0;
        w   = 0.0;
        p   = 1.0;
        m0_init_S(run,crnt,rho,u,v,w,p);
      } else {
        rho =  8.0;
        //u   =  8.25*cos(PI/6.0);
        //v   = -8.25*sin(PI/6.0);
        w   =  0.0;
        p   =  116.5;
        m0_init_S(run,crnt,rho,u,v,w,p);
      }

    }
    else if (CASE==203) { // Model 1, NS, boundary layer case for viscosity

        rho = 1.0;
        u   = 68;
        v   = 0.0;
        w   = 0.0;
        p   = 101325.0;
        m0_init_S(run,crnt,rho,u,v,w,p);
     
    } 
    else {
      printf("No case initialization found: Exiting \n");
      exit(0);
    }
  
    crnt=crnt->lnxt;
  }

}