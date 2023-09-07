#include "strdata.h"

void m0_facevalues_MUSCL(struct RUN *run,struct BRANCH * crnt,int ifc,int iside,int bc_temp, double * S_temp) {

  int iv,ing,ifc_opp;
  int ineig,Nneig;
  struct BRANCH * oppo;

  if (iside==0) {     

    for(iv=0;iv<NEQ;iv++){
      S_temp[iv] = crnt->el->S[iv];
    }
    
  } else if ((iside==1)&&(bc_temp==0)) {     // Neigbor

    if (crnt->lfc_neigs[ifc]==1) {        // 2:1 or 1:1 connectivity
     
      ing   = 0;   
      oppo=crnt->neigtr[ifc][ing];
      ifc_opp=crnt->neigfc[ifc];
      for(iv=0;iv<NEQ;iv++) {  
        S_temp[iv] = oppo->el->S[iv];
      }
    }

    if (crnt->lfc_neigs[ifc]!=1) {          // 1:2

      for (iv=0;iv<NEQ;iv++) {
        S_temp[iv]=0.0;
      }

      Nneig=crnt->lfc_neigs[ifc];         
      
      for (ineig=0;ineig<Nneig;ineig++) {
        ing  = ineig;
        oppo = crnt->neigtr[ifc][ing];
        for(iv=0;iv<NEQ;iv++){  
          S_temp[iv]  += oppo->el->S[iv]/( (double) Nneig);
        }
      }

    }    
    
  } else if ((iside==1)&&(bc_temp!=0)) { // BC found

    ineig=0;
    m0_facevalues(run,crnt,ineig,ifc,bc_temp,1);

    for(iv=0;iv<NEQ;iv++){
      S_temp[iv] = run->sol_R->vec[iv];
    }

  }
    
}