#include "strdata.h"

void m5_facevalues_MUSCL(struct RUN *run,struct BRANCH * crnt,int ifc,int iside,int bc_temp, double * S_temp) {

  int iv, ineig, Nneig, ifc_opp;
  double avf[10];
  struct BRANCH * oppo;

  if (iside==0) {      //  run->AUX_S from current element   
    for (iv = 0; iv < NEQ; iv++) {
      S_temp[iv] = crnt->el->S[iv];
    }
  } else if (iside == 1) {     
    if (bc_temp == 0) { 
      if (crnt->lfc_neigs[ifc]==1) {        // 2:1 or 1:1 connectivity
        oppo=crnt->neigtr[ifc][0];
        ifc_opp=crnt->neigfc[ifc];
        for (iv = 0; iv < NEQ; iv++) { 
          S_temp[iv] = oppo->el->S[iv];
        }
      } else { // 1:2
        Nneig=crnt->lfc_neigs[ifc];  
        for (iv=0;iv<NEQ;iv++) {
          S_temp[iv]=0.0;
        }       
        for (ineig = 0; ineig < Nneig; ineig++) {
          oppo = crnt->neigtr[ifc][ineig];
          for (iv = 0; iv < NEQ; iv++) {
            S_temp[iv] += oppo->el->S[iv] / ((double)Nneig);
          }
        }
      }    
    } else { // BC found
      ineig=0;
      m5_facevalues(run, crnt, ineig, ifc ,bc_temp, 1);
      for (iv = 0; iv < NEQ; iv++) {
        S_temp[iv] = run->sol_R->vec[iv];
      }
    } 
  }

  for (iv=0;iv<eqtypn[3];++iv) { avf[iv] = S_temp[iv + eqtypi[3]]; }
  avf_limit(avf, eqtypn[3], 100.0);
  for (iv=0;iv<eqtypn[3];++iv) { S_temp[iv+eqtypi[3]] = avf[iv]; }

}