#include "strdata.h"

void face_interpolation_WG (struct RUN* run) {

  int f,iv;
  double dx,dy,dz;
  struct BRANCH * crnt;

  crnt=run->topo->locl;     
  while (crnt!=NULL){

    for(f=0; f<(crnt->nlfc); ++f) {

      //dx =
      //dy =
      //dz =

      for(iv=0;iv<NEQ;iv++){
        crnt->el->SF[f][iv] = crnt->el->WG[iv][0]*dx + crnt->el->WG[iv][1]*dy + crnt->el->WG[iv][2]*dz;
      }


    }
      
    crnt=crnt->lnxt;
  }  

}