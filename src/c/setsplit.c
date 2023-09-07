#include "strdata.h"

void setsplit(RUN *run){
struct BRANCH * crnt;

  crnt=run->topo->locl; 
  while (crnt!=NULL){
    setsplitcalc(crnt,run);
    crnt=crnt->lnxt;
  }
}


void setsplitcalc(BRANCH * crnt,RUN *run){
  
  int ifc,fcqs;
  double maxprod,prod;

  maxprod=-10000.0;
  if (SPLITMODE==4){
    for (ifc=0;ifc<crnt->nlfc;ifc++){
      prod=crnt->cl->nx[ifc]*NXQUAD+crnt->cl->ny[ifc]*NYQUAD+crnt->cl->nz[ifc]*NZQUAD;
      if(maxprod<prod){
        crnt->fcqsplit=ifc;
        maxprod=prod;
      }
    }
    crnt->nsplit=4;
    if (crnt->type==2&&crnt->cl->fc[crnt->fcqsplit].type==0){
      crnt->nsplit=2;
    }
  } else if (SPLITMODE==8){
    crnt->nsplit=8;
    crnt->fcqsplit=-1;
  }
  
}