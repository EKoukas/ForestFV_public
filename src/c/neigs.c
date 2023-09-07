#include "strdata.h"

void neigs (RUN * run) {

  int ifc;
  struct BRANCH * crnt;

  crnt=run->topo->locl; 
  while (crnt!=NULL){
    for (ifc=0;ifc<crnt->nlfc;ifc++){
      elconn(run,crnt,ifc);
      printf(" neigs %d %d %d ",crnt->root,run->con->rank,crnt->neigtr[ifc][0]->root);
      if (crnt->neigtr[ifc][0]->root!=-1){
	      crnt->cl->fc[ifc].bc=0;
      }
    }
    crnt=crnt->lnxt;
  }

}