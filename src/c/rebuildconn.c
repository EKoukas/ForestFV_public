#include "strdata.h"

void rebuildconn(RUN * run) {

  struct BRANCH * crnt;
  struct BRANCH * buff;
  int iel,ifc,tag2,iel2,ipts,ibuf; 
  int ilvl,nlvl;
  int * adr;
  adr=malloc(10*sizeof(int));


  // associate neigbours
  crnt=run->topo->locl;
  while (crnt!=NULL){
    iel=crnt->root;
    for (ifc=0;ifc<crnt->nlfc;ifc++){
      elconn(run,crnt,ifc);
      elconnfin(crnt,ifc);
      crnt->nsfc[ifc]=max(1,crnt->lfc_neigs[ifc]);
    }
    crnt=crnt->lnxt;
  }

  // assign connectivity from buffer to prox
  for (ipts=0;ipts<run->con->size;ipts++){
    for (ibuf=0;ibuf<run->topo->nbuff[ipts];ibuf++){
      iel2=run->topo->buffdrys[ipts][ibuf];
      tag2=run->topo->bufftags[ipts][ibuf];

      tag2adr(adr,tag2,&nlvl);
      buff=run->topo->drys[iel2]->brch;
      ilvl=0;
      while (buff->nkids!=0){
        buff=buff->kids[adr[ilvl]-1];ilvl++;
      }
      ifc=run->topo->bufflfc1[ipts][ibuf];
      elconn(run,buff,ifc);
      elconnfin(buff,ifc);
      buff->nsfc[ifc]=max(1,buff->lfc_neigs[ifc]);
    }
  }
  free(adr);

}