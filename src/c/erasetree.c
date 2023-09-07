#include "strdata.h"

void erasetree (struct RUN * run,struct TREE * tree, int iel) {

  int ind,ifc,ilvl,ikid,term;
  int indx[10];
  struct BRANCH * crnt;
  struct BRANCH * prnt;
  struct BRANCH * brch;
  struct BRANCH * frst;
  struct BRANCH * last;

  for (ilvl=0;ilvl<10;ilvl++){
    indx[ilvl]=0;
  }
  ilvl=1;
  brch=tree->brch;
  
  // Restore lists
  if (brch!=NULL){
    if (brch->hangn==0){
    crnt=brch;
    frst=crnt;
    last=crnt;
    while (crnt->nkids!=0){
      crnt=crnt->kids[0];
    }
    frst=crnt; //bottom
    while(crnt!=NULL&&crnt->root==iel){
    last=crnt; //last
    crnt=crnt->lnxt;}
     
    if (run->topo->locl==frst){
      run->topo->locl=last->lnxt;
      if (last->lnxt!=NULL){last->lnxt->lprv=NULL;}
    }
    else{
      frst->lprv->lnxt=last->lnxt;
 //     if (last->lprv!=NULL){frst->lprv->lnxt=last->lnxt;}
      if (last->lnxt!=NULL){last->lnxt->lprv=frst->lprv;}
    }
    }
    // zero indx
    for (ilvl=0;ilvl<10;ilvl++){
      indx[ilvl]=0;
    }
    if (brch->nkids==0){
    }
    ilvl=1;
    term=0;
    crnt=brch;
    if (brch->nkids>0){
 //    printf("erasing split  %d %d %d %d \n",iel,ilvl,brch->nkids,ilvl);
  //     while (indx[1]<=brch->nkids&&ilvl>1){
       while (indx[1]<=brch->nkids&&term==0){
      // decent
        ilvl=crnt->level;
        if (crnt->nkids!=0){
          if (indx[ilvl]<crnt->nkids){
            crnt=crnt->kids[indx[ilvl]];
            indx[ilvl]++;
          }
          else{
	  //  printf("erasing a %d %d %d \n",crnt->level,crnt->root,ilvl);
            prnt=crnt->prnt;

            if (crnt!=NULL){
		    if (crnt->el!=NULL){leafdeallocation(run,crnt);}
	    celldeallocation(crnt);
            destroybrch(crnt);
	    }
	    if(ilvl==1){term=1;}
            indx[ilvl]=0;
  	  crnt=prnt;
          }
        }
        else{
          prnt=crnt->prnt;
  	//if (indx[ilvl]==crnt->nkids-1){prnt->nkids=0;}
	//  printf("erasing b %d %d \n",crnt->level,crnt->root);
	  leafdeallocation(run,crnt);
	  celldeallocation(crnt);
          destroybrch(crnt);
          crnt=prnt;
        }
      }
    } 
    else{leafdeallocation(run,brch);
    celldeallocation(brch);
    destroybrch(brch);
    }
    tree->brch=NULL;
	//  printf("erasing c %d \n",iel);
  }
}
