#include "strdata.h"

void rebuildlist (RUN * run) {

  int iel,ilvl;
  int istart;
  struct TREE * tree;
  struct BRANCH * brch;
  struct BRANCH * crnt;
  struct BRANCH * prvs;
  int indx[10];
  
  // Build lists
  
  prvs=NULL;
  istart=0;
  for (iel=0;iel<run->topo->ntrees;iel++) { // Loop elements
    tree=run->topo->drys[iel];
    brch=tree->brch;
    
    if (tree->part==run->con->rank){
      
      // zero indx
      for (ilvl=0;ilvl<10;ilvl++){
	      indx[ilvl]=0;
      }
      crnt=brch;
      ilvl=0;
      if (brch->nkids==0){

	      if (istart==1){
	        brch->lnxt=NULL;
     	    brch->lprv=prvs;
	        prvs->lnxt=brch;
	        prvs=brch;
	      }
	      if (istart==0){
          brch->lnxt=NULL;
       	  brch->lprv=prvs;
          prvs=brch;
	        istart=1;
	        run->topo->locl=brch; 
	      }

      } else {

        while (crnt->root!=-1){   // decent
          ilvl=crnt->level;
          if (crnt->nkids!=0){
            if (indx[ilvl]<crnt->nkids){  // still looping at ilvl
              crnt=crnt->kids[indx[ilvl]];
              indx[ilvl]++;
            } else {                         // has kids but finished looping at ilvl
              crnt=crnt->prnt;
              indx[ilvl]=0;               // go up and reset counter for ilvl
            }
          } else {    	// hit bottom	  
            if (istart==1){
              crnt->lnxt=NULL;
              crnt->lprv=prvs;
              prvs->lnxt=crnt;
              prvs=crnt;
            }
            if (istart==0){
              run->topo->locl=crnt; 
              crnt->lnxt=NULL;
              crnt->lprv=prvs;
              prvs=crnt;
              istart=1;
            }
            crnt=crnt->prnt;
          }
        }

      }

    }
  }  

}