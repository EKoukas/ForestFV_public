#include "strdata.h"

struct BRANCH * createbranch (int type, int root, int nkids, int nlnd,int nlvs) { 

  struct BRANCH * brch;
  int i;
  brch=malloc(sizeof(BRANCH));
  brch->merge=0;
  brch->split=0;
  brch->hangn=0;

  if (type==-1){

    brch->root=root;
    brch->type=type;
    brch->nkeens=0;
    brch->nkids=0;
    brch->nlnd=nlnd;
    brch->nlfc=0;
    brch->tot_neigs=0;
  
  } else if (type==0){
    
    brch->root=root;
    brch->type=type;
    brch->nkeens=6;
    brch->nkids=0;
    brch->nlnd=8;
    brch->nlfc=6;
    brch->tot_neigs=6;
  
  } else if (type==2){
  
    brch->root=root;
    brch->type=type;
    brch->nkeens=5;
    brch->nkids=0;
    brch->nlnd=6;
    brch->nlfc=5;
    brch->tot_neigs=5;
  
  } else if (type==1){
  
    brch->root=root;
    brch->type=type;
    brch->nkeens=4;
    brch->nkids=0;
    brch->nlnd=4;
    brch->nlfc=4;
    brch->tot_neigs=4;
  
  }

  brch->keen          = malloc(brch->nkeens*sizeof(BRANCH *));
  brch->keenfc        = malloc(brch->nkeens*sizeof(int));
  brch->keenfcangle   = malloc(brch->nkeens*sizeof(int));

  brch->lfc_neigs     = malloc(brch->nkeens*sizeof(int));
  brch->nsfc          = malloc(brch->nkeens*sizeof(int));
  brch->neigtr        = malloc(brch->nkeens*sizeof(BRANCH **));
  brch->neignd        = malloc(brch->nkeens*sizeof(int *));

  for (i=0; i < brch->nkeens; i++){
    brch->neigtr[i]   = malloc(4*sizeof(BRANCH *)); 
    brch->neignd[i]   = malloc(2*sizeof(int )); 
    brch->nsfc[i]=1;
    brch->lfc_neigs[i]=1;
  }

  brch->ipartbound    = malloc(brch->nkeens*sizeof(int *));

  for (i=0; i < brch->nkeens; i++){
    brch->ipartbound[i]=malloc(4*sizeof(int)); 
  }

  brch->neigfc        = malloc(brch->nkeens*sizeof(int));
  brch->neigag        = malloc(brch->nkeens*sizeof(int));
  brch->adrs          = malloc(nlvs*sizeof(int));

  brch->kids=NULL;
  brch->lnxt=NULL;
  brch->lprv=NULL;
  
  return brch;
  
}