#include "strdata.h"

void merge (int linked, RUN* run, BRANCH * brch) {

  int ikid,iv;
  
  brch->part=brch->kids[0]->part;
  if (linked==1){
    if (run->topo->locl==brch->kids[0]){
      run->topo->locl = brch;
      brch->lnxt = brch->kids[brch->nkids-1]->lnxt;
      brch->lprv = brch->kids[0]->lprv;
    } else {
      brch->lnxt=brch->kids[brch->nkids-1]->lnxt;
      brch->lprv=brch->kids[0]->lprv;
    }
  }
  
  leafallocation(run,brch);
  brch->split=0;
  brch->el->crith=0.0;
  for (ikid=0;ikid<brch->nsplit;ikid++){
    brch->split=max(brch->split,brch->kids[ikid]->split);
    brch->el->crith=max(brch->el->crith,brch->kids[ikid]->el->crith);
  }
  normalvectorcalc(run,brch);
  tagvectorcalc(run,brch);
  volmproperties(brch,run);
   
  run->topo->nleavestmp++;
  run->topo->lleaves[brch->part]-=1;

  for (iv=0;iv<2;iv++) { brch->el->SI[iv]=0.0; }

  for(iv=0;iv<NEQ;iv++){
    brch->el->S[iv]=0.0;
    
    //brch->el->SG[iv][0]=0.0;
    //brch->el->SG[iv][1]=0.0;
    //brch->el->SG[iv][2]=0.0;

    for (ikid=0;ikid<brch->nsplit;ikid++){
      //brch->el->S[iv]+=brch->kids[ikid]->el->S[iv]*brch->kids[ikid]->cl->Vol/brch->cl->Vol;
      brch->el->S[iv]+=brch->kids[ikid]->el->S[iv]/((double) brch->nsplit);

      //brch->el->SG[iv][0]+=brch->kids[ikid]->el->SG[iv][0]*brch->kids[ikid]->cl->Vol/brch->cl->Vol;
      //brch->el->SG[iv][1]+=brch->kids[ikid]->el->SG[iv][1]*brch->kids[ikid]->cl->Vol/brch->cl->Vol;
      //brch->el->SG[iv][2]+=brch->kids[ikid]->el->SG[iv][2]*brch->kids[ikid]->cl->Vol/brch->cl->Vol;
    }
  }

  if ((MODEL==4)||(MODEL==5)){
    for(iv=0;iv<NEQ;++iv){ run->Vec_temp_R[iv] = brch->el->S[iv]; }
    if      (MODEL==4){ m4_relaxationvec(run,brch); }
    else if (MODEL==5){ m5_relaxationvec(run,brch); }
    for(iv=0;iv<NEQ;++iv){ brch->el->S[iv] = run->Vec_temp_R[iv]; }
  }


  if (linked==1){
    if (brch->lnxt!=NULL) {brch->lnxt->lprv=brch;}
    if (brch->lprv!=NULL) {brch->lprv->lnxt=brch;}
  }
  for (ikid=0;ikid<brch->nsplit;ikid++){
    if (brch->kids[ikid]->nkids!=0){printf ("Merge: has kids \n");}

    leafdeallocation(run,brch->kids[ikid]);
    celldeallocation(brch->kids[ikid]);
    destroybrch(brch->kids[ikid]);
  } 
  free(brch->kids);
  brch->kids=NULL;
  brch->nkids=0;

}