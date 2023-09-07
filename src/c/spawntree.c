#include "strdata.h"

// Tag is the ikid plus 1 (tag cannot contain zeros exept the leading ones, 
// tag is base-10 integer with digits smaller or equal to 8 (1<=d<=8)
// example I
// tag=24221                   // t=24221          // l=0  ->   tag>0
// n[0]=(24221%10)/1=1         // t-=1*1=24220     // l=1  tag>9
// n[1]=24220%100/10=20/10=2   // t-=2*10=24200    // l=2 tag>99
// n[2]=24200%1000/100=2       // t-=2*100=24000   // l=3 tag>999
// n[3]=24000%10000/1000=4     // t-==4*1000=20000 // l=4 tag>9999
// n[4]=20000%100000/10000=2   // t-=2*10000=0     // l=5 tag<99999 stop // Result: l=5  n=1 2 2 4 2 

// example II
// tag=0 // l=0 n[0]=0  thus spawn the tree only to a single branch => nlvl=1

struct BRANCH * spawntree (RUN * run,TREE * tree, int tag, int iel, int icon) {
  
  struct BRANCH * brch;
  struct BRANCH * locl;
  struct BRANCH * last;
  struct BRANCH * crnt;
  struct BRANCH * keen;
  int n[10];
  int l,t,iv;
  int ilvl,nlvl;
  int ifc,ind,iel2;
  int match;
  l=0;
  for (ilvl=0;ilvl<10;ilvl++){
    n[ilvl]=0;
  }
  t=tag;

  while (tag>(pow(10,l)-1)){
    n[l]=t%((int) pow(10,l+1))/((int) pow(10,l));
    t-=n[l] * ((int) pow(10,l));
    l++;
  }
  brch=tree->brch;
  if (brch==NULL){
    //   printf ("in spawn tag %d iel %d rank %d \n",tag,iel,run->con->rank);
    brch=createbranch(run->topo->drys[iel]->type,iel,0,0,2);
    brch->prnt=run->topo->glob; 
    brch->root=iel;
    brch->level=1;
    brch->adrs[0]=-1;
    brch->adrs[1]=iel;
    brch->tag=tagaddress(brch->adrs,brch->level);
    brch->part=run->con->rank;
    brch->nsplit=SPLITMODE;
    if (icon==1){
      crnt=run->topo->locl;
      last=crnt;
      while (crnt!=NULL){
        last=crnt;
        crnt=crnt->lnxt;
      }

      brch->lprv=last;
      if (run->topo->locl!=NULL){
        last->lnxt=brch;
        brch->lnxt=NULL;
      }
      if (run->topo->locl==NULL){
        run->topo->locl=brch;
        brch->lnxt=NULL;
        brch->lprv=last;
      }
    } else {
      brch->lnxt=NULL;
      brch->lprv=NULL;
    }

    // allocate leaf for final only they will be deallocated in split
    cellallocation(brch);
    leafallocation(run,brch);
    createlfc(brch);
    for (iv=0;iv<NEQ;iv++){
	    brch->el->S[iv]=0.0;
    }
    for (ind=0;ind<brch->nlnd;ind++){
      brch->cl->nd[ind].num=run->topo->drys[iel]->vrtx[ind];
    }
    for (ind=0;ind<brch->nlnd;ind++){
      brch->cl->nd[ind].x=run->topo->vertex[run->topo->drys[iel]->vrtx[ind]];
      brch->cl->nd[ind].y=run->topo->vertey[run->topo->drys[iel]->vrtx[ind]];
      brch->cl->nd[ind].z=run->topo->vertez[run->topo->drys[iel]->vrtx[ind]];
    }
    normalvectorcalc(run,brch);
    tagvectorcalc(run,brch);
    volmproperties(brch,run);
    //   printf("spawning null first branch tag %d nlvl %d adrs %d %d \n",tag,nlvl,brch->adrs[0],brch->adrs[1]);
    brch->part=run->topo->drys[iel]->part;
    run->topo->drys[iel]->brch=brch;
  }
  nlvl=l+1;
  ilvl=0;
  //  printf("spawning %d tag %d nlvl %d n %d \n",iel,tag,nlvl,n[ilvl]);
  for (ilvl=0;ilvl<nlvl-1;ilvl++){
    //   printf("spawning %d tag %d nlvl %d n %d \n",iel,tag,nlvl,n[ilvl]);
    if (n[ilvl]>0){
      //     printf("to split %d tag %d nkids %d nlvl %d ilvl %d n %d %d %d \n",iel,tag,brch->nkids,nlvl,ilvl,n[0],n[1],n[2]);
      if (brch->nkids==0){split(icon,run,brch);} 
      brch=brch->kids[n[ilvl]-1];
    }
  }
  return brch;
}
