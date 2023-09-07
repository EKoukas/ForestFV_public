#include "strdata.h"

void spawntopo (RUN * run) {
  
  int type,typg,nlfc,nlnd;
  int iel,ielp1,istart;
  int ind;
  int node;

  struct BRANCH * crnt;
  struct BRANCH * rcnt;


  istart=0;
  crnt=run->topo->glob;
  if((VERBOSE==1)&&(run->con->rank==0)){ printf ("ForestFV: MAIN: TREES SPAWNED TO BRANCHES 0 \n");}
  for (iel=0;iel<run->topo->nleaves;iel++){ // Loop elements
    if (run->topo->drys[iel]->part==run->con->rank){
      rcnt=run->topo->drys[iel]->brch;	  
      if (istart==0){
        rcnt->lnxt=NULL;
        rcnt->lprv=run->topo->glob;
	      run->topo->locl=rcnt;
	      istart=1;
      } else {
        rcnt->lprv=crnt;
        crnt->lnxt=rcnt;
        rcnt->lnxt=NULL;
      }
      crnt=rcnt;
      crnt->prnt=run->topo->glob;  //sure?
      crnt->root=iel;
      crnt->level=1;
      crnt->adrs[0]=0;
      crnt->adrs[1]=iel;
      crnt->tag=tagaddress (crnt->adrs,crnt->level);

      for (ind=0;ind<crnt->nlnd;ind++){
        crnt->cl->nd[ind].num=run->topo->drys[iel]->vrtx[ind];
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if((VERBOSE==1)&&(run->con->rank==0)){ printf ("ForestFV: MAIN: TREES SPAWNED TO BRANCHES I \n");}
  for (iel=0;iel<run->topo->nleaves;iel++){ // Loop elements
    if (run->topo->drys[iel]->part==run->con->rank){
      crnt=run->topo->drys[iel]->brch;
      for (ind=0;ind<crnt->nlnd;ind++){
        crnt->cl->nd[ind].x=run->topo->vertex[run->topo->drys[iel]->vrtx[ind]];
        crnt->cl->nd[ind].y=run->topo->vertey[run->topo->drys[iel]->vrtx[ind]];
        crnt->cl->nd[ind].z=run->topo->vertez[run->topo->drys[iel]->vrtx[ind]];
      }
    }
    //createlfc(crnt);
  }
	
  MPI_Barrier(MPI_COMM_WORLD);
  if((VERBOSE==1)&&(run->con->rank==0)){ printf ("ForestFV: MAIN: TREES SPAWNED TO BRANCHES II \n");}
}
