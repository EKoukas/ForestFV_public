//-------------------------------------------------------
//
// Read data from .nd (nodes) and .el (elements) and creates topology
//
//-------------------------------------------------------

#include "strdata.h"

void createtopo (struct RUN * run) {
  
  int type,typg,nlfc,nlnd;
  int iel,ielp1;
  int ind;
  int node;
  double ndx,ndy,ndz;
  FILE * el;
  FILE * nd;

  struct BRANCH * crnt;
  struct BRANCH * rcnt;


  el=fopen(run->con->fileel,"r");
  nd=fopen(run->con->filend,"r");

  fscanf(el,"%d\n",&run->topo->ntrees);
  fscanf(nd,"%d\n",&run->topo->ngnd);
  
  run->topo->ndepth=1;
  run->topo->nleaves=run->topo->ntrees;

  run->topo->drys=malloc(run->topo->nleaves*sizeof(struct TREE *));
  // Local leaves number per partition
  run->topo->lleaves=malloc(run->con->size*sizeof(int));
  for (iel=0;iel<run->topo->ntrees;iel++){              // Loop elements
    run->topo->drys[iel]=malloc(sizeof(struct TREE));

    if (fscanf(el,"%d\t%d",&ielp1,&typg) != 2) {
      printf("Problem in createtopo: exiting \n");
      exit(0);
    }

    if      (typg==4) {type=0;nlnd=8;nlfc=6;} // hex 
    else if (typg==5) {type=2;nlnd=6;nlfc=5;} // pri
    else if (typg==6) {type=1;nlnd=4;nlfc=4;} // tet
    run->topo->drys[iel]->type=type;
    run->topo->drys[iel]->conf=malloc(nlfc*sizeof(int));
    run->topo->drys[iel]->vrtx=malloc(nlnd*sizeof(int));
    run->topo->drys[iel]->brch=NULL;

    for (ind=0;ind<nlnd;ind++){
      fscanf(el,"\t%d",&node);
      run->topo->drys[iel]->vrtx[ind]=node-1;
    }
  }
  run->topo->vertex=malloc(run->topo->ngnd*sizeof(double ));
  run->topo->vertey=malloc(run->topo->ngnd*sizeof(double ));
  run->topo->vertez=malloc(run->topo->ngnd*sizeof(double ));

  for (ind=0;(ind<run->topo->ngnd);ind++){ // Loop nodes
    fscanf(nd,"%d\t%le\t%le\t%le\n ",&node,&ndx,&ndy,&ndz);
    run->topo->vertex[ind]=ndx*MSCALE;
    run->topo->vertey[ind]=ndy*MSCALE;
    run->topo->vertez[ind]=ndz*MSCALE;
  }
  fclose(el);
  fclose(nd);


  run->topo->glob=createbranch(-1,-1,0,0,1);

  run->topo->glob->adrs[0]=0;
  run->topo->glob->nkids=0;
  run->topo->glob->prnt=run->topo->glob;
  run->topo->glob->root=-1;
  run->topo->glob->keen=NULL;
  run->topo->glob->level=0;
  run->topo->glob->type=-1;
  run->topo->glob->nsplit=0;
  run->topo->glob->hangn=0;

}