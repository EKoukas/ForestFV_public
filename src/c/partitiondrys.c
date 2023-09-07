#include "strdata.h"

void partitiondrys (RUN * run) {

int ifc,nlfc;
int iel,iedge,ineig;
int nedge;
int i;

idx_t *optioNS;
idx_t nvert, npart,ncon,edgecut;
idx_t * XADJ;
idx_t * ADJY;
idx_t * VWGT;
idx_t * part;
int   * ipart;
idx_t numflag=0;
idx_t wgtflag=0;
idx_t objval=0;

XADJ=malloc((1+run->topo->nleaves)*sizeof(idx_t));
VWGT=malloc(run->topo->nleaves*sizeof(idx_t));

iedge=0;
for (iel=0;iel<run->topo->nleaves;iel++){
  nlfc=numlfc(run->topo->drys[iel]->type);
  for (ifc=0;ifc<nlfc;ifc++){
    if (run->topo->drys[iel]->conf[ifc]>=0) {iedge++;}
  }
}
nedge=iedge;
ADJY=malloc(nedge*sizeof(idx_t));

iedge=0;
ineig=0;
for (iel=0;iel<run->topo->nleaves;iel++){
  VWGT[iel]=1;
  nlfc=numlfc(run->topo->drys[iel]->type);
  if (iel==0) {XADJ[iel]=0;}
  if (iel!=0) {XADJ[iel]=XADJ[iel-1]+ineig;}
  ineig=0;
  for (ifc=0;ifc<nlfc;ifc++){
    if (run->topo->drys[iel]->conf[ifc]>=0) {
      ADJY[iedge]=run->topo->drys[iel]->conf[ifc];
      iedge++;
      ineig++;
    }
  }
}
XADJ[run->topo->nleaves]=XADJ[run->topo->nleaves-1]+ineig;

//for (i=0;i<run->topo->nleaves;i++){ printf ("%d adjx %d :",i,XADJ[i]);}// for(j=0;j<XADJ[i+1]-XADJ[i];j++){printf (" %d %d ",XADJ[i]+j,ADJY[XADJ[i]+j]) ;} printf ("\n");}
//for (i=0;i<10;i++){ printf ("%d adjx %d : \n",i,XADJ[i]);}// for(j=0;j<XADJ[i+1]-XADJ[i];j++){printf (" %d %d ",XADJ[i]+j,ADJY[XADJ[i]+j]) ;} printf ("\n");}

nvert=run->topo->nleaves;
npart=run->con->size;

ncon=1;
part=malloc(nvert*sizeof(idx_t));
ipart=malloc(nvert*sizeof(int));
METIS_PartGraphRecursive(&nvert,&ncon,XADJ,ADJY,VWGT,NULL,NULL,&npart,NULL,NULL,NULL,&objval,part);


for (iel=0;iel<run->topo->nleaves;iel++){
  run->topo->drys[iel]->brch=NULL;
  run->topo->drys[iel]->part=part[iel];
  // printf(" partition drys %d %d %d \n",iel,part[iel],run->topo->drys[iel]->part);
  if (part[iel]==run->con->rank){
    //   printf(" partition drys made it in %d %d %d %d \n",iel,part[iel],run->topo->drys[iel]->part,run->con->rank);
    run->topo->drys[iel]->brch=createbranch(run->topo->drys[iel]->type,iel,0,0,2);
    leafallocation(run,run->topo->drys[iel]->brch);
    cellallocation(run->topo->drys[iel]->brch);
    createlfc(run->topo->drys[iel]->brch);
    run->topo->drys[iel]->brch->part=run->con->rank;
  }
  MPI_Barrier(MPI_COMM_WORLD);
}


free(VWGT);
free(XADJ);
free(ADJY);
free(part);
}
