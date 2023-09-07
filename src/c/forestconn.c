#include "strdata.h"

void forestconn (RUN* run) {
  
  int ignd;
  int ilnd;
  int nnds;
  int iel,iel1,iel2,ifc,ind,ifc2,iang,indarranged,indr;
  static int match;
  struct BRANCH *** NDCON;
  struct BRANCH * crnt;
  struct BRANCH * keen;
  int *  NDELS;
  int *  ELLIST;

  for (iel=0;iel<run->topo->nleaves;iel++){ // Loop elements
    if (run->topo->drys[iel]->part==run->con->rank){
      for (ifc=0;ifc<run->topo->drys[iel]->brch->nlfc;ifc++){
	      if (run->topo->drys[iel]->conf[ifc]!=-1){
	        if (run->topo->drys[run->topo->drys[iel]->conf[ifc]]->part!=run->con->rank){
            keen=createbranch(run->topo->drys[run->topo->drys[iel]->conf[ifc]]->type,run->topo->drys[iel]->conf[ifc],0,0,2);
            keen->lnxt=NULL;
            keen->lprv=NULL;
            keen->prnt=run->topo->glob;
            keen->root=run->topo->drys[iel]->conf[ifc];
            keen->level=1;
            keen->adrs[0]=-1;
            keen->adrs[1]=run->topo->drys[iel]->conf[ifc];
            keen->tag=tagaddress(keen->adrs,1);
            keen->hangn=1;

            leafallocation(run,keen);
            cellallocation(keen);
            keen->part=run->con->rank;
            for (ind=0;ind<keen->nlnd;ind++){
              keen->cl->nd[ind].num=run->topo->drys[keen->root]->vrtx[ind];
              keen->cl->nd[ind].x=run->topo->vertex[run->topo->drys[keen->root]->vrtx[ind]];
              keen->cl->nd[ind].y=run->topo->vertey[run->topo->drys[keen->root]->vrtx[ind]];
              keen->cl->nd[ind].z=run->topo->vertez[run->topo->drys[keen->root]->vrtx[ind]];
	          }
            createlfc(keen);
            run->topo->drys[iel]->brch->keen[ifc]=keen;
   
	        } else {
  	        run->topo->drys[iel]->brch->keen[ifc]=run->topo->drys[run->topo->drys[iel]->conf[ifc]]->brch;
	        }
	      } else {
          run->topo->drys[iel]->brch->keen[ifc]=run->topo->glob;
        }
      }
    }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  crnt=run->topo->locl; 
  while (crnt!=NULL){
    for (ifc=0;ifc<crnt->nlfc;ifc++){
      nnds=crnt->cl->fc[ifc].nnds;
      if (crnt->keen[ifc]->root!=-1){
        for(ifc2=0;ifc2<crnt->keen[ifc]->nlfc;ifc2++){
          if (crnt->keen[ifc]->cl->fc[ifc2].type==crnt->cl->fc[ifc].type){
            for(iang=0;iang<crnt->keen[ifc]->cl->fc[ifc2].nnds;iang++){
              match=0;
              for(ind=0;ind<crnt->keen[ifc]->cl->fc[ifc2].nnds;ind++){
                indr=fcndrot(ind,iang,crnt->keen[ifc]->cl->fc[ifc2].type);

                if(crnt->cl->nd[fcnd2elnd(ifc,ind,crnt->type)].num==crnt->keen[ifc]->cl->nd[fcrefnd2elnd(ifc2,indr,crnt->keen[ifc]->type)].num){match++;}
              }
              if (match==nnds){ 
                crnt->keenfc[ifc]=ifc2;
                crnt->keenfcangle[ifc]=iang;
              }
            }
          }
        }
      } else {
        crnt->keenfc[ifc]=0;
        crnt->keenfcangle[ifc]=0;
      }
    }
    crnt=crnt->lnxt;
  }


}