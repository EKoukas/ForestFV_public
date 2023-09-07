#include "strdata.h"

void repartition(RUN * run) {

  int ifc,ifc2,iv,nlfc;
  int iel,iedge,ineig,iel2;
  int nedge;
  int match,iang,ind,indr,found;
  int i,istart;
  struct TREE * tree;
  struct BRANCH * keen;
  struct BRANCH * brch;
  struct BRANCH * crnt;
  struct BRANCH * prvs;
  struct BRANCH * prnt;
  int nleaves;
  int ileaves;
  int ilvl,ikid;
  int * tags;
  double * buff;
  int ierr;
  int indx[10];
  MPI_Status status;
  idx_t *optioNS;
  idx_t nvert, npart,ncon,edgecut;
  idx_t * XADJ;
  idx_t * ADJY;
  idx_t * VWGT;
  idx_t * VWGTTOT;
  idx_t * part;
  long * nlvs;
  long * nlvstot;
  int * oprt;
  int * ipart;
  idx_t numflag=0;
  idx_t wgtflag=0;
  idx_t objval=0;
  
  nvert=run->topo->ntrees;
  npart=run->con->size;
  ncon=1;
  part=malloc(nvert*sizeof(idx_t));
  oprt=malloc(nvert*sizeof(int));

  ipart=malloc(nvert*sizeof(int));
  XADJ=malloc((1+run->topo->ntrees)*sizeof(idx_t));
  VWGT=malloc(run->topo->ntrees*sizeof(idx_t));
  VWGTTOT=malloc(run->topo->ntrees*sizeof(idx_t));
  nlvs=malloc(run->topo->ntrees*sizeof(long));
  nlvstot=malloc(run->topo->ntrees*sizeof(long));
  
  iedge=0;
  for (iel=0;iel<run->topo->ntrees;iel++){
    nlfc=numlfc(run->topo->drys[iel]->type);
    for (ifc=0;ifc<nlfc;ifc++){
       if (run->topo->drys[iel]->conf[ifc]>=0) {iedge++;}
    }
  }
  nedge=iedge;
  ADJY=malloc(nedge*sizeof(idx_t));
  
  iedge=0;
  ineig=0;
  for (iel=0;iel<run->topo->ntrees;iel++){
    VWGT[iel]=0;
    MPI_Barrier(MPI_COMM_WORLD);
    oprt[iel]=run->topo->drys[iel]->part;
    tree=run->topo->drys[iel];
    brch=tree->brch;
    if (run->topo->drys[iel]->part==run->con->rank){
      oprt[iel]=run->con->rank;
      while (brch->nkids!=0){
        brch=brch->kids[0];
      }
      while(brch!=NULL&&brch->root==iel){
      VWGT[iel]++;
      brch=brch->lnxt;}
    }
    nlvs[iel]=VWGT[iel];
  }  
  MPI_Allreduce(nlvs,nlvstot,run->topo->ntrees,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
  
  
  for (iel=0;iel<run->topo->ntrees;iel++){
    VWGT[iel]=nlvstot[iel];
    nlvs[iel]=nlvstot[iel];
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
  XADJ[run->topo->ntrees]=XADJ[run->topo->ntrees-1]+ineig;
  
  METIS_PartGraphRecursive(&nvert,&ncon,XADJ,ADJY,VWGT,NULL,NULL,&npart,NULL,NULL,NULL,&objval,part);
  
  for (iel=0;iel<run->topo->ntrees;iel++){
    ipart[iel]=part[iel];
  }


// partition changes, trees must be erased and branches are migrated to new trees that are spawn
// hanging kiNS are created, or kiNS are assigned in new trees
// hanging kiNS are destroyed and assigned trees are deaccosiated in old trees
// for non-affected trees kiNS can change from assigned to hanging and vis versa


// Destroy all hanging keeNS

  for (iel=0;iel<run->topo->ntrees;iel++){
    tree=run->topo->drys[iel];
    brch=tree->brch;
    ileaves=0;
    tags=malloc((nlvs[iel])*sizeof(int));
    buff=malloc(NEQ*(nlvs[iel])*sizeof(double));
    run->topo->drys[iel]->part=ipart[iel];
    if (oprt[iel]==run->con->rank){   // In in old decomposition destroy hanging branches (essentially hanging kiNS)
      for (ifc=0;ifc<numlfc(run->topo->drys[iel]->type);ifc++){
        if (brch->keen!=NULL&&brch->keen[ifc]->hangn==1){
	  leafdeallocation(run,brch->keen[ifc]);
	  celldeallocation(brch->keen[ifc]);
          destroybrch(brch->keen[ifc]);
        }
      }
    }
  }

// erase buffer trees (neigs)
  for (iel=0;iel<run->topo->ntrees;iel++){
    tree=run->topo->drys[iel];
    if (oprt[iel]!=run->con->rank){
      erasetree(run,tree,iel);
    }
  }


// Spawn trees to level one only
// Do not erase any trees

  for (iel=0;iel<run->topo->ntrees;iel++){
    tree=run->topo->drys[iel];
    brch=tree->brch;
    ileaves=0;
    tags=malloc((nlvs[iel])*sizeof(int));
    buff=malloc(NEQ*(nlvs[iel])*sizeof(double));
    run->topo->drys[iel]->part=ipart[iel];


    if (oprt[iel]==run->con->rank&&ipart[iel]==run->con->rank){  } // No change
    MPI_Barrier(MPI_COMM_WORLD);
    if (oprt[iel]==run->con->rank&&ipart[iel]!=run->con->rank){ // Tree lost
      while (brch->nkids!=0){
        brch=brch->kids[0];
      }
      while(brch!=NULL&&brch->root==iel){
        brch->tag=tagaddress(brch->adrs,brch->level);
        tags[ileaves]=brch->tag;
        for (iv=0;iv<NEQ;iv++){
          buff[iv+NEQ*ileaves]=brch->el->S[iv];
        }
        ileaves++;
      brch=brch->lnxt;}
      nleaves=ileaves;
      ierr=MPI_Send(&nleaves,1,MPI_INT,ipart[iel],0+3*iel,MPI_COMM_WORLD);
      if (nleaves==1){
        ierr=MPI_Send(tags,nleaves,MPI_INT,ipart[iel],1+3*iel,MPI_COMM_WORLD);
        ierr=MPI_Send(buff,NEQ*nleaves,MPI_DOUBLE,ipart[iel],2+3*iel,MPI_COMM_WORLD);
        //erasetree(run,tree,iel);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (oprt[iel]!=run->con->rank&&ipart[iel]==run->con->rank){ // Tree gained
      ierr=MPI_Recv(&nleaves,1,MPI_INT,oprt[iel],0+3*iel,MPI_COMM_WORLD,&status);
      brch=spawntree(run,tree,0,iel,1);
      for (ifc=0;ifc<brch->nlfc;ifc++){
	      if (tree->conf[ifc]==-1){brch->cl->fc[ifc].bc=1;}
	      if (tree->conf[ifc]!=-1){brch->cl->fc[ifc].bc=0;}
      }
      //leafallocation(run,brch);
      brch->el->crith=0.0;
      brch->split=0;
      setsplitcalc(brch,run);
      if (nleaves==1){
        ierr=MPI_Recv(tags,nleaves,MPI_INT,oprt[iel],1+3*iel,MPI_COMM_WORLD,&status);
        ierr=MPI_Recv(buff,NEQ*nleaves,MPI_DOUBLE,oprt[iel],2+3*iel,MPI_COMM_WORLD,&status);
        for (ileaves=0;ileaves<nleaves;ileaves++){
          for (iv=0;iv<NEQ;iv++){
            brch->el->S[iv]=buff[iv+NEQ*ileaves];
          }
        }
      }
    }
    free(tags);
    free(buff);
    MPI_Barrier(MPI_COMM_WORLD);
  }

 // Spawn keeNS

  for (iel=0;iel<run->topo->ntrees;iel++){
    tree=run->topo->drys[iel];
    brch=tree->brch;
    ileaves=0;
    if(ipart[iel]==run->con->rank){
      for (ifc=0;ifc<numlfc(run->topo->drys[iel]->type);ifc++){

	  found=0;
        iel2=run->topo->drys[iel]->conf[ifc];
        if (iel2==-1){
          brch->keen[ifc]=run->topo->glob;
          brch->keenfc[ifc]=0;
          brch->keenfcangle[ifc]=0;
	  found=3;
        }
        else{
          if(ipart[iel2]==run->con->rank){
            brch->keen[ifc]=run->topo->drys[iel2]->brch;
	    found=2;
          }
          else{
            keen=createbranch(run->topo->drys[iel2]->type,iel2,0,0,2);
            keen->prnt=run->topo->glob;
            keen->root=iel2;
            keen->level=1;
            keen->adrs[0]=-1;
            keen->adrs[1]=iel2;
            keen->tag=tagaddress(keen->adrs,brch->level);
            keen->nsplit=SPLITMODE;
	    keen->hangn=1;
            brch->keen[ifc]=keen;

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
            brch->keen[ifc]=keen;
	    

	        /*
    for (ifc=0;ifc<numlfc(run->topo->drys[iel]->type);ifc++){
      iel2=run->topo->drys[iel]->conf[ifc];
      printf("iel iel2, %d %d \n",iel,iel2);
      if (iel2==-1){
        brch->keen[ifc]=run->topo->glob;
      }
      else{
        if(run->topo->drys[iel2]->part==run->con->rank){
          brch->keen[ifc]=run->topo->drys[iel2]->brch;
          for(ifc2=0;ifc2<brch->keen[ifc]->nlfc;ifc2++){
            if (brch->keen[ifc]->cl->fc[ifc2].type==brch->cl->fc[ifc].type){
              for(iang=0;iang<brch->keen[ifc]->cl->fc[ifc2].nnds;iang++){
                match=0;
                for(ind=0;ind<brch->keen[ifc]->cl->fc[ifc2].nnds;ind++){
                  indr=fcndrot(ind,iang,brch->keen[ifc]->cl->fc[ifc2].type);
                  if(brch->cl->nd[fcnd2elnd(ifc,ind,brch->type)].num==brch->keen[ifc]->cl->nd[fcrefnd2elnd(ifc2,indr,brch->keen[ifc]->type)].num){match++;}
                }
                if (match==brch->cl->fc[ifc].nnds){
                  brch->keenfc[ifc]=ifc2;
                  brch->keenfcangle[ifc]=iang;
                }
              }
            }
          }
        }
        else{
          keen=createbranch(run->topo->drys[iel2]->type,iel2,0,0,2);
          keen->prnt=run->topo->glob;
          keen->root=iel2;
          keen->level=1;
          keen->adrs[0]=-1;
          keen->adrs[1]=iel2;
          keen->tag=tagaddress(keen->adrs,brch->level);
          keen->num=iel2;
          keen->nsplit=SPLITMODE;
          brch->keen[ifc]=keen;
        }
      }
    }
    */


          }
          for(ifc2=0;ifc2<brch->keen[ifc]->nlfc;ifc2++){
            if (brch->keen[ifc]->cl->fc[ifc2].type==brch->cl->fc[ifc].type){
              for(iang=0;iang<brch->keen[ifc]->cl->fc[ifc2].nnds;iang++){
                match=0;
                for(ind=0;ind<brch->keen[ifc]->cl->fc[ifc2].nnds;ind++){
                  indr=fcndrot(ind,iang,brch->keen[ifc]->cl->fc[ifc2].type);
                  if(brch->cl->nd[fcnd2elnd(ifc,ind,brch->type)].num==brch->keen[ifc]->cl->nd[fcrefnd2elnd(ifc2,indr,brch->keen[ifc]->type)].num){match++;}
                }
                if (match==brch->cl->fc[ifc].nnds){
                  brch->keenfc[ifc]=ifc2;
                  brch->keenfcangle[ifc]=iang;
		              found=1;
                }
              }
            }
          }
        }
 //       if (iel==11&&ifc==1){printf(" keen assignment iel %d ifc %d root %d ang %d found %d \n",iel,ifc,brch->keen[ifc]->root,brch->keenfcangle[ifc],found);}
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

// Spawn branches to their proper height

  for (iel=0;iel<run->topo->ntrees;iel++){
    tree=run->topo->drys[iel];
    brch=tree->brch;
    ileaves=0;
    tags=malloc((nlvs[iel])*sizeof(int));
    buff=malloc(NEQ*(nlvs[iel])*sizeof(double));
    run->topo->drys[iel]->part=ipart[iel];

    if (oprt[iel]==run->con->rank&&ipart[iel]==run->con->rank){  } // No change
    MPI_Barrier(MPI_COMM_WORLD);
    if (oprt[iel]==run->con->rank&&ipart[iel]!=run->con->rank){ // Tree lost
        while (brch->nkids!=0){
          brch=brch->kids[0];
        }
        while(brch!=NULL&&brch->root==iel){
          brch->tag=tagaddress(brch->adrs,brch->level);
          tags[ileaves]=brch->tag;
          for (iv=0;iv<NEQ;iv++){
            buff[iv+NEQ*ileaves]=brch->el->S[iv];
          }
          ileaves++;
        brch=brch->lnxt;}
        nleaves=ileaves;
        ierr=MPI_Send(&nleaves,1,MPI_INT,ipart[iel],0+3*iel,MPI_COMM_WORLD);
        if (nleaves>1){
          ierr=MPI_Send(tags,nleaves,MPI_INT,ipart[iel],1+3*iel,MPI_COMM_WORLD);
          ierr=MPI_Send(buff,NEQ*nleaves,MPI_DOUBLE,ipart[iel],2+3*iel,MPI_COMM_WORLD);
	      }
        erasetree(run,tree,iel);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (oprt[iel]!=run->con->rank&&ipart[iel]==run->con->rank){ // Tree gained
      ierr=MPI_Recv(&nleaves,1,MPI_INT,(int) oprt[iel],0+3*iel,MPI_COMM_WORLD,&status);
      if (nleaves>1){
        ierr=MPI_Recv(tags,nleaves,MPI_INT,(int) oprt[iel],1+3*iel,MPI_COMM_WORLD,&status);
        ierr=MPI_Recv(buff,NEQ*nleaves,MPI_DOUBLE,(int) oprt[iel],2+3*iel,MPI_COMM_WORLD,&status);
        //	printf("split tree repartined %d",iel);
        for (ileaves=0;ileaves<nleaves;ileaves++){
          brch=spawntree(run,tree,tags[ileaves],iel,1);
          //  leafallocation(run,brch);
          for (iv=0;iv<NEQ;iv++){
            brch->el->S[iv]=buff[iv+NEQ*ileaves];
          }
        }
      }
    }
    free(tags);
    free(buff);
    MPI_Barrier(MPI_COMM_WORLD);
  }


  free(VWGTTOT);
  free(VWGT);
  free(nlvstot);
  free(nlvs);
  free(oprt);
  free(XADJ);
  free(ADJY);
  free(part);
  free(ipart);
}
