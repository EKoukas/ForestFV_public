#include "strdata.h"

void partitiontopo (RUN * run) {

  int ifc,iv,nlfc;
  int iel,iedge,ineig;
  int nedge;
  int i,istart;
  struct TREE * tree;
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
  long * oprt;
  long * oprttot;
  int   * ipart;
  idx_t numflag=0;
  idx_t wgtflag=0;
  idx_t objval=0;
  
  nvert=run->topo->ntrees;
  npart=run->con->size;
  ncon=1;
  part=malloc(nvert*sizeof(idx_t));
  oprt=malloc(nvert*sizeof(long));
  oprttot=malloc(nvert*sizeof(long));

  ipart=malloc(nvert*sizeof(int));
  XADJ=malloc((1+run->topo->ntrees)*sizeof(idx_t));
  VWGT=malloc(run->topo->ntrees*sizeof(idx_t));
  VWGTTOT=malloc(run->topo->ntrees*sizeof(idx_t));
  nlvs=malloc(run->topo->ntrees*sizeof(long));
  nlvstot=malloc(run->topo->ntrees*sizeof(long));
  
  printf("malloc c \n");
  MPI_Barrier(MPI_COMM_WORLD);

  iedge=0;
  for (iel=0;iel<run->topo->ntrees;iel++){
    nlfc=numlfc(run->topo->drys[iel]->type);
    for (ifc=0;ifc<nlfc;ifc++){
       if (run->topo->drys[iel]->conf[ifc]>=0) {iedge++;}
    }
  }
  nedge=iedge;
  ADJY=malloc(nedge*sizeof(idx_t));
  
  printf("malloc d \n");
  MPI_Barrier(MPI_COMM_WORLD);
  
  iedge=0;
  ineig=0;
  for (iel=0;iel<run->topo->ntrees;iel++){
    VWGT[iel]=0;
    oprt[iel]=0;
    tree=run->topo->drys[iel];
    brch=tree->brch;
    if (brch!=NULL){
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
  printf("malloc e \n");
  MPI_Allreduce(nlvs,nlvstot,run->topo->ntrees,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(oprt,oprttot,run->topo->ntrees,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
  
  /*
  for (iel=0;iel<run->topo->ntrees;iel++){
    printf("reduced arrays iel %d vwgt %d oprt %d \n",iel,nlvstot[iel],oprttot[iel],run->topo->ntrees);
  }
  */

  for (iel=0;iel<run->topo->ntrees;iel++){
    VWGT[iel]=nlvstot[iel];
    nlvs[iel]=nlvstot[iel];
    oprt[iel]=oprttot[iel];
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
    //if (part[iel]!=oprt[iel]){   printf("repartition el %d part %d oprt %d proc %d \n",iel,part[iel],oprt[iel],run->con->rank);}
    tree=run->topo->drys[iel];
    brch=tree->brch;
    ileaves=0;
    tags=malloc((nlvs[iel])*sizeof(int));
    buff=malloc(NEQ*(nlvs[iel])*sizeof(double));
    //printf("local first el %d locl %d part %d  \n",iel,run->topo->locl->root,run->con->rank);
    if (brch!=NULL&&part[iel]==run->con->rank){ // No change
    }
    if (brch!=NULL&&part[iel]!=run->con->rank){ // Tree lost
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
        printf(" tree lost el %d ileaves %d tot %d tag %d level %d adrs %d %d \n",iel,ileaves,VWGT[iel],brch->tag,brch->level,brch->adrs[0],brch->adrs[1]);
      brch=brch->lnxt;}
      nleaves=ileaves;
      //if (iel==10){printf("el %d part %d to send %d sz %d tag %d %lf   %lf   %lf   %lf   %lf   %lf   %lf   \n",iel,run->con->rank,part[iel],nleaves,tags[0],buff[0],buff[1],buff[2],buff[3],buff[4],buff[5],buff[6]);}
      ierr=MPI_Send(&nleaves,1,MPI_INT,(int) part[iel],0+3*iel,MPI_COMM_WORLD);
      ierr=MPI_Send(tags,nleaves,MPI_INT,(int) part[iel],1+3*iel,MPI_COMM_WORLD);
      ierr=MPI_Send(buff,NEQ*nleaves,MPI_DOUBLE,(int) part[iel],2+3*iel,MPI_COMM_WORLD);
      erasetree(run,tree,iel);
    }
    else if (brch==NULL&&part[iel]==run->con->rank){ // Tree gained
      ierr=MPI_Recv(&nleaves,1,MPI_INT,(int) oprt[iel],0+3*iel,MPI_COMM_WORLD,&status);
      ierr=MPI_Recv(tags,nleaves,MPI_INT,(int) oprt[iel],1+3*iel,MPI_COMM_WORLD,&status);
      ierr=MPI_Recv(buff,NEQ*nleaves,MPI_DOUBLE,(int) oprt[iel],2+3*iel,MPI_COMM_WORLD,&status);
      //if (iel==10){printf("el %d part %d to recieve %d sz %d tag %d %lf  %lf   %lf   %lf   %lf   %lf   %lf  \n",iel,run->con->rank,oprt[iel],nleaves,tags[0],buff[0],buff[1],buff[2],buff[3],buff[4],buff[5],buff[6]);}
      for (ileaves=0;ileaves<nleaves;ileaves++){
	      
        printf("spawn tree of el %d lv %d tag %d sz \n",iel,ileaves,tags[ileaves]);
        brch=spawntree(run,tree,tags[ileaves],iel,1);
	// keep keen connectivity
        for (iv=0;iv<NEQ;iv++){
          brch->el->S[iv]=buff[iv+NEQ*ileaves];
        }
      }
 //     printf("el %d part %d to recieved %d \n",iel,run->con->rank,oprt[iel]);
 
    }
    free(tags);
    free(buff);
  
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
  }
  free(VWGTTOT);
  free(VWGT);
  free(nlvstot);
  free(nlvs);
  free(oprt);
  free(oprttot);
  free(XADJ);
  free(ADJY);
  free(part);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
     

//  crnt=run->topo->locl;
//  while (crnt!=NULL){
// //   if (run->con->rank==0){ printf("befor a list iel %d tag %d  \n",crnt->root,crnt->tag);}
//    crnt=crnt->lnxt;}
//    MPI_Barrier(MPI_COMM_WORLD);
//    MPI_Barrier(MPI_COMM_WORLD);
//  
//  crnt=run->topo->locl;
//  while (crnt!=NULL){
// //   if (run->con->rank==1){ printf("befor b list iel %d tag %d  \n",crnt->root,crnt->tag);}
//   crnt=crnt->lnxt;}
//
//    MPI_Barrier(MPI_COMM_WORLD);
//    MPI_Barrier(MPI_COMM_WORLD);
//

// Build lists
  for (iel=0;iel<run->topo->ntrees;iel++){ // Loop elements
    tree=run->topo->drys[iel];
    brch=tree->brch;
    if (brch!=NULL){
    }
  }
  prvs=run->topo->glob;
  istart=0;
  for (iel=0;iel<run->topo->ntrees;iel++){ // Loop elements
    tree=run->topo->drys[iel];
    brch=tree->brch;
    if (brch!=NULL){
      // zero indx
      for (ilvl=0;ilvl<10;ilvl++){
	      indx[ilvl]=0;
      }
      if (brch->nkids==0){
	  brch->lnxt=NULL;
	  brch->lprv=prvs;
	  prvs->lnxt=brch;
	  prvs=brch;
	  if (istart==0){run->topo->locl=brch; istart=1;}
      }
      ilvl=0;
      crnt=brch;
      while (indx[1]<brch->nkids){
	// decent
	ilvl=crnt->level;
        if (crnt->nkids!=0){
		if (indx[ilvl]<crnt->nkids){
		   crnt=crnt->kids[indx[ilvl]];
		   indx[ilvl]++;
		}
		else{
			crnt=crnt->prnt;
		        indx[ilvl]=0;
		}
	}
	// hit bottom
	else{
	  crnt->lnxt=NULL;
	  crnt->lprv=prvs;
	  prvs->lnxt=crnt;
	  prvs=crnt;
	  if (istart==0){run->topo->locl=crnt; istart=1;}
          crnt=crnt->prnt;
 //         printf("advancing on leaves %d part %d level %d adrs %d \n",iel,run->con->rank,crnt->level,crnt->adrs[ilvl]);
        }
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  for (iel=0;iel<run->topo->ntrees;iel++){
    tree=run->topo->drys[iel];
    brch=tree->brch;
    if (brch!=NULL){
      printf("check drys iel %d part %d nkids %d \n",iel,run->con->rank,brch->nkids);
      while (brch->nkids!=0){
        brch=brch->kids[0];
      }
      while(brch!=NULL&&brch->root==iel){
       printf("check list iel %d level %d tags %d adrs %d %d  last %d \n",iel,brch->level,brch->tag,brch->adrs[0],brch->adrs[1],brch->adrs[brch->level]);

      brch=brch->lnxt;}
    }
  }


// Build buffers

// Build neigbours
//
//
}
