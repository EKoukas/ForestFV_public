#include "strdata.h"

void repartition_parmetis(RUN * run)
{
int ifc,ifc2,iv,nlfc,ipt,ipt2,totlvs;
int iel,ineig,lel;
int nedge;
int ltrees;
int ltrees_sum;
int ntrees;
int lvert;
int ledge;
int iedge;
int * lpart;
int match,iang,ind,indr,found;
int i,istart;
struct TREE * tree;
struct BRANCH * keen;
struct BRANCH * brch;
struct BRANCH * crnt;
struct BRANCH * prvs;
struct BRANCH * prnt;
int totleaves;
int nleaves;
int ileaves;
int lleaves;
int ilvl,ikid;
int ibuf,tag;
int * tags;
double * buff;
int * adr;
int nlvl;
int ierr;
int indx[10];
int partlvs[1000];
MPI_Status status;
  idx_t edgecut;
  idx_t * VWGTTOT;
  int nvert;

  idx_t * VTXDIST;
  idx_t * XADJ;
  idx_t * ADJY;
  idx_t * VWGT;
  idx_t * ADJWGT;
  idx_t * VSIZE;
  idx_t wgtflag;
  idx_t numflag;
  idx_t ncon;
  idx_t nparts;
  real_t * TPWGTS;
  real_t * UBVEC;
  real_t itr;
  idx_t * OPTIONS;
  idx_t * LPRT;
  int aa,ba,ca,da;


  int * nlvs;
  int * nlvstot;
  int * oprt;
  int * part;
  int * ipart;
  int * iparttot;
  int tpp;
  int iel1;
  int iel2;




  //watch(run,2,0); // RPRT
  ltrees=0;
  tpp=((int) (((float) run->topo->ntrees)/((float) run->con->size)));
  iel1=tpp*run->con->rank;
  iel2=tpp*(run->con->rank+1);
  if (run->con->rank==run->con->size-1){ iel2=run->topo->ntrees; }
  ltrees=iel2-iel1;


 //// printf("Entering repartion ltrees %d rank %d \n",ltrees,run->con->rank);
  ntrees=run->topo->ntrees;

  nvert=run->topo->ntrees+1;
  lvert=ltrees+1;
  part=malloc(ntrees*sizeof(int));
  oprt=malloc(ntrees*sizeof(int));
  ipart=malloc(ntrees*sizeof(int));
  iparttot=malloc(ntrees*sizeof(int));
  LPRT=malloc(ltrees*sizeof(idx_t));
  XADJ=malloc(lvert*sizeof(idx_t));
  VWGT=malloc(ltrees*sizeof(idx_t));
  VSIZE=malloc(ltrees*sizeof(idx_t));
  OPTIONS=malloc(3*sizeof(idx_t));
  VTXDIST=malloc((run->con->size+1)*sizeof(idx_t));
  TPWGTS=malloc(run->con->size*sizeof(real_t));
  UBVEC=malloc(sizeof(real_t));
//  VWGTTOT=malloc(run->topo->ntrees*sizeof(idx_t));
  nlvs=malloc(ntrees*sizeof(int));
  nlvstot=malloc(ntrees*sizeof(int));



  ledge=0;
  for (iel=0;iel<run->topo->ntrees;iel++){
    nlvs[iel]=0;
    oprt[iel]=run->topo->drys[iel]->part;
  }

  for (iel=iel1;iel<iel2;iel++){
    nlfc=numlfc(run->topo->drys[iel]->type);
    for (ifc=0;ifc<nlfc;ifc++){
      if (run->topo->drys[iel]->conf[ifc]>=0) {ledge++;}
    }
  }
  ADJY=malloc(ledge*sizeof(idx_t));
  //ADJWGT=malloc(ledge*sizeof(idx_t));
  lleaves=0;
  for (iel=0;iel<run->topo->ntrees;iel++){
    nlvs[iel]=0;
    tree=run->topo->drys[iel];
    brch=tree->brch;
    if (run->topo->drys[iel]->part==run->con->rank){
      while (brch->nkids!=0){
        brch=brch->kids[0];
      }
      while(brch!=NULL&&brch->root==iel){
      nlvs[iel]++;
      lleaves++;
      brch=brch->lnxt;}
    }
  }

  //printf ("number of local leaves: %d %d \n",lleaves,run->con->rank);
  run->topo->pleaves=lleaves;
  MPI_Allreduce(nlvs,nlvstot,run->topo->ntrees,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  for (iel=0;iel<run->topo->ntrees;iel++){
    nlvs[iel]=nlvstot[iel];
  }
  MPI_Allreduce(&lleaves,&totleaves,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  run->topo->tleaves=totleaves;


  lel=0;
  for (iel=iel1;iel<iel2;iel++){
    tree=run->topo->drys[iel];
    brch=tree->brch;
    VWGT[lel]=nlvs[iel];
    VSIZE[lel]=VWGT[lel];
    lel++;
  }


 

  lel=0;
  iedge=0;
  ineig=0;
  for (iel=iel1;iel<iel2;iel++){
    if (lel==0) {XADJ[lel]=0;}
    if (lel!=0) {XADJ[lel]=XADJ[lel-1]+ineig;}
    LPRT[lel]=run->topo->drys[iel]->part;
    ineig=0;
    nlfc=numlfc(run->topo->drys[iel]->type);
    for (ifc=0;ifc<nlfc;ifc++){
      if (run->topo->drys[iel]->conf[ifc]>=0){
        ADJY[iedge]=run->topo->drys[iel]->conf[ifc];
        //ADJWGT[iedge]=nlvs[iel];
        iedge++;
        ineig++;
      }
    }
    lel++;
  } 
  XADJ[ltrees]=XADJ[ltrees-1]+ineig;
  wgtflag=2;
  numflag=0; 
  ncon=1;
  nparts=run->con->size;
  UBVEC[0]=1.02;
  ltrees_sum=0;
  VTXDIST[0]=0;
  for(ipt=0;ipt<run->con->size;ipt++){
    TPWGTS[ipt]=1.0/((double) nparts);
 //   TPWGTS[ipt]=1.0;
 //   nparts=1;
    VTXDIST[ipt+1]=0;
    if (ipt==run->con->rank){ltrees_sum+=ltrees;}
    MPI_Bcast(&(ltrees_sum),1,MPI_INT,ipt,MPI_COMM_WORLD);
    VTXDIST[ipt+1]=ltrees_sum;
 ////   printf("vtxdis ipt %d dist %d ltrees %d part %d \n",ipt,(int) VTXDIST[ipt],ltrees,run->con->rank);
  }
  ////printf("vtxdis ipt %d dist %d ltrees %d part %d \n",ipt,(int) VTXDIST[ipt],ltrees,run->con->rank);
  //ADJWGT=NULL;
  OPTIONS[0]=0;
  OPTIONS[1]=0;
  OPTIONS[2]=0;

  for(ipt=0;ipt<run->con->size;ipt++){
    if (ipt==run->con->rank){
//     for (iel=0;iel<ltrees+1;iel++){ printf(" i %3d",iel);} printf("\n "); 
//     for (iel=0;iel<ltrees+1;iel++){ printf(" x %3d",XADJ[iel]); } printf("\n "); 
    }
  }
  for(ipt=0;ipt<run->con->size;ipt++){
    if (ipt==run->con->rank){
//     for (iel=0;iel<ledge;iel++){ printf(" i %3d",iel);} printf("\n "); 
//     for (iel=0;iel<ledge;iel++){ printf(" a %3d",ADJY[iel]); } printf("\n "); 
    }
  }
//  for (iel=0;iel<ledge;iel++){
//  }
//  MPI_Comm comm_val=MPI_COMM_WORLD;
//  MPI_Comm * comm=&comm_val;
//  VSIZE=malloc(5*sizeof(idx_t));
//  TPWGTS=malloc(5*sizeof(real_t));
//  UBVEC=malloc(sizeof(real_t));
//  LPRT=malloc(5*sizeof(idx_t));
//  OPTIONS=malloc(3*sizeof(idx_t));
//  OPTIONS[0]=0;
//  OPTIONS[1]=0;
//  OPTIONS[2]=0;
//  wgtflag=0;
//  numflag=0;
//  ncon=1;
//  nparts=3;
//  VSIZE[0]=1;
//  VSIZE[1]=1;
//  VSIZE[2]=1;
//  VSIZE[3]=1;
//  VSIZE[4]=1;
//  TPWGTS[0]=1.0/3.0;
//  TPWGTS[1]=1.0/3.0;
//  TPWGTS[2]=1.0/3.0;
//  UBVEC[0]=1.25;
  itr=1000.0;

//  if (run->con->rank==0){
//    XADJ=malloc(6*sizeof(idx_t));
//    ADJY=malloc(13*sizeof(idx_t));
//    XADJ[0]=0; XADJ[1]=2; XADJ[2]=5; XADJ[3]=8; XADJ[4]=11; XADJ[5]=13;
//    ADJY[0]=0; ADJY[1]=5; ADJY[2]=1; ADJY[3]=2; ADJY[4]=6; ADJY[5]=1; ADJY[6]=3; ADJY[7]=7; ADJY[8]=2; ADJY[9]=4; ADJY[10]=8; ADJY[11]=3; ADJY[12]=9;
//  }
//  if (run->con->rank==1){
//    XADJ=malloc(6*sizeof(idx_t));
//    ADJY=malloc(18*sizeof(idx_t));
//    XADJ[0]=0; XADJ[1]=3; XADJ[2]=7; XADJ[3]=11; XADJ[4]=15; XADJ[5]=18;
//    ADJY[0]=0; ADJY[1]=6; ADJY[2]=10; ADJY[3]=1; ADJY[4]=5; ADJY[5]=7; ADJY[6]=11; ADJY[7]=2; ADJY[8]=6; ADJY[9]=8; ADJY[10]=12; ADJY[11]=3; ADJY[12]=7; ADJY[13]=9; ADJY[14]=13; ADJY[15]=4; ADJY[16]=8; ADJY[17]=14;
//  }
//  if (run->con->rank==2){
//    XADJ=malloc(6*sizeof(idx_t));
//    ADJY=malloc(13*sizeof(idx_t));
//    XADJ[0]=0; XADJ[1]=2; XADJ[2]=5; XADJ[3]=8; XADJ[4]=11; XADJ[5]=13;
//    ADJY[0]=5; ADJY[1]=11; ADJY[2]=6; ADJY[3]=10; ADJY[4]=12; ADJY[5]=7; ADJY[6]=11; ADJY[7]=13; ADJY[8]=8; ADJY[9]=12; ADJY[10]=14; ADJY[11]=9; ADJY[12]=13; 
//  }
//  LPRT[0]=run->con->rank;
//  LPRT[1]=run->con->rank;
//  LPRT[2]=run->con->rank;
//  LPRT[3]=run->con->rank;
//  LPRT[4]=run->con->rank;
////  VTXDIST=malloc((run->con->size+1)*sizeof(idx_t));
//  VTXDIST[0]=0;
//  VTXDIST[1]=5;
//  VTXDIST[2]=10;
//  VTXDIST[3]=15;

 ParMETIS_V3_AdaptiveRepart(VTXDIST,XADJ,ADJY,VWGT,VSIZE,NULL,&wgtflag,&numflag,&ncon,&nparts,TPWGTS,UBVEC,&itr,OPTIONS,&edgecut,LPRT,&(run->con->comm));
//wgtflag=0;
 // ParMETIS_V3_AdaptiveRepart(VTXDIST,XADJ,ADJY,NULL,NULL,NULL,&wgtflag,&numflag,&ncon,&nparts,TPWGTS,UBVEC,&itr,OPTIONS,&edgecut,LPRT,&(run->con->comm));

  lel=0;
  for (iel=0;iel<run->topo->ntrees;iel++){
    ipart[iel]=0;
  }
  for (iel=iel1;iel<iel2;iel++){
    ipart[iel]=((int) LPRT[lel]);
    lel++;
  }
  MPI_Allreduce(ipart,iparttot,run->topo->ntrees,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

  for (iel=0;iel<run->topo->ntrees;iel++){
    if (g_istep%2000==0)  { 
      ipart[iel]=iparttot[iel]; 
    } else {
      ipart[iel]=oprt[iel];
    }

	  part[iel]=iparttot[iel];
	  nlvs[iel]=nlvstot[iel];
    //if (run->con->geo_adapth==0) { ipart[iel]=oprt[iel];}
  }

 for (ipt=0;ipt<run->con->size;ipt++){
	 partlvs[ipt]=0;
	 totlvs=0;
 }
 for (iel=0;iel<run->topo->ntrees;iel++){
	 partlvs[ipart[iel]]+=nlvs[iel];
	 totlvs+=nlvs[iel];
 //   if (run->con->geo_adapth==0) { ipart[iel]=oprt[iel];}
  }
  run->tot_leaves=totlvs;
  
  //if (run->con->istep%5000==0)  { 
    for (ipt=0;ipt<1;ipt++){
      if (run->con->rank==ipt){
          printf("Repartion rank %d leaves: ",run->con->rank);
          for (ipt2=0;ipt2<run->con->size;ipt2++){
            printf("%d ",partlvs[ipt2]);
          }
        printf("total %d \n",totlvs);
        run->tot_leaves=totlvs;
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  //}
 //printf("Finished repartion 0 pct %d oct %d rank %d \n",partpct[run->con->rank],partoct[run->con->rank],run->con->rank);
//  printf("Finished repartion 0 pct %d oct %d rank %d \n",partpct[0],partoct[0],run->con->rank);
//  printf("Finished repartion 1 pct %d oct %d rank %d \n",partpct[1],partoct[1],run->con->rank);
// partition changes, trees must be erased and branches are migrated to new trees that are spawn
// hanging kins are created, or kins are assigned in new trees
// hanging kins are destroyed and assigned trees are deaccosiated in old trees
// for non-affected trees kins can change from assigned to hanging and vis versa


// Destroy all hanging keens

  for (iel=0;iel<run->topo->ntrees;iel++){
    tree=run->topo->drys[iel];
    brch=tree->brch;
    ileaves=0;
    tags=malloc((nlvs[iel])*sizeof(int));
    buff=malloc(NEQ*(nlvs[iel])*sizeof(double));
    run->topo->drys[iel]->part=ipart[iel];
    if (oprt[iel]==run->con->rank){   // In in old decomposition destroy hanging branches (essentially hanging kins)
      for (ifc=0;ifc<numlfc(run->topo->drys[iel]->type);ifc++){
        if (brch->keen!=NULL&&brch->keen[ifc]->hangn==1){
	  leafdeallocation(run,brch->keen[ifc]);
	  celldeallocation(brch->keen[ifc]);
          destroybrch(brch->keen[ifc]);
        }
      }
    }
    free(tags);
    free(buff);
  }

//  adr=malloc(10*sizeof(int));
//  for (ipt=0;ipt<run->con->size;ipt++){
//    if (run->topo->nbuff[ipt]>0){
//    for (ibuf=0;ibuf<run->topo->nbuff[ipt];ibuf++){
//      tag=run->topo->bufftags[ipt][ibuf];
//      iel=run->topo->buffdrys[ipt][ibuf];
//      tag2adr(adr,tag,&nlvl);
//      crnt=run->topo->drys[iel]->brch;
//      ilvl=0;
//      while (crnt->nkids!=0){
//        crnt=crnt->kids[adr[ilvl]-1];
//        ilvl++;
//      }
//      leafdeallocation(run,crnt);
//      celldeallocation(crnt);
//      destroybrch(crnt);
//      while (crnt->level>0){
//	for (ikid=0;ikid<crnt->nkids;ikid++){
//          crnt=crnt->prnt;
//          leafdeallocation(run,crnt);
//          celldeallocation(crnt);
//          destroybrch(crnt);
//	}
//      }
//    }
//    }
//  }
//  free(adr);


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
  }

 // Spawn keens

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
  }

// Spawn branches to their proper height

  for (iel=0;iel<run->topo->ntrees;iel++){
    tree=run->topo->drys[iel];
    brch=tree->brch;
    ileaves=0;
    run->topo->drys[iel]->part=ipart[iel];

    if (oprt[iel]==run->con->rank&&ipart[iel]!=run->con->rank){ // Tree lost
      tags=malloc((nlvs[iel])*sizeof(int));
      buff=malloc(NEQ*(nlvs[iel])*sizeof(double));
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
      free(tags);
      free(buff);
    }
    if (oprt[iel]!=run->con->rank&&ipart[iel]==run->con->rank){ // Tree gained
      tags=malloc((nlvs[iel])*sizeof(int));
      buff=malloc(NEQ*(nlvs[iel])*sizeof(double));
      ierr=MPI_Recv(&nleaves,1,MPI_INT,(int) oprt[iel],0+3*iel,MPI_COMM_WORLD,&status);
      if (nleaves>1){
        ierr=MPI_Recv(tags,nleaves,MPI_INT,(int) oprt[iel],1+3*iel,MPI_COMM_WORLD,&status);
        ierr=MPI_Recv(buff,NEQ*nleaves,MPI_DOUBLE,(int) oprt[iel],2+3*iel,MPI_COMM_WORLD,&status);
//	printf("split tree repartined %d",iel);
        for (ileaves=0;ileaves<nleaves;ileaves++){
          brch=spawntree(run,tree,tags[ileaves],iel,1);
    //      leafallocation(run,brch);
          for (iv=0;iv<NEQ;iv++){
            brch->el->S[iv]=buff[iv+NEQ*ileaves];
          }
        }
      }
      free(tags);
      free(buff);
    }
    // Was slow without it try now with it and its slow again
    // MPI_Barrier(MPI_COMM_WORLD);
  }

  lleaves=0;
  for (iel=0;iel<run->topo->ntrees;iel++){
    tree=run->topo->drys[iel];
    brch=tree->brch;
    if (run->topo->drys[iel]->part==run->con->rank){
      while (brch->nkids!=0){
        brch=brch->kids[0];
      }
      while(brch!=NULL&&brch->root==iel){
      lleaves++;
      brch=brch->lnxt;}
    }
  }

  run->topo->pleaves=lleaves;
  MPI_Allreduce(&lleaves,&totleaves,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  run->topo->tleaves=totleaves;


  free(part);
  free(oprt);
  free(LPRT);
  free(ipart);
  free(iparttot);
  free(XADJ);
  free(ADJY);
  free(VWGT);
  free(VSIZE);
  free(OPTIONS);
  free(TPWGTS);
  free(UBVEC);
  free(nlvs);
  free(nlvstot);
  free(VTXDIST);
  //watch(run,2,1); // RPRT
}
