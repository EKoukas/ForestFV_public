#include "strdata.h"

void rebuildprox(RUN * run) {
  
  struct BRANCH * crnt;
  struct BRANCH * prox;
  struct BRANCH * buff;
  struct BRANCH * hkin;

  int ipts;
  int iptr;
  int * iprox;
  int ierr;
  int iel,tag,ielneig,ifc;
  int ibuf,iprx,ifc1,ifc2,iel1,iel2,tag1,tag2,angle1,angle2;
  int iang,match,ind,indr,found;
  int aa,ba,ca,da;

  MPI_Status status;


  // start flag for buffer lists one per partition for each partition
  if (run->topo->partitioned==1){
    for (ipts=0;ipts<run->con->size;ipts++){
      if (run->topo->nprox[ipts]>0){
        free(run->topo->proxtags[ipts]);
        free(run->topo->proxdrys[ipts]);
        free(run->topo->proxlfc1[ipts]);
        free(run->topo->proxlfc2[ipts]);
      }
      if (run->topo->nbuff[ipts]>0){
        free(run->topo->bufftags[ipts]);
        free(run->topo->buffdrys[ipts]);
        free(run->topo->bufflfc1[ipts]);
        free(run->topo->bufflfc2[ipts]);
      }
    }
    free(run->topo->proxtags);
    free(run->topo->bufftags);
    free(run->topo->proxdrys);
    free(run->topo->buffdrys);
    free(run->topo->proxlfc1);
    free(run->topo->proxlfc2);
    free(run->topo->bufflfc1);
    free(run->topo->bufflfc2);
    free(run->topo->nprox);
    free(run->topo->nbuff);
  }

  //message("MMFV: RPRX: MALC \n",run);
  iprox=malloc(run->con->size*sizeof(int));

  run->topo->nprox=malloc(run->con->size*sizeof(int));
  run->topo->proxtags=malloc(run->con->size*sizeof(int *));
  run->topo->proxdrys=malloc(run->con->size*sizeof(int *));
  run->topo->proxlfc1=malloc(run->con->size*sizeof(int *));
  run->topo->proxlfc2=malloc(run->con->size*sizeof(int *));

  run->topo->nbuff=malloc(run->con->size*sizeof(int));
  run->topo->buffdrys=malloc(run->con->size*sizeof(int *));
  run->topo->bufftags=malloc(run->con->size*sizeof(int *));
  run->topo->bufflfc1=malloc(run->con->size*sizeof(int *));
  run->topo->bufflfc2=malloc(run->con->size*sizeof(int *));

  //message("MMFV: RPRX: ZERO \n",run);
  //  zero counter
  for (ipts=0;ipts<run->con->size;ipts++){
    run->topo->nprox[ipts]=0;
    run->topo->nbuff[ipts]=0;
  }

  // count prox branch faces tags
  crnt=run->topo->locl;
  while (crnt!=NULL){
    iel=crnt->root;
    for (ifc=0; ifc<crnt->nlfc;ifc++){
      ielneig=ielconn(crnt,ifc);
      if (ielneig!=-1){
        if (run->topo->drys[iel]->part!=run->topo->drys[ielneig]->part){
          run->topo->nprox[run->topo->drys[ielneig]->part]++;
        }
      }
    }
    crnt=crnt->lnxt;
  }
  // alloc prox branch faces tags

  //message("MMFV: RPRX: MALC \n",run);
  for (ipts=0;ipts<run->con->size;ipts++){
    run->topo->proxtags[ipts]=malloc(max(run->topo->nprox[ipts],0)*sizeof(int));
    run->topo->proxdrys[ipts]=malloc(max(run->topo->nprox[ipts],0)*sizeof(int));
    run->topo->proxlfc1[ipts]=malloc(max(run->topo->nprox[ipts],0)*sizeof(int));
    run->topo->proxlfc2[ipts]=malloc(max(run->topo->nprox[ipts],0)*sizeof(int));
  }

  //message("MMFV: RPRX: CONT \n",run);
  // zerp counter
  for (ipts=0;ipts<run->con->size;ipts++){
    iprox[ipts]=0;
  }

  // popul prox branch faces tags
  crnt=run->topo->locl;
  while (crnt!=NULL){
    iel=crnt->root;
    for (ifc=0; ifc<crnt->nlfc; ifc++){
      ielneig=ielconn(crnt,ifc);
      if (ielneig!=-1){
        if (run->topo->drys[iel]->part!=run->topo->drys[ielneig]->part){
          run->topo->proxtags[run->topo->drys[ielneig]->part][iprox[run->topo->drys[ielneig]->part]]=crnt->tag;
          run->topo->proxdrys[run->topo->drys[ielneig]->part][iprox[run->topo->drys[ielneig]->part]]=iel;
          run->topo->proxlfc1[run->topo->drys[ielneig]->part][iprox[run->topo->drys[ielneig]->part]]=ifc;
          run->topo->proxlfc2[run->topo->drys[ielneig]->part][iprox[run->topo->drys[ielneig]->part]]=0;
	        iprox[run->topo->drys[ielneig]->part]++;
        }
      }
    } 
    crnt=crnt->lnxt;
  }

  // communicate proxes
  for (ipts=0;ipts<run->con->size;ipts++){
    if (ipts!=run->con->rank){
      // send prox to ipts
      ierr=MPI_Send(&(run->topo->nprox[ipts]),1,MPI_INT,ipts,5*run->con->rank+0,MPI_COMM_WORLD);
      if (run->topo->nprox[ipts]>0){
        ierr=MPI_Send(run->topo->proxtags[ipts],run->topo->nprox[ipts],MPI_INT,ipts,5*run->con->rank+1,MPI_COMM_WORLD);
        ierr=MPI_Send(run->topo->proxdrys[ipts],run->topo->nprox[ipts],MPI_INT,ipts,5*run->con->rank+2,MPI_COMM_WORLD);
        ierr=MPI_Send(run->topo->proxlfc1[ipts],run->topo->nprox[ipts],MPI_INT,ipts,5*run->con->rank+3,MPI_COMM_WORLD);
        ierr=MPI_Send(run->topo->proxlfc2[ipts],run->topo->nprox[ipts],MPI_INT,ipts,5*run->con->rank+4,MPI_COMM_WORLD);
      }
    }
    // ipts to recv buff 
    if (ipts==run->con->rank){
      for (iptr=0;iptr<run->con->size;iptr++){
        if (iptr!=run->con->rank){
          ierr=MPI_Recv(&(run->topo->nbuff[iptr]),1,MPI_INT,iptr,5*iptr+0,MPI_COMM_WORLD,&status);
          run->topo->bufftags[iptr]=malloc(max(run->topo->nbuff[iptr],0)*sizeof(int));
          run->topo->buffdrys[iptr]=malloc(max(run->topo->nbuff[iptr],0)*sizeof(int));
          run->topo->bufflfc1[iptr]=malloc(max(run->topo->nbuff[iptr],0)*sizeof(int));
          run->topo->bufflfc2[iptr]=malloc(max(run->topo->nbuff[iptr],0)*sizeof(int));
          if (run->topo->nbuff[iptr]>0){
            ierr=MPI_Recv(run->topo->bufftags[iptr],run->topo->nbuff[iptr],MPI_INT,iptr,5*iptr+1,MPI_COMM_WORLD,&status);
            ierr=MPI_Recv(run->topo->buffdrys[iptr],run->topo->nbuff[iptr],MPI_INT,iptr,5*iptr+2,MPI_COMM_WORLD,&status);
            ierr=MPI_Recv(run->topo->bufflfc1[iptr],run->topo->nbuff[iptr],MPI_INT,iptr,5*iptr+3,MPI_COMM_WORLD,&status);
            ierr=MPI_Recv(run->topo->bufflfc2[iptr],run->topo->nbuff[iptr],MPI_INT,iptr,5*iptr+4,MPI_COMM_WORLD,&status);
          }
	      }
      }
    }
  }

  // spawn proxes as buffers to neigbouring branches one LEVEL only
  for (ipts=0;ipts<run->con->size;ipts++){
    for (ibuf=0;ibuf<run->topo->nbuff[ipts];ibuf++){
      tag=run->topo->bufftags[ipts][ibuf];
      iel=run->topo->buffdrys[ipts][ibuf];
      crnt=spawntree(run,run->topo->drys[iel],0,iel,0);
      crnt->hangn=2;  // 0->hanging to own tree 1-> no hanging keen branch 2-> hanging to neig tree
      crnt->el->crith=0.0;
      setsplitcalc(crnt,run);
      crnt->part=run->con->rank;
      // assign a kin to the buffer on buff face
      ifc=run->topo->bufflfc1[ipts][ibuf];
      iel2=run->topo->drys[iel]->conf[ifc];
      if (iel2==-1){
        crnt->keen[ifc]=run->topo->glob;
        crnt->keenfc[ifc]=0;
        crnt->keenfcangle[ifc]=0;
      }
      else{
        crnt->keen[ifc]=run->topo->drys[iel2]->brch;
        crnt->keenfc[ifc]=0;
        crnt->keenfcangle[ifc]=0;
      }
    }
  }
  
  // hook up the hanging keen to the buffer

  for (ipts=0;ipts<run->con->size;ipts++){
    for (iprx=0;iprx<run->topo->nprox[ipts];iprx++){
      tag1=0;
      iel1=run->topo->proxdrys[ipts][iprx];
      ifc1=run->topo->proxlfc1[ipts][iprx];
      prox=run->topo->drys[iel1]->brch;
      if (prox->keen[ifc1]->root!=-1&&prox->keen[ifc1]->hangn==1){
        hkin=prox->keen[ifc1];
        ifc2=prox->keenfc[ifc1];
        iang=prox->keenfcangle[ifc1];
	      //  if (iel1==11){ printf("hanging a keen iel %d hng %d ang %d ifc %d \n",iel1,buff->keen[ifc]->nkids,hkin->keenfcangle[ifc],ifc); }
        prox->keen[ifc1]=run->topo->drys[hkin->root]->brch; // hook up keen on non-hanging branch
        prox->keenfc[ifc1]=ifc2;
        prox->keenfcangle[ifc1]=iang;
        leafdeallocation(run,hkin);
        celldeallocation(hkin);
        destroybrch(hkin);
	      // if (iel1==11){ printf("hanging b keen iel %d hng %d ang %d ifc %d \n",iel1,buff->keen[ifc]->nkids,buff->keenfcangle[ifc],ifc); }
      } else {
        // printf("oups %d hng %d ang %d ifc %d \n",iel1,buff->keen[ifc1]->nkids,buff->keenfcangle[ifc1],ifc1); 
      }
    }
  }

  // kin connectivity ang and ifc2 assignment
  for (ipts=0;ipts<run->con->size;ipts++){
    for (ibuf=0;ibuf<run->topo->nbuff[ipts];ibuf++){
      tag1=0;
      iel1=run->topo->buffdrys[ipts][ibuf];
      buff=run->topo->drys[iel1]->brch;
      ifc1=run->topo->bufflfc1[ipts][ibuf];
        if (buff->keen[ifc1]!=NULL){
          for(ifc2=0;ifc2<buff->keen[ifc1]->nlfc;ifc2++){
            if (buff->keen[ifc1]->cl->fc[ifc2].type==buff->cl->fc[ifc1].type){
              for(iang=0;iang<buff->keen[ifc1]->cl->fc[ifc2].nnds;iang++){
                match=0;
                for(ind=0;ind<buff->keen[ifc1]->cl->fc[ifc2].nnds;ind++){
                  indr=fcndrot(ind,iang,buff->keen[ifc1]->cl->fc[ifc2].type);
                  if(buff->cl->nd[fcnd2elnd(ifc1,ind,buff->type)].num==buff->keen[ifc1]->cl->nd[fcrefnd2elnd(ifc2,indr,buff->keen[ifc1]->type)].num){match++;}
                }
                if (match==buff->cl->fc[ifc1].nnds){
                  buff->keenfc[ifc1]=ifc2;
                  buff->keenfcangle[ifc1]=iang;
                  found=1;
                }
              }
            }
          }
      }
    }
  }

  // build proxes as buffers to neigbouring branches full height
  for (ipts=0;ipts<run->con->size;ipts++){
    for (ibuf=0;ibuf<run->topo->nbuff[ipts];ibuf++){
      tag=run->topo->bufftags[ipts][ibuf];
      iel=run->topo->buffdrys[ipts][ibuf];
      crnt=spawntree(run,run->topo->drys[iel],tag,iel,0);
      crnt->hangn=2;  // 0->hanging to own tree 1-> no hanging keen branch 2-> hanging to neig tree
      crnt->el->crith=0.0;
      setsplitcalc(crnt,run);
      crnt->part=run->con->rank;
    }
  }

  free(iprox);
}