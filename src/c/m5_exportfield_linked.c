#include "strdata.h"

void m5_exportfield_linked (struct RUN * run) {
  
  int k,iv,in,isp,iel,ing;
  int ind;
  int ifc;
  int el2;
  int nd2;
  int ipt;
  int iang;
  int ifc2;
  int ngeltot;
  int ngndtot;
  int igndtot;
  int llinks;
  int nlinks,nelems;
  int nlinkstot,nelemstot;
  int ignd;
  int ngnd;
  int ngel;
  int type;
  int owrt;
  int imatv,jmatv;
  int TPNDS[8][8];
  int INX[8],INY[8],INZ[8];
  int INXP[8],INYP[8],INZP[8];
  int INXT[8],INYT[8],INZT[8];
  int ielem,inode,ilink;
  int elem1,node1,link1;
  int elem2,node2,link2;
  MPI_Status status;
  MPI_Request request;

  int il;
  int ill;
  int elcnt;
  int ** ELNDLK;
  int ** LKEL;
  int ** LKND;
  int * LKNN;
  int ELNDLKOF;
  double x,y,z,r;
  double dx,dy,dz;
  double Y[80];
  double var[80];
  double press[10];
  struct BRANCH * brch;
  struct BRANCH ** BRANCHES;
  struct TREE * tree;
  char  fstring[200]="                                             ";
  FILE * tp;

  int * glob_nodes;
  int * exprt_ingr;
  double * exprt_real;
  int nll;
  int nsz;
  int nllinks;

  double S1,S2,S3,S_tens,S_comp,rho_temp,rhos_temp;
  double as;
  int i,j,eq_0_3;
  double Stress[3][3];
  double mu_temp,gamma_temp,pinf_temp,rho0_temp,avf,avf_sum,p_hydro,rho;

  double max_si,min_si,T_max;

  nll=5+(NEQ+(NMATERIALS+1)+9+4);
  nsz=run->con->size;

  TPNDS[0][0]=0;
  TPNDS[1][0]=1;
  TPNDS[2][0]=3;
  TPNDS[3][0]=2;
  TPNDS[4][0]=4;
  TPNDS[5][0]=5;
  TPNDS[6][0]=7;
  TPNDS[7][0]=6;


  TPNDS[0][1]=1;
  TPNDS[1][1]=2;
  TPNDS[2][1]=2;
  TPNDS[3][1]=0;
  TPNDS[4][1]=3;
  TPNDS[5][1]=3;
  TPNDS[6][1]=3;
  TPNDS[7][1]=3;

  TPNDS[0][2]=2;
  TPNDS[1][2]=5;
  TPNDS[2][2]=3;
  TPNDS[3][2]=0;
  TPNDS[4][2]=1;
  TPNDS[5][2]=4;
  TPNDS[6][2]=4;
  TPNDS[7][2]=1;

  INX[0]=1;INY[0]=0;INZ[0]=0;
  INX[1]=1;INY[1]=1;INZ[1]=0;
  INX[2]=0;INY[2]=1;INZ[2]=0;
  INX[3]=0;INY[3]=0;INZ[3]=0;
  INX[4]=1;INY[4]=0;INZ[4]=1;
  INX[5]=1;INY[5]=1;INZ[5]=1;
  INX[6]=0;INY[6]=1;INZ[6]=1;
  INX[7]=0;INY[7]=0;INZ[7]=1;

  INXP[0]=1;INYP[0]=0;INZP[0]=0;
  INXP[1]=1;INYP[1]=1;INZP[1]=0;
  INXP[2]=0;INYP[2]=1;INZP[2]=0;
  INXP[3]=0;INYP[3]=0;INZP[3]=0;
  INXP[4]=0;INYP[4]=0;INZP[4]=1;
  INXP[5]=0;INYP[5]=1;INZP[5]=1;
  INXP[6]=0;INYP[6]=1;INZP[6]=1;
  INXP[7]=0;INYP[7]=0;INZP[7]=1;

  INXT[0]=1;INYT[0]=0;INZT[0]=0;
  INXT[1]=0;INYT[1]=1;INZT[1]=0;
  INXT[2]=0;INYT[2]=1;INZT[2]=0;
  INXT[3]=0;INYT[3]=0;INZT[3]=0;
  INXT[4]=0;INYT[4]=0;INZT[4]=1;
  INXT[5]=0;INYT[5]=0;INZT[5]=1;
  INXT[6]=0;INYT[6]=0;INZT[6]=1;
  INXT[7]=0;INYT[7]=0;INZT[7]=1;

  ELNDLK=malloc(10*run->topo->pleaves*sizeof(int *));
  BRANCHES=malloc(10*run->topo->pleaves*sizeof(struct BRANCH *));
  for(ielem=0;ielem<10*run->topo->pleaves;ielem++){
    ELNDLK[ielem]=malloc(8*sizeof(int));
    for(inode=0;inode<8;inode++){
      ELNDLK[ielem][inode]=-1;
    }
  }
  
  //message("MMFV: EXPL: STA1 \n",run);
  
  LKEL=malloc(8*run->topo->pleaves*sizeof(int *));
  LKND=malloc(8*run->topo->pleaves*sizeof(int *));
  LKNN=malloc(8*run->topo->pleaves*sizeof(int));
  for(ilink=0;ilink<8*run->topo->pleaves;ilink++){
    LKEL[ilink]=malloc(16*sizeof(int));
    LKND[ilink]=malloc(16*sizeof(int));
    LKNN[ilink]=0;
    for(ielem=0;ielem<16;ielem++){
      LKEL[ilink][ielem]=0;
      LKND[ilink][ielem]=0;
    }
  }
  
  //message("MMFV: EXPL: STA2 \n",run);
  // numnd=0;
  // numel=0;
  elcnt=0;
  for (iel=0;iel<run->topo->ntrees;iel++){ // Loop elements
    tree=run->topo->drys[iel];
    if (tree->part==run->con->rank){
    brch=tree->brch;
      while (brch->nkids!=0){
        brch=brch->kids[0];
      }
      while(brch!=NULL&&brch->root==iel){
        brch->el->gnum=elcnt;
        BRANCHES[elcnt]=brch;
        elcnt++;
        brch=brch->lnxt;
      }
    }
  }


  int adr[10];
  int itag;
  int ilvl;
  int tag;
  int nlvl;
  int indn;
  for (ipt=0;ipt<run->con->size;ipt++){ 
    for (itag=0;itag<run->topo->nbuff[ipt];itag++){
      iel=run->topo->buffdrys[ipt][itag];
      tag=run->topo->bufftags[ipt][itag];
      tag2adr(adr,tag,&nlvl);
      brch=run->topo->drys[iel]->brch;
      ilvl=0;
      while (brch->nkids!=0){
        brch=brch->kids[adr[ilvl]-1];ilvl++;
      }
      brch->el->gnum=elcnt;
      BRANCHES[elcnt]=brch;
      elcnt++;
    }
  }



  //message("MMFV: EXPL: STA3 \n",run);
  // numnd=0;
  // numel=0;
  elcnt=0;
  ilink=0;
  for (ipt=0;ipt<run->con->size;ipt++){ // Loop partitions
    if (ipt==run->con->rank){
      for (iel=0;iel<run->topo->ntrees;iel++){ // Loop elements
        tree=run->topo->drys[iel];
        if (tree->part==run->con->rank){
          brch=tree->brch;
          while (brch->nkids!=0){
            brch=brch->kids[0];
          }

          while(brch!=NULL&&brch->root==iel){
            elcnt++;
            for (ifc=0;ifc<brch->nlfc;ifc++){
              for (ind=0;ind<brch->cl->fc[ifc].nnds;ind++){
                node1=fcnd2elnd(ifc,ind,brch->type);
                elem1=brch->el->gnum;
                link1=ELNDLK[elem1][node1];
                //printf ("export ilink %d iel %d ifc %d ind %d el1 %d el2 %d lk1 %d lk2 %d nd1 %d nd2 %d \n",ilink,iel,ifc,ind,elem1,elem2,link1,link2,node1,node2);
                //if (brch->cl->fc[ifc].bc==0&&run->topo->drys[brch->neigtr[ifc][0]->root]->part==ipt&&brch->neigtr[ifc][0]->level==brch->level){
                indn=0;
                if (brch->neigtr[ifc][0]->level>brch->level){
                  if (SPLITMODE==8){indn=ind;}
                  if (SPLITMODE==4){
                    indn=-1;
                    for (ing=0;ing<brch->lfc_neigs[ifc];ing++){
                      if(brch->neignd[ifc][ing]==ind){
                        if (ing==0){indn=0;}
                        if (ing==1){indn=1;}
                        if (ing==2){indn=2;}
                        if (ing==3){indn=3;}
                      }
                    }
                  }
                }

                if (brch->cl->fc[ifc].bc==0&&brch->neigtr[ifc][0]->level==brch->level){     // 2D
                //if (brch->cl->fc[ifc].bc==0&&brch->neigtr[ifc][0]->level>=brch->level&&indn!=-1){   // 3D
                  ifc2=brch->neigfc[ifc];
                  iang=brch->neigag[ifc];
		  //indn=0;
		  //if (brch->neigtr[ifc][0]->level>brch->level){indn=ind;} // 3D
                  elem2=brch->neigtr[ifc][indn]->el->gnum;
                  node2=fcrefnd2elnd(ifc2,fcndrot(ind,iang,brch->cl->fc[ifc].type),brch->neigtr[ifc][indn]->type);
                  link2=ELNDLK[elem2][node2];
                  if (link1==-1&&link2==-1){ // elem1 node1 elem2 node2 get linked
                    ELNDLK[elem1][node1]=ilink;
                    ELNDLK[elem2][node2]=ilink;
                    LKEL[ilink][0]=elem1;
                    LKEL[ilink][1]=elem2;
                    LKND[ilink][0]=node1;
                    LKND[ilink][1]=node2;
                    LKNN[ilink]=2;
                    ilink++;
                  }
                  if (link1!=-1&&link2==-1){  // elem2 node2 inherits the link1
                    ELNDLK[elem2][node2]=link1;
                    LKEL[link1][LKNN[link1]]=elem2;
                    LKND[link1][LKNN[link1]]=node2;
                    LKNN[link1]++;
                  }
                  if (link1==-1&&link2!=-1){ // elem1 node1 inherits the link2
                    ELNDLK[elem1][node1]=link2;
                    LKEL[link2][LKNN[link2]]=elem1;
                    LKND[link2][LKNN[link2]]=node1;
                    LKNN[link2]++;
                  }
                  if (link1!=-1&&link2!=-1){ // elem1 node1 merge links with elem2 node2
                    if (link1!=link2){       // if they have different links
                      for(il=0;il<LKNN[link2];il++){
                        el2=LKEL[link2][il];
                        nd2=LKND[link2][il];
                        ELNDLK[el2][nd2]=link1;
                        LKEL[link1][LKNN[link1]+il]=el2;
                        LKND[link1][LKNN[link1]+il]=nd2;
                        LKEL[link2][il]=-1;
                        LKND[link2][il]=-1;
                      }
                      LKNN[link1]+=LKNN[link2];
                      LKNN[link2]=0;
                    }
                  }
     	       }
                else{
                  if (link1==-1){ // elem1 node1 elem2 node2 get linked
                    ELNDLK[elem1][node1]=ilink;
                    LKEL[ilink][0]=elem1;
                    LKND[ilink][0]=node1;
                    LKNN[ilink]=1;
                    ilink++;
                  }
                }
              }
            }
            //printf ("export il %d iel %d ifc %d ind %d lknn %d %d el1 %d el2 %d lk1 %d lk2 %d nd1 %d nd2 %d el2  %d nd2 %d \n",il,iel,ifc,ind,LKNN[link1],LKNN[link2],elem1,elem2,link1,link2,node1,node2,el2,nd2);
            brch=brch->lnxt;
          }
        }
      }
    } 
  }
  nlinks=0;
  //message("MMFV: EXPL: STA5 \n",run);
  for (il=0;il<8*run->topo->pleaves;il++){ // compact links
    //   for (iel=0;iel<8;iel++){ // Loop elements
    //     printf ("el  %d ",LKEL[il][iel]);
    //     printf ("nd  %d ",LKND[il][iel]);
    //   }
    //   printf (" \n ");
    if (LKNN[il]>0){
      for(in=0;in<LKNN[il];in++){
        LKEL[nlinks][in]=LKEL[il][in];
        LKND[nlinks][in]=LKND[il][in];
        ELNDLK[LKEL[il][in]][LKND[il][in]]=nlinks;
      }
      LKNN[nlinks]=LKNN[il];
      nlinks++;
    }
  }
  for (il=nlinks;il<8*run->topo->pleaves;il++){ // delete tail
    for(in=0;in<LKNN[il];in++){
      LKEL[il][in]=0;
      LKND[il][in]=0;
    }
    LKNN[il]=0;
  }
  
  nelems=run->topo->pleaves;
  
  MPI_Allreduce((&nlinks),(&nlinkstot),1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce((&nelems),(&nelemstot),1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  
  
  nllinks=nlinks;
  nelems=nelemstot;
  nlinks=nlinkstot;
  owrt=1;
  
  sprintf(fstring,"field_linked%08dprt%02d.dat", g_istep_tot,0); // Rank always 0


  if (run->con->rank==0){
       if (owrt==1) {tp=fopen(fstring, "w");}
       if (owrt==0) {tp=fopen(fstring, "a");}
       // Write premable
       //fprintf(tp," VARIABLES =  \"X\" \"Y\" \"Z\" \"iel\" \"tag\" \"gnum\" \"part\" \"lvl\" \"vol\" \"rlx\" \"crith\" \"rho\" ");
       fprintf(tp," VARIABLES =  \"X\" \"Y\" \"Z\" \"iel\" \"part\" \"rho\" ");
       for(imatv=0;imatv<eqtypn[0];imatv++){
         fprintf(tp," \"R%d\" ",imatv);
       }
       for(imatv=0;imatv<eqtypn[1];imatv++){
         fprintf(tp," \"U%d\" ",imatv);
       }
       fprintf(tp," \"E\" ");
       for(imatv=0;imatv<eqtypn[3];imatv++){
         fprintf(tp," \"a%d\" ",imatv);
       }
       
       for(imatv=0;imatv<eqtypn[4];imatv++){
         fprintf(tp," \"e%d\" ",imatv);
       }

       for(imatv=0;imatv<3;imatv++){
          for(jmatv=0;jmatv<3;jmatv++){
            fprintf(tp," \"Amat%d%d\" ",imatv,jmatv);
          }
        }
 
       for(imatv=0;imatv<eqtypn[0];imatv++){
         fprintf(tp," \"p%d\" ",imatv);
       }
       fprintf(tp," \"peq\" ");
       for(i=0;i<3;i++){ for(j=0;j<3;j++){     fprintf(tp," \"S%d%d\" ",i,j); }    }
       fprintf(tp," \"S_tens\" \"S_comp\" \"T_max\" \"SI\" \n ");
       fprintf(tp," \n ");

       // Write zone data
       fprintf(tp, "ZONE T=\"%d\",N=%d,E=%d,F=FEPOINT,ET=BRICK,DT=(double,double,double)\n",1,nlinks,nelems);
       fclose(tp);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  owrt=0;
  if (run->con->rank==0) {tp=fopen(fstring,"a");}
  nlinks=0;
  ELNDLKOF=0;
  for (ipt=0;ipt<run->con->size;ipt++){ // Loop partitions
    llinks=ELNDLKOF;
    ill=0;
    if (ipt==run->con->rank){

      //printf("export A %d %d %d %d %d \n",ipt,run->con->rank,nllinks,ill,NEQ);
      exprt_real=malloc(nll*nllinks*sizeof(double));
      exprt_ingr=malloc(6*nllinks*sizeof(int));
      ill=0;
      for (il=0;il<8*run->topo->pleaves;il++){ // Loop links
        if (LKNN[il]>0){
          for(iv=0;iv<NMATERIALS+1;iv++){press[iv]=0.0;}
          r=0.0; for (isp=0;isp<6;isp++){Y[isp]=0;}   // Zero
    	    for(iv=0;iv<(NEQ+(NMATERIALS+1)+9+4);iv++){var[iv]=0.0;}
          for (inode=0;inode<LKNN[il];inode++){ // Loop link nodes
            iel=LKEL[il][inode];
            brch=BRANCHES[iel];
            // Write results
            for(iv=0;iv<NEQ;iv++){var[iv]+=brch->el->S[iv]/LKNN[il];}   // Vars
            
            //=========================================================================================
            m5_facevalues(run,brch,0,0,0,0);
            
            press[NMATERIALS] = 0.0;
            for (iv=0;iv<NMATERIALS;++iv){
              press[iv] = run->sol_L->avf[iv]*run->sol_L->p_hydro[iv];
              press[NMATERIALS] += press[iv];
            }

            for(iv=0;iv<NMATERIALS+1; iv++){
              var[iv+NEQ] += press[iv]/LKNN[il];
            }

            for (i=0;i<3;++i){
              for (j=0;j<3;++j){
                Stress[i][j] = run->sol_L->st[i][j]; 
              }
            }
            
            S1 = ((Stress[0][0] + Stress[1][1]) / 2.0) + sqrt( pow(((Stress[0][0]-Stress[1][1])/2.0),2.0) + pow(Stress[0][1],2.0) );
            S2 = ((Stress[0][0] + Stress[1][1]) / 2.0) - sqrt( pow(((Stress[0][0]-Stress[1][1])/2.0),2.0) + pow(Stress[0][1],2.0) );
            S3 = Stress[2][2];

            S_tens =  max(S1,S3);
            S_comp = -min(S2,S3);

            max_si = max(max(S1,S3),max(S1,S2));
            min_si = min(min(S1,S3),min(S1,S2));

            T_max  = ( max_si - min_si ) / 2.0;

            iv=NEQ+NMATERIALS+1;
            for (i=0;i<3;++i){
              for (j=0;j<3;++j){
                var[iv]+=Stress[i][j]/LKNN[il];
                iv++;
              }
            } 

            var[NEQ+NMATERIALS+1+9+0] += S_tens/LKNN[il];
            var[NEQ+NMATERIALS+1+9+1] += S_comp/LKNN[il];
            var[NEQ+NMATERIALS+1+9+2] += T_max/LKNN[il];
            var[NEQ+NMATERIALS+1+9+3] += brch->el->SI[0]/LKNN[il];

            //=========================================================================================                                               
                                                                                     
          }
          inode=0;
          iel=LKEL[il][inode];
          in=LKND[il][inode];
          brch=BRANCHES[iel];
          tree=run->topo->drys[brch->root];
          type=brch->type;
          x=brch->cl->nd[LKND[il][inode]].x; y=brch->cl->nd[LKND[il][inode]].y; z=brch->cl->nd[LKND[il][inode]].z; // Coordinates
          exprt_real[ill*nll+0]=x;
          exprt_real[ill*nll+1]=y;
          exprt_real[ill*nll+2]=z;
          exprt_real[ill*nll+3]=brch->cl->Vol;
          exprt_real[ill*nll+4]=brch->el->crith;
          for(iv=0;iv<(NEQ+(NMATERIALS+1)+9+4);iv++){exprt_real[ill*nll+5+iv]=var[iv];}
          exprt_ingr[ill*6+0]=brch->root;    
          exprt_ingr[ill*6+1]=brch->tag;    
          exprt_ingr[ill*6+2]=brch->el->gnum;    
          exprt_ingr[ill*6+3]=brch->part;    
          exprt_ingr[ill*6+4]=brch->split;    
          exprt_ingr[ill*6+5]=0;
          llinks++;
          ill++;
        }
      }
      if (ipt!=0){ 
	//printf("export B %d %d %d %d \n",ipt,run->con->rank,nllinks,ill*NEQ);
        MPI_Send(&(ill),1,MPI_INT,0,0*nsz+ipt,MPI_COMM_WORLD);
        MPI_Send(exprt_real,ill*(5+(NEQ+(NMATERIALS+1)+9+4)),MPI_DOUBLE,0,1*nsz+ipt,MPI_COMM_WORLD);
        MPI_Send(exprt_ingr,ill*6,MPI_INT,0,2*nsz+ipt,MPI_COMM_WORLD);

        free(exprt_real);
        free(exprt_ingr);
      } 

    }
    if (run->con->rank==0){
      //printf("export C %d %d %d %d %d \n",ipt,run->con->rank,nllinks,ill,NEQ);
      if (ipt!=0){
        MPI_Recv(&ill,1,MPI_INT,ipt,0*nsz+ipt,MPI_COMM_WORLD,&status);
        //printf("export D %d %d %d %d %d \n",ipt,run->con->rank,nllinks,ill,NEQ);
        exprt_real=malloc(nll*ill*sizeof(double));
        exprt_ingr=malloc(6*ill*sizeof(int));
        MPI_Recv(exprt_real,ill*(5+(NEQ+(NMATERIALS+1)+9+4)),MPI_DOUBLE,ipt,1*nsz+ipt,MPI_COMM_WORLD,&status);
        MPI_Recv(exprt_ingr,ill*6,MPI_INT,ipt,2*nsz+ipt,MPI_COMM_WORLD,&status);
      }
      //printf("export E %d %d %d %d %d \n",ipt,run->con->rank,nllinks,ill,NEQ);
      for (il=0;il<ill;il++){ // Loop links
        
        //                           x                   ,y                   ,z                   ,iel               ,part
        fprintf(tp,"%e %e %e %d %d ",exprt_real[il*nll+0],exprt_real[il*nll+1],exprt_real[il*nll+2],exprt_ingr[il*6+0],exprt_ingr[il*6+3]); 

        rho_temp=0.0;
        for(iv=0;iv<(NMATERIALS);iv++) { rho_temp += exprt_real[nll*il+5+iv]; }
        fprintf(tp," %e ",rho_temp); 

        for(iv=eqtypi[1];iv<eqtypi[2];iv++) { exprt_real[nll*il+5+iv] = exprt_real[nll*il+5+iv]/rho_temp; }
  
        for(iv=0;iv<NEQ;iv++)                                           { fprintf(tp," %e ",exprt_real[nll*il+5+iv] ); }        
        for(iv=NEQ;iv<(NEQ+(NMATERIALS+1)+9+4);iv++)  { fprintf(tp," %e ",exprt_real[nll*il+5+iv] ); }

        //for(iv=0;iv<(NEQ+(NMATERIALS+1)+9+3+2);iv++) { fprintf(tp," %e ",exprt_real[nll*il+5+iv]); }

        fprintf(tp," \n");
        
        //for(iv=0;iv<NEQ;iv++) { fprintf(tp," %e ",exprt_real[nll*il+5+iv]); } for(iv=0;iv<NMATERIALS+1;iv++){ fprintf(tp," %e ",press[iv]); } fprintf(tp," \n");
      }
      free(exprt_real);
      free(exprt_ingr);
    }
    if ((run->con->rank==ipt)&&(ipt+1<run->con->size)){MPI_Send(&llinks,1,MPI_INT,ipt+1,3*nsz+ipt,MPI_COMM_WORLD);}
    if ((run->con->rank==ipt+1)&&(ipt+1<run->con->size)){MPI_Recv(&ELNDLKOF,1,MPI_INT,ipt,3*nsz+ipt,MPI_COMM_WORLD,&status);}
    MPI_Barrier(MPI_COMM_WORLD);
  }

  if (run->con->rank==0) {fclose(tp);}
  owrt=0;
  ignd=0;
  
  MPI_Barrier(MPI_COMM_WORLD);
  if (run->con->rank==0) {tp=fopen(fstring, "a");}
  ////message("MMFV: EXPL: STA6 \n",run);
  for (iel=0;iel<run->topo->ntrees;iel++){ // Loop elements
    //printf("nodes A %d %d %d %d \n",run->con->rank,iel,ill,NEQ);
    tree=run->topo->drys[iel];
    ill=0;
    if (tree->part==run->con->rank){
      brch=tree->brch;
      type=brch->type;
      while (brch->nkids!=0){
        brch=brch->kids[0];
      }
      ill=0;
      while(brch!=NULL&&brch->root==iel){
        for (ind=0;ind<8;ind++){
        }
        ill++;
        brch=brch->lnxt;
      }
      //printf("nodes B %d %d %d %d \n",run->con->rank,iel,ill,NEQ);
      glob_nodes=malloc(8*ill*sizeof(int));

      brch=tree->brch;
      type=brch->type;
      while (brch->nkids!=0){
        brch=brch->kids[0];
      }
      ill=0;
      while(brch!=NULL&&brch->root==iel){
        for (ind=0;ind<8;ind++){
          glob_nodes[8*ill+ind]=1+ELNDLKOF+ELNDLK[brch->el->gnum][TPNDS[ind][type]];
        }
        ill++;
        brch=brch->lnxt;
      }
    //printf("nodes C %d %d %d %d \n",run->con->rank,iel,ill,NEQ);
      if (tree->part!=0){
        MPI_Send(&ill,1,MPI_INT,0,nsz*10+iel,MPI_COMM_WORLD);
        MPI_Send(glob_nodes,ill*8,MPI_INT,0,nsz*10+run->topo->ntrees+iel,MPI_COMM_WORLD);
      }
    //printf("nodes D %d %d %d %d \n",run->con->rank,iel,ill,NEQ);
      if (run->con->rank!=0){
        free(glob_nodes);
      }
    }
    //printf("nodes E %d %d %d %d \n",run->con->rank,iel,ill,NEQ);
    if (run->con->rank==0){
      if (tree->part!=0){
        MPI_Recv(&ill,1,MPI_INT,tree->part,nsz*10+iel,MPI_COMM_WORLD,&status);
        glob_nodes=malloc(8*ill*sizeof(int));
        MPI_Recv(glob_nodes,ill*8,MPI_INT,tree->part,nsz*10+run->topo->ntrees+iel,MPI_COMM_WORLD,&status);
      }
    //printf("nodes F %d %d %d %d \n",run->con->rank,iel,ill,NEQ);

      for (il=0;il<ill;il++){
        for (ind=0;ind<8;ind++){
          fprintf(tp,"%d ",glob_nodes[8*il+ind]);
        }
        fprintf(tp,"\n");
      }
      free(glob_nodes);

    }
    //printf("nodes G %d %d %d %d \n",run->con->rank,iel,ill,NEQ);

    MPI_Barrier(MPI_COMM_WORLD);
  }
  if (run->con->rank==0) {fclose(tp);}

  //message("MMFV: EXPL: STA7 \n",run);
  for(ielem=0;ielem<10*run->topo->pleaves;ielem++){
    free(ELNDLK[ielem]);
  }
  free(ELNDLK);
  free(BRANCHES);
  
  
  for(ilink=0;ilink<8*run->topo->pleaves;ilink++){
    free(LKEL[ilink]);
    free(LKND[ilink]);
  }
  free(LKEL);
  free(LKND);
  free(LKNN);
}
