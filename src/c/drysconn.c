#include "strdata.h"

void drysconn (FOREST * topo) {
  
  int ignd,ilnd,iel,iel1,iel2,ifc,ind;
  int ** NDCON;
  int *  NDELS;
  int *  ELLIST;

  NDCON=malloc(topo->ngnd*sizeof(int *));
  NDELS=malloc(topo->ngnd*sizeof(int));
  ELLIST=malloc(topo->ngnd*sizeof(int));

  for (ignd=0;ignd<topo->ngnd;ignd++){
    NDELS[ignd]=0;
  }
  for (iel=0;iel<topo->nleaves;iel++){
    for (ilnd=0;ilnd<numlnd(topo->drys[iel]->type);ilnd++){
      NDELS[topo->drys[iel]->vrtx[ilnd]]++; 
    } 
  }

  for (ignd=0;ignd<topo->ngnd;ignd++){ 
    NDCON[ignd]=malloc(NDELS[ignd]*sizeof(int));
  }
  for (ignd=0;ignd<topo->ngnd;ignd++){
    NDELS[ignd]=0;
  }

  for (iel=0;iel<topo->nleaves;iel++){
    for (ilnd=0;ilnd<numlnd(topo->drys[iel]->type);ilnd++){
      ignd=topo->drys[iel]->vrtx[ilnd];
      NDCON[ignd][NDELS[ignd]]=iel;
      NDELS[ignd]++;
    }
  }

  for (iel=0;iel<topo->nleaves;iel++){
    MPI_Barrier(MPI_COMM_WORLD);
    for (ifc=0;ifc<numlfc(topo->drys[iel]->type);ifc++){
      topo->drys[iel]->conf[ifc]=-1;

      for (iel1=0;iel1<NDELS[topo->drys[iel]->vrtx[fcnd2elnd(ifc,0,topo->drys[iel]->type)]];iel1++){
        ELLIST[iel1]=1;
        for (ind=1;ind<numfcnd(ifc,topo->drys[iel]->type);ind++){
          for (iel2=0;iel2<NDELS[topo->drys[iel]->vrtx[fcnd2elnd(ifc,ind,topo->drys[iel]->type)]];iel2++){
            if (NDCON[topo->drys[iel]->vrtx[fcnd2elnd(ifc,ind,topo->drys[iel]->type)]][iel2]==NDCON[topo->drys[iel]->vrtx[fcnd2elnd(ifc,0,topo->drys[iel]->type)]][iel1]){ELLIST[iel1]++;}
          } 
        }
      }

      for (iel1=0;iel1<NDELS[topo->drys[iel]->vrtx[fcnd2elnd(ifc,0,topo->drys[iel]->type)]];iel1++){
        if ((ELLIST[iel1]==numfcnd(ifc,topo->drys[iel]->type))&&(NDCON[topo->drys[iel]->vrtx[fcnd2elnd(ifc,0,topo->drys[iel]->type)]][iel1]!=iel)){ 
          topo->drys[iel]->conf[ifc]=NDCON[topo->drys[iel]->vrtx[fcnd2elnd(ifc,0,topo->drys[iel]->type)]][iel1]; 
        }
      }
    }
  }

  for (ignd=0;ignd<topo->ngnd;ignd++){
    free(NDCON[ignd]);
  }

  free (NDCON);
  free (NDELS);
  free (ELLIST);

}