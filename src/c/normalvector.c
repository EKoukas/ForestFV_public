#include "strdata.h"

void normalvector(RUN *run){
  
  struct BRANCH * crnt;

  crnt=run->topo->locl; 
  while (crnt!=NULL){
    normalvectorcalc(run,crnt);
    crnt=crnt->lnxt;
  }

}

void normalvectorcalc(RUN *run,struct BRANCH * crnt){

  int f,type;
  double norm,xf;

  type=crnt->type;
  for(f=0;f<crnt->nlfc;f++){
    
    norm = sqrt(pow(( (crnt->cl->nd[fcnd2elnd(f,1,type)].y - crnt->cl->nd[fcnd2elnd(f,0,type)].y) * (crnt->cl->nd[fcnd2elnd(f,2,type)].z - crnt->cl->nd[fcnd2elnd(f,0,type)].z) - (crnt->cl->nd[fcnd2elnd(f,2,type)].y - crnt->cl->nd[fcnd2elnd(f,0,type)].y) * (crnt->cl->nd[fcnd2elnd(f,1,type)].z - crnt->cl->nd[fcnd2elnd(f,0,type)].z)),2.0) + 
                pow(( (crnt->cl->nd[fcnd2elnd(f,1,type)].x - crnt->cl->nd[fcnd2elnd(f,0,type)].x) * (crnt->cl->nd[fcnd2elnd(f,2,type)].z - crnt->cl->nd[fcnd2elnd(f,0,type)].z) - (crnt->cl->nd[fcnd2elnd(f,2,type)].x - crnt->cl->nd[fcnd2elnd(f,0,type)].x) * (crnt->cl->nd[fcnd2elnd(f,1,type)].z - crnt->cl->nd[fcnd2elnd(f,0,type)].z)),2.0) + 
                pow(( (crnt->cl->nd[fcnd2elnd(f,1,type)].x - crnt->cl->nd[fcnd2elnd(f,0,type)].x) * (crnt->cl->nd[fcnd2elnd(f,2,type)].y - crnt->cl->nd[fcnd2elnd(f,0,type)].y) - (crnt->cl->nd[fcnd2elnd(f,2,type)].x - crnt->cl->nd[fcnd2elnd(f,0,type)].x) * (crnt->cl->nd[fcnd2elnd(f,1,type)].y - crnt->cl->nd[fcnd2elnd(f,0,type)].y)),2.0));
    
    crnt->cl->nx[f] =  ((crnt->cl->nd[fcnd2elnd(f,1,type)].y - crnt->cl->nd[fcnd2elnd(f,0,type)].y) * (crnt->cl->nd[fcnd2elnd(f,2,type)].z - crnt->cl->nd[fcnd2elnd(f,0,type)].z) - (crnt->cl->nd[fcnd2elnd(f,2,type)].y - crnt->cl->nd[fcnd2elnd(f,0,type)].y) * (crnt->cl->nd[fcnd2elnd(f,1,type)].z - crnt->cl->nd[fcnd2elnd(f,0,type)].z))/norm;
    crnt->cl->ny[f] = -((crnt->cl->nd[fcnd2elnd(f,1,type)].x - crnt->cl->nd[fcnd2elnd(f,0,type)].x) * (crnt->cl->nd[fcnd2elnd(f,2,type)].z - crnt->cl->nd[fcnd2elnd(f,0,type)].z) - (crnt->cl->nd[fcnd2elnd(f,2,type)].x - crnt->cl->nd[fcnd2elnd(f,0,type)].x) * (crnt->cl->nd[fcnd2elnd(f,1,type)].z - crnt->cl->nd[fcnd2elnd(f,0,type)].z))/norm;
    crnt->cl->nz[f] =  ((crnt->cl->nd[fcnd2elnd(f,1,type)].x - crnt->cl->nd[fcnd2elnd(f,0,type)].x) * (crnt->cl->nd[fcnd2elnd(f,2,type)].y - crnt->cl->nd[fcnd2elnd(f,0,type)].y) - (crnt->cl->nd[fcnd2elnd(f,2,type)].x - crnt->cl->nd[fcnd2elnd(f,0,type)].x) * (crnt->cl->nd[fcnd2elnd(f,1,type)].y - crnt->cl->nd[fcnd2elnd(f,0,type)].y))/norm;

    if (crnt->cl->fc[f].bc!=0){

      if (BCOVERWRT[0]!=0){ if (crnt->cl->nx[f]<-0.6){crnt->cl->fc[f].bc=BCOVERWRT[0];} }
      if (BCOVERWRT[1]!=0){ if (crnt->cl->nx[f]> 0.6){crnt->cl->fc[f].bc=BCOVERWRT[1];} }

      if (BCOVERWRT[2]!=0){ if (crnt->cl->ny[f]<-0.6){crnt->cl->fc[f].bc=BCOVERWRT[2];} }
      if (BCOVERWRT[3]!=0){ if (crnt->cl->ny[f]> 0.6){crnt->cl->fc[f].bc=BCOVERWRT[3];} }
    
      if (BCOVERWRT[4]!=0){ if (crnt->cl->nz[f]<-0.6){crnt->cl->fc[f].bc=BCOVERWRT[4];} }
      if (BCOVERWRT[5]!=0){ if (crnt->cl->nz[f]> 0.6){crnt->cl->fc[f].bc=BCOVERWRT[5];} }     

    }
  }


}