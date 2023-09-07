#include "strdata.h"

void tagvector(RUN *run){

  struct BRANCH * crnt;

  crnt=run->topo->locl; 
  while (crnt!=NULL){
    tagvectorcalc(run,crnt);
    crnt=crnt->lnxt;
  }
}


void tagvectorcalc(RUN *run,struct BRANCH * crnt){

  int f,type;
  double norm,norm1;

  for(f=0;f<crnt->nlfc;f++){//face

    type=crnt->type;

    norm = sqrt( pow(((crnt->cl->nd[fcnd2elnd(f,1,type)].y - crnt->cl->nd[fcnd2elnd(f,0,type)].y)),2.0) + 
                 pow(((crnt->cl->nd[fcnd2elnd(f,1,type)].x - crnt->cl->nd[fcnd2elnd(f,0,type)].x)),2.0) + 
                 pow(((crnt->cl->nd[fcnd2elnd(f,1,type)].z - crnt->cl->nd[fcnd2elnd(f,0,type)].z)),2.0));

    crnt->cl->nxt1[f] = (crnt->cl->nd[fcnd2elnd(f,1,type)].x-crnt->cl->nd[fcnd2elnd(f,0,type)].x)/norm;
    crnt->cl->nyt1[f] = (crnt->cl->nd[fcnd2elnd(f,1,type)].y-crnt->cl->nd[fcnd2elnd(f,0,type)].y)/norm;
    crnt->cl->nzt1[f] = (crnt->cl->nd[fcnd2elnd(f,1,type)].z-crnt->cl->nd[fcnd2elnd(f,0,type)].z)/norm;

    norm1 = sqrt( pow((crnt->cl->ny[f]*crnt->cl->nzt1[f])-(crnt->cl->nyt1[f]*crnt->cl->nz[f]),2.0) + 
                  pow((crnt->cl->nx[f]*crnt->cl->nzt1[f])-(crnt->cl->nxt1[f]*crnt->cl->nz[f]),2.0) + 
                  pow((crnt->cl->nx[f]*crnt->cl->nyt1[f])-(crnt->cl->nxt1[f]*crnt->cl->ny[f]),2.0));

    crnt->cl->nxt2[f] =  (((crnt->cl->ny[f])*(crnt->cl->nzt1[f]))-((crnt->cl->nyt1[f])*(crnt->cl->nz[f])))/norm1;
    crnt->cl->nyt2[f] = -(((crnt->cl->nx[f])*(crnt->cl->nzt1[f]))-((crnt->cl->nxt1[f])*(crnt->cl->nz[f])))/norm1;
    crnt->cl->nzt2[f] =  (((crnt->cl->nx[f])*(crnt->cl->nyt1[f]))-((crnt->cl->nxt1[f])*(crnt->cl->ny[f])))/norm1;

  }

}