#include "strdata.h"

void celldeallocation(struct BRANCH *crnt){

  int ifc,ind;

  free(crnt->cl->nx);
  free(crnt->cl->ny);
  free(crnt->cl->nz);
  free(crnt->cl->nxt1);
  free(crnt->cl->nyt1);
  free(crnt->cl->nzt1);
  free(crnt->cl->nxt2);
  free(crnt->cl->nyt2);
  free(crnt->cl->nzt2);

  free(crnt->cl->dist_cf);
  free(crnt->cl->dist_cc);

  for(ifc=0;ifc<18;ifc++){
    free(crnt->cl->Area[ifc]);
  }
  free(crnt->cl->Area);

  for(ifc=0;ifc<crnt->nlfc;ifc++){
    free(crnt->cl->vec_cf[ifc]);
  }
  free(crnt->cl->vec_cf);

  free(crnt->cl->fc);
  crnt->cl->fc=NULL;

  free(crnt->cl->nd);
  crnt->cl->nd=NULL;

  free(crnt->cl);
  crnt->cl=NULL;

}