#include "strdata.h"

void cellallocation (struct BRANCH * brch) { 

  struct CELL * cl_temp;
  int ifc,ing,i;

  brch->cl=malloc(sizeof(CELL));

  cl_temp = brch->cl;

  cl_temp->nx = malloc(brch->nlfc*sizeof(double));
  cl_temp->ny = malloc(brch->nlfc*sizeof(double));
  cl_temp->nz = malloc(brch->nlfc*sizeof(double));

  cl_temp->nxt1 = malloc(brch->nlfc*sizeof(double));
  cl_temp->nyt1 = malloc(brch->nlfc*sizeof(double));
  cl_temp->nzt1 = malloc(brch->nlfc*sizeof(double));

  cl_temp->nxt2 = malloc(brch->nlfc*sizeof(double));
  cl_temp->nyt2 = malloc(brch->nlfc*sizeof(double));
  cl_temp->nzt2 = malloc(brch->nlfc*sizeof(double));

  cl_temp->fc  = malloc(brch->nlfc*sizeof(FC));
  cl_temp->nd  = malloc(brch->nlnd*sizeof(ND));

  cl_temp->dist_cf = malloc(brch->nlfc*sizeof(double));
  cl_temp->dist_cc = malloc(brch->nlfc*sizeof(double));
  
  cl_temp->vec_cf = malloc(brch->nlfc*sizeof(double *));
  for(i=0;i<brch->nlfc;i++){
    cl_temp->vec_cf[i] = malloc(3*sizeof(double));
  }

  cl_temp->Vol  = 0.0;
  cl_temp->Area = malloc(18*sizeof(double* ));
  for(ifc=0;ifc<18;ifc++){
    cl_temp->Area[ifc] = malloc(4*sizeof(double));
    for(ing=0;ing<4;ing++){
      cl_temp->Area[ifc][ing]=0.0;
    }
  }
         
}

