#include "strdata.h"

void m1_wave_speeds(struct RUN *run, struct BRANCH * crnt,int ifc,double* S_S,double* SL, double* SR,
                    double* uL, double* vL, double* wL,double* s11L,double* uR, double* vR, double* wR,double* s11R) {
                                                                              

  int iv;
  double ruL,ruR;
  double cL,cR;
  double u,v,w;
  double nx,mx,lx,ny,my,ly,nz,mz,lz;
  double * tensor_temp; 
  
  tensor_temp = malloc(9*sizeof(double));

  nx=crnt->cl->nx[ifc]; 
  mx=crnt->cl->nxt1[ifc]; 
  lx=crnt->cl->nxt2[ifc];

  ny=crnt->cl->ny[ifc]; 
  my=crnt->cl->nyt1[ifc]; 
  ly=crnt->cl->nyt2[ifc];
  
  nz=crnt->cl->nz[ifc]; 
  mz=crnt->cl->nzt1[ifc]; 
  lz=crnt->cl->nzt2[ifc];


  u=run->sol_L->u;
  v=run->sol_L->v;
  w=run->sol_L->w;

  *uL = u*nx + v*ny + w*nz;
  *vL = u*mx + v*my + w*mz;
  *wL = u*lx + v*ly + w*lz;

  u=run->sol_R->u;
  v=run->sol_R->v;
  w=run->sol_R->w;

  *uR = u*nx + v*ny + w*nz;
  *vR = u*mx + v*my + w*mz;
  *wR = u*lx + v*ly + w*lz;
	

  ruL = run->sol_L->r* (*uL);
  ruR = run->sol_R->r* (*uR);
    
  // Wave speeds, ok 
  cL = run->sol_L->c;
  cR = run->sol_R->c;

  *SR = max( (*uL+cL) , (*uR+cR) );   // Right wave speed
  *SL = min( (*uL-cL) , (*uR-cR) );   // Left  wave speed 
    
  // Rotation of st tensors 
  rotate_tensor(tensor_temp,run->sol_L->st,crnt->cl->nx[ifc]  ,crnt->cl->ny[ifc]  ,crnt->cl->nz[ifc],
                                         crnt->cl->nxt1[ifc],crnt->cl->nyt1[ifc],crnt->cl->nzt1[ifc],
                                         crnt->cl->nxt2[ifc],crnt->cl->nyt2[ifc],crnt->cl->nzt2[ifc],
                                         VERBOSE,run->con->rank);

  *s11L = tensor_temp[0];

  rotate_tensor(tensor_temp,run->sol_R->st,crnt->cl->nx[ifc]  ,crnt->cl->ny[ifc]  ,crnt->cl->nz[ifc],
                                          crnt->cl->nxt1[ifc],crnt->cl->nyt1[ifc],crnt->cl->nzt1[ifc],
                                          crnt->cl->nxt2[ifc],crnt->cl->nyt2[ifc],crnt->cl->nzt2[ifc],
                                          VERBOSE,run->con->rank);
  *s11R = tensor_temp[0];

  *S_S=( (run->sol_L->r*pow((*uL),2.0) + (*s11L)) - (run->sol_R->r*pow((*uR),2.0) + (*s11R)) - ( (*SL)*ruL) + ( (*SR)*ruR)) / 
       (ruL - ruR - ( (*SL)*run->sol_L->r) + ( (*SR)*run->sol_R->r));  

  if (isnan(*S_S)==1) {
    printf("S_S nan | %e | %e %e | %e %e | %e %e | %e %e | %e %e \n",*S_S,*SL,*SR,cL,cR,(*uL),(*uR),(*s11L),(*s11R),run->sol_L->r,run->sol_R->r);
    exit(0);
  }

  free(tensor_temp);
}