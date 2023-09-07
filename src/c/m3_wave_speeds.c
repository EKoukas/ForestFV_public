#include "strdata.h"

void m3_wave_speeds(struct RUN *run, struct BRANCH * crnt,int ifc,double* u_S,double* SL, double* SR,
                         double* uL, double* vL, double* wL,double* s11L,double* s12L,double* s13L,
												 double* uR, double* vR, double* wR,double* s11R,double* s12R,double* s13R) {
                                                                              

  int iv;
  double ruL,ruR;
  double c2YL,c2YR;
  double cL,cR;
  double u,v,w;
  double nx,mx,lx,ny,my,ly,nz,mz,lz;
  
  double * tensor_temp; 
  
  tensor_temp = malloc(9*sizeof(double));

  nx=crnt->cl->nx[ifc]; mx=crnt->cl->nxt1[ifc]; lx=crnt->cl->nxt2[ifc];
  ny=crnt->cl->ny[ifc]; my=crnt->cl->nyt1[ifc]; ly=crnt->cl->nyt2[ifc];
  nz=crnt->cl->nz[ifc]; mz=crnt->cl->nzt1[ifc]; lz=crnt->cl->nzt2[ifc];


  // Left, Right U_hats vectors states 
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
  c2YL = 0.0;
  c2YR = 0.0;
  for (iv=0; iv<NMATERIALS; ++iv) {
    c2YL += run->sol_L->c2Y[iv];       //  magnitude of sound speed, element
    c2YR += run->sol_R->c2Y[iv];      //  magnitude of sound speed, neighbor
  }
  cL = sqrt(c2YL);
  cR = sqrt(c2YR);

  *SR = max( (*uL+cL) , (*uR+cR) );   // Right wave speed
  *SL = min( (*uL-cL) , (*uR-cR) );   // Left  wave speed 
    
  // Rotation of st tensors 
  rotate_tensor(tensor_temp,run->sol_L->st,crnt->cl->nx[ifc]  ,crnt->cl->ny[ifc]  ,crnt->cl->nz[ifc],
                                         crnt->cl->nxt1[ifc],crnt->cl->nyt1[ifc],crnt->cl->nzt1[ifc],
                                         crnt->cl->nxt2[ifc],crnt->cl->nyt2[ifc],crnt->cl->nzt2[ifc],
                                         VERBOSE,run->con->rank);

  *s11L = tensor_temp[0];
  *s12L = tensor_temp[1];
  *s13L = tensor_temp[2];

  rotate_tensor(tensor_temp,run->sol_R->st,crnt->cl->nx[ifc]  ,crnt->cl->ny[ifc]  ,crnt->cl->nz[ifc],
                                          crnt->cl->nxt1[ifc],crnt->cl->nyt1[ifc],crnt->cl->nzt1[ifc],
                                          crnt->cl->nxt2[ifc],crnt->cl->nyt2[ifc],crnt->cl->nzt2[ifc],
                                          VERBOSE,run->con->rank);
  *s11R = tensor_temp[0];
  *s12R = tensor_temp[1];
  *s13R = tensor_temp[2];
    
  // Intermediate wave speed
  *u_S=((run->sol_L->r*pow ((*uL),2.0)- (*s11L)) - (run->sol_R->r*pow( (*uR),2.0)- (*s11R)) - ( (*SL)*ruL) + ( (*SR)*ruR)) / 
       (ruL - ruR - ( (*SL)*run->sol_L->r) + ( (*SR)*run->sol_R->r));  
  
  free(tensor_temp);
}