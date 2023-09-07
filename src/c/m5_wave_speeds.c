#include "strdata.h"

void m5_wave_speeds(struct RUN *run, struct BRANCH * crnt,int ifc,double* u_S,double* SL, double* SR,
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


  if ( (isnan(*u_S)==1) || (isnan(*SR)==1) || (isnan(*SL)==1) ) {
    
    int iv,i,j,v,ind;
    double rx,ry,rz;
    
    rx=0.0;ry=0.0;rz=0.0;
    for (ind=0;ind<crnt->nlnd;ind++){
      rx+=crnt->cl->nd[ind].x/((double)crnt->nlnd);
      ry+=crnt->cl->nd[ind].y/((double)crnt->nlnd);
      rz+=crnt->cl->nd[ind].z/((double)crnt->nlnd);
    }

    printf("\n");
    printf("m5_wave_speeds: el %d | %d %d | %f %f \n",g_istep,crnt->root,ifc,run->sol_L->phi[0],run->sol_R->phi[0]);
    //printf("normal: %e %e %e \n",crnt->cl->nx[f],crnt->cl->ny[f],crnt->cl->nz[f]);
    printf("coordinates: %e %e %e \n",rx,ry,rz);
    printf("\n");

    for(iv=0;iv<NEQ;iv++){
      printf("vec: %d | %e %e \n",iv,run->sol_L->vec[iv],run->sol_R->vec[iv]);
    }
    printf("\n");

    printf("u_S SL SR %e %e %e \n",*u_S,*SL,*SR);
    
    for (iv=0; iv<NMATERIALS; ++iv) {
      printf("avf sol_L,sol_R %d | %e %e \n",iv,run->sol_L->avf[iv],run->sol_R->avf[iv]);
    }

    for (iv=0; iv<NMATERIALS; ++iv) {
      printf("c2Y sol_L,sol_R %d | %e %e | %e %e \n",  iv,run->sol_L->c2Y[iv],run->sol_R->c2Y[iv],run->sol_L->p_hydro[iv],run->sol_R->p_hydro[iv]);    
    }

    printf("cL,cR %e %e \n",cL,cR);
    printf("ruL,ruR %e %e \n",ruL,ruR);

    for (i=0;i<3;++i){
      for (j=0;j<3;++j){
        printf("st %d %d | %e %e | \n", i,j,run->sol_L->st[i][j],run->sol_R->st[i][j]);

      }
    }
    printf("sol_L->r,sol_R->r %e %e \n",run->sol_L->r,run->sol_R->r);
    
    printf("*s11L,*s12L,*s13L %e %e %e \n",*s11L,*s12L,*s13L);
    printf("*s11R,*s12R,*s13R %e %e %e \n",*s11R,*s12R,*s13R);

    printf("uL,vL,wL  %e %e %e | %e %e %e \n",*uL,*vL,*wL,run->sol_L->u,run->sol_L->v,run->sol_L->w);
    printf("uR,vR,wR  %e %e %e | %e %e %e \n",*uR,*vR,*wR,run->sol_R->u,run->sol_R->v,run->sol_R->w);
    printf("\n");
    exit(0);
  }
  
  free(tensor_temp);
}