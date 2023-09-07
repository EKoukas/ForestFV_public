#include "strdata.h"

void m5_wave_speeds_gfm(struct RUN *run, struct BRANCH * crnt,int ifc,double* u_S,double* SL, double* SR,
                         double* uL, double* vL, double* wL,double* s11L,double* s12L,double* s13L,
												 double* uR, double* vR, double* wR,double* s11R,double* s12R,double* s13R) {
                                                                              

  int iv,i,j;
  double ruL,ruR;
  double c2YL,c2YR;
  double cL,cR;
  double * tensor_temp; 

  tensor_temp = malloc(9*sizeof(double)); 
  
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
    double are_elastic_sol,are_elastic_sol1;
    double Gmag_temp,Gmag_sol,Gmag_sol1;
    double rx,ry,rz;
    double Amat[3][3];
    
    rx=0.0;ry=0.0;rz=0.0;
    for (ind=0;ind<crnt->nlnd;ind++){
      rx+=crnt->cl->nd[ind].x/((double)crnt->nlnd);
      ry+=crnt->cl->nd[ind].y/((double)crnt->nlnd);
      rz+=crnt->cl->nd[ind].z/((double)crnt->nlnd);
    }

    printf("\n");
    printf("compute_wave_speeds V2: el %d | %d %d | %f %f \n",g_istep,crnt->root,ifc,run->sol_L->phi[0],run->sol_R->phi[0]);
    //printf("normal: %e %e %e \n",crnt->cl->nx[f],crnt->cl->ny[f],crnt->cl->nz[f]);
    printf("coordinates: %e %e %e \n",rx,ry,rz);
    printf("\n");

    for(iv=0;iv<NEQ;iv++){
      printf("vec: %d | %e %e \n",iv,run->sol_L->vec[iv],run->sol_R->vec[iv]);
    }
    printf("\n");

    v=0;
    for (i=0;i<3;++i){
      for (j=0;j<3;++j){
        Amat[i][j] = run->sol_L->vec[eqtypi[4]+v];
        v++;
      }
    }
    for (i=0;i<3;++i){
      for (j=0;j<3;++j){
        run->Gmat[i][j] = Amat[i][0]*Amat[j][0] + 
                          Amat[i][1]*Amat[j][1] + 
                          Amat[i][2]*Amat[j][2];
      }
    }

    Gmag_temp = G_det(run->Gmat); 
    printf("Gmag_temp: %e \n",Gmag_temp);

    printf("run->Gmat: \n");
    for (i=0;i<3;i++) {
      for (j=0;j<3;j++) {
        printf(" %e ", run->Gmat[i][j]);
      }
      printf(" \n");
    }
    printf(" \n");

    printf("Amat: \n");
    for (i=0;i<3;i++) {
      for (j=0;j<3;j++) {
        printf(" %e ", Amat[i][j]);
      }
      printf(" \n");
    }
    printf(" \n");

    are_elastic_sol=m5_elastic_energy_solid(run->Gmat,MATERMUSH[0]);
    //are_elastic_sol=are_elastic_sol*run->sol_L->ra[0];

    v=0;
    for (i=0;i<3;++i){
      for (j=0;j<3;++j){
        Amat[i][j] = run->sol_R->vec[eqtypi[4]+v];
        v++;
      }
    }
    for (i=0;i<3;++i){
      for (j=0;j<3;++j){
        run->Gmat[i][j] = Amat[i][0]*Amat[j][0] + 
                          Amat[i][1]*Amat[j][1] + 
                          Amat[i][2]*Amat[j][2];
      }
    }

    Gmag_temp = G_det(run->Gmat); 
    printf("Gmag_temp: %e \n",Gmag_temp);

    printf("run->Gmat: \n");
    for (i=0;i<3;i++) {
      for (j=0;j<3;j++) {
        printf(" %e ", run->Gmat[i][j]);
      }
      printf(" \n");
    }
    printf(" \n");

    printf("Amat: \n");
    for (i=0;i<3;i++) {
      for (j=0;j<3;j++) {
        printf(" %e ", Amat[i][j]);
      }
      printf(" \n");
    }
    printf(" \n");

    are_elastic_sol1=m5_elastic_energy_solid(run->Gmat,MATERMUSH[0]);
    //are_elastic_sol1=are_elastic_sol1*run->sol_R->ra[0];

    printf("are_elastic: %e %e \n",are_elastic_sol,are_elastic_sol1);
    printf("\n");

    printf("u_S SL SR %e %e %e \n",*u_S,*SL,*SR);
    
    for (iv=0; iv<NMATERIALS; ++iv) {
      printf("avf sol,sol1 %d | %e %e \n",iv,run->sol_L->avf[iv],run->sol_R->avf[iv]);
    }

    for (iv=0; iv<NMATERIALS; ++iv) {
      printf("c2Y sol,sol1 %d | %e %e | %e %e \n",  iv,run->sol_L->c2Y[iv],run->sol_R->c2Y[iv],run->sol_L->p_hydro[iv],run->sol_R->p_hydro[iv]);    
    }

    printf("cL,cR %e %e \n",cL,cR);
    printf("ruL,ruR %e %e \n",ruL,ruR);

    for (i=0;i<3;++i){
      for (j=0;j<3;++j){
        printf("st %d %d | %e %e | \n", i,j,run->sol_L->st[i][j],run->sol_R->st[i][j]);

      }
    }
    printf("sol->r,sol1->r %e %e \n",run->sol_L->r,run->sol_R->r);
    //printf(" t1     %e %e %e \n",crnt->el->nxt1[ifc],crnt->el->nyt1[ifc],crnt->el->nzt1[ifc]);
    //printf(" t2     %e %e %e \n",crnt->el->nxt2[ifc],crnt->el->nyt2[ifc],crnt->el->nzt2[ifc]);

    printf("*s11L,*s12L,*s13L %e %e %e \n",*s11L,*s12L,*s13L);
    printf("*s11R,*s12R,*s13R %e %e %e \n",*s11R,*s12R,*s13R);

    printf("uL,vL,wL  %e %e %e | %e %e %e \n",*uL,*vL,*wL,run->sol_L->u,run->sol_L->v,run->sol_L->w);
    printf("uR,vR,wR  %e %e %e | %e %e %e \n",*uR,*vR,*wR,run->sol_R->u,run->sol_R->v,run->sol_R->w);
    printf("\n");
    exit(0);
  }
  
  free(tensor_temp);
}