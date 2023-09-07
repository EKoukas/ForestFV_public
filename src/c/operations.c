#include "strdata.h"

void Init_1D(double *a,int Nx) { //Set zero the vector a of dimeNSion dim

  int i;
  for(i=0;i<Nx;i++){
		a[i]=0.0;
	}
  return;

}

void Init_2D(double **a,int Nx,int Ny) { //Set zero all elements of matrix a with dimeNSion m X n

  int i,j;
  for(i=0;i<Nx;i++){
		for(j=0;j<Ny;j++){
			a[i][j]=0.0;
		}
	}
  return;

}

void Init_3D(double ***a,int Nx,int Ny,int Nz) { //Set zero all elements of matrix a with dimeNSion m X n X k

  int i,j,k;
  for(i=0;i<Nx;i++){
		for(j=0;j<Ny;j++){
      for(k=0;k<Nz;k++){
			  a[i][j][k]=0.0;
      }
		}
	}
  return;

}