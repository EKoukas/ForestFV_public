#include "strdata.h"

void crammer(int n,double ** a, double ** b,double ** x,double ** c) {
	
	int i,ie;
	double det,det0,det1,det2;
	
	det=determinant(a);

	for (ie=0;ie<n;ie++) { 

		for (i=0;i<3;i++){ 
			b[i][0] = c[i][ie]; 
		 	b[i][1] = a[i][1];  	
			b[i][2] = a[i][2];  
		}	
		det0=determinant(b);

		for (i=0;i<3;i++){ 
			b[i][0] = a[i][0];  
			b[i][1] = c[i][ie]; 	
			b[i][2] = a[i][2];  
		}	
		det1=determinant(b);
		
		for (i=0;i<3;i++){ 
			b[i][0] = a[i][0];  	
			b[i][1] = a[i][1];  	
			b[i][2] = c[i][ie]; 
		}	
		det2=determinant(b);
		
		x[0][ie] = det0/det;
		x[1][ie] = det1/det;
		x[2][ie] = det2/det;
	}

}

double determinant(double **a) {
	
	double det;
  
	det =  a[0][0]*a[1][1]*a[2][2] - a[0][2]*a[1][1]*a[2][0]
      	+a[0][1]*a[1][2]*a[2][0] - a[2][1]*a[1][2]*a[0][0]
      	+a[0][2]*a[1][0]*a[2][1] - a[2][2]*a[1][0]*a[0][1];
  
	return det;

}