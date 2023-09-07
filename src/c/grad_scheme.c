#include "strdata.h"

void grad_scheme (struct RUN* run) {

  switch(GRAD_SCHEME) {

    case 0: // structured reconstruction
      if (ORDER==2) { face_interpolation_V0(run);    } // SF: i-1, i , i+1, old slow method
      if (NS==1)    { green_gauss_grad(run);         } // WG: Green-Gauss cell based, 'linear' approximation
    break;

    case 1: // structured reconstruction
      if (ORDER==2) { face_interpolation_V1(run);    } // SF: i-1, i , i+1, fast
      if (NS==1)    { green_gauss_grad(run);         } // WG: Green-Gauss cell based, 'linear' approximation
    break;

    case 2: // structured reconstruction
      if (ORDER==2) { face_interpolation_V2(run);    } // SF: i-1, i , i+1, fast with linear interp for dist 
      if (NS==1)    { green_gauss_grad(run);         } // WG: Green-Gauss cell based, 'linear' approximation
    break;

    case 3:  // unstructured reconstruction
      if (ORDER==2) { face_interpolation_unstr(run); } // SF: i-1, i , i+1 , projection method
      if (NS==1)    { green_gauss_grad(run);         } // WG: Green-Gauss cell based, 'linear' approximation
    break;

    case 4: // Green-Gauss with weighted linear interpolation
      green_gauss_grad(run);   // WG: Green-Gauss cell based 
      if (ORDER==2) { face_interpolation_WG(run); }
    break;

    case 5: // Least-squares with weighted distances
      least_squares_grad(run);  
      if (ORDER==2) { face_interpolation_WG(run); }   
    break;

  }

  g_mesh_change_grad_scheme=0;

}



void grad_scheme_nb (struct RUN* run) {

  switch(GRAD_SCHEME) {

    case 0: // structured reconstruction
      if (ORDER==2) { face_interpolation_V0_nb(run);   } // SF: i-1, i , i+1, dimensional splitting, old slow method
      if (NS==1)    { green_gauss_grad(run);           } // WG: Green-Gauss cell based, 'linear' approximation
    break;

    case 1: // structured reconstruction
      if (ORDER==2) { face_interpolation_V1_nb(run);   } // SF: i-1, i , i+1, dimensional splitting, fast
      if (NS==1)    { green_gauss_grad(run);           } // WG: Green-Gauss cell based, 'linear' approximation
    break;

    case 2: // structured reconstruction
      if (ORDER==2) { face_interpolation_V2_nb(run);   } // SF: i-1, i , i+1, dimensional splitting, fast with linear interp for dist
      if (NS==1)    { green_gauss_grad(run);           } // WG: Green-Gauss cell based, 'linear' approximation
    break;

    case 3:  // unstructured reconstruction
      if (ORDER==2) { face_interpolation_unstr_nb(run);    } // SF: i-1, i , i+1 , dimensional splitting with projection
      if (NS==1)    { green_gauss_grad(run);            } // WG: Green-Gauss cell based, 'linear' approximation
    break;


    case 4: // Green-Gauss with weighted linear interpolation
      //if ((ORDER==2)||(NS==1)) { green_gauss_grad(run);  }  // WG: Green-Gauss cell based -> used for SF
    break;


    case 5: // Least-squares with weighted distances
      //if ((ORDER==2)||(NS==1)) { least_squares(run);     }  // WG: Least-squares -> used for SF
    break;

  }

  g_mesh_change_grad_scheme=0;

}