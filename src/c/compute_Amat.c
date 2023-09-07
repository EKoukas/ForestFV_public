#include "strdata.h"

void compute_Amat (double * Amat,double asinit) {
  
  if (asinit>10.0*AMIN) {
    Amat[0] = 1.0; // A1
    Amat[1] = 0.0; // A2
    Amat[2] = 0.0; // A3

    Amat[3] = 0.0; // B1
    Amat[4] = 1.0; // B2
    Amat[5] = 0.0; // B3

    Amat[6] = 0.0; // C1
    Amat[7] = 0.0; // C2
    Amat[8] = 1.0; // C3

  } else if (asinit<10.0*AMIN) {
    Amat[0] = 0.0; // A1
    Amat[1] = 0.0; // A2
    Amat[2] = 0.0; // A3

    Amat[3] = 0.0; // B1
    Amat[4] = 0.0; // B2
    Amat[5] = 0.0; // B3

    Amat[6] = 0.0; // C1
    Amat[7] = 0.0; // C2
    Amat[8] = 0.0; // C3
  }
  
}