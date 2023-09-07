#include "strdata.h"
 
double timecpu(double Ttemp, int flag) {
  
  if (flag==0){
    Ttemp = MPI_Wtime(); 
  } else {
    Ttemp = MPI_Wtime()-Ttemp; 
  }

  return Ttemp;

}