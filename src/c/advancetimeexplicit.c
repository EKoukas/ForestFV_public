#include "strdata.h"

void advancetimeexplicit(struct RUN * run) {

  int iter_RK;    

  run->Tadv=timecpu(run->Tadv,0); 

  RK(0,run);   // Initialize Runge-Kutta
  
  for (iter_RK=1;iter_RK<RK_STEPS+1;iter_RK++){  // Loop over the number of Runge-Kutta iterations

    grad_scheme(run); 
    RHS(run);         

    if   (MODEL==5){ m5_RK(iter_RK,run); } 
    else {              RK(iter_RK,run); }

    if      (MODEL==4){ m4_relaxation_fluid(run); }
    else if (MODEL==5){ m5_relaxation_solid(run); }
    
      run->Tcomm=timecpu(run->Tcomm,0);   
    communicate_S(run,0); // Perform communication
      run->Tcomm=timecpu(run->Tcomm,1);  
  }
  
  run->Tadv=timecpu(run->Tadv,1); 

}