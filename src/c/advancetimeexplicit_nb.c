#include "strdata.h"

void advancetimeexplicit_nb(struct RUN * run) {

  int iter_RK;  
  struct BRANCH * crnt;    

  run->Tavd_S=MPI_Wtime(); 
  run->Tadv=timecpu(run->Tadv,0); 

    run->Trk0_S=MPI_Wtime(); 
  RK(0,run);   // Initialize Runge-Kutta
    run->Trk0_E=MPI_Wtime();

  for (iter_RK=1;iter_RK<RK_STEPS+1;iter_RK++){  // Loop over the number of Runge-Kutta iterations

    grad_scheme_nb(run); 

      run->Trhs1_S=MPI_Wtime(); 
    RHS_nb(run,1);     // rhs for inner elements 
      run->Trhs1_E=MPI_Wtime();

    if (ORDER==2){
      communicate_nb_S(run,1,1); // Non blocking recv for SF  (II)
    } else if (((g_istep)%EXPORT_STEP!=0) || (g_istep==0)) {
        run->Tsrecv_S=MPI_Wtime();
      communicate_nb_S(run,0,1);  // Non blocking recv for S  (A) (B)
        run->Tsrecv_E=MPI_Wtime();
    }

      run->Trhs0_S=MPI_Wtime(); 
    RHS_nb(run,0);     // rhs for outer elements 
      run->Trhs0_E=MPI_Wtime();

      run->TrkN_S=MPI_Wtime();
    if   (MODEL==5){m5_RK(iter_RK,run); } 
    else {             RK(iter_RK,run); }
      run->TrkN_E=MPI_Wtime();

    if      (MODEL==4){m4_relaxation_fluid(run);}
    else if (MODEL==5){m5_relaxation_solid(run);}

      run->Tssend_S=MPI_Wtime();
    communicate_nb_S(run,0,0);  // Non blocking send for S  (B)  
      run->Tssend_E=MPI_Wtime();
  }

  if (((g_istep+1)%EXPORT_STEP==0) && (g_istep!=0)) {
    communicate_nb_S(run,0,1);  // Non blocking recv for S (B)
  }
  
  run->Tadv=timecpu(run->Tadv,1); 
  run->Tavd_E=MPI_Wtime();

}