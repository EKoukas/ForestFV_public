#include "strdata.h"

void loop (struct RUN * run) {
     
  run->cputime=0.0;  

  if (NBCOMMS==1) { communicate_nb_S(run,0,0); }  // Non blocking send for S (A)

  g_mesh_change_grad_scheme=1;
  MPI_Barrier(MPI_COMM_WORLD);

  for (g_istep=0;g_istep<(NSTEP+1);g_istep++) {

    run->Ttot_S=MPI_Wtime();
    init_clocks(run);
    
    run->Ttot=timecpu(run->Ttot,0);
   
    if ( ((g_istep+1)%ADAPTHSTEP==0) && (ADAPTH==1) ) { meshadaptation(run); }
    
    if      (NBCOMMS==0) { advancetimeexplicit(run);    }
    else if (NBCOMMS==1) { advancetimeexplicit_nb(run); }
    
    if ((CRT_ADAPT_INTER_0!=0) && (ADAPTH==1) ) {interface_location(run);}
    
    if (((g_istep+1)%EXPORT_STEP==0) && (g_istep!=0)) { exportfield(run);  }
    if ((g_istep+1)%RESTART_STEP==0)                  { restart_save(run); }
    
    run->Ttot=timecpu(run->Ttot,1);
  
    sumtime(run);
    
    run->Ttot_E=MPI_Wtime();
    //exp_timeline(run);

    if (g_istep%100==0) { screen_out(run); }

  }

  MPI_Barrier(MPI_COMM_WORLD);

}