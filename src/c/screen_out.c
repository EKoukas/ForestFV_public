#include "strdata.h"
 
void screen_out(struct RUN * run) {

  double memory,memory_tot,check,time_remaining;
    
  int currRealMem,ierr; 
  int peakRealMem;
  int currVirtMem; 
  int peakVirtMem; 
  int cur_step,remaining_steps;

  MPI_Barrier(MPI_COMM_WORLD);
  // stores each word in status file
  char buffer[1024] = "";

  // linux file contaiNS this-process info
  FILE* file = fopen("/proc/self/status", "r");

  // read the entire file
  while (fscanf(file, " %1023s", buffer) == 1) {

      if (strcmp(buffer, "VmRSS:") == 0)  { ierr=fscanf(file, " %d", &currRealMem); }
      if (strcmp(buffer, "VmHWM:") == 0)  { ierr=fscanf(file, " %d", &peakRealMem); }
      if (strcmp(buffer, "VmSize:") == 0) { ierr=fscanf(file, " %d", &currVirtMem); }
      if (strcmp(buffer, "VmPeak:") == 0) { ierr=fscanf(file, " %d", &peakVirtMem); }
  }
  fclose(file);
  
  memory = ((double) currVirtMem)*0.001;

  
  MPI_Allreduce(&memory, &memory_tot, 1, MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);


  if (g_istep%1000==0) {

    if (g_istep==0) {
      time_remaining=0.0;
    } else {
      remaining_steps= NSTEP-g_istep;
      time_remaining=(((double)remaining_steps)*(run->cputime/g_istep))/3600.0;
    }

    if(run->con->rank==0){  printf ("\n");
      printf (">|________________________________________________________________________________________________________________________|\n"       );
      printf (">|__________________________________________________ ForestFV: CASE: %04d ________________________________________________|\n",CASE  );
      printf (">|________________________________________________________________________________________________________________________|\n"       );
      printf (">|________________   Istart|      Ifin|  Nexport|   Nadapt|        DT|  ORDER| GRAD SCHEME| EXPORT TYPE| MODEL| NBCOMMS|  |\n"       );
      printf (">|________________ %8d|  %8d| %8d| %8d|  %.2e|      %d|           %d|           %d|     %d|       %d|  |\n",
      g_istep_start,NSTEP+g_istep_start,EXPORT_STEP,ADAPTHSTEP,DT,ORDER,GRAD_SCHEME,EXPORT_TYPE,MODEL,NBCOMMS);
    
      printf (">|________________________________________________________________________________________________________________________|\n");
      printf (">|________________   Ntrees|   NLeaves|   NLEVEL| Nprocess| Memory usage (GB)| Expected time remaining (h)|               | \n");    
      printf (">|________________ %8d|  %8d| %8d| %8d|          %f|                       %5.2f|               |\n",
      run->topo->ntrees,run->tot_leaves,LEVEL,run->con->size,memory_tot/1000.0,time_remaining);
      printf (">|________________________________________________________________________________________________________________________|\n");
     // printf (">|Iiter  |Istep |Time     |Tcpu m| T/itr |RHS (%)|AMR (%)|WRO (%)| \n");
      printf (">|Iiter  |Istep |Sim time(s)|RT   (m)|T/100 (s)|AMR (s)|ADV (s)|2ND (s)|2COM (s)|RHS (s)|REL (s)|COM (s)|EXP (s)|MEM  (GB)| \n");
    }
    
  }

  if (g_istep==0) {

    if (run->con->rank==0) {

      printf (">I|%6d|%6d|  %.3e| %07.2f|  %07.3f| %06.2f| %06.2f| %06.2f|  %06.2f| %06.2f| %06.2f| %06.2f| %06.2f|%09.5f|\n",
      g_istep,g_istep_tot,g_time,
                                                            run->cputime/60.0,run->stepcputime,
                                                            run->Tmeshadapt,
                                                            run->Tadv,
                                                            run->T2nd,
                                                            run->T2ndcom,
                                                            run->Trhs,
                                                            run->Trelax,
                                                            run->Tcomm,
                                                            run->Tepx,memory_tot/1000.0);
    }
  } else {
    if (run->con->rank==0) {
      printf (">T|%6d|%6d|  %.3e| %07.2f|  %07.3f| %06.2f| %06.2f| %06.2f|  %06.2f| %06.2f| %06.2f| %06.2f| %06.2f|%09.5f|    \n",
      g_istep,g_istep_tot,g_time,
                                                            run->cputime/60.0,run->stepcputime_sum,
                                                            run->Tmeshadapt_sum,
                                                            run->Tadv_sum,
                                                            run->T2nd_sum,
                                                            run->T2ndcom_sum,
                                                            run->Trhs_sum,
                                                            run->Trelax_sum,
                                                            run->Tcomm_sum,
                                                            run->Tepx_sum,memory_tot/1000.0);


    }

  }
  
  
  run->stepcputime_sum = 0.0;
  run->Tepx_sum        = 0.0;
  run->Tmeshadapt_sum  = 0.0;
  
  run->Ttot_sum        = 0.0;
  
  run->Tadv_sum        = 0.0;
  run->T2nd_sum        = 0.0;
  run->T2ndcom_sum     = 0.0;
  run->Trhs_sum        = 0.0;
  run->Trelax_sum      = 0.0;
  run->Tcomm_sum       = 0.0;
  run->Tinterfloc_sum  = 0.0;
  

  run->Tdist_sum     = 0.0;
  run->Tdelta_sum    = 0.0;
  run->Trec_sum      = 0.0;
  run->T2ndrelax_sum = 0.0;
  run->T2ndface_sum  = 0.0;
  

}