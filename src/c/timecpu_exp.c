#include "strdata.h"
 
void timecpu_exp(struct RUN * run) {
  
  FILE * cpu;
  char  fstring[200]="                                             ";
  
  if (g_istep_tot==1) {
    sprintf(fstring,"Tcpu%02d.dat",run->con->rank); 
    cpu=fopen(fstring, "w");
  } else {
    sprintf(fstring,"Tcpu%02d.dat",run->con->rank); 
    cpu=fopen(fstring, "a");
  }
  
  
  fprintf(cpu,"%d %d %f %f %f %f %f %f %f %f  \n",run->con->rank,g_istep,
                                                                 run->Ttot_sum,
                                                                 run->Tmeshadapt_sum,
                                                                 run->Tadv_sum,  
                                                                  run->T2nd_sum,
                                                                  run->Trhs_sum,
                                                                  run->Trelax_sum,
                                                                  run->Tcomm_sum,
                                                                 run->Tepx_sum);                                                                  
  fclose(cpu);
  
    
}