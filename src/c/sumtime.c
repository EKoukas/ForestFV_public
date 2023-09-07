#include "strdata.h"
 
void sumtime(struct RUN * run) {

  run->stepcputime  = MPI_Wtime()-run->stepcputime;   // Time of iteration
  run->cputime     += run->stepcputime;               // Totall time of simulation

  run->stepcputime_sum += run->stepcputime;
  run->Ttot_sum        += run->Ttot;
  run->Tmeshadapt_sum  += run->Tmeshadapt;
  run->Tadv_sum        += run->Tadv;
  run->T2nd_sum        += run->T2nd; 
  run->T2ndcom_sum     += run->T2ndcom;
  run->Trhs_sum        += run->Trhs;
  run->Trelax_sum      += run->Trelax;
  run->Tcomm_sum       += run->Tcomm;
  run->Tinterfloc_sum  += run->Tinterfloc;
  run->Tepx_sum        += run->Tepx;
    
}