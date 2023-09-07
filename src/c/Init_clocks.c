#include "strdata.h"

void init_clocks (struct RUN * run) {

  run->stepcputime=MPI_Wtime(); // Time of iteration, start clock

  run->Tmeshadapt = 0.0;
  run->T2nd       = 0.0;
  run->T2ndcom    = 0.0;
  run->Tepx       = 0.0;
  run->Tcomm      = 0.0;
  run->Tinterfloc = 0.0;
  run->Tadv       = 0.0;
  run->Trhs       = 0.0;
  run->Trelax     = 0.0;
  run->Ttot       = 0.0;

  if (g_istep==0) {

    run->stepcputime_sum = 0.0;
    run->Ttot_sum        = 0.0;
    run->Tmeshadapt_sum  = 0.0;
    run->Tadv_sum        = 0.0;
    run->T2nd_sum        = 0.0;
    run->T2ndcom_sum     = 0.0;
    run->Trhs_sum        = 0.0;
    run->Trelax_sum      = 0.0;
    run->Tcomm_sum       = 0.0;
    run->Tinterfloc_sum  = 0.0;
    run->Tepx_sum        = 0.0;

    run->Tdist_sum     = 0.0;
    run->Tdelta_sum    = 0.0;
    run->Trec_sum      = 0.0;
    run->T2ndrelax_sum = 0.0;
    run->T2ndface_sum  = 0.0;
  }

  g_time+=DT;    // Totall time of simulation
  g_istep_tot++; // time steping with added restart time

}