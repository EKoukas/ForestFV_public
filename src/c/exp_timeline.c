#include "strdata.h"

void exp_timeline(struct RUN * run) {

  char  fstring[200];
  
  FILE * Res;


  sprintf(fstring,"time_%04d.dat",run->con->rank);   
  
  if (g_istep==0) { Res=fopen(fstring, "w"); }
  else            { Res=fopen(fstring, "a"); }
 
  
  

  fprintf(Res,"%e %f \n",run->Ttot_S,(double)run->con->rank);
  fprintf(Res,"%e %f \n",run->Ttot_E,(double)run->con->rank);
  fprintf(Res,"\n");
  fprintf(Res,"%e %f \n",run->Ttot_S,(double)run->con->rank+0.4);
  fprintf(Res,"%e %f \n",run->Ttot_S,(double)run->con->rank-0.4);
  fprintf(Res,"\n");
  fprintf(Res,"%e %f \n",run->Ttot_E,(double)run->con->rank+0.4);
  fprintf(Res,"%e %f \n",run->Ttot_E,(double)run->con->rank-0.4);
  fprintf(Res,"\n");

  /*
  fprintf(Res,"%e %f \n",run->Tavd_S,(double)run->con->rank+0.3);
  fprintf(Res,"%e %f \n",run->Tavd_S,(double)run->con->rank-0.3);
  fprintf(Res,"\n");
  fprintf(Res,"%e %f \n",run->Tavd_E,(double)run->con->rank+0.3);
  fprintf(Res,"%e %f \n",run->Tavd_E,(double)run->con->rank-0.3);
  fprintf(Res,"\n");
  fprintf(Res,"%e %f \n",run->Tavd_S,(double)run->con->rank);
  fprintf(Res,"%e %f \n",run->Tavd_E,(double)run->con->rank);
  fprintf(Res,"\n");
  */

  fprintf(Res,"%e %f \n",run->Trk0_S,(double)run->con->rank+0.3);
  fprintf(Res,"%e %f \n",run->Trk0_S,(double)run->con->rank-0.3);
  fprintf(Res,"\n");
  fprintf(Res,"%e %f \n",run->Trk0_E,(double)run->con->rank+0.3);
  fprintf(Res,"%e %f \n",run->Trk0_E,(double)run->con->rank-0.3);
  fprintf(Res,"\n");
  fprintf(Res,"%e %f \n",run->Trk0_S,(double)run->con->rank);
  fprintf(Res,"%e %f \n",run->Trk0_E,(double)run->con->rank);
  fprintf(Res,"\n");

  
  fprintf(Res,"%f %f \n",run->Trhs1_S,(double)run->con->rank+0.25);
  fprintf(Res,"%f %f \n",run->Trhs1_S,(double)run->con->rank-0.25);
  fprintf(Res,"\n");
  fprintf(Res,"%f %f \n",run->Trhs1_E,(double)run->con->rank+0.25);
  fprintf(Res,"%f %f \n",run->Trhs1_E,(double)run->con->rank-0.25);
  fprintf(Res,"\n");
  fprintf(Res,"%f %f \n",run->Trhs1_S,(double)run->con->rank);
  fprintf(Res,"%f %f \n",run->Trhs1_E,(double)run->con->rank);
  fprintf(Res,"\n");

  
  fprintf(Res,"%f %f \n",run->Tsrecv_S,(double)run->con->rank+0.1);
  fprintf(Res,"%f %f \n",run->Tsrecv_S,(double)run->con->rank-0.1);
  fprintf(Res,"\n");
  fprintf(Res,"%f %f \n",run->Tsrecv_E,(double)run->con->rank+0.1);
  fprintf(Res,"%f %f \n",run->Tsrecv_E,(double)run->con->rank-0.1);
  fprintf(Res,"\n");
  fprintf(Res,"%f %f \n",run->Tsrecv_S,(double)run->con->rank);
  fprintf(Res,"%f %f \n",run->Tsrecv_E,(double)run->con->rank);
  fprintf(Res,"\n");

  fprintf(Res,"%f %f \n",run->Trhs0_S,(double)run->con->rank+0.25);
  fprintf(Res,"%f %f \n",run->Trhs0_S,(double)run->con->rank-0.25);
  fprintf(Res,"\n");
  fprintf(Res,"%f %f \n",run->Trhs0_E,(double)run->con->rank+0.25);
  fprintf(Res,"%f %f \n",run->Trhs0_E,(double)run->con->rank-0.25);
  fprintf(Res,"\n");
  fprintf(Res,"%f %f \n",run->Trhs0_S,(double)run->con->rank);
  fprintf(Res,"%f %f \n",run->Trhs0_E,(double)run->con->rank);
  fprintf(Res,"\n");

  
  fprintf(Res,"%f %f \n",run->TrkN_S,(double)run->con->rank+0.3);
  fprintf(Res,"%f %f \n",run->TrkN_S,(double)run->con->rank-0.3);
  fprintf(Res,"\n");
  fprintf(Res,"%f %f \n",run->TrkN_E,(double)run->con->rank+0.3);
  fprintf(Res,"%f %f \n",run->TrkN_E,(double)run->con->rank-0.3);
  fprintf(Res,"\n");
  fprintf(Res,"%f %f \n",run->TrkN_S,(double)run->con->rank);
  fprintf(Res,"%f %f \n",run->TrkN_E,(double)run->con->rank);
  fprintf(Res,"\n");

  
  fprintf(Res,"%f %f \n",run->Tssend_S,(double)run->con->rank+0.05);
  fprintf(Res,"%f %f \n",run->Tssend_S,(double)run->con->rank-0.05);
  fprintf(Res,"\n");
  fprintf(Res,"%f %f \n",run->Tssend_E,(double)run->con->rank+0.05);
  fprintf(Res,"%f %f \n",run->Tssend_E,(double)run->con->rank-0.05);
  fprintf(Res,"\n");
  fprintf(Res,"%f %f \n",run->Tssend_S,(double)run->con->rank);
  fprintf(Res,"%f %f \n",run->Tssend_E,(double)run->con->rank);
  fprintf(Res,"\n");
 

  /*
  if (g_istep==0) {
    fprintf(Res,"%f %f \n",0.0          ,(double)run->con->rank);
  }
  fprintf(Res,"%f %f \n",run->Ttot_sum,(double)run->con->rank);
  fprintf(Res,"\n");
  fprintf(Res,"%f %f \n",run->Ttot_sum,(double)run->con->rank+0.4);
  fprintf(Res,"%f %f \n",run->Ttot_sum,(double)run->con->rank-0.4);
  fprintf(Res,"\n");
  fprintf(Res,"%f %f \n",run->Trhs_sum,(double)run->con->rank+0.2);
  fprintf(Res,"%f %f \n",run->Trhs_sum,(double)run->con->rank-0.2);
  fprintf(Res,"\n");
  fprintf(Res,"%f %f \n",run->Ttot_sum,(double)run->con->rank);
  */

  fclose(Res);

}