#include "strdata.h"

void meshadaptation (struct RUN * run) {

  run->Tmeshadapt=timecpu(run->Tmeshadapt,0);

  int ipass;
  FILE *fp;
  char  fstring[50];
  
  g_mesh_change_grad_scheme=1;
  
  for(ipass=0;ipass<(LEVEL+1);ipass++) {  // Maybe not be needed

    criterionsplt(run);     // Sets target level for refinement for all cells (brch->split)
    criterionsmth(run);     // Smooth levels (communication required)
    refineh(run,ipass%2);   // Split & Merge cells
    
    restructtopo(run);

    meshproperties(run);
   
    communicate_S(run,0);
    
    if (CRT_ADAPT_INTER_0!=0) {interface_location(run);}

  }

  if (NBCOMMS==1) {
    inner_el(run); 
    buffer_size_free(run); 
    buffer_size_malloc(run);
  }
  
  if(run->con->rank==0){
    sprintf(fstring,"tot_cells.dat");
    fp=fopen(fstring,"a"); 
    fprintf(fp," %d %le %d %d \n",g_istep_tot,g_time,run->topo->ntrees,run->tot_leaves); 
    fclose(fp);
  }
 
  run->Tmeshadapt=timecpu(run->Tmeshadapt,1);
  
}