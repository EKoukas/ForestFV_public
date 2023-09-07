#include "strdata.h"

void leafdeallocation(RUN * run,struct BRANCH *crnt){

  register int i;
  int m,x,y,f,j,k,isp;

  free(crnt->el->SI);
  free(crnt->el->SN);
  free (crnt->el->S);
  
  
  for(m=0;m<6;m++){
    free(crnt->el->SF[m]);
    free(crnt->el->flux_face[m]);
  }
  for(m=0;m<3;m++){
    free(crnt->el->WG[m]);
  }
  
  free(crnt->el->SF);
  free(crnt->el->flux_face);
  free(crnt->el->WG);


  for (m=0;m<6;m++) {
    free(crnt->el->u_face[m]);
  }
  free(crnt->el->u_face);
  free(crnt->el->flux_flag);
  
  free(crnt->el->RHS);
  free(crnt->el);
	crnt->el=NULL;

}
