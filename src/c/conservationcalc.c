#include "strdata.h"

void conservationcalc(RUN *run){
char  fstring[50]="                    ";
int iv;
double coNS[50];
double coNStot[50];
FILE *fp;
struct BRANCH * crnt;

  for(iv=0;iv<50;iv++) {
    coNS[iv]=0.0;
  }
 
  crnt=run->topo->locl;
  while (crnt!=NULL){ 
    for(iv=0;iv<NEQ;iv++){
      coNS[iv]+=crnt->el->S[iv]*crnt->cl->Vol;    
    }
    coNS[NEQ]+=crnt->cl->Vol;
    crnt=crnt->lnxt;
  }

  MPI_Allreduce((&coNS),(&coNStot),50,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  if(run->con->rank==0){
    sprintf(fstring,"coNServation.dat");
    fp=fopen(fstring,"a");
    fprintf(fp," %d %le   ",g_istep_tot,g_time);
    for(iv=0;iv<NEQ+1;iv++){
      fprintf(fp," %le   ",coNStot[iv]);
    }
    fprintf(fp," \n");
    fclose(fp);
  }

}