#include "strdata.h"

void restart_load(RUN *run){

  int ierr,ileaf,n,istep_restart;
  int typ,rot,tag;
  int nleaves,ntrees,nv;
  FILE *rs;
  FILE *speciesfile;
  double buf,time_restart;
  struct TREE * tree;
  struct BRANCH * crnt;

  rs=fopen("Init","r");

  ierr=fscanf(rs,"%d  %lf %d %d %d \n",&istep_restart,&time_restart,&nleaves,&ntrees,&nv);
  g_istep_tot   = istep_restart;
  g_istep_start = istep_restart;
  g_time        = time_restart;
  
  for (ileaf=0;ileaf<nleaves;ileaf++){

    ierr=fscanf(rs,"%d %d ",&rot,&tag);
    tree=run->topo->drys[rot];

    if(tree->part==run->con->rank) {

      crnt=run->topo->drys[rot]->brch;
      typ=crnt->type;
      crnt=spawntree(run,tree,tag,rot,1);
      
      for(n=0;n<nv;n++){
        ierr=fscanf(rs,"%lf ",&crnt->el->S[n]);  
      }
    
    } else {

      for(n=0;n<nv;n++){
        ierr=fscanf(rs,"%lf ",&buf); 
      }

    }
    ierr=fscanf(rs,"\n");
    
 }
 fclose(rs);

}