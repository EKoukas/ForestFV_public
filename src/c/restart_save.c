#include "strdata.h"

void restart_save (struct RUN * run) {

  int n,iel;
  int owrt;
  struct BRANCH * brch;
  struct TREE * tree;
  char  fstring[100]="                      ";
  FILE * rs;

  sprintf(fstring,"restart%08d.dat", g_istep_tot);

  rs=fopen(fstring, "w");
  MPI_Barrier(MPI_COMM_WORLD);
 if(run->con->rank==0){
   fprintf(rs,"%d %15.10e %d %d %d \n",g_istep_tot,g_time,run->topo->tleaves,run->topo->ntrees,NEQ);
 }
 MPI_Barrier(MPI_COMM_WORLD);
  fclose(rs);
  

 owrt=0;
 for (iel=0;iel<run->topo->ntrees;iel++) { // Loop elements

   tree=run->topo->drys[iel];
   brch=tree->brch;

   if (tree->part==run->con->rank) {

     while (brch->nkids!=0){
       brch=brch->kids[0];
     }

     
     rs=fopen(fstring, "a");

     while(brch!=NULL&&brch->root==iel){
       fprintf(rs,"%d %d ",brch->root,brch->tag);
        for(n=0;n<NEQ;n++){ 
          fprintf(rs," %15.10e ",brch->el->S[n]);
        }
        fprintf(rs,"\n");
      brch=brch->lnxt;
     }
     fclose(rs);
   }

   MPI_Barrier(MPI_COMM_WORLD);

 }



}
