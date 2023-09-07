#include "strdata.h"

void m4_exportfield_1D(struct RUN * run) {
  
  int iel,iv;
  int type,owrt;
  struct BRANCH * brch;
  struct TREE * tree;
  char  fstring2[200]="                                             ";
  
  FILE * Res;

  owrt=1;
  sprintf(fstring2,"Res_gnuplot.dat");   
  
  if (run->con->rank==0) {
    if (owrt==1) {Res=fopen(fstring2, "w");}
    if (owrt==0) {Res=fopen(fstring2, "a");}
  }

  MPI_Barrier(MPI_COMM_WORLD);
  owrt=0;
  for (iel=0;iel<run->topo->ntrees;iel++){ // Loop elements
    tree=run->topo->drys[iel];
    brch=tree->brch;
    if (tree->part==run->con->rank){
      while (brch->nkids!=0){
        brch=brch->kids[0];
      }

      if (owrt==1) {Res=fopen(fstring2, "w");}
      if (owrt==0) {Res=fopen(fstring2, "a");}
      
      while(brch!=NULL&&brch->root==iel){

        type=brch->type;
        m4_facevalues(run,brch,0,0,0,0);        

        fprintf(Res,"%d ",1000-brch->root);
        
        for (iv=0;iv<eqtypn[0];++iv) { fprintf(Res,"%e ",      run->sol_L->ra[iv]                       ); } // ra
                                       fprintf(Res,"%e %e %e ",run->sol_L->u,run->sol_L->v,run->sol_L->w);   // u,v,w
                                       fprintf(Res,"%e ",      run->sol_L->e                            );   // E
        for (iv=0;iv<eqtypn[3];++iv) { fprintf(Res,"%e ",      run->sol_L->avf[iv]                      ); } // avf
        for (iv=0;iv<eqtypn[4];++iv) { fprintf(Res,"%e ",      run->sol_L->are[iv]                      ); } // are
        for (iv=0;iv<eqtypn[4];++iv) { fprintf(Res,"%e ",      run->sol_L->p_hydro[iv]                  ); } // p_hydro
        

        //for (iv=0;iv<NEQ;++iv) { fprintf(Res,"%e ",brch->el->S[iv]); } 

        fprintf(Res," \n");
 
        brch=brch->lnxt;
       
      } 
      fclose(Res); 
    } 
    MPI_Barrier(MPI_COMM_WORLD);
  }
 
}