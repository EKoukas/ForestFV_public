#include "strdata.h"

/*
 * Updates the boundary conditions of the branches.
 *
 * This function reads boundary conditions from four different files
 * and updates the boundary conditions.
 *
 */

void boundary(struct RUN *run) {
  
  int i,ifc,fc,ibc,el,iel,nfcbc,ierr;
  struct BRANCH * crnt;
  double ** bctrees;
  FILE * fp;

  bctrees=malloc(run->topo->ntrees*sizeof(double));
  for (iel=0;iel<run->topo->ntrees;iel++) { 
	   bctrees[iel]=malloc(6*sizeof(double));
      for (ifc=0;ifc<6;ifc++){
		   bctrees[iel][ifc]=0;
	   }
  }
  
  for (ibc=1;ibc<5;ibc++) {
    if      (ibc==1) {fp=fopen(run->con->filebr1,"r");}
    else if (ibc==2) {fp=fopen(run->con->filebr2,"r");}
    else if (ibc==3) {fp=fopen(run->con->filebr3,"r");}
    else if (ibc==4) {fp=fopen(run->con->filebr4,"r");}

    ierr=fscanf(fp,"%d \n",&nfcbc);
    for (i=0;i<nfcbc;i++){
      ierr=fscanf(fp,"%d\t%d \n",&el,&fc);
      bctrees[el-1][fc-1]=ibc;
    }
    rewind(fp);
    fclose(fp);

  }

  crnt=run->topo->locl; 
  while (crnt!=NULL) {
    for (ifc=0;ifc<crnt->nlfc;ifc++){
      crnt->cl->fc[ifc].bc=0;
      ibc=bctrees[crnt->root][ifc];
      if      (ibc==1){crnt->cl->fc[ifc].bc=1;}
      else if (ibc==2){crnt->cl->fc[ifc].bc=2;}
      else if (ibc==3){crnt->cl->fc[ifc].bc=3;}
      else if (ibc==4){crnt->cl->fc[ifc].bc=4;}
    }
    crnt=crnt->lnxt; 
  }

  for (iel=0;iel<run->topo->ntrees;iel++) { 
	   free(bctrees[iel]);
  }
  free(bctrees);

}