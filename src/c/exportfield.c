#include "strdata.h"

void exportfield(struct RUN * run) {

  run->Tepx=timecpu(run->Tepx,0);

  switch(EXPORT_TYPE){
    case 0: break;
    case 1: exportfield_1D(run);      break;
    case 2: exportfield_linked(run);  break;
    case 3: export_ultra_fast(run);   break;
  }

  run->Tepx=timecpu(run->Tepx,1);
  
}

void exportfield_1D(struct RUN *run) {

  switch(MODEL){
    case 0: m0_exportfield_1D(run);  break;
    case 1: m1_exportfield_1D(run);  break;
    case 2: m2_exportfield_1D(run);  break;
    case 3: m3_exportfield_1D(run);  break;
    case 4: m4_exportfield_1D(run);  break;
    case 5: m5_exportfield_1D(run);  break;
  }

   
}

void exportfield_linked(struct RUN *run) {

  switch(MODEL){
    case 0: m0_exportfield_linked(run);  break;
    //case 1: m1_exportfield_linked(run);  break;
    //case 2: m2_exportfield_linked(run);  break;
    //case 3: m3_exportfield_linked(run);  break;
    case 4: m4_exportfield_linked(run);  break;
    case 5: m5_exportfield_linked(run);  break;
  }
 
}

void export_ultra_fast (struct RUN * run) {

  int n,iel;
  int owrt;
  struct BRANCH * brch;
  struct TREE * tree;
  char  fstring[100]="                      ";
	char  command[100]="                      ";
  FILE * rs;

	
	if (run->con->rank==0) {
		sprintf(command,"mkdir F_%08d",g_istep_tot);
		system(command);
	}
	MPI_Barrier(MPI_COMM_WORLD);


  sprintf(fstring,"F_%08d/field_%04d_%08d.dat",g_istep_tot,run->con->rank,g_istep_tot);
  rs=fopen(fstring, "w");

  fprintf(rs,"%d %15.10e %d %d %d \n",g_istep_tot,g_time,run->topo->tleaves,run->topo->ntrees,NEQ);
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
	
	}

}