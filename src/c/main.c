/*
 *
 * ForestFV of oct-Trees
 *
 *  \/\/\/\/
 *   \/ / /
 *    \/ /
 *     \/
 * 
 * 	Models:
 * 		-0: Euler EquatioNS or NS      (  gas/liquid          )
 * 		-1: Barotropic model           (      liquid-vapor    )
 * 		-2: Two-fluid Barotropic model (  gas-liquid-vapor    )
 * 		-3: Kapilla's 5-equation model (  gas-liquid          )
 * 		-4: Saurel's model 6-equation  (N:gas-N:liquid        )
 * 		-5: Saurel's model 6-equation  (N:gas-N:liquid-N:solid)
 * 		-6 Barton's  model ()
 * 
*/

#include "strdata.h"

int main (int input, char **inputc) {

	int it,in,ing;

	int irank;
	int ipass;
	struct RUN  runv;
	struct RUN * run;
	struct BRANCH * crnt;
	int ifc;

	MPI_Init(&input,&inputc); 		
	PetscInitialize(&input,&inputc,0,0);

	run=&runv; 	
	createrunstruct(run);

	run->con->comm=MPI_COMM_WORLD;
	
	run->con->casename=inputc[1];
	run->con->runname=inputc[2];

	MPI_Comm_rank(MPI_COMM_WORLD,&run->con->rank);  
	MPI_Comm_size(MPI_COMM_WORLD,&run->con->size);  

	filenames(run); 	    
	MPI_Barrier(MPI_COMM_WORLD);

	global_variables_def(run);

	inputset(run); if((VERBOSE==1)&&(run->con->rank==0)){ printf ("ForestFV: MAIN: INPUT SET \n");}
	
	if(run->con->rank==0){
		if(TRANSLATOR==1) {translator(run); 	  printf ("ForestFV: MAIN: NEUTRAL FILE TRANSLATED. (GAMBIT) \n");}  
		if(TRANSLATOR==2) {translator_pw(run);  printf ("ForestFV: MAIN: NEUTRAL FILE TRANSLATED. (POINTWISE)\n");}        
	}
	MPI_Barrier(MPI_COMM_WORLD);

	for (irank=0;irank<run->con->size;irank++){	
	  if (run->con->rank==irank){
	    createtopo(run); 	// reads mesh files
	  }
	  MPI_Barrier(MPI_COMM_WORLD);
	}
	
	drysconn(run->topo);
	partitiondrys(run);
	spawntopo(run);
	forestconn(run);

	MPI_Barrier(MPI_COMM_WORLD);
	run->topo->partitioned=1;

	g_istep=0;
		
	for (irank=0;irank<run->con->size;irank++){
    if (run->con->rank==irank){
	    boundary(run);
	  }
	  MPI_Barrier(MPI_COMM_WORLD);
	}

	createutility(run);
	memallocation(run);

	g_time=0.0;

	meshproperties(run);

	if((VERBOSE==1)&&(run->con->rank==0)){ printf ("ForestFV: MAIN: restructtopo 0\n");}
	run->topo->partitioned=0;
	restructtopo(run);
	run->topo->partitioned=1;
	if((VERBOSE==1)&&(run->con->rank==0)){ printf ("ForestFV: MAIN: restructtopo 1\n");}

	meshproperties(run);

	if (RESTART==1){
		restart_load(run);
		restructtopo(run);   
    meshproperties(run);
		communicate_S(run,0);
	} else {

		if (GEO_ADAPTH==1){
			
			for(ipass=0;ipass<LEVEL;ipass++){
				
				if (run->con->rank==0){ printf("ForestFV: MAIN: Geo pass %d \n",ipass); }

				init_domain(run);
				communicate_S(run,0);

				if (CRT_ADAPT_INTER_0!=0) {interface_location(run);}			
				
				criterionsplt(run);

				//communicate_S(run,3);

				criterionsmth(run);
				refineh(run,ipass%2);
				restructtopo(run);				
				meshproperties(run);

			}

			init_domain(run);	
			communicate_S(run,0);
		} else {  
			init_domain(run);
			communicate_S(run,0);
		}

	}
	
	if (NBCOMMS          ==1) {inner_el(run);	buffer_size_malloc(run);}
	if (CRT_ADAPT_INTER_0!=0) {interface_location(run);}

  if (run->con->rank==0){ printf("ForestFV: MAIN: Export 0 started \n"); }

	exportfield(run);	
	
  MPI_Barrier(MPI_COMM_WORLD);

	if (run->con->rank==0){ printf("ForestFV: MAIN: Export 0 finished \n"); }
	if (run->con->rank==0){ printf("ForestFV: MAIN: Time stepping starting.... \n"); }
	
	//exit(0);

	loop(run);

	if (run->con->rank==0){ printf("ForestFV: MAIN: Execution successful \n"); }
}