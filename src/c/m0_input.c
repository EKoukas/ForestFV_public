#include "strdata.h"

void m0_input (struct RUN* run) {

  FILE *fp;
  char *A;
  char temp[100];
  int i,iv,line,ierr;

  
  A  = malloc(1000*sizeof(char));
  fp = fopen(run->con->filecon,"r");

  line=31;
  for(i=0;i<line;i++){
    if (fgets(temp,100,fp) == NULL) {
        printf("Error skipping lines in input file \n");
        exit(0);
    }
  }

  line=29;
  
  if (fscanf(fp,"%s\t%lf\n",A,&MATERPINF[0])==2){check_db (MATERPINF[0], 0.0,10.0e15, ++line);if(run->con->rank==0){printf("31:MATERPINF.....:%lf\n",MATERPINF[0]);}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
  if (fscanf(fp,"%s\t%lf\n",A,&MATERGAMA[0])==2){check_db (MATERGAMA[0], 0.0,200.0,   ++line);if(run->con->rank==0){printf("32:MATERGAMA.....:%lf\n",MATERGAMA[0]);}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
  if (fscanf(fp,"%s\t%lf\n",A,&MATERVISC[0])==2){check_db (MATERVISC[0], 0.0,10.0e15, ++line);if(run->con->rank==0){printf("33:MATERVISC.....:%lf\n",MATERVISC[0]);}}else{printf("Problem with input file line: %d \n",++line);exit(0);}

  rewind(fp);
  fclose(fp);


  NUEQ      = 3;
  NEQMASS   = 1;
  NEQMOMENT = 3;
  NEQENERGY = 1;

  NEQ       = NEQMASS + NEQMOMENT + NEQENERGY;      
  NPRIMITIV = NEQMASS + NEQMOMENT + NEQENERGY;   

  if (PRIMTV==0) { NEQ_TEMP = NEQ;       }
  if (PRIMTV==1) { NEQ_TEMP = NPRIMITIV; }

  NMATERIALS = 1;

  eqtypn    = malloc(NUEQ*sizeof(int));
  eqtypn[0] = NEQMASS;      
  eqtypn[1] = NEQMOMENT;    
  eqtypn[2] = NEQENERGY;      

  eqtypi    = malloc(NUEQ*sizeof(int));
  eqtypi[0] = 0;
  eqtypi[1] = eqtypi[0] + eqtypn[0]; 
  eqtypi[2] = eqtypi[1] + eqtypn[1]; 

  if (ORDER==2) {
    g_limitedvars=malloc(NEQ_TEMP*sizeof(double));
    for (iv=0;iv<NEQ_TEMP;++iv){
      g_limitedvars[iv]=1.0;
    }       
  }

  return;

}