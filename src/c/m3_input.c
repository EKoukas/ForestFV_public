#include "strdata.h"

void m3_input (struct RUN* run) {

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

  if (fscanf(fp,"%s\t%lf\n",A,&AVF_LIM) ==2){check_db (AVF_LIM, 0.0,1.0, ++line);if(run->con->rank==0){printf("31:AVF_LIM....:%lf\n",AVF_LIM );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
  if (fscanf(fp,"%s\t%lf\n",A,&AMIN)    ==2){check_db (AMIN,    0.0,1.0, ++line);if(run->con->rank==0){printf("32:AMIN.......:%lf\n",AMIN    );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
  
  fscanf(fp,"%s",A); 
  fscanf(fp,"\t%d ",&CRT_ADAPT_INTER_0);     check_int(CRT_ADAPT_INTER_0    ,0,1 ,++line); 
  fscanf(fp,"\t%d ",&CRT_LVL_ADAPT_INTER_0); check_int(CRT_LVL_ADAPT_INTER_0,0,10,  line);
  fscanf(fp,"\t%d ",&CRT_LRS_ADAPT_INTER_0); check_int(CRT_LRS_ADAPT_INTER_0,0,100, line);
  if(run->con->rank==0){printf("38:ADAPT_INTER_0.: %d | %d | %d\n" ,CRT_ADAPT_INTER_0,CRT_LVL_ADAPT_INTER_0,CRT_LRS_ADAPT_INTER_0);}

  fscanf(fp,"%s",A); 
  fscanf(fp,"\t%d ",&CRT_ADAPT_INTER_1);     check_int(CRT_ADAPT_INTER_1    ,0,1 ,++line); 
  fscanf(fp,"\t%d ",&CRT_LVL_ADAPT_INTER_1); check_int(CRT_LVL_ADAPT_INTER_1,0,10,  line);
  fscanf(fp,"\t%d ",&CRT_LRS_ADAPT_INTER_1); check_int(CRT_LRS_ADAPT_INTER_1,0,100, line);
  if(run->con->rank==0){printf("39:ADAPT_INTER_1.: %d | %d | %d\n" ,CRT_ADAPT_INTER_1,CRT_LVL_ADAPT_INTER_1,CRT_LRS_ADAPT_INTER_1);}

  if ((CRT_ADAPT_INTER_0==0) && (CRT_ADAPT_INTER_1==1)) {printf("Criteria for interface adaptation wrong, lines 35/36, exiting");exit(0);}
  LEVEL = max(max(LEVEL,CRT_LVL_ADAPT_INTER_0),CRT_LVL_ADAPT_INTER_1);

  ADAPT_INTER_MAX=0;
  if ((CRT_ADAPT_INTER_0==1) && (CRT_ADAPT_INTER_1==1)) {  ADAPT_INTER_MAX = CRT_LRS_ADAPT_INTER_0 + CRT_LRS_ADAPT_INTER_1; }
  if ((CRT_ADAPT_INTER_0==1) && (CRT_ADAPT_INTER_1==0)) {  ADAPT_INTER_MAX = CRT_LRS_ADAPT_INTER_0;                         }
  
  if(run->con->rank==0){ printf("==============================================================================\n");}
  if(run->con->rank==0){ printf("gamma:   ");}
  fscanf(fp,"%s",A);
  for (i=0;i<2;i++){
    fscanf(fp,"\t%lf ",&MATERGAMA[i]); check_db (MATERGAMA[i], 0.0,100.0, ++line);  if(run->con->rank==0){  printf(" %.2e | ",MATERGAMA[i]   ); }
  }
  fscanf(fp,"\n"); if(run->con->rank==0){  printf("\n");}

 
  if(run->con->rank==0){ printf("p_inf:   ");}
  fscanf(fp,"%s",A);
  for (i=0;i<2;i++){
    fscanf(fp,"\t%lf ",&MATERPINF[i]);  check_db (MATERPINF[i], 0.0,10.0e15, ++line); if(run->con->rank==0){  printf(" %.2e | ",MATERPINF[i]   ); }
  }
  fscanf(fp,"\n"); if(run->con->rank==0){  printf("\n");}

  
  if(run->con->rank==0){ printf("viscocity:  ");}
  fscanf(fp,"%s",A);
  for (i=0;i<2;i++){
      fscanf(fp,"\t%lf ",&MATERVISC[i]); check_db (MATERVISC[i], 0.0,10.0e15, ++line); if(run->con->rank==0){  printf(" %.2e | ",MATERVISC[i]   ); }
  }
  fscanf(fp,"\n"); if(run->con->rank==0){  printf("\n");}
  if(run->con->rank==0){ printf("==============================================================================\n");}

  rewind(fp);
  fclose(fp);

  NMATERIALS = 2;

  NUEQ      = 4;
  NEQMASS   = NMATERIALS;
  NEQMOMENT = 3;
  NEQENERGY = 1;
  NEQVF     = NMATERIALS-1;

  NEQ       = NEQMASS + NEQMOMENT + NEQENERGY + NEQVF;    
  NPRIMITIV = NEQMASS + NEQMOMENT + NEQENERGY + NEQVF;  
  NCONSEQ   = NEQMASS + NEQMOMENT + NEQENERGY;    

  if (PRIMTV==0) { NEQ_TEMP = NEQ;       }
  if (PRIMTV==1) { NEQ_TEMP = NPRIMITIV; }

  eqtypn    = malloc(NUEQ*sizeof(int));
  eqtypn[0] = NEQMASS;      
  eqtypn[1] = NEQMOMENT;    
  eqtypn[2] = NEQENERGY;
  eqtypn[3] = NEQVF;      

  eqtypi    = malloc(NUEQ*sizeof(int));
  eqtypi[0] = 0;
  eqtypi[1] = eqtypi[0] + eqtypn[0]; 
  eqtypi[2] = eqtypi[1] + eqtypn[1]; 
  eqtypi[3] = eqtypi[2] + eqtypn[2]; 
  
    
  if (ORDER==2) {
    g_limitedvars=malloc(NEQ_TEMP*sizeof(double));
    for (iv=0;iv<NEQ_TEMP;++iv){
      g_limitedvars[iv]=1.0;
    }       
  }

 
  return;

}