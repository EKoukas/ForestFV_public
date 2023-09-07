#include "strdata.h"

void m4_input (struct RUN* run) {
  
  FILE *fp;
  char *A;
  char temp[100];
  int i,iv,line,ierr;
  

  A  = malloc(1000*sizeof(char));
  fp = fopen(run->con->filecon,"r");

  line=35;

  if (((CASE>=2100)&&(CASE<=2199)) || 
      ((CASE>=3100)&&(CASE<=3199))) {
        ++line;
        ++line;
        ++line;          
  }
  


  for(i=0;i<line;i++){
    if (fgets(temp,100,fp) == NULL) {
        printf("Error skipping lines in input file \n");
        exit(0);
    }
  }
   
  if (fscanf(fp,"%s\t%lf\n",A,&AVF_LIM)    ==2){check_db (AVF_LIM,    0.0,1.0, ++line);if(run->con->rank==0){printf("31:AVF_LIM.......: %lf\n",AVF_LIM   );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
  if (fscanf(fp,"%s\t%d \n",A,&PRELAX)     ==2){check_int(PRELAX,     0,3,     ++line);if(run->con->rank==0){printf("32:PRELAX........: %d\n" ,PRELAX    );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
  if (fscanf(fp,"%s\t%d \n",A,&NMATERIALS) ==2){check_int(NMATERIALS, 0,10,    ++line);if(run->con->rank==0){printf("33:NMATERIALS....: %d\n" ,NMATERIALS);}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
  if (fscanf(fp,"%s\t%lf\n",A,&AMIN)       ==2){check_db (AMIN,       0.0,1.0, ++line);if(run->con->rank==0){printf("34:AMIN..........: %le\n",AMIN      );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
  
  fscanf(fp,"%s",A); 
  fscanf(fp,"\t%d ",&CRT_ADAPT_INTER_0);     check_int(CRT_ADAPT_INTER_0    ,0,1 ,++line); 
  fscanf(fp,"\t%d ",&CRT_LVL_ADAPT_INTER_0); check_int(CRT_LVL_ADAPT_INTER_0,0,10,  line);
  fscanf(fp,"\t%d ",&CRT_LRS_ADAPT_INTER_0); check_int(CRT_LRS_ADAPT_INTER_0,0,100, line);
  if(run->con->rank==0){printf("35:ADAPT_INTER_0...: %d | %d | %d\n" ,CRT_ADAPT_INTER_0,CRT_LVL_ADAPT_INTER_0,CRT_LRS_ADAPT_INTER_0);}

  fscanf(fp,"%s",A); 
  fscanf(fp,"\t%d ",&CRT_ADAPT_INTER_1);     check_int(CRT_ADAPT_INTER_1    ,0,1 ,++line); 
  fscanf(fp,"\t%d ",&CRT_LVL_ADAPT_INTER_1); check_int(CRT_LVL_ADAPT_INTER_1,0,10,  line);
  fscanf(fp,"\t%d ",&CRT_LRS_ADAPT_INTER_1); check_int(CRT_LRS_ADAPT_INTER_1,0,100, line);
  if(run->con->rank==0){printf("36:ADAPT_INTER_1...: %d | %d | %d\n" ,CRT_ADAPT_INTER_1,CRT_LVL_ADAPT_INTER_1,CRT_LRS_ADAPT_INTER_1);}

  if ((CRT_ADAPT_INTER_0==0) && (CRT_ADAPT_INTER_1==1)) {printf("Criteria for interface adaptation wrong, lines 35/36, exiting");exit(0);}
  LEVEL = max(max(LEVEL,CRT_LVL_ADAPT_INTER_0),CRT_LVL_ADAPT_INTER_1);

  ADAPT_INTER_MAX=0;
  if ((CRT_ADAPT_INTER_0==1) && (CRT_ADAPT_INTER_1==1)) {  ADAPT_INTER_MAX = CRT_LRS_ADAPT_INTER_0 + CRT_LRS_ADAPT_INTER_1; }
  if ((CRT_ADAPT_INTER_0==1) && (CRT_ADAPT_INTER_1==0)) {  ADAPT_INTER_MAX = CRT_LRS_ADAPT_INTER_0;                         }
  
  if(run->con->rank==0){ printf("=============================================================================\n");}
  if(run->con->rank==0){ printf("p_inf............:");}
  fscanf(fp,"%s",A);
  for (i=0;i<NMATERIALS;i++){
    fscanf(fp,"\t%lf ",&MATERPINF[i]);  check_db (MATERPINF[i], 0.0,10.0e15, ++line); if(run->con->rank==0){  printf(" %.2e |",MATERPINF[i]   ); }
  }
  fscanf(fp,"\n"); if(run->con->rank==0){  printf("\n");}


  if(run->con->rank==0){ printf("gamma............:");}
  fscanf(fp,"%s",A);
  for (i=0;i<NMATERIALS;i++){
    fscanf(fp,"\t%lf ",&MATERGAMA[i]); check_db (MATERGAMA[i], 0.0,100.0, ++line);  if(run->con->rank==0){  printf(" %.2e |",MATERGAMA[i]   ); }
  }
  fscanf(fp,"\n"); if(run->con->rank==0){  printf("\n");}

  
  if(run->con->rank==0){ printf("viscocity........:");}
  fscanf(fp,"%s",A);
  for (i=0;i<NMATERIALS;i++){
    fscanf(fp,"\t%lf ",&MATERVISC[i]); check_db (MATERVISC[i], 0.0,10.0e15, ++line); if(run->con->rank==0){  printf(" %.2e |",MATERVISC[i]   ); }
  }
  fscanf(fp,"\n"); if(run->con->rank==0){  printf("\n");}

  if(run->con->rank==0){ printf("rho init.........:");}
  fscanf(fp,"%s",A);
  for (i=0;i<NMATERIALS;i++){
    fscanf(fp,"\t%lf ",&MATERRINI[i]); check_db (MATERRINI[i], 0.0,100000.0, ++line); if(run->con->rank==0){  printf(" %.2e |",MATERRINI[i]   ); }
  }
  fscanf(fp,"\n"); if(run->con->rank==0){  printf("\n");}

  if(run->con->rank==0){ printf("pr. init.........:");}
  fscanf(fp,"%s",A);
  for (i=0;i<NMATERIALS;i++){
    fscanf(fp,"\t%lf ",&MATERPINI[i]); check_db (MATERPINI[i], 0.0,10.0e15, ++line); if(run->con->rank==0){  printf(" %.2e |",MATERPINI[i]   ); }
  }
  fscanf(fp,"\n"); if(run->con->rank==0){  printf("\n");}

  if(run->con->rank==0){ printf("u_init...........:");}
  fscanf(fp,"%s",A);
  for (i=0;i<NMATERIALS;i++){
    fscanf(fp,"\t%lf ",&u_init[i]); check_db (u_init[i], 0.0,1000.0, ++line); if(run->con->rank==0){  printf(" %.2e |",u_init[i]   ); }
  }
  fscanf(fp,"\n"); if(run->con->rank==0){  printf("\n");}


  if(run->con->rank==0){ printf("v_init...........:");}
  fscanf(fp,"%s",A);
  for (i=0;i<NMATERIALS;i++){
    fscanf(fp,"\t%lf ",&v_init[i]); check_db (v_init[i], 0.0,1000.0, ++line); if(run->con->rank==0){  printf(" %.2e |",v_init[i]   ); }
  }
  fscanf(fp,"\n"); if(run->con->rank==0){  printf("\n");}


  if(run->con->rank==0){ printf("w_init...........:");}
  fscanf(fp,"%s",A);
  for (i=0;i<NMATERIALS;i++){
    fscanf(fp,"\t%lf ",&w_init[i]); check_db (w_init[i], 0.0,1000.0, ++line); if(run->con->rank==0){  printf(" %.2e |",w_init[i]   ); }
  }
  fscanf(fp,"\n"); if(run->con->rank==0){  printf("\n");}
  if(run->con->rank==0){ printf("=============================================================================\n");}

  rewind(fp);
  fclose(fp);

  NUEQ       = 5;
  NEQMASS    = NMATERIALS;
  NEQMOMENT  = 3;
  NEQENERGY  = 1;
  NEQVF      = NMATERIALS-1;
  NEQSENERGY = NMATERIALS;

  NEQ       = NEQMASS + NEQMOMENT + NEQENERGY + NEQVF + NEQSENERGY;      
  NCONSEQ   = NEQMASS + NEQMOMENT + NEQENERGY;  
  NPRIMITIV = NEQMASS + NEQMOMENT + NEQSENERGY + NEQVF;    

  if      (PRIMTV==0) { NEQ_TEMP = NEQ;       }
  else if (PRIMTV==1) { NEQ_TEMP = NPRIMITIV; }

  if(run->con->rank==0){
    printf("\n");
    if(run->con->rank==0){ printf("==============================================================================\n");}
    printf("NEQ..............: %d \n",NEQ);
    printf("NEQMASS..........: %d \n",NEQMASS);
    printf("NEQMOMENT........: %d \n",NEQMOMENT);
    printf("NEQENERGY........: %d \n",NEQENERGY);
    printf("NEQVF............: %d \n",NEQVF);
    printf("NEQSENERGY.......: %d \n",NEQSENERGY);
    if(run->con->rank==0){ printf("==============================================================================\n");}
  }

  eqtypn    = malloc(NUEQ*sizeof(int));
  eqtypn[0] = NEQMASS;      
  eqtypn[1] = NEQMOMENT;    
  eqtypn[2] = NEQENERGY;
  eqtypn[3] = NEQVF; 
  eqtypn[4] = NEQSENERGY;     

  eqtypi    = malloc(NUEQ*sizeof(int));
  eqtypi[0] = 0;
  eqtypi[1] = eqtypi[0] + eqtypn[0]; 
  eqtypi[2] = eqtypi[1] + eqtypn[1]; 
  eqtypi[3] = eqtypi[2] + eqtypn[2];
  eqtypi[4] = eqtypi[3] + eqtypn[3]; 
  
  if (ORDER==2) {
    g_limitedvars=malloc(NEQ_TEMP*sizeof(double));
    for (iv=0;iv<NEQ_TEMP;++iv){
      g_limitedvars[iv]=1.0;
    }       
  }

  return;

}