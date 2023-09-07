#include "strdata.h"

void m5_input (struct RUN* run) {
  
  FILE *fp;
  char *A;
  char temp[100];
  int i,iv,line,ierr,i_Amat_primitv;

  A  = malloc(1000*sizeof(char));
  fp = fopen(run->con->filecon,"r");

  line=35;

  if (CASE==2370){++line;}
  if (((CASE>=2100)&&(CASE<=2199)) || ((CASE>=3100)&&(CASE<=3199))) {
      ++line;
      ++line;
      ++line;  
      ++line;
  }
  if (((CASE>=2150)&&(CASE<=2159)) || ((CASE>=3150)&&(CASE<=3159))) { ++line; } // 2 solids
  if (((CASE>=2160)&&(CASE<=2169)) || ((CASE>=3160)&&(CASE<=3169))) { ++line; } // 3 solids
  if (((CASE>=2170)&&(CASE<=2179)) || ((CASE>=3170)&&(CASE<=3179))) { ++line; } // 4 solids

  for(i=0;i<line;i++){
    if (fgets(temp,100,fp) == NULL) {
        printf("Error skipping lines in input file \n");
        exit(0);
    }
  }
  
  if (fscanf(fp,"%s\t%lf\n",A,&AVF_LIM)    ==2){check_db (AVF_LIM,    0.0,1.0, ++line);if(run->con->rank==0){printf("35:AVF_LIM.......: %le\n",AVF_LIM   );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
  if (fscanf(fp,"%s\t%d \n",A,&PRELAX)     ==2){check_int(PRELAX,     0,3,     ++line);if(run->con->rank==0){printf("36:PRELAX........: %d\n" ,PRELAX    );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
  if (fscanf(fp,"%s\t%d \n",A,&NMATERIALS) ==2){check_int(NMATERIALS, 0,10,    ++line);if(run->con->rank==0){printf("37:NMATERIALS....: %d\n" ,NMATERIALS);}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
  if (fscanf(fp,"%s\t%lf\n",A,&AMIN)       ==2){check_db (AMIN,       0.0,1.0, ++line);if(run->con->rank==0){printf("38:AMIN..........: %le\n",AMIN      );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
  
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

  CRT_LRS_ADAPT_INTER_1 = CRT_LRS_ADAPT_INTER_0 + CRT_LRS_ADAPT_INTER_1;

  if ((CRT_ADAPT_INTER_0==0) && (CRT_ADAPT_INTER_1==1)) {printf("Criteria for interface adaptation wrong, lines 35/36, exiting");exit(0);}
  LEVEL = max(max(LEVEL,CRT_LVL_ADAPT_INTER_0),CRT_LVL_ADAPT_INTER_1);

  ADAPT_INTER_MAX=0;
  if ((CRT_ADAPT_INTER_0==1) && (CRT_ADAPT_INTER_1==1)) {  ADAPT_INTER_MAX = CRT_LRS_ADAPT_INTER_0 + CRT_LRS_ADAPT_INTER_1; }
  if ((CRT_ADAPT_INTER_0==1) && (CRT_ADAPT_INTER_1==0)) {  ADAPT_INTER_MAX = CRT_LRS_ADAPT_INTER_0;                         }

  if(run->con->rank==0){ printf("==============================================================================\n");}
  if(run->con->rank==0){ printf("p_inf............:");}
  fscanf(fp,"%s",A);
  for (i=0;i<NMATERIALS;i++){
    fscanf(fp,"\t%lf ",&MATERPINF[i]); check_db (MATERPINF[i], 0.0,10.0e15, ++line); if(run->con->rank==0){  printf(" %.2e |",MATERPINF[i]   ); }
  }
  fscanf(fp,"\n"); if(run->con->rank==0){  printf("\n");}

  if(run->con->rank==0){ printf("gamma............:");}
  fscanf(fp,"%s",A);
  for (i=0;i<NMATERIALS;i++){
    fscanf(fp,"\t%lf ",&MATERGAMA[i]); check_db (MATERGAMA[i], 0.0,100.0, ++line);  if(run->con->rank==0){  printf(" %.2e |",MATERGAMA[i]   ); }
  }
  fscanf(fp,"\n"); if(run->con->rank==0){  printf("\n");}

  if(run->con->rank==0){ printf("mu...............:");}
  fscanf(fp,"%s",A);
  for (i=0;i<NMATERIALS;i++){
    fscanf(fp,"\t%lf ",&MATERMUSH[i]); check_db (MATERMUSH[i], 0.0,10.0e15, ++line);         if(run->con->rank==0){  printf(" %.2e |",MATERMUSH[i]   ); }
    if (MATERMUSH[i]!=0.0) {
      n_solids++;
    }
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
    fscanf(fp,"\t%lf ",&u_init[i]); check_db (u_init[i], 0.0,10000.0, ++line); if(run->con->rank==0){  printf(" %.2e |",u_init[i]   ); }
  }
  fscanf(fp,"\n"); if(run->con->rank==0){  printf("\n");}


  if(run->con->rank==0){ printf("v_init...........:");}
  fscanf(fp,"%s",A);
  for (i=0;i<NMATERIALS;i++){
    fscanf(fp,"\t%lf ",&v_init[i]); check_db (v_init[i], 0.0,10000.0, ++line); if(run->con->rank==0){  printf(" %.2e |",v_init[i]   ); }
  }
  fscanf(fp,"\n"); if(run->con->rank==0){  printf("\n");}


  if(run->con->rank==0){ printf("w_init...........:");}
  fscanf(fp,"%s",A);
  for (i=0;i<NMATERIALS;i++){
    fscanf(fp,"\t%lf ",&w_init[i]); check_db (w_init[i], 0.0,10000.0, ++line); if(run->con->rank==0){  printf(" %.2e |",w_init[i]   ); }
  }
  fscanf(fp,"\n"); if(run->con->rank==0){  printf("\n");}
  if(run->con->rank==0){ printf("==============================================================================\n");}

  rewind(fp);
  fclose(fp);

  NUEQ       = 6;
  NEQMASS    = NMATERIALS;
  NEQMOMENT  = 3;
  NEQENERGY  = 1;
  NEQVF      = NMATERIALS-1;
  NEQSENERGY = NMATERIALS;
  NEQAMAT    = 9;

  NEQ       = NEQMASS + NEQMOMENT + NEQENERGY + NEQVF + NEQSENERGY + NEQAMAT;      
  NCONSEQ   = NEQMASS + NEQMOMENT + NEQENERGY;  
  NPRIMITIV = NEQMASS + NEQMOMENT + NEQSENERGY + NEQVF + NEQAMAT;     

  if(run->con->rank==0){
    printf("\n");
    if(run->con->rank==0){ printf("==============================================================================\n");}
    printf("NEQ..............: %d \n",NEQ);
    printf("NEQMASS..........: %d \n",NEQMASS);
    printf("NEQMOMENT........: %d \n",NEQMOMENT);
    printf("NEQENERGY........: %d \n",NEQENERGY);
    printf("NEQVF............: %d \n",NEQVF);
    printf("NEQSENERGY.......: %d \n",NEQSENERGY);
    printf("NEQAMAT..........: %d \n",NEQAMAT);
    if(run->con->rank==0){ printf("==============================================================================\n");}
  }

  if      (PRIMTV==0) { NEQ_TEMP = NEQ;       }
  else if (PRIMTV==1) { NEQ_TEMP = NPRIMITIV; }

  eqtypn    = malloc(NUEQ*sizeof(int));
  eqtypn[0] = NEQMASS;      
  eqtypn[1] = NEQMOMENT;    
  eqtypn[2] = NEQENERGY;
  eqtypn[3] = NEQVF; 
  eqtypn[4] = NEQSENERGY;     
  eqtypn[5] = NEQAMAT;

  eqtypi    = malloc(NUEQ*sizeof(int));
  eqtypi[0] = 0;
  eqtypi[1] = eqtypi[0] + eqtypn[0]; 
  eqtypi[2] = eqtypi[1] + eqtypn[1]; 
  eqtypi[3] = eqtypi[2] + eqtypn[2];
  eqtypi[4] = eqtypi[3] + eqtypn[3]; 
  eqtypi[5] = eqtypi[4] + eqtypn[4];
  

  if (ORDER==2){

    g_limitedvars=malloc(NEQ_TEMP*sizeof(double));

    if (PRIMTV==1){

      for (iv=0;iv<NEQ_TEMP;++iv){
        g_limitedvars[iv]=1.0;
      } 
      i_Amat_primitv = eqtypn[3] + 
                       eqtypn[0] + 
                       eqtypn[1] + 
                       eqtypn[4];
      for (iv=i_Amat_primitv;iv<i_Amat_primitv+9;++iv){
        g_limitedvars[iv]=0.0;
      } 
      
    } else {

      for (iv=0;iv<NEQ_TEMP;++iv){
        g_limitedvars[iv]=1.0;
      }
      for (iv=eqtypi[4];iv<eqtypi[5];++iv){    
        g_limitedvars[iv]=0.0;
      } 

    }

  }


  return;

}