#include "strdata.h"

void inputset(struct RUN *run) {

  FILE *fp;
  char *A;
  int i,line,irank;
 

  A  = malloc(500*sizeof(char));
  fp = fopen(run->con->filecon,"r");

  if (fp==NULL) {
    printf("ForestFV: inputset: input file not found exiting \n ");
    exit(0);
  }

  if(run->con->rank==0){ printf("=============================================================================\n");}
  if(run->con->rank==0){ printf("    ___________                               __  ______________   ___       \n");}
  if(run->con->rank==0){ printf("    \\_   _____/  ____  _______  ____   ______/  |_\\_   _____/   \\ /  /    \n");}
  if(run->con->rank==0){ printf("      |  ___)   / __ \\ \\_  __ \\/ __ \\ /  ___/   __\\ |  ___)  \\   \\  / \n");}
  if(run->con->rank==0){ printf("      |  \\__   (  \\_\\ ) |  | \\/  ___/_\\___ \\ |  |   |  \\__    \\    / \n");}
  if(run->con->rank==0){ printf("     /___  /    \\____/  |__|   \\___  /____  \\|__|  /___  /     \\  /      \n");}
  if(run->con->rank==0){ printf("         \\/                        \\/     \\/           \\/       \\/      \n");}
  if(run->con->rank==0){ printf("=============================================================================\n");}

  for (irank=0;irank<run->con->size;irank++){
    if (run->con->rank==irank){

      line = 0;
      if (fscanf(fp,"%s\t%d \n",A,&MODEL)       ==2){check_int(MODEL,       0,5,     ++line);if(run->con->rank==0){printf("1:MODEL..........: %d\n" ,MODEL       );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
      if (fscanf(fp,"%s\t%d \n",A,&TYPE_INIT)   ==2){check_int(TYPE_INIT,   0,1,     ++line);if(run->con->rank==0){printf("2:TYPE_INIT......: %d\n" ,TYPE_INIT   );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
      if (fscanf(fp,"%s\t%d \n",A,&VERBOSE)     ==2){check_int(VERBOSE,     0,1,     ++line);if(run->con->rank==0){printf("3:VERBOSE........: %d\n" ,VERBOSE     );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
      if (fscanf(fp,"%s\t%d \n",A,&CASE)        ==2){check_int(CASE,        0,5000,  ++line);if(run->con->rank==0){printf("4:CASE...........: %d\n" ,CASE        );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}

      if ((CASE==2370)&&(MODEL==5)){
      if (fscanf(fp,"%s\t%le \n",A,&S_CENTER)   ==2){check_db(S_CENTER, -1.0,100.0,  line);if(run->con->rank==0){printf("X:S_CENTER.......: %le\n" ,S_CENTER    );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}   
      }

      if ((MODEL==4)||(MODEL==5)) {
        if (((CASE>=2100)&&(CASE<=2199)) || 
            ((CASE>=3100)&&(CASE<=3199))) {
        if (fscanf(fp,"%s\t%le \n",A,&NDFI_R0_BB)      ==2){check_db(NDFI_R0_BB,      0.0,1000.0,  line);if(run->con->rank==0){printf("X:NDFI_R0_BB.....: %le\n" ,NDFI_R0_BB     );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}   
        if (fscanf(fp,"%s\t%le \n",A,&NDFI_R0_MNC)     ==2){check_db(NDFI_R0_MNC,     0.0,1000.0,  line);if(run->con->rank==0){printf("X:NDFI_R0_MNC....: %le\n" ,NDFI_R0_MNC    );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
        if (fscanf(fp,"%s\t%le \n",A,&NDFI_DIST_MNC)   ==2){check_db(NDFI_DIST_MNC,   0.0,1000.0,  line);if(run->con->rank==0){printf("X:NDFI_DIST_MNC..: %le\n" ,NDFI_DIST_MNC  );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
        if (MODEL==5) {
          if (fscanf(fp,"%s\t%le \n",A,&NDFI_DIST_SOLID_0) ==2){check_db(NDFI_DIST_SOLID_0, 0.0,10.0,  line);if(run->con->rank==0){printf("X:NDFI_DIST_SOLID_0: %le\n" ,NDFI_DIST_SOLID_0);}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
        
          if (((CASE>=2150)&&(CASE<=2159)) || ((CASE>=3150)&&(CASE<=3159))) { // 2 solids
            if (fscanf(fp,"%s\t%le \n",A,&NDFI_DIST_SOLID_1) ==2){check_db(NDFI_DIST_SOLID_0, 0.0,10.0,  line);if(run->con->rank==0){printf("X:NDFI_DIST_SOLID_1: %le\n" ,NDFI_DIST_SOLID_1);}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
          }
          if (((CASE>=2160)&&(CASE<=2169)) || ((CASE>=3160)&&(CASE<=3169))) { // 3 solids
            if (fscanf(fp,"%s\t%le \n",A,&NDFI_DIST_SOLID_2) ==2){check_db(NDFI_DIST_SOLID_0, 0.0,10.0,  line);if(run->con->rank==0){printf("X:NDFI_DIST_SOLID_2: %le\n" ,NDFI_DIST_SOLID_2);}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
          }
          if (((CASE>=2170)&&(CASE<=2179)) || ((CASE>=3170)&&(CASE<=3179))) { // 4 solids
            if (fscanf(fp,"%s\t%le \n",A,&NDFI_DIST_SOLID_3) ==2){check_db(NDFI_DIST_SOLID_0, 0.0,10.0,  line);if(run->con->rank==0){printf("X:NDFI_DIST_SOLID_3: %le\n" ,NDFI_DIST_SOLID_3);}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
          }
        }
        }
      }
    
      if (fscanf(fp,"%s\t%d \n",A,&NS)          ==2){check_int(NS,          0,1,     ++line);if(run->con->rank==0){printf("5:NS.............: %d\n" ,NS          );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
      if (fscanf(fp,"%s\t%d \n",A,&NBCOMMS)     ==2){check_int(NBCOMMS,     0,1,     ++line);if(run->con->rank==0){printf("6:NBCOMMS........: %d\n" ,NBCOMMS     );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
      if (fscanf(fp,"%s\t%d \n",A,&ORDER)       ==2){check_int(ORDER,       0,2,     ++line);if(run->con->rank==0){printf("7:ORDER..........: %d\n" ,ORDER       );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
      if (fscanf(fp,"%s\t%d \n",A,&GRAD_SCHEME) ==2){check_int(GRAD_SCHEME, 0,4,     ++line);if(run->con->rank==0){printf("8:GRAD_SCHEME....: %d\n" ,GRAD_SCHEME );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
      if (fscanf(fp,"%s\t%d \n",A,&PRIMTV)      ==2){check_int(PRIMTV,      0,1,     ++line);if(run->con->rank==0){printf("9:PRIMTV.........: %d\n" ,PRIMTV      );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
      if (fscanf(fp,"%s\t%d \n",A,&LIMITER)     ==2){check_int(LIMITER,     0,1,     ++line);if(run->con->rank==0){printf("10:LIMITER.......: %d\n" ,LIMITER     );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
      if (fscanf(fp,"%s\t%d \n",A,&RK_STEPS)    ==2){check_int(RK_STEPS,    0,4,     ++line);if(run->con->rank==0){printf("11:RK_STEPS......: %d\n" ,RK_STEPS    );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
      if (fscanf(fp,"%s\t%le\n",A,&DT)          ==2){check_db (DT,          0.0,1.0, ++line);if(run->con->rank==0){printf("12:DT............: %le\n",DT          );}}else{printf("Problem with input file line: %d \n",++line);exit(0);} 
      if (fscanf(fp,"%s\t%d \n",A,&NSTEP)       ==2){check_int(NSTEP,       0,1e6,   ++line);if(run->con->rank==0){printf("13:NSTEP.........: %d\n" ,NSTEP       );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
      if (fscanf(fp,"%s\t%d \n",A,&EXPORT_TYPE) ==2){check_int(EXPORT_TYPE, 0,3,     ++line);if(run->con->rank==0){printf("14:EXPORT_TYPE...: %d\n" ,EXPORT_TYPE );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
      if (fscanf(fp,"%s\t%d \n",A,&EXPORT_STEP) ==2){check_int(EXPORT_STEP, 0,1e6,   ++line);if(run->con->rank==0){printf("15:EXPORT_STEP...: %d\n" ,EXPORT_STEP );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
      if (fscanf(fp,"%s\t%d \n",A,&RESTART_STEP)==2){check_int(RESTART_STEP,0,1e6,   ++line);if(run->con->rank==0){printf("16:RESTART_STEP..: %d\n" ,RESTART_STEP);}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
      if (fscanf(fp,"%s\t%d \n",A,&ADAPTH)      ==2){check_int(ADAPTH,      0,1,     ++line);if(run->con->rank==0){printf("17:ADAPTH........: %d\n" ,ADAPTH      );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
      if (fscanf(fp,"%s\t%d \n",A,&ADAPTHSTEP)  ==2){check_int(ADAPTHSTEP,  0,5000,  ++line);if(run->con->rank==0){printf("18:ADAPTSTEP.....: %d\n" ,ADAPTHSTEP  );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
      if (fscanf(fp,"%s\t%d \n",A,&GEO_ADAPTH)  ==2){check_int(GEO_ADAPTH,  0,1,     ++line);if(run->con->rank==0){printf("19:GEO_ADAPTH....: %d\n" ,GEO_ADAPTH  );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}

      fscanf(fp,"%s",A); 
      fscanf(fp,"\t%d ",&CRT_UNIFORM);     check_int(CRT_UNIFORM    ,0,1 ,++line); 
      fscanf(fp,"\t%d ",&CRT_LVL_UNIFORM); check_int(CRT_LVL_UNIFORM,0,10,  line);
      if(run->con->rank==0){printf("20:CRT_UNIFORM...: %d | %d\n" ,CRT_UNIFORM,CRT_LVL_UNIFORM);}

      fscanf(fp,"%s",A); 
      fscanf(fp,"\t%d ",&CRT_GEO_AREA);     check_int(CRT_GEO_AREA    ,0,1 ,++line); 
      fscanf(fp,"\t%d ",&CRT_LVL_GEO_AREA); check_int(CRT_LVL_GEO_AREA,0,10,  line);
      if(run->con->rank==0){printf("21:CRT_GEO_AREA..: %d | %d\n" ,CRT_GEO_AREA,CRT_LVL_GEO_AREA);}

      fscanf(fp,"%s",A); 
      fscanf(fp,"\t%d " ,&CRT_GRAD_PRS);     check_int(CRT_GRAD_PRS    ,0  ,1  ,++line); 
      fscanf(fp,"\t%d " ,&CRT_LVL_GRAD_PRS); check_int(CRT_LVL_GRAD_PRS,0  ,10 ,  line);
      fscanf(fp,"\t%lf ",&CRT_GRAD_PRS_THR); check_db(CRT_GRAD_PRS_THR ,0.0,1.0,  line);
      if(run->con->rank==0){printf("22:CRT_GRAD_PRS..: %d | %d | %f \n" ,CRT_GRAD_PRS,CRT_LVL_GRAD_PRS,CRT_GRAD_PRS_THR);}

      fscanf(fp,"%s",A); 
      fscanf(fp,"\t%d " ,&CRT_GRAD_RHO);     check_int(CRT_GRAD_RHO    ,0  ,1  ,++line); 
      fscanf(fp,"\t%d " ,&CRT_LVL_GRAD_RHO); check_int(CRT_LVL_GRAD_RHO,0  ,10 ,  line);
      fscanf(fp,"\t%lf ",&CRT_GRAD_RHO_THR); check_db(CRT_GRAD_RHO_THR ,0.0,1.0,  line);
      if(run->con->rank==0){printf("23:CRT_GRAD_RHO..: %d | %d | %f \n" ,CRT_GRAD_RHO,CRT_LVL_GRAD_RHO,CRT_GRAD_RHO_THR);}

      fscanf(fp,"%s",A); 
      fscanf(fp,"\t%d " ,&CRT_GRAD_VEL);     check_int(CRT_GRAD_VEL    ,0  ,1  ,++line); 
      fscanf(fp,"\t%d " ,&CRT_LVL_GRAD_VEL); check_int(CRT_LVL_GRAD_VEL,0  ,10 ,  line);
      fscanf(fp,"\t%lf ",&CRT_GRAD_VEL_THR); check_db(CRT_GRAD_VEL_THR ,0.0,1.0,  line);
      if(run->con->rank==0){printf("24:CRT_GRAD_VEL..: %d | %d | %f \n" ,CRT_GRAD_VEL,CRT_LVL_GRAD_VEL,CRT_GRAD_VEL_THR);}

      LEVEL = max(max(max(max(CRT_LVL_UNIFORM,CRT_LVL_GEO_AREA),CRT_LVL_GRAD_PRS),CRT_LVL_GRAD_RHO),CRT_LVL_GRAD_VEL);
      
      if (fscanf(fp,"%s\t%d \n",A,&RESTART)     ==2){check_int(RESTART,     0,1,     ++line);if(run->con->rank==0){printf("25:RESTART.......: %d\n" ,RESTART     );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
      if (fscanf(fp,"%s\t%d \n",A,&TRANSLATOR)  ==2){check_int(TRANSLATOR,  0,2,     ++line);if(run->con->rank==0){printf("26:TRANSLATOR....: %d\n" ,TRANSLATOR  );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
      if (fscanf(fp,"%s\t%lf\n",A,&MSCALE)      ==2){check_db (MSCALE,      0.0,50.0,++line);if(run->con->rank==0){printf("27:MSCALE........: %lf\n",MSCALE      );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
      if (fscanf(fp,"%s\t%d \n",A,&SPLITMODE)   ==2){check_int(SPLITMODE,   0,8,     ++line);if(run->con->rank==0){printf("28:SPLITMODE.....: %d\n" ,SPLITMODE   );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
      if (fscanf(fp,"%s\t%lf\n",A,&NXQUAD)      ==2){check_db (NXQUAD,      0.0,1.0, ++line);if(run->con->rank==0){printf("29:NXQUAD........: %lf\n",NXQUAD      );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
      if (fscanf(fp,"%s\t%lf\n",A,&NYQUAD)      ==2){check_db (NYQUAD,      0.0,1.0, ++line);if(run->con->rank==0){printf("30:NYQUAD........: %lf\n",NYQUAD      );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
      if (fscanf(fp,"%s\t%lf\n",A,&NZQUAD)      ==2){check_db (NZQUAD,      0.0,1.0, ++line);if(run->con->rank==0){printf("31:NZQUAD........: %lf\n",NZQUAD      );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
      if (fscanf(fp,"%s\t%d \n",A,&ZEROU0)      ==2){check_int(ZEROU0,      0,1,     ++line);if(run->con->rank==0){printf("32:ZEROU0........: %d\n" ,ZEROU0      );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
      if (fscanf(fp,"%s\t%d \n",A,&ZEROU1)      ==2){check_int(ZEROU1,      0,1,     ++line);if(run->con->rank==0){printf("33:ZEROU1........: %d\n" ,ZEROU1      );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
      if (fscanf(fp,"%s\t%d \n",A,&ZEROU2)      ==2){check_int(ZEROU2,      0,1,     ++line);if(run->con->rank==0){printf("34:ZEROU2........: %d\n" ,ZEROU2      );}}else{printf("Problem with input file line: %d \n",++line);exit(0);}
      
      if(run->con->rank==0){  printf("bc | -x | +x | -y | +y | -z | +z | \n");printf("   ");}
      fscanf(fp,"%s",A);
      for (i=0;i<6;i++){
        fscanf(fp,"\t%d ",&BCOVERWRT[i]);  if(run->con->rank==0){  printf("|  %d ",BCOVERWRT[i]   );}
        check_int(BCOVERWRT[i],0,5,++line);
      }
      fscanf(fp,"\n");                                       if(run->con->rank==0){  printf("|\n");}

      rewind(fp);
      fclose(fp);

      switch(MODEL){
        case 0: m0_input(run); break;  // Euler EquatioNS or NS      (  gas/liquid          )
        case 1: m1_input(run); break;  // Barotropic model           (      liquid-vapor    )
        case 2: m2_input(run); break;  // Two-fluid Barotropic model (  gas-liquid-vapor    )
        case 3: m3_input(run); break;  // Kapilla's 5-equation model (  gas-liquid          )
        case 4: m4_input(run); break;  // Saurel's model 6-equation  (N:gas-N:liquid        )
        case 5: m5_input(run); break;  // Saurel's model 6-equation  (N:gas-N:liquid-N:solid)
      }

      //if (TYPE_INIT==1) { TYPE_INIT_1(run); } // Geometric initilazation NOT in init_domain
      
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  g_istep_tot=0;
  g_istep_start=0;
  free(A);

  return;

}