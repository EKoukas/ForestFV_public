#include "strdata.h"

void face_interpolation_V0(struct RUN * run){
  
  run->T2nd=timecpu(run->T2nd,0);

  switch(MODEL) {
    case 0: m0_faceinterp_V0(run); break;
    case 1: m1_faceinterp_V0(run); break;
    case 2: m2_faceinterp_V0(run); break;
    case 3: m3_faceinterp_V0(run); break;
    case 4: m4_faceinterp_V0(run); break;
    case 5: m5_faceinterp_V0(run); break;
  }

  run->T2nd=timecpu(run->T2nd,1);

    run->T2ndcom=timecpu(run->T2ndcom,0);   
  communicate_S(run,1);     
    run->T2ndcom=timecpu(run->T2ndcom,1);
    
}

void face_interpolation_V1(struct RUN * run){
  
  run->T2nd=timecpu(run->T2nd,0);

  switch(MODEL) {
    case 0: m0_faceinterp_V1(run); break;
    case 1: m1_faceinterp_V1(run); break;
    case 2: m2_faceinterp_V1(run); break;
    case 3: m3_faceinterp_V1(run); break;
    case 4: m4_faceinterp_V1(run); break;
    case 5: m5_faceinterp_V1(run); break;
  }

  run->T2nd=timecpu(run->T2nd,1);

    run->T2ndcom=timecpu(run->T2ndcom,0);   
  communicate_S(run,1);     
    run->T2ndcom=timecpu(run->T2ndcom,1);
    
}

void face_interpolation_V2(struct RUN * run){
  
  run->T2nd=timecpu(run->T2nd,0);

  switch(MODEL) {
    case 0: m0_faceinterp_V2(run); break;
    //case 1: m1_faceinterp_V2(run); break;
    //case 2: m2_faceinterp_V2(run); break;
    //case 3: m3_faceinterp_V2(run); break;
    //case 4: m4_faceinterp_V2(run); break;
    //case 5: m5_faceinterp_V2(run); break;
  }

  run->T2nd=timecpu(run->T2nd,1);

    run->T2ndcom=timecpu(run->T2ndcom,0);   
  communicate_S(run,1);     
    run->T2ndcom=timecpu(run->T2ndcom,1);
    
}

void face_interpolation_unstr(struct RUN * run){
  
  run->T2nd=timecpu(run->T2nd,0);

  switch(MODEL){
    case 0: m0_faceinterp_unstr(run,-1); break;
    case 1: m1_faceinterp_unstr(run,-1); break;
    case 2: m2_faceinterp_unstr(run,-1); break;
    case 3: m3_faceinterp_unstr(run,-1); break;
    case 4: m4_faceinterp_unstr(run,-1); break;
    case 5: m5_faceinterp_unstr(run,-1); break;
  }

  run->T2nd=timecpu(run->T2nd,1);

    run->T2ndcom=timecpu(run->T2ndcom,0);   
  communicate_S(run,1);     
    run->T2ndcom=timecpu(run->T2ndcom,1);
    
}

void face_interpolation_V0_nb(struct RUN * run){

  run->T2nd=timecpu(run->T2nd,0);

  switch(MODEL) {
    // Face interp for inner cells
    case 0: m0_faceinterp_V0_nb(run,1); break;
    case 1: m1_faceinterp_V0_nb(run,1); break;
    case 2: m2_faceinterp_V0_nb(run,1); break;
    case 3: m3_faceinterp_V0_nb(run,1); break;
    case 4: m4_faceinterp_V0_nb(run,1); break;
    case 5: m5_faceinterp_V0_nb(run,1); break;
  }

  if (((g_istep)%EXPORT_STEP!=0) || (g_istep==0)) {
    communicate_nb_S(run,0,1);  // Non blocking recv for S  (A) (B)
  }
 
  switch(MODEL) {
    // Face interp for outer cells
    case 0: m0_faceinterp_V0_nb(run,0); break;
    case 1: m1_faceinterp_V0_nb(run,0); break;
    case 2: m2_faceinterp_V0_nb(run,0); break;
    case 3: m3_faceinterp_V0_nb(run,0); break;
    case 4: m4_faceinterp_V0_nb(run,0); break;
    case 5: m5_faceinterp_V0_nb(run,0); break;
  }

  run->T2nd=timecpu(run->T2nd,0);
 
    run->T2ndcom=timecpu(run->T2ndcom,0);   
  communicate_nb_S(run,1,0);  // Non blocking send for SF (I)
    run->T2ndcom=timecpu(run->T2ndcom,1);
  
}

void face_interpolation_V1_nb(struct RUN * run){

  run->T2nd=timecpu(run->T2nd,0);

  switch(MODEL){
    // Face interp for inner cells
    case 0: m0_faceinterp_V1_nb(run,1); break;
    case 1: m1_faceinterp_V1_nb(run,1); break;
    case 2: m2_faceinterp_V1_nb(run,1); break;
    case 3: m3_faceinterp_V1_nb(run,1); break;
    case 4: m4_faceinterp_V1_nb(run,1); break;
    case 5: m5_faceinterp_V1_nb(run,1); break;
  }

  if (((g_istep)%EXPORT_STEP!=0) || (g_istep==0)) {
    communicate_nb_S(run,0,1);  // Non blocking recv for S  (A) (B)
  }
 
  switch(MODEL){
    // Face interp for outer cells
    case 0: m0_faceinterp_V1_nb(run,0); break;
    case 1: m1_faceinterp_V1_nb(run,0); break;
    case 2: m2_faceinterp_V1_nb(run,0); break;
    case 3: m3_faceinterp_V1_nb(run,0); break;
    case 4: m4_faceinterp_V1_nb(run,0); break;
    case 5: m5_faceinterp_V1_nb(run,0); break;
  }

  run->T2nd=timecpu(run->T2nd,1);
 
    run->T2ndcom=timecpu(run->T2ndcom,0);   
  communicate_nb_S(run,1,0);  // Non blocking send for SF (I)
    run->T2ndcom=timecpu(run->T2ndcom,1);
  
}

void face_interpolation_V2_nb(struct RUN * run){

  /*
  run->T2nd=timecpu(run->T2nd,0);

  switch(MODEL){
    // Face interp for inner cells
    case 0: m0_faceinterp_V2_nb(run,1); break;
    case 1: m1_faceinterp_V2_nb(run,1); break;
    case 2: m2_faceinterp_V2_nb(run,1); break;
    case 3: m3_faceinterp_V2_nb(run,1); break;
    case 4: m4_faceinterp_V2_nb(run,1); break;
    case 5: m5_faceinterp_V2_nb(run,1); break;
  }

  if (((g_istep)%EXPORT_STEP!=0) || (g_istep==0)) {
    communicate_nb_S(run,0,1);  // Non blocking recv for S  (A) (B)
  }
 
  switch(MODEL){
    // Face interp for outer cells
    case 0: m0_faceinterp_V2_nb(run,0); break;
    case 1: m1_faceinterp_V2_nb(run,0); break;
    case 2: m2_faceinterp_V2_nb(run,0); break;
    case 3: m3_faceinterp_V2_nb(run,0); break;
    case 4: m4_faceinterp_V2_nb(run,0); break;
    case 5: m5_faceinterp_V2_nb(run,0); break;
  }

  run->T2nd=timecpu(run->T2nd,1);
 
    run->T2ndcom=timecpu(run->T2ndcom,0);   
  communicate_nb_S(run,1,0);  // Non blocking send for SF (I)
    run->T2ndcom=timecpu(run->T2ndcom,1);
  */

}

void face_interpolation_unstr_nb(struct RUN * run){
  
  run->T2nd=timecpu(run->T2nd,0);

  switch(MODEL){
    case 0: m0_faceinterp_unstr(run,1); break;
    case 1: m1_faceinterp_unstr(run,1); break;
    case 2: m2_faceinterp_unstr(run,1); break;
    case 3: m3_faceinterp_unstr(run,1); break;
    case 4: m4_faceinterp_unstr(run,1); break;
    case 5: m5_faceinterp_unstr(run,1); break;
  }

  if (((g_istep)%EXPORT_STEP!=0) || (g_istep==0)) {
    communicate_nb_S(run,0,1);  // Non blocking recv for S  (A) (B)
  }

  switch(MODEL){
    case 0: m0_faceinterp_unstr(run,0); break;
    case 1: m1_faceinterp_unstr(run,0); break;
    case 2: m2_faceinterp_unstr(run,0); break;
    case 3: m3_faceinterp_unstr(run,0); break;
    case 4: m4_faceinterp_unstr(run,0); break;
    case 5: m5_faceinterp_unstr(run,0); break;
  }

  run->T2nd=timecpu(run->T2nd,1);

    run->T2ndcom=timecpu(run->T2ndcom,0);   
  communicate_nb_S(run,1,0);     
    run->T2ndcom=timecpu(run->T2ndcom,1);
    
}