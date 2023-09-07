#include "strdata.h"

void m4_init_domain(struct RUN *run){

  int i,iv,ieq,j;
  double xw,yw,zw;
  double AMIN_temp,R0,dist_wall,x_solid,x_start,S_center;
  double temp_p,temp_rho;
  struct BRANCH * crnt;
  
  crnt=run->topo->locl;
  while (crnt!=NULL){
    
    xw=crnt->cl->xc;
    yw=crnt->cl->yc;
    zw=crnt->cl->zc;

    for(iv=0;iv<NEQ;iv++){
      crnt->el->S[iv] = 0.0;
    }

    if (CASE==0) {            // Interface advection [Liquid Gas] 

      if((xw)>0.5) {
        m4_init_S(run,crnt,1); // air
      } else {
        m4_init_S(run,crnt,0); // water
      }        

    } else if (CASE==1) {       // Shock tube case [Liquid Gas] 

      if((xw)>0.75) {
        m4_init_S(run,crnt,1); // air
      } else {
        m4_init_S(run,crnt,0);
      }

    } else if (CASE==3) {       // Cavitation Test [Liquid Gas]

      if((xw)>0.5){ 
        u_init[0] = 100.0;
        m4_init_S(run,crnt,0); 
      } else {
        u_init[0] = -100.0;
        m4_init_S(run,crnt,0);
      }

    } else if (CASE==200) {     // Solid - Gas advection 2D

      if(((xw-0.5)*(xw-0.5))+(((yw-0.5)*(yw-0.5)))<pow((0.3/2.0),2.0)){    // Initialize circle in the middle
        m4_init_S(run,crnt,0); 
      } else {
        m4_init_S(run,crnt,1);
      }

    } else if (CASE==201) {     // Shock induced Bubble collapse

      double amin_temp;
      amin_temp = AMIN;
      if (((xw-0.05)*(xw-0.05))+((yw-0.0)*(yw-0.0))<=pow((0.05/2.0),2.0)){ // Bubble
        AMIN = 0.05;
        u_init[0]    = 0.0;
        u_init[1]    = 0.0;

        MATERRINI[0] = 0.158;
        MATERRINI[1] = 1.204;

        MATERPINI[0] = 101325.0;
        MATERPINI[1] = 101325.0;

        m4_init_S(run,crnt,0);
      } else if (xw>0.1) {
        AMIN = pow(10.0,-8);
        
        u_init[0]    =    0.0;
        u_init[1]    = -114.49;
        
        MATERRINI[0] = 0.158;
        MATERRINI[1] = 1.608;

        MATERPINI[0] = 159060.0;
        MATERPINI[1] = 159060.0;

        m4_init_S(run,crnt,1);
      } else {                  // Pre-shock
        AMIN = pow(10.0,-8);
        u_init[0]    = 0.0;
        u_init[1]    = 0.0;

        MATERRINI[0] = 0.158;
        MATERRINI[1] = 1.204;

        MATERPINI[0] = 101325.0;
        MATERPINI[1] = 101325.0;

        m4_init_S(run,crnt,1);
      }
      AMIN = AMIN_temp;

    } else if (CASE==202) {     // RMI

      if((((xw-1.8)*(xw-1.8))+((yw-0.0)*(yw-0.0))>=pow(0.6,2.0)) && (xw>1.8)) {
        m4_init_S(run,crnt,0);
      } else {
        m4_init_S(run,crnt,1);
      }
        
    } else if (CASE==2100) {    // needle free injection kyrizis paper
        
      R0 = 0.35/1000.0;
      R0 = R0/5.0;
      MATERPINI[1]=5.0*pow(10.0,7.0); // high pressure
      MATERRINI[1]=30.0;              // high pressure
      if ( (pow((xw),2.0) + pow((yw),2.0)) <= pow(R0,2.0)) {     // high pressure bubble
        m4_init_S(run,crnt,1);  // bubble
      } else {
        m4_init_S(run,crnt,0);  // Water
      }

      dist_wall = 1.7/1000.0;
      R0 = 0.26/1000.0;
      MATERPINI[1]=101325.0;
      MATERRINI[1]=1.225;     // low pressure
      if ( (pow((xw-dist_wall),2.0) + pow((yw),2.0)) <= pow(R0,2.0)) {     // atm pressure gas
        m4_init_S(run,crnt,1);       // gas
      }
      if (xw>=dist_wall ) {
        m4_init_S(run,crnt,1);       // gas
      }

    } else if ((CASE==2110)||(CASE==2111)) { // needle free injection Twente     2D/axi [ok]
        
      R0 = 60.0/pow(10.0,6.0);
      temp_p   = MATERPINI[1];
      temp_rho = MATERRINI[1];
      MATERPINI[1]=temp_p;      // high pressure
      MATERRINI[1]=temp_rho;    // high pressure
      if ( (pow((xw),2.0) + pow((yw),2.0)) <= pow(R0,2.0)) {     // high pressure bubble
        m4_init_S(run,crnt,1);  // bubble
      } else {
        m4_init_S(run,crnt,0);  // Water
      }

      
      dist_wall = 1350.0/pow(10.0,6.0);
      R0 = 200.0/pow(10.0,6.0);
      MATERPINI[1]=101325.0;
      MATERRINI[1]=1.225;     // low pressure
      if ( (pow((xw-dist_wall),2.0) + pow((yw),2.0)) <= pow(R0,2.0)) {     // atm pressure gas
        m4_init_S(run,crnt,1);       // gas
      }
      if (xw>=dist_wall ) {
        m4_init_S(run,crnt,1);       // gas
      }
      
      MATERPINI[1] = temp_p;
      MATERRINI[1] = temp_rho;
    
    } else if ((CASE==2120)||(CASE==2121)) { // needle free injection twente wba 2D/axi
        
      //R0 = 0.12/1000.0;
      R0 = 0.06/1000.0;
      MATERPINI[1]=2.7*pow(10.0,6.0); // high pressure
      MATERRINI[1]=32.0;              // high pressure
      if ( (pow((xw),2.0) + pow((yw),2.0)) <= pow(R0,2.0)) {     // high pressure bubble
        m4_init_S(run,crnt,1);  // bubble
      } else {
        m4_init_S(run,crnt,0);  // Water
      }

      dist_wall = 0.895/1000.0;
      R0 = 0.282/1000.0;
      MATERPINI[1]=101325.0;
      MATERRINI[1]=1.225;     // low pressure
      if ( (pow((xw-dist_wall),2.0) + pow((yw),2.0)) <= pow(R0,2.0)) {     // atm pressure gas
        m4_init_S(run,crnt,1);       // gas
      }

      if (xw>=dist_wall ) {
        m4_init_S(run,crnt,1);       // gas
      }
    } else if (CASE==2360) {    // KS bubble colonius wang

      shock_prof = 0;

      R0        = 0.05/1000.0;
      S_center  = 0.1/1000.0;
      dist_wall = 4.0/1000.0 - S_center;
      x_start   = dist_wall - 1.25*R0;

      if ( (pow((xw-dist_wall),2.0) + pow((yw-0.0),2.0)) <= pow(R0,2.0)) {     
        m4_init_S(run,crnt,1);              // bubble
      } else {
        m4_init_S(run,crnt,0);       // Water
      }

      if ( xw<=x_start ) {
        m4_init_S_SW(run,crnt,0,x_start,xw);       // Water
      }
    } else if (CASE==2370) {    // KS bubble attachment / paper_0 / No solid

      shock_prof = 0;

      R0 = MSCALE*10.0;   
      x_solid   = 16.0*R0;  

      dist_wall = x_solid   - S_center*R0;
      x_start   = dist_wall - 1.25*R0;
                 
      if ( (pow((xw-dist_wall),2.0) + pow((yw-0.0),2.0)) <= pow(R0,2.0)) {     
        m4_init_S(run,crnt,1);       // bubble
      } else {
        m4_init_S(run,crnt,0);       // Water
      }

      if ( xw<=x_start ) {
        m4_init_S_SW(run,crnt,0,x_start,xw);       // Water
      }
    
    

    } else if (CASE==2500) {    // bb col / experiment
        
      R0 = 0.05/1000.0;  
      //R0 = (MSCALE*10)/2.0;  
      dist_wall = S_center*R0;
      
      if ( (pow((xw-dist_wall),2.0) + pow((yw),2.0)) <= pow(R0,2.0)) {     
        m4_init_S(run,crnt,1);             // bubble
      } else {
        m4_init_S(run,crnt,0);             // Water
      }
    } else if (CASE==3000) {    // Shock bubble / blood vessel


    } else if (CASE==3010) {    // blood vessel / expanding bubble

    } else if (CASE==3100) {    // needle free injection 3D curved
       
       m4_init_S(run,crnt,0);  // Water

    } else if (CASE==3110) {    // needle free injection 3D  [ok]
        
      R0 = 40.0/pow(10.0,6.0);
      temp_p   = MATERPINI[1];
      temp_rho = MATERRINI[1];
      MATERPINI[1]=temp_p;      // high pressure
      MATERRINI[1]=temp_rho;    // high pressure
      if ( (pow((xw),2.0) + pow((yw),2.0) + pow((zw-(50.0/pow(10.0,6.0))),2.0) ) <= pow(R0,2.0)) {     // high pressure bubble
        m4_init_S(run,crnt,1);  // bubble
      } else {
        m4_init_S(run,crnt,0);  // Water
      }

      dist_wall = 1350.0/pow(10.0,6.0);
      R0 = 200.0/pow(10.0,6.0);
      MATERPINI[1]=101325.0;
      MATERRINI[1]=1.225;     // low pressure
      if ( (pow((xw-dist_wall),2.0) + pow((yw),2.0)) <= pow(R0,2.0)) {     // atm pressure gas
        m4_init_S(run,crnt,1);       // gas
      }
      if (xw>=dist_wall ) {
        m4_init_S(run,crnt,1);       // gas
      }
      
      MATERPINI[1] = temp_p;
      MATERRINI[1] = temp_rho;

    } else if (CASE==3120) {    // needle free injection 3D half [ok]
        
      R0 = 40.0/pow(10.0,6.0);
      temp_p   = MATERPINI[1];
      temp_rho = MATERRINI[1];
      MATERPINI[1]=temp_p;      // high pressure
      MATERRINI[1]=temp_rho;    // high pressure
      if ( (pow((xw),2.0) + pow((yw),2.0) + pow((zw),2.0) ) <= pow(R0,2.0)) {     // high pressure bubble
        m4_init_S(run,crnt,1);  // bubble
      } else {
        m4_init_S(run,crnt,0);  // Water
      }
      
      dist_wall = 1350.0/pow(10.0,6.0);
      R0 = 200.0/pow(10.0,6.0);
      MATERPINI[1]=101325.0;
      MATERRINI[1]=1.225;     // low pressure
      if ( (pow((xw-dist_wall),2.0) + pow((yw),2.0)) <= pow(R0,2.0)) {     // atm pressure gas
        m4_init_S(run,crnt,1);       // gas
      }
      if (xw>=dist_wall ) {
        m4_init_S(run,crnt,1);       // gas
      }
      
      MATERPINI[1] = temp_p;
      MATERRINI[1] = temp_rho;
    } else if (CASE==3130) {    // needle free injection 3D half [ok]

      temp_p   = MATERPINI[1];
      temp_rho = MATERRINI[1];

      R0 = 0.12/1000.0;
      MATERPINI[1]=2.7*pow(10.0,6.0); // high pressure
      MATERRINI[1]=32.0;              // high pressure
      if ( (pow((xw),2.0) + pow((yw),2.0) + pow((zw),2.0) ) <= pow(R0,2.0)) {     // high pressure bubble
        m4_init_S(run,crnt,1);  // bubble
      } else {
        m4_init_S(run,crnt,0);  // Water
      }

      
      dist_wall = 0.895/1000.0;
      R0 = 0.282/1000.0;
      MATERPINI[1]=101325.0;
      MATERRINI[1]=1.225;     // low pressure
      if ( (pow((xw-dist_wall),2.0) + pow((yw),2.0)) <= pow(R0,2.0)) {     // atm pressure gas
        m4_init_S(run,crnt,1);       // gas
      }
      if (xw>=dist_wall ) {
        m4_init_S(run,crnt,1);       // gas
      }
      
      MATERPINI[1] = temp_p;
      MATERRINI[1] = temp_rho;

    } else if (CASE==3140) {

      R0 = 0.06/1000.0;
      MATERPINI[1]=2.7*pow(10.0,6.0); // high pressure
      MATERRINI[1]=32.0;              // high pressure
      if ( (pow((zw),2.0) + pow((xw),2.0) + pow((yw),2.0)) <= pow(R0,2.0)) {     // high pressure bubble
        m4_init_S(run,crnt,1);  // bubble
      } else {
        m4_init_S(run,crnt,0);  // Water
      }

      dist_wall = 0.895/1000.0;
      R0 = 0.282/1000.0;
      MATERPINI[1]=101325.0;
      MATERRINI[1]=1.225;     // low pressure
      if ( (pow((zw-dist_wall),2.0) + pow((xw),2.0) + pow((yw),2.0)) <= pow(R0,2.0)) {     // atm pressure gas
        m4_init_S(run,crnt,1);       // gas
      }

      if (zw>=dist_wall ) {
        m4_init_S(run,crnt,1);       // gas
      }

    } else if (CASE==3150) {

      R0 = 0.03/1000.0;
      MATERPINI[1]=2.7*pow(10.0,6.0); // high pressure
      MATERRINI[1]=32.0;              // high pressure
      if ( (pow((xw),2.0) + pow((zw),2.0) + pow((yw),2.0)) <= pow(R0,2.0)) {     // high pressure bubble
        m4_init_S(run,crnt,1);  // bubble
      } else {
        m4_init_S(run,crnt,0);  // Water
      }

      dist_wall = 0.895/1000.0;
      R0 = 0.31/1000.0;
      MATERPINI[1]=101325.0;
      MATERRINI[1]=1.225;     // low pressure
      if ( (pow((xw-dist_wall),2.0) + pow((zw),2.0) + pow((yw),2.0)) <= pow(R0,2.0)) {     // atm pressure gas
        m4_init_S(run,crnt,1);       // gas
      }

      if (xw>=dist_wall ) {
        m4_init_S(run,crnt,1);       // gas
      }

    } else if (CASE==3400) {    // Injector case 3D andreas 

      m4_init_S(run,crnt,0);
      
    
    } else if (CASE==3410) {    // Injector case 3D Sandia 

      m4_init_S(run,crnt,0);

    } else {
      printf("No case initialization found: Exiting");
      exit(0);
    }


    crnt=crnt->lnxt;
  }
}