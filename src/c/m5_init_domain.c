#include "strdata.h"

void m5_init_domain(struct RUN *run){

  int i,iv,ieq,j;
  double xw,yw,zw;
  double AMIN_temp,R0,dist_wall,x_solid,x_start,dist_solid;
  double x_solid_0,x_solid_1;
  double x_cent_V,y_cent_V,z_cent_V,x_cent_B,y_cent_B,z_cent_B,dr;
  double temp_p,temp_rho;
  double a1,a2,a3;
  double a1_0,a2_0,a3_0;
  double a1_1,a2_1,a3_1;
  double length,ratio,width;
  double x_start_sm,x_sm;
  double y_start_sm,y_sm;
  double x_start_big,x_big;
  double y_start_big,y_big;

  struct BRANCH * crnt;
  
  crnt=run->topo->locl;
  while (crnt!=NULL){
    
    xw=crnt->cl->xc;
    yw=crnt->cl->yc;
    zw=crnt->cl->zc;

    for(iv=0;iv<NEQ;iv++){
      crnt->el->S[iv] = 0.0;
    }

    if (CASE==0) {         // solid advection [Solid Gas]
      if((xw)<0.4){
        m5_init_S(run,crnt,0); 
      } else {
        m5_init_S(run,crnt,1);
      }
    } else if (CASE==1)    { // Shock tube Test [Solid Gas]
      if(xw<0.6){
        m5_init_S(run,crnt,0); 
      } else {
        m5_init_S(run,crnt,1);
      }
    } else if (CASE==200)  { // Solid - Gas advection 2D 
      if(((xw-0.5)*(xw-0.5))+(((yw-0.5)*(yw-0.5)))<pow((0.3/2.0),2.0)){    // Initialize circle in the middle
        m5_init_S(run,crnt,0); 
      } else {
        m5_init_S(run,crnt,1);
      }
    } else if (CASE==205)  { // KS - shock

      shock_prof = 0;
      
      a1      = 15.0/1000.0;
      a2      = (7.5/1000.0)+a1;
      a3      = 3.25 / (1000.0);
      x_start = 0.0145;

      if( ((xw>=a1) && (xw<=a2)) && (yw<=a3)) {
        m5_init_S(run,crnt,0);  // Kidney stone
      } else if( xw<=x_start ) {
        m5_init_S_SW(run,crnt,1,x_start,xw);  // Water
      } else {
        m5_init_S(run,crnt,1);  // Water
      } 

    } else if (CASE==2051) { // KS - shock

      shock_prof = 0;

      a1      = 16.0/1000.0;
      a2      = (14.0/1000.0)+a1;
      a3      = 7.0 / 1000.0;
      x_start = 15.5 / 1000.0;

      if( ((xw>=a1) && (xw<=a2)) && (yw<=a3)) {
        m5_init_S(run,crnt,0);  // Kidney stone
      } else if( xw<=x_start ) {
        m5_init_S_SW(run,crnt,1,x_start,xw);  // Water
      } else {
        m5_init_S(run,crnt,1);  // Water
      } 

    } else if (CASE==206)  { // KS - pulse 2 materials
      a1_0 = 15.0/1000.0;
      a2_0 = (7.5/1000.0)+a1;
      a3_0 = 3.25 / (1000.0);

      a1_1 = 15.0/1000.0+(7.5/4000.0);
      a2_1 = (7.5/1000.0)/2.0 + a1_1;
      a3_1 = a3_0 / 2.0;

      x_start=0.0145;
      if( ((xw>=a1_1) && (xw<=a2_1)) && (yw<=a3_1)) {
        m5_init_S(run,crnt,1);  // Kidney stone 1
      } else if( ((xw>=a1_0) && (xw<=a2_0)) && (yw<=a3_0)) {
        m5_init_S(run,crnt,0);  // Kidney stone 2
      } else if( xw<= x_start) {
        m5_init_S_SW(run,crnt,2,x_start,xw);  // Water
      } else {
        m5_init_S(run,crnt,2);  // Water
      }  

    } else if (CASE==220)  { // KS bubble size
    
      R0 = 0.05 /1000.0;  
      length = 2.0;
      ratio  = 1.15;

      a1 = 5.0/1000.0;
      a2 = (length/1000.0)+a1;
      a3 = (length/ratio) / (2000.0);

      dist_wall = a1 - 2.0*R0;
      x_start   = dist_wall - 1.5*R0;

      if ( (pow((xw-dist_wall-0.0),2.0) + pow((yw-0.0),2.0)) <= pow(R0,2.0)) {     
        m5_init_S(run,crnt,2);         // bubble
      } else if( ((xw>=a1) && (xw<=a2)) && (yw<=a3)) {
        m5_init_S(run,crnt,0);         // KS
      } else if ( xw<=x_start ) {
        m5_init_S_SW(run,crnt,1,x_start,xw);  // Water
      } else {
        m5_init_S(run,crnt,1);         // Water
      }
        
    } else if ((CASE>=2100)&&(CASE<=2199)) { // needle free injection + 1 solid
      
      i_liquid = 1;
      i_gas    = 2;

      if ((CASE>=2150)&&(CASE<=2159)) { i_liquid++; i_gas++; } 
      if ((CASE>=2160)&&(CASE<=2169)) { i_liquid++; i_gas++; }
      if ((CASE>=2170)&&(CASE<=2179)) { i_liquid++; i_gas++; }
      
      R0       = NDFI_R0_BB/1000.0;
      temp_p   = MATERPINI[i_gas]; 
      temp_rho = MATERRINI[i_gas];            
      if ( (pow((xw),2.0) + pow((yw),2.0)) <= pow(R0,2.0)) {     // high pressure bubble
       m5_init_S(run,crnt,i_gas);  // gas
      } else {
       m5_init_S(run,crnt,i_liquid);  // liquid
      }

      R0           = NDFI_R0_MNC/1000.0;
      dist_wall    = NDFI_DIST_MNC/1000.0;
      MATERPINI[2] = 101325.0;
      MATERRINI[2] = 1.225;    
      if ( (pow((xw-dist_wall),2.0) + pow((yw),2.0)) <= pow(R0,2.0)) {     // atm pressure gas
       m5_init_S(run,crnt,i_gas);       // gas
      }

      if (xw>=dist_wall ) {
       m5_init_S(run,crnt,i_gas);       // gas
      }

      if ( xw>=(NDFI_DIST_SOLID_0/1000.0;) ) { m5_init_S(run,crnt,0); } // solid 0
      
      if ((CASE>=2150)&&(CASE<=2159)) {
      if ( xw>=(NDFI_DIST_SOLID_1/1000.0;) ) { m5_init_S(run,crnt,1); } // solid 1
      }

      if ((CASE>=2160)&&(CASE<=2169)) {
      if ( xw>=(NDFI_DIST_SOLID_2/1000.0;) ) { m5_init_S(run,crnt,2); } // solid 2
      }

      if ((CASE>=2170)&&(CASE<=2179)) {
      if ( xw>=(NDFI_DIST_SOLID_3/1000.0;) ) { m5_init_S(run,crnt,3); } // solid 3
      }

      MATERPINI[2] = temp_p;
      MATERRINI[2] = temp_rho;   

    } else if (CASE==2370) { // KS bubble attachment / paper_0

      shock_prof = 0;

      R0 = MSCALE*10.0;   
      x_solid   = 16.0*R0;  

      dist_wall = x_solid - S_CENTER*R0;
      x_start   = dist_wall - 1.25*R0;

      if ( (pow((xw-dist_wall),2.0) + pow((yw-0.0),2.0)) <= pow(R0,2.0)) {     
        m5_init_S(run,crnt,2);             // bubble
      } else if ( xw<=x_start ) {
        m5_init_S_SW(run,crnt,1,x_start,xw);       // Water
      } else {
        m5_init_S(run,crnt,1);             // Water
      }

      if (xw>x_solid)  {
        m5_init_S(run,crnt,0);
      }
    } else if (CASE==2600) { // bubble / blood brain barrier
      
      R0 = MSCALE*10.0;   

      double width =0.4;

      dist_wall = 8.0*R0;
      
      double y_solid_st  = 1.25*R0;
      double y_solid_end = width*R0 + y_solid_st;

      if ( (pow((xw-dist_wall),2.0) + pow((yw-0.0),2.0)) <= pow(R0,2.0)) {     
        m5_init_S(run,crnt,2);             // bubble
      } else {
        m5_init_S(run,crnt,1);             // Water
      }

      if ((yw>y_solid_st)&&(yw<y_solid_end))  {
        m5_init_S(run,crnt,0);
      }

    } else if (CASE==2610) { // bubble / blood vessel
      
      double press_dist;
      double press_temp_0,press_temp_1,press_temp_2;

      R0 = MSCALE*10.0;   

      width = 0.4;

      dist_wall = 0.0*R0;
      
      double y_solid_st  = 1.25*R0;
      double y_solid_end = width*R0 + y_solid_st;

      press_temp_0 = MATERPINI[0];
      press_temp_1 = MATERPINI[1];
      press_temp_2 = MATERPINI[2];

      if ((yw>(3.5*R0))&&(yw<(6.0*R0))) {
        press_temp_0 = MATERPINI[0];
        press_dist   = press_temp_0*fabs((6.0-(yw/R0))/2.5) + 101325.0*fabs(((yw/R0)-3.5)/2.5);
        MATERPINI[0] = press_dist; 
        MATERPINI[1] = press_dist;
        MATERPINI[2] = press_dist;
      } else if (yw>(6.0*R0)) {
        press_temp_0 = MATERPINI[0];
        MATERPINI[0] = 101325.0; 
        MATERPINI[1] = 101325.0;
        MATERPINI[2] = 101325.0;
      }
      

      if ( (pow((xw-dist_wall),2.0) + pow((yw-0.0),2.0)) <= pow(R0,2.0)) {     
        m5_init_S(run,crnt,2);             // bubble
      } else {
        m5_init_S(run,crnt,1);             // Water
      }

      if ((yw>y_solid_st)&&(yw<y_solid_end))  {
        m5_init_S(run,crnt,0);
      }

      MATERPINI[0] = press_temp_0;
      MATERPINI[1] = press_temp_1;
      MATERPINI[2] = press_temp_2;
    
    } else if ((CASE>=2700)&&(CASE<=2799)) { // 2 solids test cases
    
      if ((CASE==2700)||(CASE==2701)) { // 1/2 solids projectiles [test for plastic deformation]

        x_start_sm = 0.2; x_sm  = 0.1;
        y_start_sm = 0.0; y_sm  = 0.05;         
        
        x_start_big = x_start_sm+x_sm;  x_big = 0.1;
        y_start_big = 0.0;              y_big = 0.25; 
        
        if (CASE==2700) { // 1 solid
          if (((xw>=(x_start_sm)) && (xw<=(x_start_sm+x_sm))) && ((yw>=(y_start_sm)) && (yw<=(y_start_sm+y_sm)))) {  
            u_init[0] = 800.0;
            m5_init_S(run,crnt,0);  // Copper
          } else if (((xw>=(x_start_big)) && (xw<=(x_start_big+x_big))) && ((yw>=(y_start_big)) && (yw<=(y_start_big+y_big)))) {  
            u_init[0] = 0.0;
            m5_init_S(run,crnt,0);  // Copper
          } else {
            m5_init_S(run,crnt,1);  // Air
          }
        }

        if (CASE==2701) { // 2 solids
          if (((xw>=(x_start_sm)) && (xw<=(x_start_sm+x_sm))) && ((yw>=(y_start_sm)) && (yw<=(y_start_sm+y_sm)))) {
            u_init[0] = 800.0;  
            m5_init_S(run,crnt,0);  // Copper
          } else if (((xw>=(x_start_big)) && (xw<=(x_start_big+x_big))) && ((yw>=(y_start_big)) && (yw<=(y_start_big+y_big)))) {  
            m5_init_S(run,crnt,1);  // aluminium
          } else {
            m5_init_S(run,crnt,2);  // Air
          }
        }

      } else if (CASE==2710) { // Shock induced bubble collapse with 2 solids 

        shock_prof = 0;

        R0 = MSCALE*10.0;   
        x_solid_0   = 16.0*R0;  

        dist_wall = x_solid_0 - S_CENTER*R0;
        x_start   = dist_wall - 1.25*R0;

        if ( (pow((xw-dist_wall),2.0) + pow((yw-0.0),2.0)) <= pow(R0,2.0)) {     
          m5_init_S(run,crnt,3);             // bubble
        } else if ( xw<=x_start ) {
          m5_init_S_SW(run,crnt,2,x_start,xw);       // Water
        } else {
          m5_init_S(run,crnt,2);             // Water
        }

        if (xw>x_solid_0)  {
          m5_init_S(run,crnt,0);
        }

        x_solid_1 = x_solid_0 + 0.25*R0;
        if (xw>x_solid_1)  {
          m5_init_S(run,crnt,1);
        }
        
      } else if (CASE==2711) { // Shock impact on 2 solids 

        shock_prof = 0;

        R0 = MSCALE*10.0;   
        x_solid_0   = 14.0*R0;  

        x_start   = x_solid_0 - 0.10*R0;

        if ( xw<=x_start ) {
          m5_init_S_SW(run,crnt,2,x_start,xw);       // Water
        } else {
          m5_init_S(run,crnt,2);             // Water
        }

        if (xw>x_solid_0)  {
          m5_init_S(run,crnt,0);
        }

        x_solid_1 = x_solid_0 + 0.25*R0;
        if (xw>x_solid_1)  {
          m5_init_S(run,crnt,1);
        }
        
      }

    } else if (CASE==3111) { // needle free injection 3D + solid 
        
      R0 = 40.0/pow(10.0,6.0);
      temp_p   = MATERPINI[2];
      temp_rho = MATERRINI[2];
      MATERPINI[2]=temp_p;      // high pressure
      MATERRINI[2]=temp_rho;    // high pressure
      if ( (pow((xw),2.0) + pow((yw),2.0) + pow((zw-(50.0/pow(10.0,6.0))),2.0) ) <= pow(R0,2.0)) {     // high pressure bubble
        m5_init_S(run,crnt,2);  // bubble
      } else {
        m5_init_S(run,crnt,1);  // Water
      }

      dist_wall = 1350.0/pow(10.0,6.0);
      R0 = 200.0/pow(10.0,6.0);
      MATERPINI[2]=101325.0;
      MATERRINI[2]=1.225;     // low pressure
      if ( (pow((xw-dist_wall),2.0) + pow((yw),2.0)) <= pow(R0,2.0)) {     // atm pressure gas
        m5_init_S(run,crnt,2);       // gas
      }
      if (xw>=dist_wall ) {
        m5_init_S(run,crnt,2);       // gas
      }
      
      dist_solid = 2100.0/pow(10.0,6.0);
      if (xw>=dist_solid ) {
        m5_init_S(run,crnt,0);       
      }

      MATERPINI[2] = temp_p;
      MATERRINI[2] = temp_rho;

    } else if (CASE==3121) { // needle free injection 3D half + solid [ok]
        
        R0 = 40.0/pow(10.0,6.0);
        temp_p   = MATERPINI[2];
        temp_rho = MATERRINI[2];
        MATERPINI[2]=temp_p;      // high pressure
        MATERRINI[2]=temp_rho;    // high pressure
        if ( (pow((xw),2.0) + pow((yw),2.0) + pow((zw),2.0) ) <= pow(R0,2.0)) {     // high pressure bubble
         m5_init_S(run,crnt,2);  // bubble
        } else {
         m5_init_S(run,crnt,1);  // Water
        }

        
        dist_wall = 1350.0/pow(10.0,6.0);
        R0 = 200.0/pow(10.0,6.0);
        MATERPINI[2]=101325.0;
        MATERRINI[2]=1.225;     // low pressure
        if ( (pow((xw-dist_wall),2.0) + pow((yw),2.0)) <= pow(R0,2.0)) {     // atm pressure gas
         m5_init_S(run,crnt,2);       // gas
        }
        if (xw>=dist_wall ) {
         m5_init_S(run,crnt,2);       // gas
        }
        
        dist_solid = 2100.0/pow(10.0,6.0);
        if (xw>=dist_solid ) {
         m5_init_S(run,crnt,0);       // gas
        }

        MATERPINI[2] = temp_p;
        MATERRINI[2] = temp_rho;

    } else if (CASE==3000) { // Blood vessel shock induced injury
        
      R0     = 0.00008;
      length = 12.0*R0;

      x_cent_V = 9.0*R0;
      y_cent_V = 2.0*R0;
      z_cent_V = 2.0*R0;

      x_cent_B = 9.4*R0;
      y_cent_B = 2.0*R0;
      z_cent_B = 2.0*R0;

      x_start   = 7.9*R0; 
    
      dr = 0.05*R0;
      
      shock_prof = 0;
      
      if ( (pow((xw-x_cent_B),2.0) + pow((yw-y_cent_B),2.0) + pow((zw-z_cent_B),2.0)) <= pow((0.5*R0),2.0)) { 
        m5_init_S(run,crnt,2); 
      } else if ( ((pow((xw-x_cent_V),2.0) + pow((zw-z_cent_V),2.0) ) <= pow(R0+dr,2.0)) && 
                  ((pow((xw-x_cent_V),2.0) + pow((zw-z_cent_V),2.0) ) >= pow(R0-dr,2.0)) ) { 
        m5_init_S(run,crnt,0);
      } else if (xw<x_start) {
        m5_init_S_SW(run,crnt,1,x_start,xw);
      } else {
        m5_init_S(run,crnt,1);
      }
    } else if (CASE==3010) { // Blood vessel shock induced injury half
        
      R0     = 0.00008;
      length = 12.0*R0;

      x_cent_V = 9.0*R0;
      y_cent_V = 0.0*R0;
      z_cent_V = 0.0*R0;

      x_cent_B = 9.4*R0;
      y_cent_B = 0.0*R0;
      z_cent_B = 0.0*R0;

      x_start   = 7.9*R0; 
    
      dr = 0.05*R0;
      
      if ( (pow((xw-x_cent_B),2.0) + pow((yw-y_cent_B),2.0) + pow((zw-z_cent_B),2.0)) <= pow((0.5*R0),2.0)) { 
        m5_init_S(run,crnt,2); 
      } else if ( ((pow((xw-x_cent_V),2.0) + pow((zw-z_cent_V),2.0) ) <= pow(R0+dr,2.0)) && 
                  ((pow((xw-x_cent_V),2.0) + pow((zw-z_cent_V),2.0) ) >= pow(R0-dr,2.0)) ) { 
        m5_init_S(run,crnt,0);
      } else if (xw<x_start) {
        m5_init_S_SW(run,crnt,1,x_start,xw);
      } else {
        m5_init_S(run,crnt,1);
      }
    } else {
      printf("No case initialization found: Exiting");
      exit(0);
    }


    crnt=crnt->lnxt;
  }
}