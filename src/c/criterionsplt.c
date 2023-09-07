#include "strdata.h"

void criterionsplt (struct RUN* run) {
  
  int v,level_temp;
  struct BRANCH * crnt;

  crnt=run->topo->locl; 
  while (crnt!=NULL){
    level_temp=1;

    if (CRT_UNIFORM==1)  { level_temp = max(level_temp,CRT_LVL_UNIFORM);              }
    if (CRT_GEO_AREA==1) { level_temp = max(level_temp,criterion_geo_area(run,crnt)); }
    
    if ((CRT_GRAD_PRS==1)||(CRT_GRAD_RHO==1)||(CRT_GRAD_VEL==1)) { 
      level_temp = max(level_temp,criterion_grad(run,crnt));     
    }

    if (CRT_ADAPT_INTER_0!=0) { 
      
      for (v=0;v<2;++v){
        if ((crnt->el->SI[v]>=1)&&(crnt->el->SI[v]<CRT_LRS_ADAPT_INTER_0)) {
          level_temp = max(level_temp,CRT_LVL_ADAPT_INTER_0);
        }

        if ((crnt->el->SI[v]>=CRT_LRS_ADAPT_INTER_0) &&
            (crnt->el->SI[v]< CRT_LRS_ADAPT_INTER_1)) {
          level_temp = max(level_temp,CRT_LVL_ADAPT_INTER_1);
        }
      }
             
    }

    crnt->split=level_temp;
    crnt=crnt->lnxt;
  }

}

double criterion_grad(struct RUN *run,struct BRANCH * crnt) {

  int iv,v,f,ineig,ing;
  double x_cur,x_nb;

  double p_hydro_L,p_hydro_R;
  double grad_press,grad_press_max;

  double rho_L,rho_R;
  double grad_rho,grad_rho_max;

  double vel_mag_L,vel_mag_R;
  double grad_vel,grad_vel_max;

  double level_temp;
  
  double epsilon=1.0e-8;

  level_temp=1;

  //======================================================================
  grad_press_max = 0.0;
  grad_rho_max   = 0.0;
  grad_vel_max   = 0.0;
  for(f=0; f<(crnt->nlfc); f++){
    if(crnt->cl->fc[f].bc==0){
 
      for (ing=0;ing<crnt->nsfc[f];ing++) {
        
        switch(MODEL) {
          case 1: m1_facevalues (run,crnt,ing,f,0,0); m1_facevalues (run,crnt,ing,f,0,1); break;
          case 2: m2_facevalues (run,crnt,ing,f,0,0); m2_facevalues (run,crnt,ing,f,0,1); break;
          case 3: m3_facevalues (run,crnt,ing,f,0,0); m3_facevalues (run,crnt,ing,f,0,1); break;
          case 4: m4_facevalues (run,crnt,ing,f,0,0); m4_facevalues (run,crnt,ing,f,0,1); break;
          case 5: m5_facevalues (run,crnt,ing,f,0,0); m5_facevalues (run,crnt,ing,f,0,1); break;
        }

        if (CRT_GRAD_PRS==1) {
          p_hydro_L=0.0; p_hydro_R=0.0;
          for(iv=0; iv<NMATERIALS; ++iv){
            p_hydro_L += run->sol_L->avf[iv]*run->sol_L->p_hydro[iv]; 
            p_hydro_R += run->sol_R->avf[iv]*run->sol_R->p_hydro[iv]; 
          }
          
          grad_press = 0.0;
          x_cur = p_hydro_L;  x_nb  = p_hydro_R;
          if ((fabs(x_nb)>epsilon) || (fabs(x_cur)>epsilon) ){
            grad_press = fabs(x_nb - x_cur)/max(fabs(x_nb),fabs(x_cur)); 
          }
          grad_press_max  = max(grad_press ,grad_press_max);  
        }

        if (CRT_GRAD_RHO==1) {
          rho_L = run->sol_L->r; 
          rho_R = run->sol_R->r;
          
          grad_rho = 0.0;
          x_cur = rho_L;  x_nb  = rho_R;
          if ((fabs(x_nb)>epsilon) || (fabs(x_cur)>epsilon) ){
            grad_rho = fabs(x_nb - x_cur)/max(fabs(x_nb),fabs(x_cur)); 
          }
          grad_rho_max  = max(grad_rho ,grad_rho_max); 
        }

        if (CRT_GRAD_VEL==1) {
          vel_mag_L = sqrt( pow(run->sol_L->u,2.0) + pow(run->sol_L->v,2.0) + pow(run->sol_L->w,2.0)); 
          vel_mag_R = sqrt( pow(run->sol_R->u,2.0) + pow(run->sol_R->v,2.0) + pow(run->sol_R->w,2.0));
          
          grad_vel = 0.0;
          x_cur = vel_mag_L;  x_nb  = vel_mag_R;
          if ((fabs(x_nb)>epsilon) || (fabs(x_cur)>epsilon) ){
            grad_vel = fabs(x_nb - x_cur)/max(fabs(x_nb),fabs(x_cur)); 
          }
          grad_vel_max  = max(grad_rho ,grad_vel_max);  
        }


      }
    }      
  }

  if (CRT_GRAD_PRS==1) {
    if (fabs(grad_press_max)>=CRT_GRAD_PRS_THR) {
      level_temp = CRT_LVL_GRAD_PRS; 
    }
  }

  if (CRT_GRAD_RHO==1) {
    if (fabs(grad_rho_max)>=CRT_GRAD_RHO_THR) {
      level_temp = CRT_LVL_GRAD_RHO; 
    }
  }

  if (CRT_GRAD_VEL==1) {
    if (fabs(grad_vel_max)>=CRT_GRAD_VEL_THR) {
      level_temp = CRT_LVL_GRAD_VEL; 
    }
  }

  return level_temp;
}

double criterion_geo_area(struct RUN *run,struct BRANCH * crnt) {

  int ind,level_temp;
  double rx,ry,rz;
  double R0,x_cent,y_cent,z_cent;
  double x_start_sm,x_sm;
  double y_start_sm,y_sm;
  double x_start_big,x_big;
  double y_start_big,y_big;

  level_temp=1;

  rx=0.0;ry=0.0;rz=0.0;
  for (ind=0;ind<crnt->nlnd;ind++){
    rx+=crnt->cl->nd[ind].x;
    ry+=crnt->cl->nd[ind].y;
    rz+=crnt->cl->nd[ind].z;
  }
  rx /=	((double)crnt->nlnd);
  ry /= ((double)crnt->nlnd);
  rz /= ((double)crnt->nlnd);


  if (MODEL==0) {

  } else if (MODEL==1) {

  } else if (MODEL==2) {

  } else if (MODEL==3) {

  } else if (MODEL==4) {

  } else if (MODEL==5) {

    if (CASE==2051) {

      double a1      = 15.8/1000.0;
      double a2      = (14.4/1000.0)+a1;
      double a3      = 7.2 / 1000.0;

      if( ((rx>=a1) && (rx<=a2)) && (ry<=a3)) {
        level_temp = CRT_LVL_GEO_AREA; 
      }

    } else if (CASE==2370) {
      
      R0     = MSCALE*10.0;
      x_cent = 15.0*R0;
      y_cent = 0.0*R0;
      z_cent = 0.0*R0;
      if ( (pow((rx-x_cent),2.0) + 2.5*pow((ry-y_cent),2.0) + pow((rz-z_cent),2.0)) <= pow((5.0*R0),2.0)) { 
        level_temp = CRT_LVL_GEO_AREA;             
      } 
  
    } else if (CASE==2600) {

      R0 = MSCALE*10.0;   

      double width=0.5;
      double dist_wall = 8.0*R0;
      
      double y_solid_st  = 1.25*R0;
      double y_solid_end = width*R0 + y_solid_st;

      if ( (pow((rx-dist_wall),2.0) + pow((ry-0.0),2.0)) <= pow(R0,2.0)) { 
        level_temp = CRT_LVL_GEO_AREA;
      }

    } else if ((CASE==2700)||(CASE==2701)) {

      x_start_sm = 0.2; x_sm  = 0.1;
      y_start_sm = 0.0; y_sm  = 0.05;         
      
      x_start_big = x_start_sm+x_sm;  x_big = 0.1;
      y_start_big = 0.0;              y_big = 0.25;

      if (((rx>=(0.9*x_start_sm)) && (rx<=(1.1*(x_start_sm+x_sm)))) && ((ry>=(0.9*y_start_sm)) && (ry<=(1.1*(y_start_sm+y_sm))))) {  
        level_temp = CRT_LVL_GEO_AREA;
      } else if (((rx>=(0.9*x_start_big)) && (rx<=(1.1*(x_start_big+x_big)))) && ((ry>=(0.9*y_start_big)) && (ry<=(1.1*(y_start_big+y_big))))) {  
        level_temp = CRT_LVL_GEO_AREA;
      }

    }


  }   

  return level_temp;

}