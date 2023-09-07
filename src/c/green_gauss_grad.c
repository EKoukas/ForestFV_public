/*
  Green-Gauss scheme with 'linear' interpolation for face values
  ->'bad' approximation for interpolation on different level elements 

  Non-orthogonal correctors to be added
*/

#include "strdata.h"

void green_gauss_grad (struct RUN* run) {

  switch(MODEL){

    case 0: m0_gg_gradcalc(run); break;
    case 1: m1_gg_gradcalc(run); break;
    case 2: m2_gg_gradcalc(run); break;
    case 3: m3_gg_gradcalc(run); break;
    case 4: m4_gg_gradcalc(run); break;
    case 5: m5_gg_gradcalc(run); break;

  }

  communicate_S(run,2); 
}

void m0_gg_gradcalc (struct RUN* run) {
  
  int iv,f,ing,bc_temp;
  double temp_vec[3];
  double diff,corrector,crnt_per,neig_per,temp_dist;
  struct BRANCH * crnt;

  crnt=run->topo->locl; 
  while (crnt!=NULL){

    for(iv=0;iv<NEQ;iv++){
      crnt->el->WG[iv][0]=0.0;
      crnt->el->WG[iv][1]=0.0;
      crnt->el->WG[iv][2]=0.0;
    }

    m0_facevalues(run,crnt,0,0,0,0);
    m0_cons2primtv(run->sol_L->vec,run->sol_L->wec);

    if (g_mesh_change_grad_scheme==1) {
      for(f=0; f<(crnt->nlfc); ++f) {
        //distance_cf(run,crnt,f,ing);
        for (ing=0;ing<crnt->nsfc[f];ing++) {
          //distance_cc(run,crnt,f,ing);  
        }
      } 
    }

    for(f=0;f<crnt->nlfc;f++){  // Loop for faces

      for (ing=0;ing<crnt->nsfc[f];ing++) {
        
        //vector_diff(temp_vec,temp_dist,crnt->cl->vec_cf[f],crnt->cl->vec_cc[f][ing]);

        //corrector = crnt->cl->dist_cc[f][ing]/(temp_dist+crnt->cl->dist_cf[f]);

        // Get center to face vector of neig cell
        //crnt_per = corrector*crnt->cl->dist_cf[f]/crnt->cl->dist_cc[f][ing];
        //neig_per = corrector*     temp_dist      /crnt->cl->dist_cc[f][ing];

        bc_temp = crnt->cl->fc[f].bc;
        m0_facevalues(run,crnt,ing,f,bc_temp,1);
        m0_cons2primtv(run->sol_R->vec,run->sol_R->wec);

	      for(iv=0;iv<NEQ;iv++){

          diff = neig_per*run->sol_R->wec[iv] + crnt_per*run->sol_L->wec[iv];

          crnt->el->WG[iv][0]+=(diff)*crnt->cl->nx[f]*crnt->cl->Area[f][ing]/crnt->cl->Vol;
	        crnt->el->WG[iv][1]+=(diff)*crnt->cl->ny[f]*crnt->cl->Area[f][ing]/crnt->cl->Vol;
	        crnt->el->WG[iv][2]+=(diff)*crnt->cl->nz[f]*crnt->cl->Area[f][ing]/crnt->cl->Vol;
	      }

      }
    }
    
    crnt=crnt->lnxt;
  }
  
}

void m1_gg_gradcalc (struct RUN* run) {
  
  int iv,f,ing,bc_temp;
  struct BRANCH * crnt;

  crnt=run->topo->locl; 
  while (crnt!=NULL){

    for(iv=0;iv<NEQ;iv++){
      crnt->el->WG[iv][0]=0.0;
      crnt->el->WG[iv][1]=0.0;
      crnt->el->WG[iv][2]=0.0;
    }

    for(f=0;f<crnt->nlfc;f++){  // Loop for faces
      for (ing=0;ing<crnt->nsfc[f];ing++) {
        
        bc_temp = crnt->cl->fc[f].bc;
        m1_facevalues(run,crnt,ing,f,0,      0);
        m1_facevalues(run,crnt,ing,f,bc_temp,1);

	      for(iv=0;iv<NEQ;iv++){
          crnt->el->WG[iv][0]+=0.5*(run->sol_R->vec[iv]+run->sol_L->vec[iv])*crnt->cl->nx[f]*crnt->cl->Area[f][ing]/crnt->cl->Vol;
	        crnt->el->WG[iv][1]+=0.5*(run->sol_R->vec[iv]+run->sol_L->vec[iv])*crnt->cl->ny[f]*crnt->cl->Area[f][ing]/crnt->cl->Vol;
	        crnt->el->WG[iv][2]+=0.5*(run->sol_R->vec[iv]+run->sol_L->vec[iv])*crnt->cl->nz[f]*crnt->cl->Area[f][ing]/crnt->cl->Vol;
	      }

      }
    }
    
    crnt=crnt->lnxt;
  }
  
}

void m2_gg_gradcalc (struct RUN* run) {
  
  int iv,f,ing,bc_temp;
  struct BRANCH * crnt;

  crnt=run->topo->locl; 
  while (crnt!=NULL){

    for(iv=0;iv<NEQ;iv++){
      crnt->el->WG[iv][0]=0.0;
      crnt->el->WG[iv][1]=0.0;
      crnt->el->WG[iv][2]=0.0;
    }

    for(f=0;f<crnt->nlfc;f++){  // Loop for faces
      for (ing=0;ing<crnt->nsfc[f];ing++) {
        
        bc_temp = crnt->cl->fc[f].bc;
        m2_facevalues(run,crnt,ing,f,0,      0);
        m2_facevalues(run,crnt,ing,f,bc_temp,1);

	      for(iv=0;iv<NEQ;iv++){
          crnt->el->WG[iv][0]+=0.5*(run->sol_R->vec[iv]+run->sol_L->vec[iv])*crnt->cl->nx[f]*crnt->cl->Area[f][ing]/crnt->cl->Vol;
	        crnt->el->WG[iv][1]+=0.5*(run->sol_R->vec[iv]+run->sol_L->vec[iv])*crnt->cl->ny[f]*crnt->cl->Area[f][ing]/crnt->cl->Vol;
	        crnt->el->WG[iv][2]+=0.5*(run->sol_R->vec[iv]+run->sol_L->vec[iv])*crnt->cl->nz[f]*crnt->cl->Area[f][ing]/crnt->cl->Vol;
	      }

      }
    }
    
    crnt=crnt->lnxt;
  }
  
}

void m3_gg_gradcalc (struct RUN* run) {
  
  int iv,f,ing,bc_temp;
  struct BRANCH * crnt;

  crnt=run->topo->locl; 
  while (crnt!=NULL){

    for(iv=0;iv<NEQ;iv++){
      crnt->el->WG[iv][0]=0.0;
      crnt->el->WG[iv][1]=0.0;
      crnt->el->WG[iv][2]=0.0;
    }

    for(f=0;f<crnt->nlfc;f++){  // Loop for faces
      for (ing=0;ing<crnt->nsfc[f];ing++) {
        
        bc_temp = crnt->cl->fc[f].bc;
        m3_facevalues(run,crnt,ing,f,0,      0);
        m3_facevalues(run,crnt,ing,f,bc_temp,1);

	      for(iv=0;iv<NEQ;iv++){
          crnt->el->WG[iv][0]+=0.5*(run->sol_R->vec[iv]+run->sol_L->vec[iv])*crnt->cl->nx[f]*crnt->cl->Area[f][ing]/crnt->cl->Vol;
	        crnt->el->WG[iv][1]+=0.5*(run->sol_R->vec[iv]+run->sol_L->vec[iv])*crnt->cl->ny[f]*crnt->cl->Area[f][ing]/crnt->cl->Vol;
	        crnt->el->WG[iv][2]+=0.5*(run->sol_R->vec[iv]+run->sol_L->vec[iv])*crnt->cl->nz[f]*crnt->cl->Area[f][ing]/crnt->cl->Vol;
	      }

      }
    }
    
    crnt=crnt->lnxt;
  }
  
}

void m4_gg_gradcalc (struct RUN* run) {
  
  int iv,f,ing,bc_temp;
  struct BRANCH * crnt;

  crnt=run->topo->locl; 
  while (crnt!=NULL){

    for(iv=0;iv<NEQ;iv++){
      crnt->el->WG[iv][0]=0.0;
      crnt->el->WG[iv][1]=0.0;
      crnt->el->WG[iv][2]=0.0;
    }

    for(f=0;f<crnt->nlfc;f++){  // Loop for faces
      for (ing=0;ing<crnt->nsfc[f];ing++) {
        
        bc_temp = crnt->cl->fc[f].bc;
        m4_facevalues(run,crnt,ing,f,0,      0);
        m4_facevalues(run,crnt,ing,f,bc_temp,1);

	      for(iv=0;iv<NEQ;iv++){
          crnt->el->WG[iv][0]+=0.5*(run->sol_R->vec[iv]+run->sol_L->vec[iv])*crnt->cl->nx[f]*crnt->cl->Area[f][ing]/crnt->cl->Vol;
	        crnt->el->WG[iv][1]+=0.5*(run->sol_R->vec[iv]+run->sol_L->vec[iv])*crnt->cl->ny[f]*crnt->cl->Area[f][ing]/crnt->cl->Vol;
	        crnt->el->WG[iv][2]+=0.5*(run->sol_R->vec[iv]+run->sol_L->vec[iv])*crnt->cl->nz[f]*crnt->cl->Area[f][ing]/crnt->cl->Vol;
	      }

      }
    }
    
    crnt=crnt->lnxt;
  }
  
}

void m5_gg_gradcalc (struct RUN* run) {
  
  int iv,f,ing,bc_temp;
  struct BRANCH * crnt;

  crnt=run->topo->locl; 
  while (crnt!=NULL){

    for(iv=0;iv<NEQ;iv++){
      crnt->el->WG[iv][0]=0.0;
      crnt->el->WG[iv][1]=0.0;
      crnt->el->WG[iv][2]=0.0;
    }

    for(f=0;f<crnt->nlfc;f++){  // Loop for faces
      for (ing=0;ing<crnt->nsfc[f];ing++) {
        
        bc_temp = crnt->cl->fc[f].bc;
        m5_facevalues(run,crnt,ing,f,0,      0);
        m5_facevalues(run,crnt,ing,f,bc_temp,1);

	      for(iv=0;iv<NEQ;iv++){
          crnt->el->WG[iv][0]+=0.5*(run->sol_R->vec[iv]+run->sol_L->vec[iv])*crnt->cl->nx[f]*crnt->cl->Area[f][ing]/crnt->cl->Vol;
	        crnt->el->WG[iv][1]+=0.5*(run->sol_R->vec[iv]+run->sol_L->vec[iv])*crnt->cl->ny[f]*crnt->cl->Area[f][ing]/crnt->cl->Vol;
	        crnt->el->WG[iv][2]+=0.5*(run->sol_R->vec[iv]+run->sol_L->vec[iv])*crnt->cl->nz[f]*crnt->cl->Area[f][ing]/crnt->cl->Vol;
	      }

      }
    }
    
    crnt=crnt->lnxt;
  }
  
}
