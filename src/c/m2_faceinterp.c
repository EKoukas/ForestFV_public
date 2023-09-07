#include "strdata.h"

void m2_faceinterp_V0(struct RUN *run){

	int f,fo;
  double dist_p1_cf,dist_p1_cc,dist_m1_cf,dist_m1_cc;
  double dist_cf,dist_cc;
  struct BRANCH * crnt;

  crnt=run->topo->locl;
  while (crnt!=NULL){
    
    for(f=0; f<(crnt->nlfc); ++f) {

      fo=fcopp(f,crnt->type);
       
      distance_V0(run,crnt,f ,&dist_p1_cf,&dist_p1_cc); // center-face, center-center
      distance_V0(run,crnt,fo,&dist_m1_cf,&dist_m1_cc); //  

      m2_deltavalues_V0(run,crnt,g_delta_back, fo);  // i-1 , i   -> delta_back    (1)
      m2_deltavalues_V0(run,crnt,g_delta_front,f );  // i   , i+1 -> delta_front   (0)
        
      m2_reconctruct_V0(run,crnt,f,dist_p1_cf,dist_p1_cc,dist_m1_cf,dist_m1_cc); // new grad
    }
    
    crnt=crnt->lnxt;
  }
 
}

void m2_faceinterp_V0_nb(struct RUN *run,int inner_outer){

	int f,fo;
  double dist_p1_cf,dist_p1_cc,dist_m1_cf,dist_m1_cc;
  double dist_cf,dist_cc;
  struct BRANCH * crnt;

  crnt=run->topo->locl;
  while (crnt!=NULL){
    if (crnt->el->inner==inner_outer){
      for(f=0; f<(crnt->nlfc); ++f) {

        fo=fcopp(f,crnt->type);
        
        distance_V0(run,crnt,f ,&dist_p1_cf,&dist_p1_cc); // center-face, center-center
        distance_V0(run,crnt,fo,&dist_m1_cf,&dist_m1_cc); //  

        m2_deltavalues_V0(run,crnt,g_delta_back, fo);  // i-1 , i   -> delta_back    (1)
        m2_deltavalues_V0(run,crnt,g_delta_front,f );  // i   , i+1 -> delta_front   (0)
          
        m2_reconctruct_V0(run,crnt,f,dist_p1_cf,dist_p1_cc,dist_m1_cf,dist_m1_cc); // new grad
      }
    }
    crnt=crnt->lnxt;
  }
 	
}

void m2_faceinterp_V1(struct RUN *run) {

	int iv,f,fo;
  double dist_p1_cf,dist_p1_cc;
  struct BRANCH * crnt;

  crnt=run->topo->locl;
  while (crnt!=NULL){

    m2_facevalues_MUSCL(run,crnt,0,0,0,g_S_cons_L);           // gives value to g_S_cons_L: matrix with conservative values from the center of the cell 
    if (PRIMTV==1){ m2_cons2primtv(g_S_cons_L,g_S_prim_L); }  // g_S_cons_L  to g_S_prim_L

    if (g_mesh_change_grad_scheme==1) {
      for(f=0; f<(crnt->nlfc); ++f) {   
        distance_V1(run,crnt,f ,&dist_p1_cf,&dist_p1_cc); // center-face, center-center      
        crnt->cl->dist_cf[f] = dist_p1_cf;
        crnt->cl->dist_cc[f] = dist_p1_cc;
      } 
    }

    for(f=0; f<(crnt->nlfc); ++f) {
      fo=fcopp(f,crnt->type);
      m2_deltavalues_V1(run,crnt,f,fo);  
      m2_reconctruct_V1(run,crnt,f,fo); // new grad
    }
    
    crnt=crnt->lnxt;
  }
	
}

void m2_faceinterp_V1_nb(struct RUN *run,int inner_outer) {

	int iv,f,fo;
  double dist_p1_cf,dist_p1_cc;
  struct BRANCH * crnt;

  crnt=run->topo->locl;
  while (crnt!=NULL){
    if (crnt->el->inner==inner_outer){
      m2_facevalues_MUSCL(run,crnt,0,0,0,g_S_cons_L);           // gives value to g_S_cons_L: matrix with conservative values from the center of the cell 
      if (PRIMTV==1){ m2_cons2primtv(g_S_cons_L,g_S_prim_L); }  // g_S_cons_L  to g_S_prim_L

      if (g_mesh_change_grad_scheme==1) {
        for(f=0; f<(crnt->nlfc); ++f) {       
          distance_V1(run,crnt,f ,&dist_p1_cf,&dist_p1_cc); // center-face, center-center   
          crnt->cl->dist_cf[f] = dist_p1_cf;
          crnt->cl->dist_cc[f] = dist_p1_cc;
        } 
      }

      for(f=0; f<(crnt->nlfc); ++f) {
        fo=fcopp(f,crnt->type);
        m2_deltavalues_V1(run,crnt,f,fo);  
        m2_reconctruct_V1(run,crnt,f,fo); // new grad
      }
    }

    crnt=crnt->lnxt;
  }
 
}

void m2_faceinterp_unstr(struct RUN * run,int inner_outer)  {

  double x0,y0,z0;
  double x1,y1,z1;
  double x2,y2,z2;
  double xc0,yc0,zc0;
  double xc2,yc2,zc2;
  int i,ind,f,ing,ifopp,fo,inopp,iv;
  int nopp;
  int neq_temp;
  double dist_p1_cf,dist_p1_cc,dist_m1_cf,dist_m1_cc;
  double fx,fy,fz;
  double d0,d2;
  double r0,r2;
  double vecx0[3][3];
  double xx[3];
  struct BRANCH * crnt;
  struct BRANCH * fore;
  struct BRANCH * back;

  // P0 // P1 // P2

  crnt=run->topo->locl;
  while (crnt!=NULL){

    if ((crnt->el->inner==inner_outer)||(inner_outer==-1)) {

      x1=0.0;y1=0.0;z1=0.0; // P1 point
      for (ind=0;ind<crnt->nlnd;ind++){
        x1+=crnt->cl->nd[ind].x;
        y1+=crnt->cl->nd[ind].y;
        z1+=crnt->cl->nd[ind].z;
      }

      x1 /= ((double)crnt->nlnd);
      y1 /= ((double)crnt->nlnd);
      z1 /= ((double)crnt->nlnd);    

      for(f=0; f<(crnt->nlfc); ++f) {   //  Fore-face  coordinates 

        xc2=0.0; yc2=0.0; zc2=0.0;      // F2 point;
        for (ind=0;ind<crnt->cl->fc[f].nnds;ind++) {
          xc2 += crnt->cl->nd[fcnd2elnd(f,ind,crnt->type)].x;
          yc2 += crnt->cl->nd[fcnd2elnd(f,ind,crnt->type)].y;
          zc2 += crnt->cl->nd[fcnd2elnd(f,ind,crnt->type)].z;
        }
        xc2 /= ((double)crnt->cl->fc[f].nnds);
        yc2 /= ((double)crnt->cl->fc[f].nnds);
        zc2 /= ((double)crnt->cl->fc[f].nnds);

        //  Fore-point coordinates 
        if (crnt->cl->fc[f].bc==0){
          x2=0.0; y2=0.0; z2=0.0;         // P2 point
          for (ing=0;ing<crnt->nsfc[f];ing++){
            fore=crnt->neigtr[f][ing];
            for (ind=0;ind<fore->nlnd;ind++){
              x2 += fore->cl->nd[ind].x;
              y2 += fore->cl->nd[ind].y;
              z2 += fore->cl->nd[ind].z;
            }
            x2 /= ((double)fore->nlnd);
            y2 /= ((double)fore->nlnd);
            z2 /= ((double)fore->nlnd);
          }
          x2 /= ((double) crnt->nsfc[f]);
          y2 /= ((double) crnt->nsfc[f]);
          z2 /= ((double) crnt->nsfc[f]);
        } else {
          x2=0.0;y2=0.0;z2=0.0;
        }

        d2=fabs( (xc2-x1)*crnt->cl->nx[f] + (yc2-y1)*crnt->cl->ny[f] + (zc2-z1)*crnt->cl->nz[f]); // Fore face distance

        if(crnt->cl->fc[f].bc==0) { // Fore P2 distance
          r2 = fabs((x2-x1)*crnt->cl->nx[f]+(y2-y1)*crnt->cl->ny[f]+(z2-z1)*crnt->cl->nz[f]);
        } else {
          r2 = 2.0*d2;
        } 
        //  if(crnt->cl->fc[f].bc==0){r2=pow((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1),0.5);}
    
        //  P0 evaluation point
        fx =-r2*crnt->cl->nx[f];
        fy =-r2*crnt->cl->ny[f];
        fz =-r2*crnt->cl->nz[f];

        m2_deltavalues_V0(run,crnt,g_delta_front,f);  // i-1 , i   -> delta_back    (1)
        nopp=noppface(f,crnt->type);        // Number of opposite faces based on geometry properties
        for (ifopp=0;ifopp<nopp;ifopp++) {

          fo=oppface(ifopp,f,crnt->type);
          m2_deltavalues_V0(run,crnt,run->delta_temp[ifopp],fo);

          if (crnt->cl->fc[fo].bc==0) { 
            
            x0=0.0; y0=0.0; z0=0.0;
            for (inopp=0;inopp<crnt->nsfc[fo];inopp++) {

              back=crnt->neigtr[fo][inopp];
              for (ind=0;ind<back->nlnd;ind++){
                x0 += back->cl->nd[ind].x;
                y0 += back->cl->nd[ind].y;
                z0 += back->cl->nd[ind].z;
              }

              x0 /= ((double)back->nlnd);
              y0 /= ((double)back->nlnd);
              z0 /= ((double)back->nlnd);
            }

            vecx0[ifopp][0] = x0/((double) crnt->nsfc[fo])-x1;
            vecx0[ifopp][1] = y0/((double) crnt->nsfc[fo])-y1;
            vecx0[ifopp][2] = z0/((double) crnt->nsfc[fo])-z1;

          } else {

            xc0=0.0;yc0=0.0;zc0=0.0;
            for (ind=0;ind<crnt->cl->fc[f].nnds;ind++){
              xc0 += crnt->cl->nd[fcnd2elnd(fo,ind,crnt->type)].x;
              yc0 += crnt->cl->nd[fcnd2elnd(fo,ind,crnt->type)].y;
              zc0 += crnt->cl->nd[fcnd2elnd(fo,ind,crnt->type)].z;
            } 
            
            xc0 /= ((double)crnt->cl->fc[fo].nnds);
            yc0 /= ((double)crnt->cl->fc[fo].nnds);
            zc0 /= ((double)crnt->cl->fc[fo].nnds);
            
            vecx0[ifopp][0] = 2.0*(xc0-x1);
            vecx0[ifopp][1] = 2.0*(yc0-y1);
            vecx0[ifopp][2] = 2.0*(zc0-z1);

          }
        
        } 
        
        if (nopp==1) {  // Hexes
          
          // r0=vecx0[0][0]*fx+vecx0[0][1]*fy+vecx[0][2]*fz;
          r0 = fabs(vecx0[0][0]*crnt->cl->nx[f] + 
                    vecx0[0][1]*crnt->cl->ny[f] + 
                    vecx0[0][2]*crnt->cl->nz[f]);

          //  if (crnt->root==250&&f!=82){printf("branch: %d face %d deltatemp %e %e \n",crnt->root,f,run->delta_temp[f][0],g_delta_back[0]);}
          for(iv=0;iv<NEQ_TEMP;iv++){   
            g_delta_back[iv] = run->delta_temp[0][iv];
          }

        } else if (nopp==2) { // cross product neig coor

          xx[0] = vecx0[0][1]*vecx0[1][2] - vecx0[1][1]*vecx0[0][2]; 
          xx[1] = vecx0[0][2]*vecx0[1][0] - vecx0[1][2]*vecx0[0][0];
          xx[2] = vecx0[0][0]*vecx0[1][1] - vecx0[1][0]*vecx0[0][1];

          for (i=0;i<3;i++) { 
            run->amat[0][i] = vecx0[0][i];  // coor of 1st neig
            run->amat[1][i] = vecx0[1][i];  // coor of 2nd neig
            run->amat[2][i] = xx[i];       
          }

          for(iv=0;iv<NEQ_TEMP;iv++) {    //variables
            run->cmat[0][iv] = run->delta_temp[0][iv];
            run->cmat[1][iv] = run->delta_temp[1][iv];
            run->cmat[2][iv] = 0.0;
          }

          crammer(NEQ_TEMP,run->amat,run->bmat,run->xmat,run->cmat);
          r0 = sqrt(fx*fx + fy*fy + fz*fz);
          d0 = r0/2.0;
          for(iv=0;iv<NEQ_TEMP;iv++) {   

            g_delta_back[iv] = run->xmat[0][iv]*fx + 
                              run->xmat[1][iv]*fy + 
                              run->xmat[2][iv]*fz;
          }
        
        } else if (nopp==3) {
          
          for (i=0;i<3;i++) {  
            run->amat[0][i] = vecx0[0][i];
            run->amat[1][i] = vecx0[1][i];
            run->amat[2][i] = vecx0[2][i];
          }

          for(iv=0;iv<NEQ_TEMP;iv++){    //variables
            run->cmat[0][iv] = run->delta_temp[0][iv];
            run->cmat[1][iv] = run->delta_temp[1][iv];
            run->cmat[2][iv] = run->delta_temp[2][iv];
          }

          crammer(NEQ_TEMP,run->amat,run->bmat,run->xmat,run->cmat);

          r0 = sqrt(fx*fx + fy*fy + fz*fz);
          d0 = r0/2.0;
          for(iv=0;iv<neq_temp;iv++){    //variables
            g_delta_back[iv] = run->xmat[0][iv]*fx + 
                                run->xmat[1][iv]*fy + 
                                run->xmat[2][iv]*fz;
          }
        }
        // if (crnt->root==250&&f==2){printf("branch: %d face %d reconstruct r0 %e d0 %e r2 %e d2 %e g_delta_front %e \n",crnt->root,f,r0,d0,r2,d2,g_delta_front[0]);}
        dist_p1_cf = d2;
        dist_p1_cc = r2;
        dist_m1_cf = d2;
        dist_m1_cc = r0;
        
        m2_reconctruct_V0(run,crnt,f,dist_p1_cf,dist_p1_cc,dist_m1_cf,dist_m1_cc); // new grad

      }
           
    }
    crnt=crnt->lnxt;
  }


  return;

}