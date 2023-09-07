#include "strdata.h"

void distance_V0(struct RUN *run,struct BRANCH *crnt,int f,double* dist_cf,double* dist_cc) {

  int ind;
  double dist_temp;
  double xc,yc,zc,nd_temp;
  double xc_neig,yc_neig,zc_neig;
  struct BRANCH * neig_temp;
  struct CELL * cl_temp;
  
  dist_temp=0.0;
  *dist_cf=0.0;
  *dist_cc=0.0;

  neig_temp=crnt->neigtr[f][0];

  if (neig_temp->level-crnt->level==0) { // Same level   

    cl_temp=crnt->cl;
    nd_temp=((double)crnt->nlnd);

    xc=0.0;yc=0.0;zc=0.0;
    for (ind=0;ind<crnt->nlnd;ind++){
      xc+=cl_temp->nd[ind].x/nd_temp;
      yc+=cl_temp->nd[ind].y/nd_temp;
      zc+=cl_temp->nd[ind].z/nd_temp;
    }

    cl_temp=neig_temp->cl;
    nd_temp=((double)neig_temp->nlnd);

    xc_neig=0.0;yc_neig=0.0;zc_neig=0.0;
    for (ind=0;ind<neig_temp->nlnd;ind++){
      xc_neig+=cl_temp->nd[ind].x/nd_temp;
      yc_neig+=cl_temp->nd[ind].y/nd_temp;
      zc_neig+=cl_temp->nd[ind].z/nd_temp;
    }

    dist_temp = 0.5*sqrt( pow((xc-xc_neig),2.0) + pow((yc-yc_neig),2.0) + pow((zc-zc_neig),2.0));
    *dist_cf  = 0.5*dist_temp;
    *dist_cc  = dist_temp;
  
  } else if (neig_temp->level-crnt->level==1) { // Neigboring element has been refined

    cl_temp=crnt->cl;
    nd_temp=((double)crnt->nlnd);

    xc=0.0;yc=0.0;zc=0.0;
    for (ind=0;ind<crnt->nlnd;ind++){
      xc+=cl_temp->nd[ind].x/nd_temp;
      yc+=cl_temp->nd[ind].y/nd_temp;
      zc+=cl_temp->nd[ind].z/nd_temp;
    }

    cl_temp=neig_temp->prnt->cl;
    nd_temp=((double)neig_temp->prnt->nlnd);

    xc_neig=0.0;yc_neig=0.0;zc_neig=0.0;
    for (ind=0;ind<neig_temp->prnt->nlnd;ind++){
      xc_neig+=cl_temp->nd[ind].x/nd_temp;
      yc_neig+=cl_temp->nd[ind].y/nd_temp;
      zc_neig+=cl_temp->nd[ind].z/nd_temp;
    }

    dist_temp = 0.5*sqrt( pow((xc-xc_neig),2.0) + pow((yc-yc_neig),2.0) + pow((zc-zc_neig),2.0));
    *dist_cf  = 0.5*dist_temp;
    *dist_cc  = 0.5*dist_temp+0.25*dist_temp; 

  } else if ((neig_temp->level-crnt->level==-1) && (neig_temp->level!=0)) { // Current element has been refined

    cl_temp=crnt->prnt->cl;
    nd_temp=((double)crnt->prnt->nlnd);

    xc=0.0;yc=0.0;zc=0.0;
    for (ind=0;ind<crnt->prnt->nlnd;ind++){
      xc+=cl_temp->nd[ind].x/nd_temp;
      yc+=cl_temp->nd[ind].y/nd_temp;
      zc+=cl_temp->nd[ind].z/nd_temp;
    }

    cl_temp=neig_temp->cl;
    nd_temp=((double)neig_temp->nlnd);

    xc_neig=0.0;yc_neig=0.0;zc_neig=0.0;
    for (ind=0;ind<neig_temp->nlnd;ind++){
      xc_neig+=cl_temp->nd[ind].x/nd_temp;
      yc_neig+=cl_temp->nd[ind].y/nd_temp;
      zc_neig+=cl_temp->nd[ind].z/nd_temp;
    }

    dist_temp = 0.25*sqrt( pow((xc-xc_neig),2.0) + pow((yc-yc_neig),2.0) + pow((zc-zc_neig),2.0));
    *dist_cf  = 0.25*dist_temp;
    *dist_cc  = 0.25*dist_temp+0.5*dist_temp;
  } 

  return;

}

void calculate_average_coordinates(struct BRANCH *branch, double *xc, double *yc, double *zc) {
  
  int ind;
  double nd_temp = (double) branch->nlnd;
  struct CELL *cl_temp = branch->cl;

  *xc = 0.0;
  *yc = 0.0;
  *zc = 0.0;
  for (ind = 0; ind < branch->nlnd; ind++) {
    *xc += cl_temp->nd[ind].x;
    *yc += cl_temp->nd[ind].y;
    *zc += cl_temp->nd[ind].z;
  }
  *xc /= nd_temp;
  *yc /= nd_temp;
  *zc /= nd_temp;
}

void distance_V1(struct RUN *run, struct BRANCH *crnt, int f, double *dist_cf, double *dist_cc) {

  double dist_temp;
  double xc, yc, zc, xc_neig, yc_neig, zc_neig;
  struct BRANCH *neig_temp;
  
  dist_temp = 0.0;
  *dist_cf = 0.0;
  *dist_cc = 0.0;

  neig_temp = crnt->neigtr[f][0];

  int level_diff = neig_temp->level - crnt->level;

  if (level_diff == 0) { // Same level
    
    calculate_average_coordinates(crnt, &xc, &yc, &zc);
    calculate_average_coordinates(neig_temp, &xc_neig, &yc_neig, &zc_neig);
    //dist_temp = 0.5 * sqrt(pow((xc - xc_neig), 2.0) + pow((yc - yc_neig), 2.0) + pow((zc - zc_neig), 2.0));
    dist_temp = sqrt(pow((xc - xc_neig), 2.0) + pow((yc - yc_neig), 2.0) + pow((zc - zc_neig), 2.0));
    *dist_cf = 0.5 * dist_temp;
    *dist_cc = dist_temp;

  } else if (level_diff==1) { // 1:2

    struct BRANCH *neig_ref = neig_temp->prnt;

    calculate_average_coordinates(crnt, &xc, &yc, &zc);
    calculate_average_coordinates(neig_ref, &xc_neig, &yc_neig, &zc_neig);
    //dist_temp = 0.5 * sqrt(pow((xc - xc_neig), 2.0) + pow((yc - yc_neig), 2.0) + pow((zc - zc_neig), 2.0));
    dist_temp = sqrt(pow((xc - xc_neig), 2.0) + pow((yc - yc_neig), 2.0) + pow((zc - zc_neig), 2.0));
    *dist_cf = 0.5*dist_temp;
    *dist_cc = 0.5*dist_temp+0.25*dist_temp; 
  
  } else if ((level_diff==-1) && (neig_temp->level!=0)) { // 2:1

    struct BRANCH *crnt_ref = crnt->prnt;

    calculate_average_coordinates(crnt_ref, &xc, &yc, &zc);
    calculate_average_coordinates(neig_temp, &xc_neig, &yc_neig, &zc_neig);
    //dist_temp = 0.25 * sqrt(pow((xc - xc_neig), 2.0) + pow((yc - yc_neig), 2.0) + pow((zc - zc_neig), 2.0));
    dist_temp = sqrt(pow((xc - xc_neig), 2.0) + pow((yc - yc_neig), 2.0) + pow((zc - zc_neig), 2.0));
    *dist_cf = 0.25*dist_temp;
    *dist_cc = 0.25*dist_temp+0.5*dist_temp;

  }

}
/*
void calculate_face_center(struct BRANCH *crnt, int f, double *xf, double *yf, double *zf) {

  int nnds_temp = crnt->cl->fc[f].nnds;

  xf=0.0; yf=0.0; zf=0.0;       

  for (ind=0;ind<nnds_temp;ind++) {
    xf += crnt->cl->nd[fcnd2elnd(f,ind,crnt->type)].x;
    yf += crnt->cl->nd[fcnd2elnd(f,ind,crnt->type)].y;
    zf += crnt->cl->nd[fcnd2elnd(f,ind,crnt->type)].z;
  }
  xf /= ((double)nnds_temp);
  yf /= ((double)nnds_temp);
  zf /= ((double)nnds_temp);

}

void distance_V2(struct RUN *run, struct BRANCH *crnt, int f) {

  double xc,yc,zc;
  double xf,yf,zf;
  double xc_neig,yc_neig,zc_neig;
  struct BRANCH *neig_temp;

  neig_temp = crnt->neigtr[f][0];

  int level_diff = neig_temp->level - crnt->level;

  calculate_average_coordinates(crnt, &xc, &yc, &zc);
  calculate_face_center(crnt, f, &xf, &yf, &zf);

  crnt->cl->vec_cf[f][0] = xf - xc;
  crnt->cl->vec_cf[f][1] = yf - yc;
  crnt->cl->vec_cf[f][2] = zf - zc;

  crnt->cl->dist_cf[f] = sqrt( pow(crnt->cl->vec_cf[f][0],2.0) + pow(crnt->cl->vec_cf[f][1],2.0) + pow(crnt->cl->vec_cf[f][2],2.0));

  calculate_average_coordinates(neig_temp, &xc_neig, &yc_neig, &zc_neig);
  crnt->cl->dist_cc[f][0] = sqrt(pow((xc - xc_neig), 2.0) + pow((yc - yc_neig), 2.0) + pow((zc - zc_neig), 2.0));
    
  return;
}

void distance_cf(struct RUN *run, struct BRANCH *crnt, int f) {

  double xc,yc,zc;
  double xf,yf,zf;

  calculate_average_coordinates(crnt, &xc, &yc, &zc);
  calculate_face_center(crnt, f, &xf, &yf, &zf);

  crnt->cl->vec_cf[f][0] = xf - xc;
  crnt->cl->vec_cf[f][1] = yf - yc;
  crnt->cl->vec_cf[f][2] = zf - zc;

  crnt->cl->dist_cf[f] = sqrt( pow(crnt->cl->vec_cf[f][0],2.0) + pow(crnt->cl->vec_cf[f][1],2.0) + pow(crnt->cl->vec_cf[f][2],2.0));

  return;
}

void distance_cc(struct RUN *run, struct BRANCH *crnt, int f, int ing) {

  double xc,yc,zc;
  double xc_neig,yc_neig,zc_neig;
  struct BRANCH *neig_temp;

  neig_temp = crnt->neigtr[f][ing];

  calculate_average_coordinates(crnt, &xc, &yc, &zc);
  calculate_average_coordinates(neig_temp, &xc_neig, &yc_neig, &zc_neig);

  crnt->cl->vec_cc[f][ing][0] = xc_neig - xc;
  crnt->cl->vec_cc[f][ing][1] = yc_neig - yc;
  crnt->cl->vec_cc[f][ing][2] = zc_neig - zc;

  crnt->cl->dist_cc[f][ing] = sqrt(pow((xc - xc_neig), 2.0) + pow((yc - yc_neig), 2.0) + pow((zc - zc_neig), 2.0));
    
  return;
}
*/