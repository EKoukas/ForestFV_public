#include "strdata.h"

void rotate_vector (double * vel_temp,double u,
                                      double v,
                                      double w,double nx,double ny,double nz,
                                               double mx,double my,double mz,
                                               double lx,double ly,double lz,int flag_write,int pros) {

    vel_temp[0] = u*nx + v*ny + w*nz;
    vel_temp[1] = u*mx + v*my + w*mz;
    vel_temp[2] = u*lx + v*ly + w*lz;

    //return vel_temp;
}

void rotate_vector_gtl_v2 (struct RUN *run, struct BRANCH * crnt,int *ifc, double *u,double *v,double *w) {

  run->vel_temp[0] = (*u)*crnt->cl->nx[(*ifc)]   + (*v)*crnt->cl->ny[(*ifc)]   + (*w)*crnt->cl->nz[(*ifc)];
  run->vel_temp[1] = (*u)*crnt->cl->nxt1[(*ifc)] + (*v)*crnt->cl->nyt1[(*ifc)] + (*w)*crnt->cl->nzt1[(*ifc)];
  run->vel_temp[2] = (*u)*crnt->cl->nxt2[(*ifc)] + (*v)*crnt->cl->nyt2[(*ifc)] + (*w)*crnt->cl->nzt2[(*ifc)];
     
}

void rotate_vector_ltg_v2 (struct RUN *run, struct BRANCH * crnt,int *ifc, double *u,double *v,double *w) {

  run->vel_temp[0] = (*u)*crnt->cl->nx[(*ifc)] + (*v)*crnt->cl->nxt1[(*ifc)] + (*w)*crnt->cl->nxt2[(*ifc)];
  run->vel_temp[1] = (*u)*crnt->cl->ny[(*ifc)] + (*v)*crnt->cl->nyt1[(*ifc)] + (*w)*crnt->cl->nyt2[(*ifc)];
  run->vel_temp[2] = (*u)*crnt->cl->nz[(*ifc)] + (*v)*crnt->cl->nzt1[(*ifc)] + (*w)*crnt->cl->nzt2[(*ifc)];
     
}
/*
void rotate_tensor(double* tensor_temp,
                          double** st_temp, double nx, double ny, double nz,
                          double mx, double my, double mz,
                          double lx, double ly, double lz, int flag_write, int pros) {

    int i, j, k;

    double vec[3][3] = { { nx, ny, nz },
                         { mx, my, mz },
                         { lx, ly, lz } };
                         
    double vecT[3][3] = { { nx, mx, lx },
                          { ny, my, ly },
                          { nz, mz, lz } };
                          
    double Mat_temp[3][3] = { { 0.0 } };

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            for (k = 0; k < 3; ++k) {
                Mat_temp[i][j] += vec[i][k] * st_temp[k][j];
            }
        }
    }

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            double temp = 0.0;
            for (k = 0; k < 3; ++k) {
                temp += Mat_temp[i][k] * vecT[k][j];
            }
            Mat_temp[i][j] = temp;
        }
    }

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            tensor_temp[i * 3 + j] = Mat_temp[i][j];
        }
    }

}



void rotate_tensor (double * tensor_temp,
                    double ** st_temp,    double nx,double ny,double nz,
                                          double mx,double my,double mz,
                                          double lx,double ly,double lz,int flag_write,int pros) {

    int i, j, k;

    double vec[3][3] = { { nx, ny, nz },
                         { mx, my, mz },
                         { lx, ly, lz } };
                         
    double vecT[3][3] = { { nx, mx, lx },
                          { ny, my, ly },
                          { nz, mz, lz } };
                          
    double Mat_temp_1[3][3] = { { 0.0 } };

    double Mat_temp_2[3][3] = { { 0.0 } };

    for (i=0;i<3;++i) {
      for (j=0;j<3;++j) {
        for (k=0;k<3;++k) {
          Mat_temp_1[i][j] += vec[i][k] * st_temp[k][j];
        }
      }
    }

    for (i=0;i<3;++i) {
      for (j=0;j<3;++j) {
        for (k=0;k<3;++k) {
          Mat_temp_2[i][j] += Mat_temp_1[i][k] * vecT[k][j];
        }
      }
    }

    tensor_temp[0] = Mat_temp_2[0][0];
    tensor_temp[1] = Mat_temp_2[0][1];
    tensor_temp[2] = Mat_temp_2[0][2];

    tensor_temp[3] = Mat_temp_2[1][0];
    tensor_temp[4] = Mat_temp_2[1][1];
    tensor_temp[5] = Mat_temp_2[1][2];

    tensor_temp[6] = Mat_temp_2[2][0];
    tensor_temp[7] = Mat_temp_2[2][1];
    tensor_temp[8] = Mat_temp_2[2][2];

} 
*/


void rotate_tensor (double * tensor_temp,
                    double ** st_temp,    double nx,double ny,double nz,
                                          double mx,double my,double mz,
                                          double lx,double ly,double lz,int flag_write,int pros) {

  int i,j,k;

  double vec[3][3];
  double vecT[3][3];
  double Mat_temp_1[3][3];
  double Mat_temp_2[3][3];

  vec[0][0] = nx;
  vec[0][1] = ny;
  vec[0][2] = nz;

  vec[1][0] = mx;
  vec[1][1] = my;
  vec[1][2] = mz;

  vec[2][0] = lx;
  vec[2][1] = ly;
  vec[2][2] = lz;

  for (i=0;i<3;++i) {
    for (j=0;j<3;++j) {
      vecT[j][i]       = vec[i][j];
      Mat_temp_1[i][j] = 0.0;
      Mat_temp_2[i][j] = 0.0;
    }
  }

  for (i=0;i<3;++i) {
    for (j=0;j<3;++j) {
      for (k=0;k<3;++k) {
        Mat_temp_1[i][j] += vec[i][k] * st_temp[k][j];
      }
    }
  }

  for (i=0;i<3;++i) {
    for (j=0;j<3;++j) {
      for (k=0;k<3;++k) {
        Mat_temp_2[i][j] += Mat_temp_1[i][k] * vecT[k][j];
      }
    }
  }

  tensor_temp[0] = Mat_temp_2[0][0];
  tensor_temp[1] = Mat_temp_2[0][1];
  tensor_temp[2] = Mat_temp_2[0][2];

  tensor_temp[3] = Mat_temp_2[1][0];
  tensor_temp[4] = Mat_temp_2[1][1];
  tensor_temp[5] = Mat_temp_2[1][2];

  tensor_temp[6] = Mat_temp_2[2][0];
  tensor_temp[7] = Mat_temp_2[2][1];
  tensor_temp[8] = Mat_temp_2[2][2];

} 