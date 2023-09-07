#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<stddef.h>
#include<time.h>
#include<math.h>
#include<sys/types.h>
#include<sys/resource.h>
#include<float.h>
#include<mpi.h>

#include"metis.h"
#include"parmetis.h"

#include "petscsys.h"
#include "petsc.h"
#include "petscksp.h"
#include "petscmat.h"
#include "petscsnes.h"
#include "petscvec.h"
#include "petscts.h"

#include <xmmintrin.h> // Include SSE header

/* Define macros */
#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

#ifndef sign
#define sign(a,b) ((b)<0 ? -fabs(a) : fabs(a))
#endif

#define printmat3_fmt(A)	#A " =\n%g %g %g\n%g %g %g\n%g %g %g\n"

#define printmat3(A)	( printf(printmat3_fmt(A), 					\
	((double *)(A))[3*0+0], ((double *)(A))[3*1+0], ((double *)(A))[3*2+0],		\
	((double *)(A))[3*0+1], ((double *)(A))[3*1+1], ((double *)(A))[3*2+1], 	\
	((double *)(A))[3*0+2], ((double *)(A))[3*1+2], ((double *)(A))[3*2+2]) )

// =====================================================================
// Inputs: (all caps always)
extern int MODEL;
extern int TYPE_INIT;
extern int VERBOSE;
extern int CASE;

// par for cases
extern double S_CENTER;       // c_2370

extern double NDFI_R0_BB;
extern double NDFI_R0_MNC;
extern double NDFI_DIST_MNC;
extern double NDFI_DIST_SOLID;

extern int NS;
extern int NBCOMMS;
extern int ORDER;
extern int GRAD_SCHEME;
extern int PRIMTV;
extern int LIMITER;
extern int RK_STEPS;
extern double DT;
extern int NSTEP;

extern int EXPORT_TYPE;
extern int EXPORT_STEP;
extern int RESTART_STEP;
extern int ADAPTH;
extern int ADAPTHSTEP;
extern int START;
extern int GEO_ADAPTH;

extern int CRT_UNIFORM;
extern int CRT_LVL_UNIFORM;

extern int CRT_GEO_AREA;
extern int CRT_LVL_GEO_AREA;

extern int CRT_GRAD_PRS;
extern int CRT_LVL_GRAD_PRS;
extern double CRT_GRAD_PRS_THR;

extern int CRT_GRAD_RHO;
extern int CRT_LVL_GRAD_RHO;
extern double CRT_GRAD_RHO_THR;

extern int CRT_GRAD_VEL;
extern int CRT_LVL_GRAD_VEL;
extern double CRT_GRAD_VEL_THR;

extern int CRT_ADAPT_INTER_0;
extern int CRT_LVL_ADAPT_INTER_0;
extern int CRT_LRS_ADAPT_INTER_0;

extern int CRT_ADAPT_INTER_1;
extern int CRT_LVL_ADAPT_INTER_1;
extern int CRT_LRS_ADAPT_INTER_1;

extern int ADAPT_INTER_MAX;

extern int LEVEL; // max level from criterions;

extern int RESTART;
extern int TRANSLATOR;
extern double MSCALE;
extern int SPLITMODE;
extern double NXQUAD;
extern double NYQUAD;
extern double NZQUAD;
extern int ZEROU0;
extern int ZEROU1;
extern int ZEROU2;

extern int BCOVERWRT[10];  // boundary conditions overwrite

extern double YMIN;
extern double AVF_LIM;
extern double AMIN;

extern int PRELAX;
extern int PLASTIC;
extern int NMATERIALS; // number of materials

extern int NUEQ;       // number of unique equations
extern int NEQMASS;    // number of mass continuity equations
extern int NEQMOMENT;  // number of momentum equations
extern int NEQENERGY;  // number of total energy
extern int NYMASS;     // number of mass fraction advection equations
extern int NEQVF;      // number of volume fraction advection equations
extern int NEQSENERGY; // number of internal energy equations
extern int NEQAMAT;    // number of equation for deformation matrix

extern int NEQ;        // sum of equations
extern int NPRIMITIV;  // number of primitive variables
extern int NCONSEQ;    // number of conservative equations 
extern int NEQ_TEMP;   // temporary number of equationd for second order

extern int *eqtypn;
extern int *eqtypi;



extern int ninterf;



// model 5


extern double solid;

extern int shock_prof;

extern double GFM_flag;
extern double n_solids;

// Material properties
extern double MATERPINF[10];
extern double MATERMUSH[10];
extern double MATERGAMA[10];
extern double MATERVISC[10];
extern double MATERRINI[10];
extern double MATERPINI[10];
extern double u_init[10];
extern double v_init[10];
extern double w_init[10];
    


// =====================================================================

// =====================================================================
// Global variables
extern int g_istep;            // iteration counter from start of simulation
extern int g_istep_tot;        // Total counter with added restart iterations if any 
extern int g_istep_start;      // 
extern double g_time;          // simulation time (istep*dt) (with added sim time if any)
extern double g_mesh_change_grad_scheme;

extern double * g_limitedvars;
extern double * g_avf;
extern double * g_avf_c;

// 2nd order
extern double * g_S_cons_L;
extern double * g_S_prim_L;

extern double * g_S_cons_R;
extern double * g_S_prim_R;

extern double * g_delta_back;
extern double * g_delta_front;

extern double * g_grad_limited;

extern double * g_S_prim_rec_L;
extern double * g_S_cons_rec_L;

// =====================================================================
// Flux variables

extern double ** g_flux_viscous;
extern double ** g_flux_adv;
extern double ** g_flux_sym_adv;
extern double *  g_flux_HLLC;

// =====================================================================

extern double ** g_sendbuff_S;
extern double ** g_sendbuff_SF;
extern double ** g_recvbuff_S;
extern double ** g_recvbuff_SF;

typedef struct ND{
  int num;
  double x;
  double y;
  double z;
}ND;

typedef struct FC{
  int num;
  int type;
  int nnds;
  int bc;
}FC;

typedef struct LEAF{

  int gnum;
  double *RHS;
  double *S;
  double *SN;
  double **SF;
  double **WG;
  double crith;

  // GFM
  double *gfm_face;
  int *SI;

  // Face only loop
  double **flux_face;
  double **u_face;
  int *flux_flag;

  int inner;

}LEAF;

// 76 bytes x 10M x 32 tasks = 24.32 G
typedef struct TREE{
  int type;               // 4 bytes
  int * vrtx;             // 32 bytes
  int * conf;             // 24 bytes max
  int part;             // 4  bytes max
  struct BRANCH * brch;   // 8 bytes
}TREE;

typedef struct BRANCH{
  int igroup;

  int part;
  int oldpart;
  int level;
  int type;
  int nkeens;
  int nkids;
  int nlnd;
  int nlfc;
  int tot_neigs;
  int * adrs;
  int split;
  int merge;
  int tag;
  int hangn;
  struct BRANCH * prnt;   // Pointer to parent tree
  int root;     // Pointer to Ancestor tree

  struct BRANCH * next;   // Pointer to next tree (NULL if last)
  struct BRANCH * lnxt;   // 
  struct BRANCH * lprv;   // 

  struct BRANCH *** neigtr;
  int ** neignd;
  int * neigfc;
  int * neigag;
  int * lfc_neigs;
  int * nsfc;
  int ** ipartbound;

  struct BRANCH ** keen;
  int * keenfc;
  int * keenfcangle;

  int fcqsplit;
  int nsplit;
  int ctfc;

  struct BRANCH ** kids;
  struct LEAF * el;
  struct CELL * cl;

}BRANCH;

typedef struct CELL{
    double *nx;
    double *ny;
    double *nz;
    double *nxt1;
    double *nyt1;
    double *nzt1;
    double *nxt2;
    double *nyt2;
    double *nzt2;
    double *dist_cf;
    double *dist_cc;
    double **vec_cf;  // vector from: center of cell to center of face
    double Vol;
    double **Area;
    double xc;
    double yc;
    double zc;

    struct ND * nd;
    struct FC * fc;
}CELL;

typedef struct FOREST{
    int partitioned;
    int ntrees;
    int ndepth;
    int nleaves;
    int * lleaves;
    int pleaves;
    int tleaves;
    int nleavestmp;
    int ngnd;   
    int ** proxlfc1;
    int ** proxlfc2;
    int ** proxtags;
    int ** proxdrys;
    int * nprox;
    int ** bufflfc1;
    int ** bufflfc2;
    int ** bufftags;
    int ** buffdrys;
    int * nbuff;
    struct TREE ** drys; 
    struct BRANCH * locl;
    struct BRANCH * glob;
    double * vertex;   // 3 x8 bytes x10M x32 = 7.680G
    double * vertey;   // 
    double * vertez;   // 
    
}FOREST;

typedef struct CON{

  MPI_Comm comm;  
  
  int rank;
  int size; 

  char * casename;
  char * runname;

  int cnamelength;

  char * filend;          //Pointer for filename which includes nodes info
  char * fileel;          //Pointer for filename which includes elements info
  char * filecon;         //Pointer for filename which includes control values
  char * filebr1;         //Pointer for filename which includes boundary conditions 1 info
  char * filebr2;         //Pointer for filename which includes boundary conditions 2 info
  char * filebr3;         //Pointer for filename which includes boundary conditions 3 info
  char * filebr4;         //Pointer for filename which includes boundary conditions 4 info
  
}CON;

typedef struct SOL{
  double r;        //density
  double u;        //velocity u
  double v;        //velocity v
  double w;        //velocity w
  double e;        //energy
  double p;        //pressure
  double c;
  double ymass;
  double e_internal;
  
  double u_temp;
  double v_temp;
  double w_temp;

  double *vec;
  double *wec;
  double *rho;
  double *ra;
  double *Y;
  double *avf;
  double *are;
  double *p_hydro;
  double *c2Y;
  double *c2;
  double *phi;
  
  double *Amat;
  double **st;
  
  double ***avf_stress;
  
  double u_vel;
  double v_vel;
  double w_vel;

  double ux;
  double uy;
  double uz;
  double vx;
  double vy;
  double vz;
  double wx;
  double wy;
  double wz;    
}SOL;


typedef struct RUN{
  
  struct SOL    * sol_L;
  struct SOL    * sol_R;
  struct CON    * con;
  struct PAR    * par;
  struct FOREST * topo;

  // Auxiliary variables 
  double ** Gmat;
  double * ra;
  double * avf;
  double * are;

  double *source_avf;
  double *** source_st;

  double u_nc; // Velocities for non-conservative part in the residual (gradient) and in fluxes
  double v_nc;
  double w_nc;
  double *avf_S;
  double *are_S; 

  double **Amat;
  double **sstensor;
  double **stress_solid;
  double * avf_c;

  // Auxiliary variables for hllc 
  double *rm;
  double *rm_S;
  double *e_S;
  double *Us_temp;
  double *ar_S;

  double * tensor_temp;
  double *vel_temp;       
  double **Amat_S;

  // Auxiliary variables for relaxation
  double *Vec_temp_R;
  double *rhoKs;
  double *pkinit;
  double *ar_R;
  double *are_R;
  double *avf_R;

  double *t0_m;
  double *p0_m;
  double *t1_m;
  double *r_m;

  double temp_u_S,temp_SL,temp_SR;

  double ** delta_temp;

  // face_inter_unstr
  double ** amat; // 3x3
  double ** bmat; // 3x3
  double ** cmat; // 3xNEQ
  double ** xmat; // 3xNEQ

  // Timing variables
  double cputime;
  double stepcputime;
  double Ttot;
  double Tmeshadapt;
  double Tadv;
  double T2nd;
  double Tdist;
  double Tdelta;
  double Trec;
  double T2ndrelax;
  double T2ndcom;
  double T2ndface;

  double Tdist_sum;
  double Tdelta_sum;
  double Trec_sum;
  double T2ndrelax_sum;
  double T2ndcom_sum;
  double T2ndface_sum;
    
  double Trhs;
  double Trelax;
  double Tcomm;
  double Tinterfloc;
  double Tepx;

  double stepcputime_sum;
  double Ttot_sum;
  double Tmeshadapt_sum;
  double Tadv_sum;
  double T2nd_sum;
  double Trhs_sum;
  double Trelax_sum;
  double Tcomm_sum;
  double Tinterfloc_sum;
  double Tepx_sum;
  int tot_leaves;

  double Ttot_S;    double Ttot_E;
  double Tavd_S;    double Tavd_E;
  double Trk0_S;    double Trk0_E;
  double Tfi_S;     double Tfi_E;
  double Tfi1_S;    double Tfi1_E;
  double Tsrecv_S;  double Tsrecv_E;
  double Tfi0_S;    double Tfi0_E;
  double Tsfsend_S; double Tsfsend_E;
  double Trhs1_S;   double Trhs1_E;
  
  double Tsfrecv_S; double Tsfrecv_E;
  double Trhs0_S;   double Trhs0_E;
  double TrkN_S;    double TrkN_E;
  double Tssend_S;  double Tssend_E;


}RUN;


long int treetag (struct TREE * tree);

struct TREE * createtree (int type, int num, int nkids, int nlnd,int nlvs);

int main (int input, char **inputc);
void exp_timeline(struct RUN * run);
// ==================================================================================
//                            Init subroutines
void global_variables_def(struct RUN *run);
void inputset(struct RUN *run);
void m0_input(struct RUN* run);
void m1_input(struct RUN* run);
void m2_input(struct RUN* run);
void m3_input(struct RUN* run);
void m4_input(struct RUN* run);
void m5_input(struct RUN* run);

void init_domain (struct RUN * run);
void m0_init_domain(struct RUN *run);
void m1_init_domain(struct RUN *run);
void m2_init_domain(struct RUN *run);
void m3_init_domain(struct RUN *run);
void m4_init_domain(struct RUN *run);
void m5_init_domain(struct RUN *run);

void m0_init_S(struct RUN * run,struct BRANCH * crnt,double rho, double u,double v,double w,double p);
void m1_init_S(struct RUN * run,struct BRANCH * crnt,double rho, double u,double v,double w);
void m2_init_S(struct RUN * run,struct BRANCH * crnt,double rho, double u,double v,double w,double ymass);
void m3_init_S(struct RUN * run, struct BRANCH * crnt,int mat);
void m4_init_S(struct RUN * run, struct BRANCH * crnt,int mat);
void m5_init_S(struct RUN * run, struct BRANCH * crnt,int mat);

void m4_init_S_SW(struct RUN * run, struct BRANCH * crnt,int mat,double x_start, double xw);
void m5_init_S_SW(struct RUN * run, struct BRANCH * crnt,int mat,double x_start, double xw);

void check_int(int i_temp, int i_bound_low, int i_bound_high, int line);
void check_db(double temp, double bound_low, double bound_high, int line);

void inner_el(struct RUN * run);
void interface_location(struct RUN * run);
// ==================================================================================

// ==================================================================================
//                            Solution subroutines
void loop(struct RUN * run);
void advancetimeexplicit(struct RUN * run);
void advancetimeexplicit_nb(struct RUN * run);

void face_interpolation_V0(struct RUN * run);
void face_interpolation_V1(struct RUN * run);
void face_interpolation_V2(struct RUN * run);
void face_interpolation_WG (struct RUN* run);
void face_interpolation_V0_nb(struct RUN * run);
void face_interpolation_V1_nb(struct RUN * run);
void face_interpolation_V2_nb(struct RUN * run);
void face_interpolation_unstr(struct RUN * run);
void face_interpolation_unstr_nb(struct RUN * run);

void m0_residual(struct RUN * run, struct BRANCH * crnt);
void m1_residual(struct RUN * run, struct BRANCH * crnt);
void m2_residual(struct RUN * run, struct BRANCH * crnt);
void m3_residual(struct RUN * run, struct BRANCH * crnt);
void m4_residual(struct RUN * run, struct BRANCH * crnt);
void m5_residual(struct RUN * run, struct BRANCH * crnt);

void m0_facevalues(struct RUN *run,struct BRANCH * crnt,int in,int ifc,int bc_temp,int imode);
void m1_facevalues(struct RUN *run,struct BRANCH * crnt,int in,int ifc,int bc_temp,int imode);
void m2_facevalues(struct RUN *run,struct BRANCH * crnt,int in,int ifc,int bc_temp,int imode);
void m3_facevalues(struct RUN *run,struct BRANCH * crnt,int in,int ifc,int bc_temp,int imode);
void m4_facevalues(struct RUN *run,struct BRANCH * crnt,int in,int ifc,int bc_temp,int imode);
void m5_facevalues(struct RUN *run,struct BRANCH * crnt,int in,int ifc,int bc_temp,int imode);

void m0_facercnstr(struct RUN *run,struct BRANCH * crnt,int in,int ifc,int bc_temp,int imode);
void m1_facercnstr(struct RUN *run,struct BRANCH * crnt,int in,int ifc,int bc_temp,int imode);
void m2_facercnstr(struct RUN *run,struct BRANCH * crnt,int in,int ifc,int bc_temp,int imode);
void m3_facercnstr(struct RUN *run,struct BRANCH * crnt,int in,int ifc,int bc_temp,int imode);
void m4_facercnstr(struct RUN *run,struct BRANCH * crnt,int in,int ifc,int bc_temp,int imode);
void m5_facercnstr(struct RUN *run,struct BRANCH * crnt,int in,int ifc,int bc_temp,int imode);

void m0_faceinterp_V0(struct RUN *run);  void m0_faceinterp_V0_nb(struct RUN *run,int inner_outer);
void m1_faceinterp_V0(struct RUN *run);  void m1_faceinterp_V0_nb(struct RUN *run,int inner_outer);
void m2_faceinterp_V0(struct RUN *run);  void m2_faceinterp_V0_nb(struct RUN *run,int inner_outer);
void m3_faceinterp_V0(struct RUN *run);  void m3_faceinterp_V0_nb(struct RUN *run,int inner_outer);
void m4_faceinterp_V0(struct RUN *run);  void m4_faceinterp_V0_nb(struct RUN *run,int inner_outer);
void m5_faceinterp_V0(struct RUN *run);  void m5_faceinterp_V0_nb(struct RUN *run,int inner_outer);

void m0_faceinterp_V1(struct RUN *run); void m0_faceinterp_V1_nb(struct RUN *run,int inner_outer);
void m1_faceinterp_V1(struct RUN *run); void m1_faceinterp_V1_nb(struct RUN *run,int inner_outer);
void m2_faceinterp_V1(struct RUN *run); void m2_faceinterp_V1_nb(struct RUN *run,int inner_outer);
void m3_faceinterp_V1(struct RUN *run); void m3_faceinterp_V1_nb(struct RUN *run,int inner_outer);
void m4_faceinterp_V1(struct RUN *run); void m4_faceinterp_V1_nb(struct RUN *run,int inner_outer);
void m5_faceinterp_V1(struct RUN *run); void m5_faceinterp_V1_nb(struct RUN *run,int inner_outer);

void m0_faceinterp_V2(struct RUN *run); void m0_faceinterp_V2_nb(struct RUN *run,int inner_outer);
void m1_faceinterp_V2(struct RUN *run); void m1_faceinterp_V2_nb(struct RUN *run,int inner_outer);
void m2_faceinterp_V2(struct RUN *run); void m2_faceinterp_V2_nb(struct RUN *run,int inner_outer);
void m3_faceinterp_V2(struct RUN *run); void m3_faceinterp_V2_nb(struct RUN *run,int inner_outer);
void m4_faceinterp_V2(struct RUN *run); void m4_faceinterp_V2_nb(struct RUN *run,int inner_outer);
void m5_faceinterp_V2(struct RUN *run); void m5_faceinterp_V2_nb(struct RUN *run,int inner_outer);

void m0_faceinterp_unstr(struct RUN * run,int inner_outer);
void m1_faceinterp_unstr(struct RUN * run,int inner_outer);
void m2_faceinterp_unstr(struct RUN * run,int inner_outer);
void m3_faceinterp_unstr(struct RUN * run,int inner_outer);
void m4_faceinterp_unstr(struct RUN * run,int inner_outer);
void m5_faceinterp_unstr(struct RUN * run,int inner_outer);

void m0_cons2primtv(double * S_cons,double * S_prim);
void m1_cons2primtv(double * S_cons,double * S_prim);
void m2_cons2primtv(double * S_cons,double * S_prim);
void m3_cons2primtv(double * S_cons,double * S_prim);
void m4_cons2primtv(double * S_cons,double * S_prim);
void m5_cons2primtv(double * S_cons,double * S_prim);

void m0_prim2conc(double * S_prim,double * S_cons);
void m1_prim2conc(double * S_prim,double * S_cons);
void m2_prim2conc(double * S_prim,double * S_cons);
void m3_prim2conc(double * S_prim,double * S_cons);
void m4_prim2conc(double * S_prim,double * S_cons);
void m5_prim2conc(double * S_prim,double * S_cons);

void grad_scheme    (struct RUN* run); 
void grad_scheme_nb (struct RUN* run);

void green_gauss_grad(struct RUN* run);
void least_squares_grad(struct RUN* run);

void m0_gg_gradcalc (struct RUN* run);
void m1_gg_gradcalc (struct RUN* run);
void m2_gg_gradcalc (struct RUN* run);
void m3_gg_gradcalc (struct RUN* run);
void m4_gg_gradcalc (struct RUN* run);
void m5_gg_gradcalc (struct RUN* run);

void face_interpolation_SG (struct RUN* run);

void m0_deltavalues_V0(struct RUN *run,struct BRANCH *crnt,double * delta_temp,int f_temp);
void m1_deltavalues_V0(struct RUN *run,struct BRANCH *crnt,double * delta_temp,int f_temp);
void m2_deltavalues_V0(struct RUN *run,struct BRANCH *crnt,double * delta_temp,int f_temp);
void m3_deltavalues_V0(struct RUN *run,struct BRANCH *crnt,double * delta_temp,int f_temp);
void m4_deltavalues_V0(struct RUN *run,struct BRANCH *crnt,double * delta_temp,int f_temp);
void m5_deltavalues_V0(struct RUN *run,struct BRANCH *crnt,double * delta_temp,int f_temp);

void m0_deltavalues_V1(struct RUN *run,struct BRANCH *crnt,int f,int fo);
void m1_deltavalues_V1(struct RUN *run,struct BRANCH *crnt,int f,int fo);
void m2_deltavalues_V1(struct RUN *run,struct BRANCH *crnt,int f,int fo);
void m3_deltavalues_V1(struct RUN *run,struct BRANCH *crnt,int f,int fo);
void m4_deltavalues_V1(struct RUN *run,struct BRANCH *crnt,int f,int fo);
void m5_deltavalues_V1(struct RUN *run,struct BRANCH *crnt,int f,int fo);

void m0_facevalues_MUSCL(struct RUN *run,struct BRANCH * crnt,int ifc,int iside,int bc_temp, double * S_temp);
void m1_facevalues_MUSCL(struct RUN *run,struct BRANCH * crnt,int ifc,int iside,int bc_temp, double * S_temp);
void m2_facevalues_MUSCL(struct RUN *run,struct BRANCH * crnt,int ifc,int iside,int bc_temp, double * S_temp);
void m3_facevalues_MUSCL(struct RUN *run,struct BRANCH * crnt,int ifc,int iside,int bc_temp, double * S_temp);
void m4_facevalues_MUSCL(struct RUN *run,struct BRANCH * crnt,int ifc,int iside,int bc_temp, double * S_temp);
void m5_facevalues_MUSCL(struct RUN *run,struct BRANCH * crnt,int ifc,int iside,int bc_temp, double * S_temp);

void m0_reconctruct_V0(struct RUN *run,struct BRANCH *crnt,int f,double dist_p1_cf,double dist_p1_cc,double dist_m1_cf,double dist_m1_cc);
void m1_reconctruct_V0(struct RUN *run,struct BRANCH *crnt,int f,double dist_p1_cf,double dist_p1_cc,double dist_m1_cf,double dist_m1_cc);
void m2_reconctruct_V0(struct RUN *run,struct BRANCH *crnt,int f,double dist_p1_cf,double dist_p1_cc,double dist_m1_cf,double dist_m1_cc);
void m3_reconctruct_V0(struct RUN *run,struct BRANCH *crnt,int f,double dist_p1_cf,double dist_p1_cc,double dist_m1_cf,double dist_m1_cc);
void m4_reconctruct_V0(struct RUN *run,struct BRANCH *crnt,int f,double dist_p1_cf,double dist_p1_cc,double dist_m1_cf,double dist_m1_cc);
void m5_reconctruct_V0(struct RUN *run,struct BRANCH *crnt,int f,double dist_p1_cf,double dist_p1_cc,double dist_m1_cf,double dist_m1_cc);

void m0_reconctruct_V1(struct RUN *run,struct BRANCH *crnt,int f, int fo);
void m1_reconctruct_V1(struct RUN *run,struct BRANCH *crnt,int f, int fo);
void m2_reconctruct_V1(struct RUN *run,struct BRANCH *crnt,int f, int fo);
void m3_reconctruct_V1(struct RUN *run,struct BRANCH *crnt,int f, int fo);
void m4_reconctruct_V1(struct RUN *run,struct BRANCH *crnt,int f, int fo);
void m5_reconctruct_V1(struct RUN *run,struct BRANCH *crnt,int f, int fo);

void m0_hllc(struct RUN *run,struct BRANCH * crnt,int in,int ifc);
void m1_hllc(struct RUN *run,struct BRANCH * crnt,int in,int ifc);
void m2_hllc(struct RUN *run,struct BRANCH * crnt,int in,int ifc);
void m3_hllc(struct RUN *run,struct BRANCH * crnt,int in,int ifc);
void m4_hllc(struct RUN *run,struct BRANCH * crnt,int in,int ifc);
void m5_hllc(struct RUN *run,struct BRANCH * crnt,int in,int ifc);
void m5_hllc_gfm(struct RUN *run,struct BRANCH * crnt,int in,int ifc);

void m5_RK(int irk,struct RUN *run);

void m0_flux(struct SOL* sol);
void m1_flux(struct SOL* sol);
void m2_flux(struct SOL* sol);
void m3_flux(struct RUN *run,struct SOL* sol,int flag_ifs);
void m4_flux(struct RUN *run,struct SOL* sol,int flag_ifs);
void m5_flux(struct RUN *run,struct SOL* sol,int flag_ifs);
void m5_flux_gfm(struct RUN *run,struct SOL* sol,int flag_ifs);

void flux_viscous(struct RUN *run,struct BRANCH * crnt,struct SOL* sol_L,struct SOL* sol_R);

void m4_relaxation_fluid(struct RUN * run);
void m5_relaxation_solid(struct RUN * run);

void m4_relaxationvec(struct RUN * run,struct BRANCH * crnt);
void m5_relaxationvec(struct RUN * run,struct BRANCH * crnt);


void m0_wave_speeds(struct RUN *run, struct BRANCH * crnt,int ifc,double* S_S,double* SL, double* SR,
                    double* uL, double* vL, double* wL,double* s11L,
                    double* uR, double* vR, double* wR,double* s11R);

void m1_wave_speeds(struct RUN *run, struct BRANCH * crnt,int ifc,double* S_S,double* SL, double* SR,
                    double* uL, double* vL, double* wL,double* s11L,
                    double* uR, double* vR, double* wR,double* s11R);

void m2_wave_speeds(struct RUN *run, struct BRANCH * crnt,int ifc,double* S_S,double* SL, double* SR,
                    double* uL, double* vL, double* wL,double* s11L,
                    double* uR, double* vR, double* wR,double* s11R);

void m3_wave_speeds(struct RUN *run, struct BRANCH * crnt,int ifc,double* u_S,double* SL, double* SR,
                    double* uL, double* vL, double* wL,double* s11L,double* s12L,double* s13L,
										double* uR, double* vR, double* wR,double* s11R,double* s12R,double* s13R);

void m4_wave_speeds(struct RUN *run, struct BRANCH * crnt,int ifc,double* u_S,double* SL, double* SR,
                    double* uL, double* vL, double* wL,double* s11L,double* s12L,double* s13L,
										double* uR, double* vR, double* wR,double* s11R,double* s12R,double* s13R);

void m5_wave_speeds(struct RUN *run, struct BRANCH * crnt,int ifc,double* u_S,double* SL, double* SR,
                    double* uL, double* vL, double* wL,double* s11L,double* s12L,double* s13L,
										double* uR, double* vR, double* wR,double* s11R,double* s12R,double* s13R);

void m5_wave_speeds_gfm(struct RUN *run, struct BRANCH * crnt,int ifc,double* u_S,double* SL, double* SR,
                    double* uL, double* vL, double* wL,double* s11L,double* s12L,double* s13L,
										double* uR, double* vR, double* wR,double* s11R,double* s12R,double* s13R);

void m0_star_region(struct RUN *run, struct BRANCH * crnt, struct SOL * sol,
                    int ifc,int flag_ifs,double S_S,double SL,double SR,
                    double uL,double vL,double wL,double uR,double vR,double wR);

void m1_star_region(struct RUN *run, struct BRANCH * crnt, struct SOL * sol,
                    int ifc,int flag_ifs,double S_S,double SL,double SR,
                    double uL,double vL,double wL,double uR,double vR,double wR);

void m2_star_region(struct RUN *run, struct BRANCH * crnt, struct SOL * sol,
                    int ifc,int flag_ifs,double S_S,double SL,double SR,
                    double uL,double vL,double wL,double uR,double vR,double wR);

void m3_star_region(struct RUN *run, struct BRANCH * crnt, struct SOL * sol,
                            int ifc,int flag_ifs,double u_S,double SL,double SR,
                            double uL,double vL,double wL,double uR,double vR,double wR,
                            double s11L,double s12L,double s13L,double s11R,double s12R,double s13R);

void m4_star_region(struct RUN *run, struct BRANCH * crnt, struct SOL * sol,
                            int ifc,int flag_ifs,double u_S,double SL,double SR,
                            double uL,double vL,double wL,double uR,double vR,double wR,
                            double s11L,double s12L,double s13L,double s11R,double s12R,double s13R);

void m5_star_region(struct RUN *run, struct BRANCH * crnt, struct SOL * sol,
                            int ifc,int flag_ifs,double u_S,double SL,double SR,
                            double uL,double vL,double wL,double uR,double vR,double wR,
                            double s11L,double s12L,double s13L,double s11R,double s12R,double s13R);


void m5_gfm_rel(struct BRANCH * crnt,struct RUN *run);
void m5_gfm(struct RUN *run,struct BRANCH * crnt);
void m5_plastic_rlx(struct BRANCH * crnt,struct RUN * run);

double m4_eos_press(struct RUN * run,double are,double avf,double **Gmat,double mu,double gamma,double pinf);
void m4_eos_stress(struct RUN * run,double ** avf_stress, double are,double avf,double pr,double **Gmat,double Gmag,double mu,double gamma,double pinf);
 
double m4_eos_c2Y(double p, double r,double avf,double ra,double mu,double gamma,double pinf);

void m1_barotropic(struct SOL* sol_temp);
void m2_barotropic_2fluid(struct SOL* sol_temp);

void m5_eos_stress(double ** avf_stress,double avf,double pr,double **Gmat,double Gmag,double mu);
void m5_eos_stress_V2(double ** avf_stress,double avf,double pr,double *Gmat,double Gmag,double mu);
double m5_elastic_energy_solid(double ** Gmat,double mu);
double m5_elastic_energy_solid_V2(double * Gmat,double Gmag,double mu);
double m5_trace_energy(double ** Gmat,double Gmag);

void m0_bc(struct BRANCH * crnt, int p, int f, struct RUN *run);
void m1_bc(struct BRANCH * crnt, int p, int f, struct RUN *run);
void m2_bc(struct BRANCH * crnt, int p, int f, struct RUN *run);
void m3_bc(struct BRANCH * crnt, int p, int f, struct RUN *run);
void m4_bc(struct BRANCH * crnt, int p, int f, struct RUN *run);
void m5_bc(struct BRANCH * crnt, int p, int f, struct RUN *run);

void m0_bc_inlet(struct BRANCH * crnt,int f,struct RUN *run);
void m1_bc_inlet(struct BRANCH * crnt,int f,struct RUN *run);
void m2_bc_inlet(struct BRANCH * crnt,int f,struct RUN *run);
void m3_bc_inlet(struct BRANCH * crnt,int f,struct RUN *run);
void m4_bc_inlet(struct BRANCH * crnt,int f,struct RUN *run);
void m5_bc_inlet(struct BRANCH * crnt,int f,struct RUN *run);

void m4_bc_outlet(struct BRANCH * crnt,int f,struct RUN *run);

void bc_wall(struct BRANCH * crnt,int f,struct RUN *run);
void bc_symmetry(struct BRANCH * crnt,int f,struct RUN *run);
void m5_bc_symmetry(struct BRANCH * crnt,int f,struct RUN *run);
void bc_reflective(struct BRANCH * crnt,int f,struct RUN *run);

void m0_asinvalues(struct RUN *run,struct BRANCH * crnt,struct SOL* sol, double rho, double ru,double rv,double rw,double e);


void exportfield(struct RUN * run);

void exportfield_1D(struct RUN *run);
void m0_exportfield_1D(struct RUN * run);
void m1_exportfield_1D(struct RUN * run);
void m2_exportfield_1D(struct RUN * run);
void m3_exportfield_1D(struct RUN * run);
void m4_exportfield_1D(struct RUN * run);
void m5_exportfield_1D(struct RUN * run);

void exportfield_linked(struct RUN *run);
void m0_exportfield_linked(struct RUN * run);
void m1_exportfield_linked(struct RUN * run);
void m2_exportfield_linked(struct RUN * run);
void m3_exportfield_linked(struct RUN * run);
void m4_exportfield_linked(struct RUN * run);
void m5_exportfield_linked(struct RUN * run);

void export_ultra_fast (struct RUN * run);

void m0_exportfield_test (struct RUN * run);

void buffer_size_malloc(struct RUN * run);
void buffer_size_free(struct RUN * run);
// ==================================================================================
// Mesh and miscelanious
void boundary(struct RUN *run);
void calctreegm (struct RUN * run);
void statistics_calc (int id, struct TREE * tree, struct RUN * run);
void statinit (struct RUN * run);
void statsave (struct RUN * run);
void commmaster   (RUN * run);
void communicate_S (RUN * run, int ifield);
void communicate_nb_S (struct RUN * run,int ifield,int i_action);

void communicate_C (RUN * run);
void createlfc (struct BRANCH * tree);
void createrunstruct (struct RUN * run);
void createtopo (struct RUN * run);
void restructtopo (RUN * run);
struct BRANCH * createbranch (int type, int num, int nkids, int nlnd,int nlvs);
void createutility(RUN *run);
void elconn (struct RUN * run, struct BRANCH * brch,int ifc);
int elnd2fcnd (int ifc, int ind, int typ);


int fcnd2elnd (int ifc, int ifcnd,int typ);
int fcndrot (int ifcnd, int iangl,int typ);
int fcopnd2elnd (int ifc, int ifcnd,int type);
int fcopp (int ifc,int type);
int noppface (int ifc,int type);
int oppface (int i,int ifc,int type);
int fcrefnd2elnd (int ifc, int ifcnd,int typ);
int fcrefnd2Prismkid (int ifc, int ifcnd);
void filenames(struct RUN * run);

void rotate_vector(double * vel_temp,double u,double v,double w,double nx,double ny,double nz,
                                                                double mx,double my,double mz,
                                                                double lx,double ly,double lz,int flag_write,int pros);

void rotate_tensor(double * tensor_temp,double ** st_temp,double nx,double ny,double nz,
                                                          double mx,double my,double mz,
                                                          double lx,double ly,double lz,int flag_write,int pros);



void forestconn (RUN * run);


void leafallocation(struct RUN * run,struct BRANCH *crnt);
void leafdeallocation(RUN * run,struct BRANCH *crnt);


void memallocation(RUN *run);
void split (int linked, RUN * run, BRANCH * tree);
void merge (int linked, RUN * run, BRANCH * tree);

void tagvector(RUN *run);
void tagvectorcalc(RUN *run,struct BRANCH * crnt);
void textcolor(int attr, int fg, int bg);
void translator(RUN *run);
void translator_pw(RUN *run);
void brchconn (BRANCH * tree);


void normalvectorcalc(RUN *run,struct BRANCH * crnt);

void criterioncalc (struct RUN* run);

void criterionsmth (struct RUN* run);
void criterionsplt (struct RUN* run);

double criterion_grad(struct RUN *run,struct BRANCH * crnt);
double criterion_geo_area(struct RUN *run,struct BRANCH * crnt);


void meshadaptation (RUN * run);
void normalvector(struct RUN *run);

void neigsfin (struct RUN * run);
void refineh (RUN * run, int imode);
void setsplit(RUN *run);
void setsplitcalc(struct BRANCH * crnt,RUN *run);
void properties(struct RUN * run);
void RHS(struct RUN * run);
void RHS_nb(struct RUN * run,int inner_outer);

void RK(int v,RUN *run);

void compute_Amat (double * Amat,double asinit);



void restart_load (struct RUN *run);
void restart_save (struct RUN *run);




void distance_pp(struct RUN *run,struct BRANCH *crnt,int f,double* dx_cc,double* dy_cc,double* dz_cc);
void distance_pf(struct RUN *run,struct BRANCH *crnt,int f,double* dx_cf,double* dy_cf,double* dz_cf);

void faceneigfc (struct RUN * run);


void meshproperties(struct RUN *run);
void cellproperties(struct BRANCH * crnt, struct RUN * run);
void volmproperties(struct BRANCH * crnt, struct RUN * run);
void areaproperties(struct BRANCH * crnt, struct RUN * run);
void conservationcalc(RUN *run);

void crammer(int n,double ** a, double ** b,double ** x,double ** c);
double determinant(double **a);
int numlfc(int typ);
int numlnd(int typ);
int numfcnd(int ifc, int typ);
void drysconn (FOREST * topo);
void spawntopo (RUN * topo);
void partitiondrys (RUN * run);

void repartition (RUN * run);
void repartition_parmetis (RUN * run);
void rebuildlist (RUN * run);
void rebuildconn (RUN * run);
void rebuildprox (RUN * run);
void cellallocation (struct BRANCH * brch);
void celldeallocation (struct BRANCH *crnt);
void solallocation (SOL * sol, RUN *run);
struct BRANCH * spawntree (RUN * run,TREE * tree, int tag, int iel, int icon);
void erasetree (struct RUN * run, struct TREE * tree, int iel);
void partitiontopo (RUN * run);
void destroybrch (BRANCH * brch);
int tagaddress (int * adrs,int nlvl);
int ielconn (BRANCH * brch,int ifc);
int tagconn (BRANCH * brch,int ifc);
int tagconnfin (BRANCH * brch,int ifc, int ineig);
void elconnfin (BRANCH * brch,int ifc);
void tag2adr(int * n, int tag, int * nlvl);
void getMemory(struct RUN * run, int* currRealMem, int* peakRealMem, int* currVirtMem, int* peakVirtMem);
void display(RUN * run);

void init_clocks (struct RUN * run);

double G_det(double ** Gmat);


void rotate_vector_gtl_v2 (struct RUN *run, struct BRANCH * crnt,int *ifc, double *u,double *v,double *w);
void rotate_vector_ltg_v2 (struct RUN *run, struct BRANCH * crnt,int *ifc, double *u,double *v,double *w);

double F_pe(struct RUN * run,double pe,int i_paper);
double DF_pe(struct RUN * run,double pe,int i_paper);

void screen_out(struct RUN * run);
double timecpu(double Ttemp, int flag);
void sumtime(struct RUN * run);

void cons2primtv(struct SOL* sol,struct RUN * run);
void cons2primtv_V2(struct SOL* sol,struct RUN * run);


void distance_V0(struct RUN *run,struct BRANCH *crnt,int f,double * dist_cf,double * dist_cc);
void distance_V1(struct RUN *run, struct BRANCH *crnt, int f, double *dist_cf, double *dist_cc);
void distance_V2(struct RUN *run, struct BRANCH *crnt, int f);
void calculate_average_coordinates(struct BRANCH *branch, double *xc, double *yc, double *zc);
void calculate_face_center(struct RUN *run, struct BRANCH *crnt, double *xf, double *yf, double *zf);
void avf_limit(double *avf, int n, double max_limit);


// ==========================================================================================================================
// SVD
/* Computes cross product of 3D vectors x, y and stores the result in z */
static inline void cross(double * restrict z, const double * restrict x, 
	const double * restrict y);

/* Normalizes a 3D vector (with respect to L2) */
static inline void unit3(double * restrict x);

/*
 * Solves for the roots of a monic cubic polynomial with 3 coefficients 
 * ordered by degree that is assumed to have 3 real roots (D <= 0) 
 */
void solvecubic(double * restrict c);

/* Computes the LDUP decomposition in-place */
void ldu3(double * restrict A, int * restrict P);

/* Does the backward-solve step, or U*x = y */
static inline void ldubsolve3(double * restrict x, const double * restrict y, 
	const double * restrict LDU, const int * restrict P);

/* Explicitly computes the SVD of a 3x3 matrix */
void svd3(double * restrict U, double * restrict S, double * restrict V, 
	const double * restrict A);

/* Computes the matrix multiplication C = A*B */
static inline void matmul3(double * restrict C, const double * restrict A, 
	const double * restrict B);

/* Computes the matrix multiplication y = A*x */
static inline void matvec3(double * restrict y, const double * restrict A,
	const double * restrict x);

/* Computes the matrix multiplication AA = A*A^T */
static inline void aat3(double * restrict AA, const double * restrict A);

static inline void cross(double * restrict z, const double * restrict x, const double * restrict y) {
	z[0] = x[1]*y[2]-x[2]*y[1];
	z[1] = -(x[0]*y[2]-x[2]*y[0]);
	z[2] = x[0]*y[1]-x[1]*y[0];
}

static inline void unit3(double * restrict x) {
	double tmp = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	x[0] /= tmp;
	x[1] /= tmp;
	x[2] /= tmp;
}

static inline void ldubsolve3(double * restrict x, const double * restrict y, const double * restrict LDU, const int * restrict P) {
	x[P[2]] = y[2];
	x[P[1]] = y[1] - LDU[3*P[2]+1]*x[P[2]];
	x[P[0]] = y[0] - LDU[3*P[2]+0]*x[P[2]] - LDU[3*P[1]+0]*x[P[1]];
}

static inline void matmul3(double * restrict C, const double * restrict A, const double * restrict B) {
	C[3*0+0] = A[3*0+0]*B[3*0+0] + A[3*1+0]*B[3*0+1] + A[3*2+0]*B[3*0+2];
	C[3*1+0] = A[3*0+0]*B[3*1+0] + A[3*1+0]*B[3*1+1] + A[3*2+0]*B[3*1+2];
	C[3*2+0] = A[3*0+0]*B[3*2+0] + A[3*1+0]*B[3*2+1] + A[3*2+0]*B[3*2+2];

	C[3*0+1] = A[3*0+1]*B[3*0+0] + A[3*1+1]*B[3*0+1] + A[3*2+1]*B[3*0+2];
	C[3*1+1] = A[3*0+1]*B[3*1+0] + A[3*1+1]*B[3*1+1] + A[3*2+1]*B[3*1+2];
	C[3*2+1] = A[3*0+1]*B[3*2+0] + A[3*1+1]*B[3*2+1] + A[3*2+1]*B[3*2+2];

	C[3*0+2] = A[3*0+2]*B[3*0+0] + A[3*1+2]*B[3*0+1] + A[3*2+2]*B[3*0+2];
	C[3*1+2] = A[3*0+2]*B[3*1+0] + A[3*1+2]*B[3*1+1] + A[3*2+2]*B[3*1+2];
	C[3*2+2] = A[3*0+2]*B[3*2+0] + A[3*1+2]*B[3*2+1] + A[3*2+2]*B[3*2+2];
}

static inline void matvec3(double * restrict y, const double * restrict A, const double * restrict x) {
	y[0] = A[3*0+0]*x[0] + A[3*1+0]*x[1] + A[3*2+0]*x[2];
	y[1] = A[3*0+1]*x[0] + A[3*1+1]*x[1] + A[3*2+1]*x[2];
	y[2] = A[3*0+2]*x[0] + A[3*1+2]*x[1] + A[3*2+2]*x[2];
}

static inline void aat3(double * restrict AA, const double * restrict A) {
	AA[3*0+0] = A[3*0+0]*A[3*0+0] + A[3*1+0]*A[3*1+0] + A[3*2+0]*A[3*2+0];
	AA[3*1+0] = A[3*0+0]*A[3*0+1] + A[3*1+0]*A[3*1+1] + A[3*2+0]*A[3*2+1];
	AA[3*2+0] = A[3*0+0]*A[3*0+2] + A[3*1+0]*A[3*1+2] + A[3*2+0]*A[3*2+2];

	AA[3*0+1] = AA[3*1+0];
	AA[3*1+1] = A[3*0+1]*A[3*0+1] + A[3*1+1]*A[3*1+1] + A[3*2+1]*A[3*2+1];
	AA[3*2+1] = A[3*0+1]*A[3*0+2] + A[3*1+1]*A[3*1+2] + A[3*2+1]*A[3*2+2];

	AA[3*0+2] = AA[3*2+0];
	AA[3*1+2] = AA[3*2+1];
	AA[3*2+2] = A[3*0+2]*A[3*0+2] + A[3*1+2]*A[3*1+2] + A[3*2+2]*A[3*2+2];
}
// ==========================================================================================================================
