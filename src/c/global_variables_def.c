#include "strdata.h"

  // =====================================================================
  // Inputs: (all caps always)
	int MODEL;
	int TYPE_INIT;
	int VERBOSE;
	int CASE;
	double S_CENTER;

	double NDFI_R0_BB;
	double NDFI_R0_MNC;
	double NDFI_DIST_MNC;
	double NDFI_DIST_SOLID;

	int NS;
	int NBCOMMS;
	int ORDER;
	int GRAD_SCHEME;
	int PRIMTV;
	int LIMITER;
	int RK_STEPS;
	double DT;
	int NSTEP;

	int EXPORT_TYPE;
	int EXPORT_STEP;
	int RESTART_STEP;
	int ADAPTH;
	int ADAPTHSTEP;
	int START;
	int GEO_ADAPTH;

	int CRT_UNIFORM;
 	int CRT_LVL_UNIFORM;

 	int CRT_GEO_AREA;
 	int CRT_LVL_GEO_AREA;

 	int CRT_GRAD_PRS;
 	int CRT_LVL_GRAD_PRS;
 	double CRT_GRAD_PRS_THR;

	int CRT_GRAD_RHO;
	int CRT_LVL_GRAD_RHO;
	double CRT_GRAD_RHO_THR;

	int CRT_GRAD_VEL;
	int CRT_LVL_GRAD_VEL;
	double CRT_GRAD_VEL_THR;

	int CRT_ADAPT_INTER_0;
 	int CRT_LVL_ADAPT_INTER_0;
 	int CRT_LRS_ADAPT_INTER_0;

 	int CRT_ADAPT_INTER_1;
 	int CRT_LVL_ADAPT_INTER_1;
 	int CRT_LRS_ADAPT_INTER_1;

	int ADAPT_INTER_MAX;

	int LEVEL;

	int RESTART;
	int TRANSLATOR;
	double MSCALE;
	int SPLITMODE;
	double NXQUAD;
	double NYQUAD;
	double NZQUAD;
	int ZEROU0;
	int ZEROU1;
	int ZEROU2;
	
	int BCOVERWRT[10];  // boundary conditions overwrite

	double YMIN;
	double AVF_LIM;
	double AMIN;

	int PRELAX;
	int PLASTIC;
	int NMATERIALS; // number of materials

	int NUEQ;       // number of unique equations
	int NEQMASS;    // number of mass continuity equations
	int NEQMOMENT;  // number of momentum equations
	int NEQENERGY;  // number of total energy
	int NYMASS;     // number of mass fraction advection equations
	int NEQVF;      // number of volume fraction advection equations
	int NEQSENERGY; // number of internal energy equations
	int NEQAMAT;    // number of equation for deformation matrix

	int NEQ;        // sum of equations
	int NPRIMITIV;  // number of primitive variables
	int NCONSEQ;    // number of conservative equations 
	int NEQ_TEMP;   // temporary number of equationd for second order

	int *eqtypn;
	int *eqtypi;



	int ninterf;



	// model 3
	double S_center;
	double width;
	double width_1;
	double solid;

	int shock_prof;

	double GFM_flag;
	double n_solids;

	// Material properties
	double MATERPINF[10];
	double MATERMUSH[10];
	double MATERGAMA[10];
	double MATERVISC[10];
	double MATERRINI[10];
	double MATERPINI[10];
	double u_init[10];
	double v_init[10];
	double w_init[10];
		


	// =====================================================================

	// =====================================================================
	// Global variables
	int g_istep;            // iteration counter from start of simulation
	int g_istep_tot;        // Total counter with added restart iterations if any 
	int g_istep_start;
	double g_time;          // simulation time (istep*dt) (with added sim time if any)
	double g_mesh_change_grad_scheme;
	double * g_limitedvars;
	double * g_avf;
	double * g_avf_c;

	// 2nd order
	double * g_S_cons_L;
	double * g_S_prim_L;

	double * g_S_cons_R;
	double * g_S_prim_R;

	double * g_delta_back;
	double * g_delta_front;

	double * g_grad_limited;

	double * g_S_prim_rec_L;
	double * g_S_cons_rec_L;

	// =====================================================================
	// Flux variables

	double ** g_flux_viscous;
	double ** g_flux_adv;
	double ** g_flux_sym_adv;
	double *  g_flux_HLLC;

	// =====================================================================

	double ** g_sendbuff_S;
	double ** g_sendbuff_SF;
	double ** g_recvbuff_S;
	double ** g_recvbuff_SF;

void global_variables_def(struct RUN *run) {

}