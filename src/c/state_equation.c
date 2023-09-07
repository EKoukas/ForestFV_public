#include "strdata.h"

void m1_barotropic(struct SOL* sol_temp) {

  double nexp=7.15;
  double beta=293526643.0;
  double cbar=1450.0;
  double psat=2339.0;
  double rsat=998.2;

  if (sol_temp->r>=rsat) {
    sol_temp->p_hydro[0] = psat+beta*(pow((sol_temp->r/rsat),nexp)-1.0);
    sol_temp->c          = sqrt(beta*nexp*(pow((sol_temp->r/rsat),(nexp-1)))/rsat);
  } else if (sol_temp->r<rsat) {
    sol_temp->p_hydro[0] = psat+cbar*((1.0/rsat)-(1.0/sol_temp->r));
    sol_temp->c          = sqrt(cbar/((sol_temp->r)*(sol_temp->r)));
  }

}

void m2_barotropic_2fluid(struct SOL* sol_temp) {

  double R_G       = 286.0;
  double T_ref     = 293.15;
  double r_sat     = 998.1618;
  double p_sat     = 2340.0;
  double c_gas     = 343.24;
  double c_liq     = 1482.5;
  double c_mixture = 1.0;

  double c_temp,r_lm,r_gas,a_gas,a_lm;
  double a_coef,b_coef,d_coef;
  double pres_1,pres_2,check_NaN;

  //===============================================================================================================
  c_temp=c_liq;
 
  // quatradic polynomial from pressure
  a_coef = 1.0; // ok
  b_coef = pow(c_temp,2.0)*r_sat - p_sat - sol_temp->r*sol_temp->ymass*R_G*T_ref - sol_temp->r*pow(c_temp,2.0)*(1.0-sol_temp->ymass); // ok 
  d_coef = (p_sat-(r_sat*pow(c_temp,2.0)))* sol_temp->r*sol_temp->ymass*R_G*T_ref;  // ok
  
  check_NaN=sqrt(pow(b_coef,2.0)-4.0*d_coef);
  if (isnan(check_NaN)==1) {
    c_temp=c_mixture;
      
    a_coef = 1.0; // ok
    b_coef = pow(c_temp,2.0)*r_sat - p_sat - sol_temp->r*sol_temp->ymass*R_G*T_ref - sol_temp->r*pow(c_temp,2.0)*(1.0-sol_temp->ymass); // ok 
    d_coef = (p_sat-(r_sat*pow(c_temp,2.0)))* sol_temp->r*sol_temp->ymass*R_G*T_ref;  // ok
      
    check_NaN=sqrt(pow(b_coef,2.0)-4.0*d_coef);
    if (isnan(check_NaN)==1) {
      printf(" double nan \n");
      exit(0);
    } else {
      pres_1 = (-b_coef+sqrt(pow(b_coef,2.0)-4.0*d_coef))/2.0;
      pres_2 = (-b_coef-sqrt(pow(b_coef,2.0)-4.0*d_coef))/2.0;
      sol_temp->p_hydro[0] = max(pres_1,pres_2);
    }

  } else {
    pres_1 = (-b_coef+sqrt(pow(b_coef,2.0)-4.0*d_coef))/2.0;
    pres_2 = (-b_coef-sqrt(pow(b_coef,2.0)-4.0*d_coef))/2.0;
    sol_temp->p_hydro[0] = max(pres_1,pres_2);    
  }

  r_lm  = r_sat + (1.0/pow(c_temp,2.0))*(sol_temp->p_hydro[0]-p_sat);
  r_gas = sol_temp->p_hydro[0]/(R_G*T_ref);

  a_gas = ((sol_temp->ymass*sol_temp->r)/r_gas); // CHECK
  a_lm  = 1.0 - a_gas;                     // CHECK

  //sol_temp->c = c_temp;
  sol_temp->c = max(5.0,sqrt( 1.0/(sol_temp->r* ( (a_lm/(r_lm*c_temp*c_temp)) + (a_gas/(r_gas*c_gas*c_gas)) ) ) ));
 
  //===============================================================================================================

}

void m5_eos_stress(double ** avf_stress,double avf,double pr,double **Gmat,double Gmag,double mu) {
  
  if ((fabs(Gmag)!=0.0) && (mu!=0.0) ) {

    int i,j,k;
    double J1,J2;
    double Gmat2[3][3] = { 0 };
    double stress_solid[9];
    
    for (i=0; i<3; ++i) {
      for (k=0; k<3; ++k) {
        for (j=0; j<3; ++j) {
          Gmat2[i][j] += Gmat[i][k] * Gmat[k][j];
        }
      }
    }

    J1 = Gmat [0][0] + Gmat [1][1] + Gmat [2][2];
    J2 = Gmat2[0][0] + Gmat2[1][1] + Gmat2[2][2];

    //                Gmat^2      - J2/3*I -           |G|^1/3      * (Gmat      - J1/3*I) 
    stress_solid[0] = Gmat2[0][0] - J2/3.0 - pow(fabs(Gmag),1.0/3.0)*(Gmat[0][0] - J1/3.0);
    stress_solid[1] = Gmat2[0][1]          - pow(fabs(Gmag),1.0/3.0)*(Gmat[0][1]         );
    stress_solid[2] = Gmat2[0][2]          - pow(fabs(Gmag),1.0/3.0)*(Gmat[0][2]         );

    stress_solid[3] = Gmat2[1][0]          - pow(fabs(Gmag),1.0/3.0)*(Gmat[1][0]         );
    stress_solid[4] = Gmat2[1][1] - J2/3.0 - pow(fabs(Gmag),1.0/3.0)*(Gmat[1][1] - J1/3.0);
    stress_solid[5] = Gmat2[1][2]          - pow(fabs(Gmag),1.0/3.0)*(Gmat[1][2]         );

    stress_solid[6] = Gmat2[2][0]          - pow(fabs(Gmag),1.0/3.0)*(Gmat[2][0]         );
    stress_solid[7] = Gmat2[2][1]          - pow(fabs(Gmag),1.0/3.0)*(Gmat[2][1]         );
    stress_solid[8] = Gmat2[2][2] - J2/3.0 - pow(fabs(Gmag),1.0/3.0)*(Gmat[2][2] - J1/3.0);

    avf_stress[0][0] = avf*pr - (mu/pow(fabs(Gmag),1.0/6.0))*(stress_solid[0]);
    avf_stress[0][1] =        - (mu/pow(fabs(Gmag),1.0/6.0))*(stress_solid[1]);
    avf_stress[0][2] =        - (mu/pow(fabs(Gmag),1.0/6.0))*(stress_solid[2]);

    avf_stress[1][0] =        - (mu/pow(fabs(Gmag),1.0/6.0))*(stress_solid[3]);
    avf_stress[1][1] = avf*pr - (mu/pow(fabs(Gmag),1.0/6.0))*(stress_solid[4]);  
    avf_stress[1][2] =        - (mu/pow(fabs(Gmag),1.0/6.0))*(stress_solid[5]);

    avf_stress[2][0] =        - (mu/pow(fabs(Gmag),1.0/6.0))*(stress_solid[6]);  
    avf_stress[2][1] =        - (mu/pow(fabs(Gmag),1.0/6.0))*(stress_solid[7]);
    avf_stress[2][2] = avf*pr - (mu/pow(fabs(Gmag),1.0/6.0))*(stress_solid[8]);

  } else {
    
    avf_stress[0][0] = avf*pr;
    avf_stress[0][1] = 0.0;
    avf_stress[0][2] = 0.0;

    avf_stress[1][0] = 0.0;
    avf_stress[1][1] = avf*pr;  
    avf_stress[1][2] = 0.0;

    avf_stress[2][0] = 0.0;  
    avf_stress[2][1] = 0.0;
    avf_stress[2][2] = avf*pr;

  }
  
}

void m5_eos_stress_V2(double ** avf_stress,double avf,double pr,double *Gmat,double Gmag,double mu) {
  
  if ((fabs(Gmag)!=0.0) && (mu!=0.0) ) {

    int i,j,k;
    double J1,J2;
    double Gmat2[9] = { 0 };
    double stress_solid[9] = { 0 };
    
    for (i=0; i<3; ++i) {
      for (k=0; k<3; ++k) {
        for (j=0; j<3; ++j) {
          Gmat2[i*3+j] += Gmat[i*3+k] * Gmat[k*3+j];
        }
      }
    }

    J1 = Gmat[0]  +  Gmat[4]  + Gmat[8];
    J2 = Gmat2[0] +  Gmat2[4] + Gmat2[8];

    //                Gmat^2   - J2/3*I -           |G|^1/3     * (Gmat    - J1/3*I) 
    stress_solid[0] = Gmat2[0] - J2/3.0 - pow(fabs(Gmag),1.0/3.0)*(Gmat[0] - J1/3.0);
    stress_solid[1] = Gmat2[1]          - pow(fabs(Gmag),1.0/3.0)*(Gmat[1]         );
    stress_solid[2] = Gmat2[2]          - pow(fabs(Gmag),1.0/3.0)*(Gmat[2]         );

    stress_solid[3] = Gmat2[3]          - pow(fabs(Gmag),1.0/3.0)*(Gmat[3]         );
    stress_solid[4] = Gmat2[4] - J2/3.0 - pow(fabs(Gmag),1.0/3.0)*(Gmat[4] - J1/3.0);
    stress_solid[5] = Gmat2[5]          - pow(fabs(Gmag),1.0/3.0)*(Gmat[5]         );

    stress_solid[6] = Gmat2[6]          - pow(fabs(Gmag),1.0/3.0)*(Gmat[6]         );
    stress_solid[7] = Gmat2[7]          - pow(fabs(Gmag),1.0/3.0)*(Gmat[7]         );
    stress_solid[8] = Gmat2[8] - J2/3.0 - pow(fabs(Gmag),1.0/3.0)*(Gmat[8] - J1/3.0);

    avf_stress[0][0] = -avf*pr - (mu/pow(fabs(Gmag),1.0/6.0))*(stress_solid[0]);
    avf_stress[0][1] =         - (mu/pow(fabs(Gmag),1.0/6.0))*(stress_solid[1]);
    avf_stress[0][2] =         - (mu/pow(fabs(Gmag),1.0/6.0))*(stress_solid[2]);

    avf_stress[1][0] =         - (mu/pow(fabs(Gmag),1.0/6.0))*(stress_solid[3]);
    avf_stress[1][1] = -avf*pr - (mu/pow(fabs(Gmag),1.0/6.0))*(stress_solid[4]);  
    avf_stress[1][2] =         - (mu/pow(fabs(Gmag),1.0/6.0))*(stress_solid[5]);

    avf_stress[2][0] =         - (mu/pow(fabs(Gmag),1.0/6.0))*(stress_solid[6]);  
    avf_stress[2][1] =         - (mu/pow(fabs(Gmag),1.0/6.0))*(stress_solid[7]);
    avf_stress[2][2] = -avf*pr - (mu/pow(fabs(Gmag),1.0/6.0))*(stress_solid[8]);

  } else {
    
    avf_stress[0][0] = -avf*pr;
    avf_stress[0][1] = 0.0;
    avf_stress[0][2] = 0.0;

    avf_stress[1][0] = 0.0;
    avf_stress[1][1] = -avf*pr;  
    avf_stress[1][2] = 0.0;

    avf_stress[2][0] = 0.0;  
    avf_stress[2][1] = 0.0;
    avf_stress[2][2] = -avf*pr;

  }
  
}

double G_det(double ** Gmat){
    
  double Gmag;
  
  Gmag = Gmat[0][0]*(Gmat[1][1]*Gmat[2][2] - Gmat[1][2]*Gmat[2][1]) - 
         Gmat[0][1]*(Gmat[1][0]*Gmat[2][2] - Gmat[1][2]*Gmat[2][0]) +
         Gmat[0][2]*(Gmat[1][0]*Gmat[2][1] - Gmat[1][1]*Gmat[2][0]);
  
  return Gmag;
}

double m5_elastic_energy_solid(double ** Gmat,double mu) {

  double Gmag,eints;

  Gmag = G_det(Gmat);  

  if ((fabs(Gmag)!=0.0) && (mu!=0.0)) { 

    int i,j,k;
    double tr1,Gmag13;

    double Mat_temp_1[3][3];
    double Mat_temp_2[3][3] = { 0 };
    
    Gmag13 = pow(fabs(Gmag),1.0/3.0);
    if (isnan(Gmag13)==1) {
      Gmag13 = 0.0;
    } 

    // Gmat - |Gmat|^(1/3)I 
    Mat_temp_1[0][0] = Gmat[0][0] - Gmag13 ;
    Mat_temp_1[0][1] = Gmat[0][1];
    Mat_temp_1[0][2] = Gmat[0][2];

    Mat_temp_1[1][0] = Gmat[1][0];
    Mat_temp_1[1][1] = Gmat[1][1] - Gmag13 ;
    Mat_temp_1[1][2] = Gmat[1][2];

    Mat_temp_1[2][0] = Gmat[2][0];
    Mat_temp_1[2][1] = Gmat[2][1];
    Mat_temp_1[2][2] = Gmat[2][2] - Gmag13 ;

    // (Gmat - |Gmat|^(1/3)I)^2
    for (i=0;i<3;++i) {
      for (k=0;k<3;++k) {
        for (j=0;j<3;++j) { 
          Mat_temp_2[i][j] += Mat_temp_1[i][k] * Mat_temp_1[k][j];
        }
      }
    } 

    // tr(Gmat - |Gmat|^(1/3)I)^2
    tr1 = Mat_temp_2[0][0] +  Mat_temp_2[1][1] +  Mat_temp_2[2][2];

    eints = mu/(4.0*pow(fabs(Gmag),1.0/6.0))*tr1;
  } else {
    eints = 0.0;
  }

  return eints;
}

double m5_elastic_energy_solid_V2(double * Gmat,double Gmag,double mu) {

  double eints;

  if (fabs(Gmag)!=0.0) { 

    int i,j,k;
    double tr1,Gmag13;

    double Mat_temp_1[9];
    double Mat_temp_2[9] = { 0 };
    
    Gmag13 = pow(fabs(Gmag),1.0/3.0);
    if (isnan(Gmag13)==1) {
      Gmag13 = 0.0;
    } 

    // Gmat - |Gmat|^(1/3)I 
    Mat_temp_1[0] = Gmat[0] - Gmag13;
    Mat_temp_1[1] = Gmat[1];
    Mat_temp_1[2] = Gmat[2];

    Mat_temp_1[3] = Gmat[3];
    Mat_temp_1[4] = Gmat[4] - Gmag13;
    Mat_temp_1[5] = Gmat[5];

    Mat_temp_1[6] = Gmat[6];
    Mat_temp_1[7] = Gmat[7];
    Mat_temp_1[8] = Gmat[8] - Gmag13;

    // (Gmat - |Gmat|^(1/3)I)^2
    for (i=0;i<3;++i) {
      for (k=0;k<3;++k) {
        for (j=0;j<3;++j) { 
          Mat_temp_2[i*3+j] += Mat_temp_1[i*3+k] * Mat_temp_1[k*3+j];
        }
      }
    } 

    // tr(Gmat - |Gmat|^(1/3)I)^2
    tr1 = Mat_temp_2[0] +  Mat_temp_2[4] +  Mat_temp_2[8];

    eints = mu/(4.0*pow(fabs(Gmag),1.0/6.0))*tr1;
  } else {
    eints = 0.0;
  }

  return eints;
}