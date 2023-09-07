#include "strdata.h"

void m0_init_S(struct RUN * run,struct BRANCH * crnt,double rho, double u,double v,double w,double p) {

	double e_internal;
	
	e_internal = (p + MATERGAMA[0]*MATERPINF[0])/(rho*(MATERGAMA[0]-1.0)); 

	crnt->el->S[0] = rho;
	crnt->el->S[1] = rho*u;
	crnt->el->S[2] = rho*v;
	crnt->el->S[3] = rho*w;	
	crnt->el->S[4] = rho*( e_internal + 0.5*(pow(u,2.0) + pow(v,2.0) + pow(w,2.0)));
	
}

void m1_init_S(struct RUN * run,struct BRANCH * crnt,double rho, double u,double v,double w) {

	crnt->el->S[0] = rho;
	crnt->el->S[1] = rho*u;
	crnt->el->S[2] = rho*v;
	crnt->el->S[3] = rho*w;
  	
}

void m2_init_S(struct RUN * run,struct BRANCH * crnt,double rho, double u,double v,double w,double ymass) {
	
	crnt->el->S[0] = rho;
	crnt->el->S[1] = rho*u;
	crnt->el->S[2] = rho*v;
	crnt->el->S[3] = rho*w;
	crnt->el->S[4] = rho*ymass;
	
}

void m3_init_S(struct RUN * run, struct BRANCH * crnt,int mat) {
	
	int iv;
  double avf,rho,are;
  
  rho=0.0;
  crnt->el->S[eqtypi[2]] = 0.0;
  for(iv=0; iv<NMATERIALS; ++iv){ 
    if (mat==iv) {avf=1.0 - (double)(NMATERIALS-1)*AMIN; }
    if (mat!=iv) {avf=                             AMIN; }
    crnt->el->S[eqtypi[0]+iv] = avf * MATERRINI[iv];        // ar
    rho                      += avf * MATERRINI[iv];   

    crnt->el->S[eqtypi[2]]   += avf * (MATERPINI[mat] + MATERGAMA[iv]*MATERPINF[iv])/(MATERGAMA[iv]-1.0);  // rE

  }
 
  crnt->el->S[eqtypi[1]+0]  = u_init[mat] * rho; // ru
  crnt->el->S[eqtypi[1]+1]  = v_init[mat] * rho; // rv
  crnt->el->S[eqtypi[1]+2]  = w_init[mat] * rho; // rw

  crnt->el->S[eqtypi[2]] += 0.5*rho*(pow(u_init[mat],2.0) + 
                                     pow(v_init[mat],2.0) +
                                     pow(w_init[mat],2.0));  // rE
  
  for(iv=0; iv<(NMATERIALS-1); ++iv){
    if (mat==iv) {avf=1.0 - (double)(NMATERIALS-1)*AMIN; }
    if (mat!=iv) {avf=                             AMIN; }

    crnt->el->S[eqtypi[3]+iv] = avf;  // avf
  }
	
}

void m4_init_S(struct RUN * run, struct BRANCH * crnt,int mat) {

  int iv;
  double avf,rho,are;
  
  rho=0.0;
  crnt->el->S[eqtypi[2]] = 0.0;
  for(iv=0; iv<NMATERIALS; ++iv){ 
    if (mat==iv) {avf=1.0 - (double)(NMATERIALS-1)*AMIN; }
    if (mat!=iv) {avf=                             AMIN; }
    crnt->el->S[eqtypi[0]+iv] = avf * MATERRINI[iv];        // ar
    rho                      += avf * MATERRINI[iv];
    
    are = avf * (MATERPINI[mat] + MATERGAMA[iv]*MATERPINF[iv])/(MATERGAMA[iv]-1.0);    

    crnt->el->S[eqtypi[2]]   += are;  // rE
    crnt->el->S[eqtypi[4]+iv] = are;  // are

  }
 
  crnt->el->S[eqtypi[1]+0] = u_init[mat] * rho; // ru
  crnt->el->S[eqtypi[1]+1] = v_init[mat] * rho; // rv
  crnt->el->S[eqtypi[1]+2] = w_init[mat] * rho; // rw

  crnt->el->S[eqtypi[2]] += 0.5*rho*(pow(u_init[mat],2.0) + 
                                     pow(v_init[mat],2.0) +
                                     pow(w_init[mat],2.0));  // rE
  
  for(iv=0; iv<(NMATERIALS-1); ++iv){
    if (mat==iv) {avf=1.0 - (double)(NMATERIALS-1)*AMIN; }
    if (mat!=iv) {avf=                             AMIN; }

    crnt->el->S[eqtypi[3]+iv] = avf;  // avf
  }

  //printf("E tot %d | %f \n",crnt->root,crnt->el->S[eqtypi[2]]);

}

void m4_init_S_SW(struct RUN * run, struct BRANCH * crnt,int mat,double x_start, double xw) {

  int iv;
  double avf,rho,are;
  double p0,ps,alpha,omega,c0,pressure,u_vel;
  double beta,t1,kappa,p_max,sc,x_temp,x_start_temp,x_rescaled;
  double PI = 3.14159265358979323846;

  switch(shock_prof){
    case 0: // Jc
      p0    = 101000.0;
      ps    = 35.0e6;
      alpha = 1.48e6;
      omega = 1.21e6;
      c0    = 1624.943090;
      pressure = p0 + 2.0*ps*exp(-alpha*((-xw/c0)+x_start/c0))*cos(omega*((-xw/c0)+x_start/c0) + PI/3.0);
    break;

    case 1: // Wang
      p0    = 101000.0;
      ps    = 39.3e6;
      alpha = 8.0e5;
      omega = 6.4e5;
      c0    = 1624.943090;
      pressure = p0 + 2.0*ps*exp(-alpha*((-xw/c0)+x_start/c0))*cos(omega*((-xw/c0)+x_start/c0) + PI/3.0);
    break;

    case 2:  // SW-A1
      alpha        = 6.0;
      beta         = 1.0;
      t1           = 0.00217;
      kappa        = 1.09;
      p_max        = 20.0e6;
      c0           = 1624.943090;
      sc           = 10e4/c0;
      x_start_temp = -x_start*sc;
      x_rescaled   = (1.0/x_start)*xw - 1.0;
      x_temp       = (-xw*sc-x_start_temp);
      pressure = p0 + p_max*(kappa*(1.0-exp(-(x_temp)/t1))*exp(-alpha*(x_temp))*((((beta-1.0)*(x_temp)*(x_temp))-(beta*(x_temp))+1.0)));
    break;

    case 3: // SW-A2
      alpha        = 2.25;
      beta         = 4.0;
      t1           = 0.00217;
      kappa        = 1.08;
      p_max        = 20.0e6;
      c0           = 1624.943090;
      sc           = 10e4/c0;
      x_start_temp = -x_start*sc;
      x_rescaled   = (1.0/x_start)*xw - 1.0;
      x_temp       = (-xw*sc-x_start_temp);
      pressure = p0 + p_max*(kappa*(1.0-exp(-(x_temp)/t1))*exp(-alpha*(x_temp))*((((beta-1.0)*(x_temp)*(x_temp))-(beta*(x_temp))+1.0)));
    break;

  }
  
  
  if (pressure<p0) {pressure = p0;}
  
  u_vel = (pressure - p0)/(MATERRINI[mat]*c0);

  rho=0.0;
  crnt->el->S[eqtypi[2]] = 0.0;
  for(iv=0; iv<NMATERIALS; ++iv){ 
    if (mat==iv) {avf=1.0 - (double)(NMATERIALS-1)*AMIN; }
    if (mat!=iv) {avf=                             AMIN; }
    crnt->el->S[eqtypi[0]+iv] = avf * MATERRINI[iv];        // ar
    rho                      += avf * MATERRINI[iv];
    
    are = avf * (pressure + MATERGAMA[iv]*MATERPINF[iv])/(MATERGAMA[iv]-1.0);    

    crnt->el->S[eqtypi[2]]   += are;  // rE
    crnt->el->S[eqtypi[4]+iv] = are;  // are

  }
 
  crnt->el->S[eqtypi[1]+0]  = u_vel * rho; // ru
  crnt->el->S[eqtypi[1]+1]  = v_init[mat] * rho; // rv
  crnt->el->S[eqtypi[1]+2]  = w_init[mat] * rho; // rw

  crnt->el->S[eqtypi[2]] += 0.5*rho*(pow(u_vel,2.0) + 
                                     pow(v_init[mat],2.0) +
                                     pow(w_init[mat],2.0));  // rE
  
  for(iv=0; iv<(NMATERIALS-1); ++iv){
    if (mat==iv) {avf=1.0 - (double)(NMATERIALS-1)*AMIN; }
    if (mat!=iv) {avf=                             AMIN; }

    crnt->el->S[eqtypi[3]+iv] = avf;  // avf
  }

  return;

}

void m5_init_S(struct RUN * run, struct BRANCH * crnt,int mat) {

  int iv,i,j;
  double avf_s,avf,are_elast,rho,are,Gmag;
  double Gmat[9];
  double Amat[9];

  if (MATERMUSH[mat]!=0.0) {
    avf_s = 1.0 - (double)(NMATERIALS-1)*AMIN;
  } else {
    avf_s = AMIN;
  }
  
  compute_Amat(Amat,avf_s);
  
  for (i=0;i<3;++i){
    for (j=0;j<3;++j){
      Gmat[i*3+j] = Amat[i*3+0]*Amat[j*3+0] + 
                    Amat[i*3+1]*Amat[j*3+1] + 
                    Amat[i*3+2]*Amat[j*3+2]; 
    } 
  }

  Gmag = Gmat[0]*(Gmat[4]*Gmat[8] - Gmat[5]*Gmat[7]) - Gmat[1]*(Gmat[3]*Gmat[8] - Gmat[5]*Gmat[6]) + Gmat[2]*(Gmat[3]*Gmat[7] - Gmat[4]*Gmat[6]);
  
  rho=0.0;
  crnt->el->S[eqtypi[2]] = 0.0;
  for(iv=0; iv<NMATERIALS; ++iv){ 
    if (mat==iv) {avf=1.0 - (double)(NMATERIALS-1)*AMIN; }
    if (mat!=iv) {avf=                             AMIN; }
    crnt->el->S[eqtypi[0]+iv] = avf * MATERRINI[iv];        // ar
    rho                      += avf * MATERRINI[iv];
    
    are_elast=m5_elastic_energy_solid_V2(Gmat,Gmag,MATERMUSH[iv]); 
    
    are = avf * (MATERPINI[mat] + MATERGAMA[iv]*MATERPINF[iv])/(MATERGAMA[iv]-1.0) + are_elast;    

    crnt->el->S[eqtypi[2]]   += are;  // rE
    crnt->el->S[eqtypi[4]+iv] = are;  // are

  }
 
  crnt->el->S[eqtypi[1]+0]  = u_init[mat] * rho; // ru
  crnt->el->S[eqtypi[1]+1]  = v_init[mat] * rho; // rv
  crnt->el->S[eqtypi[1]+2]  = w_init[mat] * rho; // rw

  crnt->el->S[eqtypi[2]] += 0.5*rho*(pow(u_init[mat],2.0) + 
                                     pow(v_init[mat],2.0) +
                                     pow(w_init[mat],2.0));  // rE
  
  for(iv=0; iv<(NMATERIALS-1); ++iv){
    if (mat==iv) {avf=1.0 - (double)(NMATERIALS-1)*AMIN; }
    if (mat!=iv) {avf=                             AMIN; }

    crnt->el->S[eqtypi[3]+iv] = avf;  // avf
  }

  crnt->el->S[eqtypi[5]+0] = Amat[0];  // A1
  crnt->el->S[eqtypi[5]+4] = Amat[4];  // B2
  crnt->el->S[eqtypi[5]+8] = Amat[8];  // C3
  
}

void m5_init_S_SW(struct RUN * run, struct BRANCH * crnt,int mat,double x_start, double xw) {

  int iv,i,j;
  double avf_s,avf,are_elast,rho,are,e_elast,Gmag;
  double p0,ps,alpha,omega,c0,pressure,u_vel;
  double beta,t1,kappa,p_max,sc,x_temp,x_start_temp,x_rescaled;
  double PI = 3.14159265358979323846;
  double Gmat[9];
  double Amat[9];

  if ((mat==0) && (MATERMUSH[0]!=0.0)) {
    avf_s = 1.0 - (double)(NMATERIALS-1)*AMIN;
  } else {
    avf_s = AMIN;
  }
  
  compute_Amat(Amat,avf_s);
  
  for (i=0;i<3;i++){
    for (j=0;j<3;j++){
      Gmat[i*3+j] = Amat[i*3+0]*Amat[j*3+0] + 
                    Amat[i*3+1]*Amat[j*3+1] + 
                    Amat[i*3+2]*Amat[j*3+2];
    }
  }

  Gmag = Gmat[0]*(Gmat[4]*Gmat[8] - Gmat[5]*Gmat[7]) - Gmat[1]*(Gmat[3]*Gmat[8] - Gmat[5]*Gmat[6]) + Gmat[2]*(Gmat[3]*Gmat[7] - Gmat[4]*Gmat[6]);
  
  switch(shock_prof){
    case 0: // Jc
      p0    = 101000.0;
      ps    = 35.0e6;
      alpha = 1.48e6;
      omega = 1.21e6;
      c0    = 1624.943090;
      pressure = p0 + 2.0*ps*exp(-alpha*((-xw/c0)+x_start/c0))*cos(omega*((-xw/c0)+x_start/c0) + PI/3.0);
    break;

    case 1: // Wang
      p0    = 101000.0;
      ps    = 39.3e6;
      alpha = 8.0e5;
      omega = 6.4e5;
      c0    = 1624.943090;
      pressure = p0 + 2.0*ps*exp(-alpha*((-xw/c0)+x_start/c0))*cos(omega*((-xw/c0)+x_start/c0) + PI/3.0);
    break;

    case 2:  // SW-A1
      alpha        = 6.0;
      beta         = 1.0;
      t1           = 0.00217;
      kappa        = 1.09;
      p_max        = 20.0e6;
      c0           = 1624.943090;
      sc           = 10e4/c0;
      x_start_temp = -x_start*sc;
      x_rescaled   = (1.0/x_start)*xw - 1.0;
      x_temp       = (-xw*sc-x_start_temp);
      pressure = p0 + p_max*(kappa*(1.0-exp(-(x_temp)/t1))*exp(-alpha*(x_temp))*((((beta-1.0)*(x_temp)*(x_temp))-(beta*(x_temp))+1.0)));
    break;

    case 3: // SW-A2
      alpha        = 2.25;
      beta         = 4.0;
      t1           = 0.00217;
      kappa        = 1.08;
      p_max        = 20.0e6;
      c0           = 1624.943090;
      sc           = 10e4/c0;
      x_start_temp = -x_start*sc;
      x_rescaled   = (1.0/x_start)*xw - 1.0;
      x_temp       = (-xw*sc-x_start_temp);
      pressure = p0 + p_max*(kappa*(1.0-exp(-(x_temp)/t1))*exp(-alpha*(x_temp))*((((beta-1.0)*(x_temp)*(x_temp))-(beta*(x_temp))+1.0)));
    break;

  }
  
  if (pressure<p0) {pressure = p0;}
  
  u_vel = (pressure - p0)/(MATERRINI[mat]*c0);

  rho=0.0;
  crnt->el->S[eqtypi[2]] = 0.0;
  for(iv=0; iv<NMATERIALS; ++iv){ 
    if (mat==iv) {avf=1.0 - (double)(NMATERIALS-1)*AMIN; }
    if (mat!=iv) {avf=                             AMIN; }
    crnt->el->S[eqtypi[0]+iv] = avf * MATERRINI[iv];        // ar
    rho                      += avf * MATERRINI[iv];

    are_elast = 0.0;
    are_elast = m5_elastic_energy_solid_V2(Gmat,Gmag,MATERMUSH[iv]); 
    
    are = avf * (pressure + MATERGAMA[iv]*MATERPINF[iv])/(MATERGAMA[iv]-1.0) + are_elast;    

    crnt->el->S[eqtypi[2]]   += are;  // rE
    crnt->el->S[eqtypi[4]+iv] = are;  // are

  }
 
  crnt->el->S[eqtypi[1]+0]  = u_vel * rho; // ru
  crnt->el->S[eqtypi[1]+1]  = v_init[mat] * rho; // rv
  crnt->el->S[eqtypi[1]+2]  = w_init[mat] * rho; // rw

  crnt->el->S[eqtypi[2]] += 0.5*rho*(pow(u_vel,2.0) + 
                                     pow(v_init[mat],2.0) +
                                     pow(w_init[mat],2.0));  // rE
  
  for(iv=0; iv<(NMATERIALS-1); ++iv){
    if (mat==iv) {avf=1.0 - (double)(NMATERIALS-1)*AMIN; }
    if (mat!=iv) {avf=                             AMIN; }

    crnt->el->S[eqtypi[3]+iv] = avf;  // avf
  }

  crnt->el->S[eqtypi[5]+0] = Amat[0];  // A1
  crnt->el->S[eqtypi[5]+4] = Amat[4];  // B2
  crnt->el->S[eqtypi[5]+8] = Amat[8];  // C3
  
}