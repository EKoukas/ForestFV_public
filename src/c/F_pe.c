#include "strdata.h"

double F_pe(struct RUN * run,double pe,int i_paper) {
  
  int iv;
  double F_pe,gamma,pinf,p_quest;

  if (i_paper==0) { // Simple 
    F_pe=0.0;
    for (iv=0;iv<eqtypn[0];++iv){
      gamma = MATERGAMA[iv];
      pinf  = MATERPINF[iv];
      p_quest = run->p0_m[iv];
      F_pe += run->ar_R[iv]*run->t0_m[iv]*( (run->p0_m[iv] + gamma*pinf + (gamma-1.0)*p_quest ) / (pe + gamma*pinf + (gamma-1.0)*p_quest) );
    }
     F_pe += -1.0;
  }
  else if (i_paper==1) { // Favrie
    F_pe=0.0;
    for (iv=0;iv<eqtypn[0];++iv){
      gamma = MATERGAMA[iv];
      pinf  = MATERPINF[iv];
      F_pe += run->ar_R[iv]*run->t0_m[iv]*( (run->p0_m[iv] + gamma*pinf + (gamma-1.0)*pe ) / (pe + gamma*pinf + (gamma-1.0)*pe) );
    }
     F_pe += -1.0;
  } 
  else if (i_paper==2) { // Ntanou
    F_pe=0.0;
    for (iv=0;iv<eqtypn[0];++iv){
      gamma = MATERGAMA[iv];
      pinf  = MATERPINF[iv];
      F_pe += (run->avf_R[iv]/gamma)*((run->p0_m[iv] - pe)/(pe + pinf));
    }
  }
  else if (i_paper==3) { // Metayer
    F_pe=0.0;
    for (iv=0;iv<eqtypn[0];++iv){
      gamma = MATERGAMA[iv];
      pinf  = MATERPINF[iv];
      F_pe += ((run->avf_R[iv]/gamma)*((run->p0_m[iv] + pinf)/(pe + pinf))) - (run->avf_R[iv]/gamma);
    }
  }
  return F_pe;
}


double DF_pe(struct RUN * run,double pe,int i_paper) {

  int iv,i_anal;
  double A1,A2,DF_pe,pe_small;
  double gamma,pinf,p_quest,c1,c2,f,g;

  i_anal = 1;
  pe_small = pow(10.0,4.0);

  if (i_anal==0) {
      DF_pe = F_pe(run,pe+pe_small,i_paper) - F_pe(run,pe-pe_small,i_paper) / (2.0*pe_small);
  } else if (i_anal==1) {
    if (i_paper==0) {

      DF_pe=0.0;
      for (iv=0;iv<eqtypn[0];++iv){
        gamma = MATERGAMA[iv];
        pinf  = MATERPINF[iv];
        p_quest = run->p0_m[iv];

        c1 = run->ar_R[iv]*run->t0_m[iv];
        f = (run->p0_m[iv] + gamma*pinf + (gamma-1.0)*p_quest );
        g = (pe            + gamma*pinf + (gamma-1.0)*p_quest );

        DF_pe += (0.0 - c1*f) / pow(g,2.0);
        
      }

    } 
    else if (i_paper==1) {

      DF_pe=0.0;
      for (iv=0;iv<eqtypn[0];++iv){
        gamma = MATERGAMA[iv];
        pinf  = MATERPINF[iv];

        c1 = run->ar_R[iv]*run->t0_m[iv];
        f = (run->p0_m[iv] + gamma*pinf + (gamma-1.0)*pe);
        g = (pe            + gamma*pinf + (gamma-1.0)*pe);

        DF_pe += (c1*(gamma-1.0)*g - c1*f*gamma) / pow(g,2.0);
      }

    } 
    else if (i_paper==2) {
      
      DF_pe=0.0;
      for (iv=0;iv<eqtypn[0];++iv){
        gamma = MATERGAMA[iv];
        pinf  = MATERPINF[iv];

        c2 = run->avf_R[iv]/gamma;
        f = (run->p0_m[iv] - pe);
        g = (pe + pinf);

        DF_pe += (-c2*g - c2*f) / pow(g,2.0);
      }

    }
    else if (i_paper==3) {
        
      DF_pe=0.0;
      for (iv=0;iv<eqtypn[0];++iv){
        gamma = MATERGAMA[iv];
        pinf  = MATERPINF[iv];

        c2 = run->avf_R[iv]/gamma;
        f = (run->p0_m[iv] + pinf);
        g = (pe + pinf);

        DF_pe += ( 0.0 - c2*f) / pow(g,2.0);
    
      }

    }
  }
  return DF_pe;
}