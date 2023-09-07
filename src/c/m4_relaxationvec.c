#include "strdata.h"

void m4_relaxationvec(struct RUN * run,struct BRANCH * crnt) {

  int i,j,iv,v,iter,i_as;
  double A_1,A_2;
  double avfsum,avf_dif,avf_per,rho;
  double rE_total,u_vel,v_vel,w_vel;  
  double pe,p_relaxed;
  double e_hyd,num,den;

	double avf_R[10];
  double ar_R[10];
  double are_R[10];
	double t0_m[10];
	double t1_m[10];
	double p0_m[10];
	double r_m[10];
  double avf_init[10];
    
  double eps_sml = 1.0e-12;
  double eps_big = 1.0e-8;

  // =========================================================
  //                      Assign values
  // =========================================================
  rho=0.0;
  for(iv=0;iv<eqtypn[0];++iv){
		ar_R[iv]  = run->Vec_temp_R[iv+eqtypi[0]];
		are_R[iv] = run->Vec_temp_R[iv+eqtypi[4]];
		rho      += ar_R[iv];
  }

	iv = eqtypi[1];
  u_vel = run->Vec_temp_R[iv + 0]/rho;
  v_vel = run->Vec_temp_R[iv + 1]/rho;
  w_vel = run->Vec_temp_R[iv + 2]/rho;
  rE_total=run->Vec_temp_R[eqtypi[2]];

  //==================================================================
  //vf limit
  avfsum=0.0;  
  for(iv=0;iv<eqtypn[3];++iv){
    
    avf_R[iv] = run->Vec_temp_R[iv+eqtypi[3]];
    
    if (avf_R[iv]<=      AMIN/100.0)   { avf_R[iv]=     AMIN/100.0;  } 
    if (avf_R[iv]>=(1.0-(AMIN/100.0))) { avf_R[iv]=1.0-(AMIN/100.0); }   
    avfsum  += avf_R[iv]; 
  }  
    
  if (avfsum>=(1.0-(AMIN/100.0))) { 
    avf_dif = avfsum - (1.0-(AMIN/100.0));
    for(j=0;j<eqtypn[3];++j){
      avf_per       = avf_R[j]/avfsum;
      avf_R[j] = avf_R[j] - avf_per*avf_dif;
    }
    avfsum = 1.0-(AMIN/100.0);
  }
  avf_R[eqtypn[3]] = 1.0 - avfsum;
  //=====================================================================

  if (PRELAX==1) {

    for (iv=0;iv<eqtypn[0];++iv){
      t0_m[iv] = 1.0/(ar_R[iv]/avf_R[iv]);
      p0_m[iv] = (are_R[iv]*(MATERGAMA[iv]-1.0) - avf_R[iv]*MATERGAMA[iv]*MATERPINF[iv])/avf_R[iv];
    }
      
    A_1 = (avf_R[0]/MATERGAMA[0])*(p0_m[0]+MATERPINF[0]) / ( (avf_R[0]/MATERGAMA[0]) + (avf_R[1]/MATERGAMA[1]) );
    A_2 = (avf_R[1]/MATERGAMA[1])*(p0_m[1]+MATERPINF[1]) / ( (avf_R[0]/MATERGAMA[0]) + (avf_R[1]/MATERGAMA[1]) );
    
    pe  =             0.5*(A_1+A_2 - (MATERPINF[0]+MATERPINF[1])) + 
          sqrt(0.25*(pow( (A_2-A_1 - (MATERPINF[1]-MATERPINF[0])),2.0)) + A_1*A_2);

    if ((isnan(pe)==1)){
      printf(" \n");
      printf("m4_relaxationvec: 1: pe NaN %d %d %d | %f %f | %f %f \n",crnt->root,crnt->type,g_istep,A_1,A_2,p0_m[0],p0_m[1]);
      exit(0);
    }

    for (iv=0;iv<eqtypn[0];++iv){
      t1_m[iv] = t0_m[iv] * (p0_m[iv] + MATERGAMA[iv]*MATERPINF[iv] + (MATERGAMA[iv]-1.0)*pe) / 
                            (pe       + MATERGAMA[iv]*MATERPINF[iv] + (MATERGAMA[iv]-1.0)*pe);

      r_m[iv] = 1.0/t1_m[iv];
      avf_R[iv] = ar_R[iv]/r_m[iv];
    }       

    // vf limit========================================
    avfsum=0.0;  
    for(iv=0;iv<eqtypn[3];++iv){ 
      if (avf_R[iv]<=      AMIN/100.0)   { avf_R[iv]=     AMIN/100.0;  } 
      if (avf_R[iv]>=(1.0-(AMIN/100.0))) { avf_R[iv]=1.0-(AMIN/100.0); }   
      avfsum  += avf_R[iv]; 
    }  
      
    if (avfsum>=(1.0-(AMIN/100.0))) { 
      avf_dif = avfsum - (1.0-(AMIN/100.0));
      for(j=0;j<eqtypn[3];++j){
        avf_per  = avf_R[j]/avfsum;
        avf_R[j] = avf_R[j] - avf_per*avf_dif;
      }
      avfsum = 1.0-(AMIN/100.0);
    }
    avf_R[eqtypn[3]] = 1.0 - avfsum;

    //=================================================

    p_relaxed = pe;

    // Re-initialization
    e_hyd = rE_total/rho - 0.5*(pow((u_vel),2.0) + pow((v_vel),2.0) + pow((w_vel),2.0));

    num = 0.0; den = 0.0;
    for (iv=0;iv<eqtypn[0];++iv){
      num += avf_R[iv] * MATERGAMA[iv] * MATERPINF[iv] / (MATERGAMA[iv]-1.0);
      den += avf_R[iv]                                 / (MATERGAMA[iv]-1.0);
    }
          
    p_relaxed = (rho*e_hyd - num) / den;

    for (iv=0;iv<eqtypn[0];++iv){
      are_R[iv] = (avf_R[iv]*(p_relaxed + MATERGAMA[iv]*MATERPINF[iv]) / (MATERGAMA[iv]-1.0)); 
    }
 
  } 
  else if (PRELAX==2) {

		int i_paper = 1;
		int iter_NR,principal_phase;
		
		double pe_new,pe_old,pe_anal,pe_old_2,p_init;
		double max_avf,max_press,min_press;
		double cond_1,cond_2;

    pe_old=0.0;
    max_avf = - 1.0;
    max_press = 0.0;
    min_press = 1e20;
    for (iv=0;iv<eqtypn[0];++iv){
    
      t0_m[iv] = 1.0/(ar_R[iv]/avf_R[iv]);    
      p0_m[iv] = (are_R[iv]*(MATERGAMA[iv]-1.0) - avf_R[iv]*MATERGAMA[iv]*MATERPINF[iv])/avf_R[iv];
      
      if (p0_m[iv]<-(1.0-eps_big)*MATERPINF[iv]+eps_big){
        //printf("Limit pressure 1 | %d  %d | %e %e | %e %e %e | \n",iv,crnt->root,p0_m[iv],-(1.0-eps_big)*MATERPINF[iv]+eps_big,avf_R[0],avf_R[1],avf_R[2]);
        p0_m[iv]=-(1.0-eps_big)*MATERPINF[iv]+1.0;
        
      }
      
      if (max_press <p0_m[iv]) {
        max_press =p0_m[iv];
      }
      if (min_press >p0_m[iv]) {
        min_press =p0_m[iv];
      }

      pe_old_2 += avf_R[iv]*p0_m[iv];
      
      if (avf_R[iv]>max_avf) {
        max_avf = avf_R[iv];
        principal_phase  = iv;
      }

    }
    pe_old = 0.001;
    
    
    iter_NR=0;
    pe_new = pe_old; //pe_old;
    p_init = pe_old;
    
    cond_1 = fabs(pe_new-pe_old); ///fabs(pe_new);
    cond_2 = fabs(F_pe(run,pe_new,i_paper));

    //while ((fabs(F_pe(run,pe_new,i_paper))>=pow(10.0,-10.0)) && (iter_NR<=101)) {
    //while (((fabs(pe_new-pe_old)>=pow(10,-6.0)) && (iter_NR<=101)) || (iter_NR==0)) {
    while (((cond_2>=pow(10.0,-12.0)) || (cond_1>=pow(10.0,-4.0))) && (iter_NR<=101)) {  

      pe_new = pe_old - 0.8*(F_pe(run,pe_old,i_paper)/DF_pe(run,pe_old,i_paper));

      cond_1 = fabs(pe_new-pe_old); ///fabs(pe_new);
      cond_2 = fabs(F_pe(run,pe_new,i_paper));

      if ((iter_NR>=100)) { 

        if (((fabs(pe_old-pe_new)<pow(10.0,3.0))) || (isnan(pe_new)==1)) {
          printf(" NR FAILED 1: %d %d | %d %d  \n",i_paper,iter_NR,crnt->root,g_istep);
          printf(" p_init:  %e  \n",p_init);
          printf(" pe_old: %e  \n",pe_anal);
          printf(" press: %e %e \n",pe_old,pe_new);
          printf(" F_pe: %e %e \n",F_pe(run,pe_old,i_paper),DF_pe(run,pe_old,i_paper));
          printf(" \n"); 

          for (iv=0;iv<eqtypn[0];++iv){
            printf("ar_R %d | %e \n",iv,ar_R[iv]);
            printf("rho  %d | %e \n",iv,1.0/t0_m[iv]);
            printf("p0_m %d | %e \n",iv,p0_m[iv]);
            printf("avf_R %d | %e \n",iv,avf_R[iv]); 
            printf("are_R %d | %e \n",iv,are_R[iv] ); 
            
            printf(" \n");       
          }
          exit(0);
          
        }
          
      }
            
      pe_old = pe_new;
      iter_NR++;

    }

    if (i_paper==2) {
      if ((pe_new>1.0001*max_press) || (pe_new<0.99999*min_press) ) {
        printf("out of bounds %e | %e %e | %e %e %e \n",pe_new,max_press,min_press,p0_m[0],p0_m[1],p0_m[2]);
        exit(0);
      }
    }

    if ((isnan(pe_new)==1)){
      printf(" \n");
      printf("Relaxation 2: pe NaN %d %d %d | %e %e %e %e | %e %e | \n",
      crnt->root,g_istep,iter_NR,pe_new,pe_old,cond_1,cond_2,F_pe(run,pe_old,i_paper),DF_pe(run,pe_old,i_paper));
      exit(0);
    }

    if (pe_new<(-(1.0-eps_big)*MATERPINF[principal_phase]+eps_big)){
      printf("Limit pressure 2 %d %d \n",crnt->root,g_istep);
      pe_new=-(1.0-eps_big)*MATERPINF[principal_phase]+1.0;
    }
 
    pe=pe_new;
    for (iv=0;iv<eqtypn[0];++iv){

      if (i_paper==0) { // Simple 
       	t1_m[iv] = t0_m[iv] * (p0_m[iv] + MATERGAMA[iv]*MATERPINF[iv] + (MATERGAMA[iv]-1.0)*p0_m[iv]) / 
                              (pe       + MATERGAMA[iv]*MATERPINF[iv] + (MATERGAMA[iv]-1.0)*p0_m[iv]);
      } else if (i_paper==1) { // Favrie
       	t1_m[iv] = t0_m[iv] * (p0_m[iv] + MATERGAMA[iv]*MATERPINF[iv] + (MATERGAMA[iv]-1.0)*pe) / 
                              (pe       + MATERGAMA[iv]*MATERPINF[iv] + (MATERGAMA[iv]-1.0)*pe);
      } else if (i_paper==2) {
       	t1_m[iv] = t0_m[iv] * (p0_m[iv] + MATERGAMA[iv]*MATERPINF[iv] + (MATERGAMA[iv]-1.0)*pe) / 
                        			(pe 	    + MATERGAMA[iv]*MATERPINF[iv] + (MATERGAMA[iv]-1.0)*pe);
      }
      
     	r_m[iv] = 1.0/t1_m[iv];
      avf_R[iv] = ar_R[iv]/r_m[iv];
        
    } 

    // vf limit ===============================================================
    avfsum=0.0;  
    for(iv=0;iv<eqtypn[3];++iv){        
      if (avf_R[iv]<=      AMIN/100.0)   { avf_R[iv]=     AMIN/100.0;  } 
      if (avf_R[iv]>=(1.0-(AMIN/100.0))) { avf_R[iv]=1.0-(AMIN/100.0); }   
      avfsum  += avf_R[iv]; 
    }  
      
    if (avfsum>=(1.0-(AMIN/100.0))) { 
      avf_dif = avfsum - (1.0-(AMIN/100.0));
      for(j=0;j<eqtypn[3];++j){
        avf_per       = avf_R[j]/avfsum;
        avf_R[j] = avf_R[j] - avf_per*avf_dif;
      }
      avfsum  = 1.0-(AMIN/100.0);
    }
    avf_R[eqtypn[3]] = 1.0 - avfsum;   
    // vf limit ===============================================================
    
    p_relaxed = pe;
    
    e_hyd = rE_total/rho - 0.5*(pow((u_vel),2.0) + pow((v_vel),2.0) + pow((w_vel),2.0));

    num = 0.0;
    den = 0.0;

    for (iv=0;iv<eqtypn[0];++iv){
      num += avf_R[iv] * MATERGAMA[iv] * MATERPINF[iv] / (MATERGAMA[iv]-1.0);
      den += avf_R[iv]                                / (MATERGAMA[iv]-1.0);
    }
          
    p_relaxed = (rho*e_hyd - num) / den;

    if ((isnan(p_relaxed)==1)){
      printf(" \n");
      printf("Relaxation NR: p_relaxed NaN %d %d %d | %e %e %e | %e | %e %e  \n",crnt->root,crnt->type,g_istep,e_hyd,num,den,pe,avf_R[0],avf_R[1]);
      exit(0);
    }

    for (iv=0;iv<eqtypn[0];++iv){
      are_R[iv] = (avf_R[iv]*(p_relaxed + MATERGAMA[iv]*MATERPINF[iv]) / (MATERGAMA[iv]-1.0)); 
    }

  }  
  else if (PRELAX==3) {

		int irelax;
		double F,dF,Gmag;
		double rhoKs[10];
		double pkinit[10];

    for(iv=0;iv<NMATERIALS;++iv){
      if (avf_R[iv]<0.0||ar_R[iv]<0.0) {
        avf_R[iv]  = eps_sml;
        ar_R[iv]   = eps_sml;
        are_R[iv]  = eps_sml;
      }
      if (avf_R[iv] >= 1.0) {avf_R[iv] = 1.0-eps_sml;}
    }

    // Pressures relaxation procedure ===================================
    // Is the pressure relaxation procedure necessary?
    irelax=1;
    for(iv=0;iv<NMATERIALS;++iv){
      if (avf_R[iv]>1.0-eps_sml){irelax=0;}
    }
      
    if (irelax==1){
  
      // initialisation
      p_relaxed = 0.0;
      for(iv=0;iv<NMATERIALS;++iv){
        if(avf_R[iv]>eps_sml){
          
          pkinit[iv] = (are_R[iv])/avf_R[iv]*(MATERGAMA[iv]-1.0)-MATERGAMA[iv]*MATERPINF[iv];
          
          if (isnan(pkinit[iv])==1){

            Gmag = G_det(run->Gmat);

            printf(" \n");
            printf("Relaxation 3: pkinit NaN | %d %d %d | %e \n",crnt->root,g_istep,iv,Gmag);
            printf("ar_R: ");
            for(v=0;v<NMATERIALS;++v){
              printf(" %e ",ar_R[v]);
            }
            printf(" \n");

            printf("avf_R: ");
            for(v=0;v<NMATERIALS;++v){
              printf(" %e ",avf_R[v]);
            }
            printf(" \n");

            printf("are_R: ");
            for(v=0;v<NMATERIALS;++v){
              printf(" %e ",are_R[v]);
            }
            printf(" \n");

            printf("pkinit: ");
            for(v=0;v<NMATERIALS;++v){
              printf(" %e ",pkinit[v]);
            }
            printf(" \n");

          }

          if (pkinit[iv]<=-(1.0-eps_big)*MATERPINF[iv]+eps_big) {
            pkinit[iv]=-(1.0-eps_big)*MATERPINF[iv]+eps_big;
          }
        } else {
          pkinit[iv]=0.0;
        }
        p_relaxed+=avf_R[iv]*pkinit[iv];
      }

      // NR for p_relaxed
      iter = 0;
      F    = 1.0e-9;
      dF   = 1.0e9;
      for (iv=0;iv<NMATERIALS;++iv) {
        rhoKs[iv] = 0.0;
      }

      while((fabs(F)>1.0e-12)&&(iter<500)){
          
        p_relaxed = p_relaxed-F/dF;
        iter = iter + 1;
   
        for(iv=0;iv<NMATERIALS;++iv){
          if (p_relaxed<-(1.0-eps_big)*MATERPINF[iv]+eps_big){
            //printf("Limit pressure 2 \n");
            p_relaxed=-(1.0-eps_big)*MATERPINF[iv]+1.0;
          }
        }
        
        F  = -1.0;
        dF =  0.0;
        for(iv=0;iv<NMATERIALS;++iv){
          if (avf_R[iv]>eps_sml){
            rhoKs[iv]=ar_R[iv]/max(avf_R[iv],eps_sml)*pow((p_relaxed+MATERPINF[iv])/(pkinit[iv]+MATERPINF[iv]),1.0/MATERGAMA[iv]);
            F  = F+ar_R[iv]/rhoKs[iv];
            dF = dF-ar_R[iv]/(MATERGAMA[iv]*rhoKs[iv]*(p_relaxed+MATERPINF[iv]));
          }
        }
        
        if ((isnan(p_relaxed)==1) || (isnan(F)==1) || (isnan(dF)==1)){
          printf(" \n");
          printf("Relaxation 3: p_relaxed NaN | %d %d %d  \n",crnt->root,g_istep,iter);
          printf("p_relaxed F dF | %e %e %e | \n",p_relaxed,F,dF);
          printf("avf: ");
          for(iv=0;iv<NMATERIALS;++iv){
            printf(" %e ",avf_R[iv]);
          }
          printf(" \n");

          printf("ar_R: ");
          for(iv=0;iv<NMATERIALS;++iv){
            printf(" %e ",ar_R[iv]);
          }
          printf(" \n");

          printf("rhoKs: ");
          for(iv=0;iv<NMATERIALS;++iv){
            printf(" %e ",rhoKs[iv]);
          }
          printf(" \n");

          printf("pkinit: ");
          for(iv=0;iv<NMATERIALS;++iv){
            printf(" %e ",pkinit[iv]);
          }
          printf(" \n");

          exit(0);
        }

      }

      // adjust volume fraction
      avfsum = 0.0;
      for(iv=0;iv<NMATERIALS;++iv){
        if (avf_R[iv]>eps_sml){
          avf_R[iv]=ar_R[iv]/rhoKs[iv];
        }
      }

      // vf limit ===============================================================
      avfsum=0.0;  
      for(iv=0;iv<eqtypn[3];++iv){        
        if (avf_R[iv]<=      AMIN/100.0)   { avf_R[iv]=     AMIN/100.0;  } 
        if (avf_R[iv]>=(1.0-(AMIN/100.0))) { avf_R[iv]=1.0-(AMIN/100.0); }   
        avfsum  += avf_R[iv]; 
      }  
        
      if (avfsum>=(1.0-(AMIN/100.0))) { 
        avf_dif = avfsum - (1.0-(AMIN/100.0));
        for(j=0;j<eqtypn[3];++j){
          avf_per       = avf_R[j]/avfsum;
          avf_R[j] = avf_R[j] - avf_per*avf_dif;
        }
        avfsum  = 1.0-(AMIN/100.0);
      }
      avf_R[eqtypn[3]] = 1.0 - avfsum;   
      // vf limit ===============================================================

    }
    
    e_hyd = rE_total/rho - 0.5*(pow((u_vel),2.0)  + pow((v_vel),2.0)  + pow((w_vel),2.0));
    
    num = 0.0;
    den = 0.0;
    for (iv=0;iv<NMATERIALS;++iv){
      num += avf_R[iv] * MATERGAMA[iv] * MATERPINF[iv] / (MATERGAMA[iv]-1.0);
      den += avf_R[iv]                                 / (MATERGAMA[iv]-1.0);
    }

    // Reseting of internal energies
    p_relaxed = (rho*e_hyd - num) / den;

    for (iv=0;iv<NMATERIALS;++iv){
      are_R[iv] = (avf_R[iv]*(p_relaxed + MATERGAMA[iv]*MATERPINF[iv]) / (MATERGAMA[iv]-1.0)); 
    } 

  }

  // =========================================================
  //                Update run->Vec_temp_Rtor
  // =========================================================

  for(iv=0;iv<eqtypn[3];++iv){ 
    
    if ((isnan(avf_R[iv])==1)){
      printf(" \n");
      printf("Rlx_v2: FAIL at vf %d %d | %e | %e %e | \n",g_istep,crnt->root,avf_R[iv],pe,p_relaxed);
      exit(0);
    }
    
    run->Vec_temp_R[iv+eqtypi[3]] = avf_R[iv];
  }
  
  for(iv=0;iv<eqtypn[4];++iv){
    if ((isnan(are_R[iv])==1)){
      printf(" \n");
      printf("Rlx_v2: FAIL at internal energies %d %d %f %f %f \n",g_istep,crnt->root,run->Vec_temp_R[iv+eqtypi[4]],p_relaxed,avf_R[iv]);
      exit(0);
    }
    run->Vec_temp_R[iv+eqtypi[4]] = are_R[iv];
  }


}