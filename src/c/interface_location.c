//-----------------------------------------------------------
//			Find Interface location and mark adjacent cells
//-----------------------------------------------------------

#include "strdata.h"

void interface_location(struct RUN * run) {

	int f,v,flag_inter,flag_moved,SI_lvl,ing,oppo_SI;
	struct BRANCH * crnt;
  struct BRANCH * oppo;
	
  // Initial, 1
	crnt=run->topo->locl;
	while (crnt!=NULL){   // Loop 
    			
    crnt->el->SI[0] = 0;
    crnt->el->SI[1] = 0;
    flag_inter=0;
    f=0;

		while ((f<(crnt->nlfc)) && (flag_inter==0)) {
      if(crnt->cl->fc[f].bc==0){
      
        for (ing=0;ing<crnt->nsfc[f];ing++) {
          
          switch(MODEL) {
            case 1: m1_facevalues (run,crnt,ing,f,0,0); m1_facevalues (run,crnt,ing,f,0,1); break;
            case 2: m2_facevalues (run,crnt,ing,f,0,0); m2_facevalues (run,crnt,ing,f,0,1); break;
            case 3: m3_facevalues (run,crnt,ing,f,0,0); m3_facevalues (run,crnt,ing,f,0,1); break;
            case 4: m4_facevalues (run,crnt,ing,f,0,0); m4_facevalues (run,crnt,ing,f,0,1); break;
            case 5: m5_facevalues (run,crnt,ing,f,0,0); m5_facevalues (run,crnt,ing,f,0,1); break;
          }

          for(v=0;v<eqtypn[3];++v) {
            if (run->sol_L->phi[v]*run->sol_R->phi[v]<0.0) {
              if      (v==0) { crnt->el->SI[0] = 1; }
              else if (v!=0) { crnt->el->SI[1] = 1; }
              flag_inter=1;
            } 
          }
        } 
      }
      f++;
		}

		crnt=crnt->lnxt;
  }
  
  communicate_S(run,4);
 
  for (SI_lvl=2;SI_lvl<ADAPT_INTER_MAX;SI_lvl++) {

    crnt=run->topo->locl;
    while (crnt!=NULL){   // Loop 
      
      flag_inter=0;
      f=0;

      while ((f<(crnt->nlfc)) && (flag_inter==0)) {
        if(crnt->cl->fc[f].bc==0) {

          for (ing=0;ing<crnt->nsfc[f];ing++){
            oppo=crnt->neigtr[f][ing];
            for (v=0;v<2;++v){
              oppo_SI=oppo->el->SI[v];
              if ((oppo_SI==(SI_lvl-1)) && (crnt->el->SI[v]==0)) {
                crnt->el->SI[v] = SI_lvl;
                flag_inter=1;
              }  
            }     
          } 
        }
        f++;
      }
      
      crnt=crnt->lnxt;
    }
    
    communicate_S(run,4);

  }
  


}
