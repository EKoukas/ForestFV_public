#include "strdata.h"

void RHS_nb(struct RUN * run,int inner_outer){
  
  int iv;
  struct BRANCH * crnt;

  switch(MODEL){

		case 0:
      run->Trhs=timecpu(run->Trhs,0);
      crnt=run->topo->locl;
      while (crnt!=NULL){
        if (crnt->el->inner==inner_outer){
          m0_residual(run,crnt);
        }
        crnt=crnt->lnxt;
      }
      run->Trhs=timecpu(run->Trhs,1);
    break;

		case 1:
      run->Trhs=timecpu(run->Trhs,0);
      crnt=run->topo->locl;
      while (crnt!=NULL){
        if (crnt->el->inner==inner_outer){
          m1_residual(run,crnt);
        }
        crnt=crnt->lnxt;
      }
      run->Trhs=timecpu(run->Trhs,1);
    break;

    case 2:
      run->Trhs=timecpu(run->Trhs,0);
      crnt=run->topo->locl;
      while (crnt!=NULL){
        if (crnt->el->inner==inner_outer){
          m2_residual(run,crnt);
        }
        crnt=crnt->lnxt;
      }
      run->Trhs=timecpu(run->Trhs,1);
    break;

    case 3:
      run->Trhs=timecpu(run->Trhs,0);
      crnt=run->topo->locl;
      while (crnt!=NULL){
        if (crnt->el->inner==inner_outer){
          m3_residual(run,crnt);
        }
        crnt=crnt->lnxt;
      }
      run->Trhs=timecpu(run->Trhs,1);
    break;

    case 4:
      run->Trhs=timecpu(run->Trhs,0);
      crnt=run->topo->locl;
      while (crnt!=NULL){
        if (crnt->el->inner==inner_outer){
          m4_residual(run,crnt);
        }
        crnt=crnt->lnxt;
      }
      run->Trhs=timecpu(run->Trhs,1);
    break;

    case 5:
      run->Trhs=timecpu(run->Trhs,0);
      crnt=run->topo->locl;
      while (crnt!=NULL){
        if (crnt->el->inner==inner_outer){
          m5_residual(run,crnt);
        }
        crnt=crnt->lnxt;
      }
      run->Trhs=timecpu(run->Trhs,1);
    break;

    default:
		break;
	} 
  
}