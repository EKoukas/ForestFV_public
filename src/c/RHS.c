#include "strdata.h"

void RHS(struct RUN * run){
  
  struct BRANCH * crnt;

  switch(MODEL){

		case 0:
      run->Trhs=timecpu(run->Trhs,0);
      crnt=run->topo->locl;
      while (crnt!=NULL){
        m0_residual(run,crnt);
        crnt=crnt->lnxt;
      }
      run->Trhs=timecpu(run->Trhs,1);
    break;

		case 1:
      run->Trhs=timecpu(run->Trhs,0);
      crnt=run->topo->locl;
      while (crnt!=NULL){
        m1_residual(run,crnt);
        crnt=crnt->lnxt;
      }
      run->Trhs=timecpu(run->Trhs,1);
    break;

    case 2:
      run->Trhs=timecpu(run->Trhs,0);
      crnt=run->topo->locl;
      while (crnt!=NULL){
        m2_residual(run,crnt);
        crnt=crnt->lnxt;
      }
      run->Trhs=timecpu(run->Trhs,1);
    break;

    case 3:
      run->Trhs=timecpu(run->Trhs,0);
      crnt=run->topo->locl;
      while (crnt!=NULL){
        m3_residual(run,crnt);
        crnt=crnt->lnxt;
      }
      run->Trhs=timecpu(run->Trhs,1);
    break;

    case 4:
      run->Trhs=timecpu(run->Trhs,0);
      crnt=run->topo->locl;
      while (crnt!=NULL){
        m4_residual(run,crnt);
        crnt=crnt->lnxt;
      }
      run->Trhs=timecpu(run->Trhs,1);
    break;

    case 5:
      run->Trhs=timecpu(run->Trhs,0);
      crnt=run->topo->locl;
      while (crnt!=NULL){
        m5_residual(run,crnt);
        crnt=crnt->lnxt;
      }
      run->Trhs=timecpu(run->Trhs,1);
    break;

    default:
		break;
	} 
  
}