#include "strdata.h"

void  createlfc (BRANCH * crnt) {
  
	if (crnt->type==0){    
		crnt->cl->fc[0].type = 0;
		crnt->cl->fc[0].nnds = 4;
		crnt->cl->fc[0].bc   = 0;

		crnt->cl->fc[1].type = 0;
		crnt->cl->fc[1].nnds = 4;
		crnt->cl->fc[1].bc   = 0;

		crnt->cl->fc[2].type = 0;
		crnt->cl->fc[2].nnds = 4;
		crnt->cl->fc[2].bc   = 0;

		crnt->cl->fc[3].type = 0;
		crnt->cl->fc[3].nnds = 4;
		crnt->cl->fc[3].bc   = 0;

		crnt->cl->fc[4].type = 0;
		crnt->cl->fc[4].nnds = 4;
		crnt->cl->fc[4].bc   = 0;

		crnt->cl->fc[5].type = 0;
		crnt->cl->fc[5].nnds = 4;
		crnt->cl->fc[5].bc   = 0;

	} else if (crnt->type==2){
		
		crnt->cl->fc[0].type = 0;
		crnt->cl->fc[0].nnds = 4;
		crnt->cl->fc[0].bc   = 0;

		crnt->cl->fc[1].type = 0;
		crnt->cl->fc[1].nnds = 4;
		crnt->cl->fc[1].bc   = 0;

		crnt->cl->fc[2].type = 0;
		crnt->cl->fc[2].nnds = 4;
		crnt->cl->fc[2].bc   = 0;

		crnt->cl->fc[3].type = 1;
		crnt->cl->fc[3].nnds = 3;
		crnt->cl->fc[3].bc   = 0;

		crnt->cl->fc[4].type = 1;
		crnt->cl->fc[4].nnds = 3;
		crnt->cl->fc[4].bc   = 0;

	} else if (crnt->type==1){
		crnt->cl->fc[0].type = 1;
		crnt->cl->fc[0].nnds = 3;
		crnt->cl->fc[0].bc   = 0;


		crnt->cl->fc[1].type = 1;
		crnt->cl->fc[1].nnds = 3;
		crnt->cl->fc[1].bc   = 0;


		crnt->cl->fc[2].type = 1;
		crnt->cl->fc[2].nnds = 3;
		crnt->cl->fc[2].bc   = 0;


		crnt->cl->fc[3].type = 1;
		crnt->cl->fc[3].nnds = 3;
		crnt->cl->fc[3].bc   = 0;
	}

}