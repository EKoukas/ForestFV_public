#include "strdata.h"

void m2_cons2primtv(double * S_cons,double * S_prim) {

	S_prim[0] = S_cons[0];
	S_prim[1] = S_cons[1]/S_cons[0]; 
	S_prim[2] = S_cons[2]/S_cons[0];
	S_prim[3] = S_cons[3]/S_cons[0];

	S_prim[4] = S_cons[4]; ///S_cons[0];

	//if (S_prim[4]<     (YMIN/1000.0) ) { S_prim[4] =     (YMIN/1000.0); }
	//if (S_prim[4]>(1.0-(YMIN/1000.0))) { S_prim[4] = 1.0-(YMIN/1000.0); }

}