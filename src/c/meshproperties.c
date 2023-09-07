#include "strdata.h"

void meshproperties(struct RUN *run) {

  struct BRANCH * crnt;
  
  crnt = run->topo->locl;
  while (crnt!=NULL){   // Loop 
    
    normalvectorcalc(run,crnt);
    tagvectorcalc(run,crnt);
    cellproperties(crnt,run);
    setsplitcalc(crnt,run);
    
    crnt=crnt->lnxt;
  }

}

void cellproperties(struct BRANCH * crnt, struct RUN * run) {
  
  int in;
  
  areaproperties(crnt,run);
  volmproperties(crnt,run);
	
  crnt->cl->xc=0.0;
	crnt->cl->yc=0.0;
	crnt->cl->zc=0.0;
  
  for(in=0;in<crnt->nlnd;in++){  
    crnt->cl->xc+=crnt->cl->nd[in].x;
    crnt->cl->yc+=crnt->cl->nd[in].y;
    crnt->cl->zc+=crnt->cl->nd[in].z;
  }
  crnt->cl->xc /= crnt->nlnd;
  crnt->cl->yc /= crnt->nlnd;
  crnt->cl->zc /= crnt->nlnd;

}
	
void areaproperties(struct BRANCH * crnt, struct RUN * run) {

	int f,it,ing,type,fopp;
	double dx1,dy1,dz1;
	double dx2,dy2,dz2;
  struct BRANCH * oppo;


	for(f=0;f<crnt->nlfc;f++){             // Loop for faces 
	  for(ing=0;ing<crnt->nsfc[f];ing++){  // Loop for neigs 
      crnt->cl->Area[f][ing]=0.0;	
	    if (crnt->nsfc[f]==1){
	      oppo=crnt;
	      fopp=f;
	    } else {
        oppo=crnt->neigtr[f][ing];
	      fopp=crnt->neigfc[f];
	    }
      type=oppo->type;
      for(it=0;it<crnt->cl->fc[f].nnds-2;it++){  // Loop for triangles
        dx1 = oppo->cl->nd[fcnd2elnd(fopp,1+it,type)].x - oppo->cl->nd[fcnd2elnd(fopp,0,type)].x;
        dx2 = oppo->cl->nd[fcnd2elnd(fopp,2+it,type)].x - oppo->cl->nd[fcnd2elnd(fopp,0,type)].x;
        dy1 = oppo->cl->nd[fcnd2elnd(fopp,1+it,type)].y - oppo->cl->nd[fcnd2elnd(fopp,0,type)].y;
        dy2 = oppo->cl->nd[fcnd2elnd(fopp,2+it,type)].y - oppo->cl->nd[fcnd2elnd(fopp,0,type)].y;
        dz1 = oppo->cl->nd[fcnd2elnd(fopp,1+it,type)].z - oppo->cl->nd[fcnd2elnd(fopp,0,type)].z;
        dz2 = oppo->cl->nd[fcnd2elnd(fopp,2+it,type)].z - oppo->cl->nd[fcnd2elnd(fopp,0,type)].z;
        crnt->cl->Area[f][ing]+=fabs(0.5*sqrt(pow(dy1*dz2-dz1*dy2,2)+pow(dz1*dx2-dx1*dz2,2)+pow(dx1*dy2-dy1*dx2,2)));
      }
	  }
  }

}


void volmproperties(struct BRANCH * crnt, struct RUN * run) {

	int in,it,num_tetra;
  int vn[6][3];
	double dx[3],dy[3],dz[3];
	double cross[3];

	if (crnt->type==0){

    num_tetra=6; 
	  vn[0][0]=1; vn[0][1]=3; vn[0][2]=5;
	  vn[1][0]=5; vn[1][1]=3; vn[1][2]=7;
	  vn[2][0]=5; vn[2][1]=7; vn[2][2]=4;
	  vn[3][0]=3; vn[3][1]=2; vn[3][2]=7;
	  vn[4][0]=7; vn[4][1]=2; vn[4][2]=6;
	  vn[5][0]=7; vn[5][1]=6; vn[5][2]=4;
	
  } else if (crnt->type==1){

    num_tetra=1;
    vn[0][0]=1; vn[0][1]=2; vn[0][2]=3;

  } else if (crnt->type==2){

    num_tetra=3;
    vn[0][0]=1; vn[0][1]=2; vn[0][2]=5;
    vn[1][0]=1; vn[1][1]=5; vn[1][2]=4;
    vn[2][0]=4; vn[2][1]=5; vn[2][2]=3;

  } else {
    printf("volmproperties: Not supported element found, exiting");
    exit(0);
  }

  crnt->cl->Vol=0.0;
  for(it=0; it<num_tetra;it++){  

    for(in=0; in<3;in++){  
      dx[in]=(crnt->cl->nd[vn[it][in]].x-crnt->cl->nd[0].x);
      dy[in]=(crnt->cl->nd[vn[it][in]].y-crnt->cl->nd[0].y);
      dz[in]=(crnt->cl->nd[vn[it][in]].z-crnt->cl->nd[0].z);
	  }
    
    cross[0]=dy[0]*dz[1]-dz[0]*dy[1];
    cross[1]=dz[0]*dx[1]-dx[0]*dz[1];
    cross[2]=dx[0]*dy[1]-dy[0]*dx[1];
	  crnt->cl->Vol+=(cross[0]*dx[2]+cross[1]*dy[2]+cross[2]*dz[2])/6.0;
    
  }

}