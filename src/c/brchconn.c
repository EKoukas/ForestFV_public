#include "strdata.h"

/*
 * Important information for tetras:
 * a) Restriction A: inner elements 4 5 6 7 are oriented so that keen neighbouring faces show the same face
 * b) Restriction B: inner elements 4 5 6 7 are oriented so that keen neighbouring faces have the reverse orientation ie: middlepoint 0 -> big face 0 etc for midpoibts 1 and 2
 * c) For triangular faces we first rotate and then mirror thus if the opposite to 0 is 2 then the orientation angle is 1 if the opposite to 0 is 1 then the orientation angle is 2. If the oposite to 0 is 0 the orientation angle is 0.
 */ 

void brchconn (BRANCH * brch) {
  
  int ikid,ikeen,ifc,i,j,k,is,id,ifd;
  int fcs[3][2];
  int idfc[6][3];

  if (brch->type==0&brch->nsplit==8){
    fcs[0][0]=3;
    fcs[0][1]=1;
    fcs[1][0]=0;
    fcs[1][1]=2;
    fcs[2][0]=4;
    fcs[2][1]=5;
    for (i=0;i<2;i++){
      for (j=0;j<2;j++){
        for (k=0;k<2;k++){
          ikid=i+2*j+4*k;
        
          for (id=0;id<3;id++){
            for (is=0;is<2;is++){

              ifc   = fcs[id][is];
              if ((id==0&&((i+2*is-1>=0)&(i+2*is-1<=1)))||
                  (id==1&&((j+2*is-1>=0)&(j+2*is-1<=1)))||
                  (id==2&&((k+2*is-1>=0)&(k+2*is-1<=1)))){
                    ikeen = ikid + (2*is-1)*pow(2,id);
                    brch->kids[ikid]->keen[ifc]=brch->kids[ikeen];
                    brch->kids[ikid]->keenfc[ifc]=fcs[id][1-is]  ;
                    brch->kids[ikid]->keenfcangle[ifc]=0;
              } else {
                    brch->kids[ikid]->keen[ifc]=brch;
                    brch->kids[ikid]->keenfc[ifc]=ifc;
                    brch->kids[ikid]->keenfcangle[ifc]=0;
              }
            }
          }
        }
      }
    }
  }


  if (brch->type==0&&brch->nsplit==4){
    fcs[0][0]=2;
    fcs[0][1]=0;
    fcs[1][0]=3;
    fcs[1][1]=1;
    fcs[2][0]=4;
    fcs[2][1]=5;

    idfc[0][0]=1; idfc[0][1]=2; idfc[0][2]=0;
    idfc[1][0]=2; idfc[1][1]=0; idfc[1][2]=2;
    idfc[2][0]=2; idfc[2][1]=1; idfc[2][2]=0;
    idfc[3][0]=0; idfc[3][1]=2; idfc[3][2]=1;
    idfc[4][0]=1; idfc[4][1]=0; idfc[4][2]=2;
    idfc[5][0]=0; idfc[5][1]=1; idfc[5][2]=2;

    for (ikid=0;ikid<brch->nkids;ikid++){
      for (ifc=0;ifc<brch->nlfc;ifc++){
        brch->kids[ikid]->keen[ifc]=brch;
        brch->kids[ikid]->keenfc[ifc]=ifc;
        brch->kids[ikid]->keenfcangle[ifc]=0;
      }
    }


    brch->kids[0]->keen[fcs[idfc[brch->fcqsplit][0]][1]]=brch->kids[1];
    brch->kids[0]->keen[fcs[idfc[brch->fcqsplit][1]][1]]=brch->kids[3];

    brch->kids[1]->keen[fcs[idfc[brch->fcqsplit][0]][0]]=brch->kids[0];
    brch->kids[1]->keen[fcs[idfc[brch->fcqsplit][1]][1]]=brch->kids[2];

    brch->kids[2]->keen[fcs[idfc[brch->fcqsplit][0]][0]]=brch->kids[3];
    brch->kids[2]->keen[fcs[idfc[brch->fcqsplit][1]][0]]=brch->kids[1];

    brch->kids[3]->keen[fcs[idfc[brch->fcqsplit][0]][1]]=brch->kids[2];
    brch->kids[3]->keen[fcs[idfc[brch->fcqsplit][1]][0]]=brch->kids[0];

    brch->kids[0]->keenfc[fcs[idfc[brch->fcqsplit][0]][1]]=fcs[idfc[brch->fcqsplit][0]][0];
    brch->kids[0]->keenfc[fcs[idfc[brch->fcqsplit][1]][1]]=fcs[idfc[brch->fcqsplit][1]][0];

    brch->kids[1]->keenfc[fcs[idfc[brch->fcqsplit][0]][0]]=fcs[idfc[brch->fcqsplit][0]][1];
    brch->kids[1]->keenfc[fcs[idfc[brch->fcqsplit][1]][1]]=fcs[idfc[brch->fcqsplit][1]][0];

    brch->kids[2]->keenfc[fcs[idfc[brch->fcqsplit][0]][0]]=fcs[idfc[brch->fcqsplit][0]][1];
    brch->kids[2]->keenfc[fcs[idfc[brch->fcqsplit][1]][0]]=fcs[idfc[brch->fcqsplit][1]][1];

    brch->kids[3]->keenfc[fcs[idfc[brch->fcqsplit][0]][1]]=fcs[idfc[brch->fcqsplit][0]][0];
    brch->kids[3]->keenfc[fcs[idfc[brch->fcqsplit][1]][0]]=fcs[idfc[brch->fcqsplit][1]][1];

    brch->kids[0]->keenfcangle[fcs[idfc[brch->fcqsplit][0]][1]]=0;
    brch->kids[0]->keenfcangle[fcs[idfc[brch->fcqsplit][1]][1]]=0;

    brch->kids[1]->keenfcangle[fcs[idfc[brch->fcqsplit][0]][0]]=0;
    brch->kids[1]->keenfcangle[fcs[idfc[brch->fcqsplit][1]][1]]=0;

    brch->kids[2]->keenfcangle[fcs[idfc[brch->fcqsplit][0]][0]]=0;
    brch->kids[2]->keenfcangle[fcs[idfc[brch->fcqsplit][1]][0]]=0;

    brch->kids[3]->keenfcangle[fcs[idfc[brch->fcqsplit][0]][1]]=0;
    brch->kids[3]->keenfcangle[fcs[idfc[brch->fcqsplit][1]][0]]=0;

  }

  if (brch->type==2&&brch->nsplit==8){
    for (ikid=0;ikid<brch->nkids;ikid++){
      for (ifc=0;ifc<brch->nlfc;ifc++){
        brch->kids[ikid]->keen[ifc]=brch;
        brch->kids[ikid]->keenfc[ifc]=ifc;
        brch->kids[ikid]->keenfcangle[ifc]=0;
      }
    }
    brch->kids[0]->keen[1]=brch->kids[6];
    brch->kids[1]->keen[2]=brch->kids[6];
    brch->kids[2]->keen[0]=brch->kids[6];

    brch->kids[0]->keenfc[1]=1;
    brch->kids[1]->keenfc[2]=2;
    brch->kids[2]->keenfc[0]=0;

    brch->kids[0]->keenfcangle[1]=1;
    brch->kids[1]->keenfcangle[2]=1;
    brch->kids[2]->keenfcangle[0]=3;

    brch->kids[0]->keen[4]=brch->kids[3];
    brch->kids[1]->keen[4]=brch->kids[4];
    brch->kids[2]->keen[4]=brch->kids[5];

    brch->kids[0]->keenfc[4]=3;
    brch->kids[1]->keenfc[4]=3;
    brch->kids[2]->keenfc[4]=3;

    brch->kids[0]->keenfcangle[4]=0;
    brch->kids[1]->keenfcangle[4]=0;
    brch->kids[2]->keenfcangle[4]=0;



    brch->kids[3]->keen[1]=brch->kids[7];
    brch->kids[4]->keen[2]=brch->kids[7];
    brch->kids[5]->keen[0]=brch->kids[7];

    brch->kids[3]->keenfc[1]=1;
    brch->kids[4]->keenfc[2]=2;
    brch->kids[5]->keenfc[0]=0;

    brch->kids[3]->keenfcangle[1]=1;
    brch->kids[4]->keenfcangle[2]=1;
    brch->kids[5]->keenfcangle[0]=3;

    brch->kids[3]->keen[3]=brch->kids[0];
    brch->kids[4]->keen[3]=brch->kids[1];
    brch->kids[5]->keen[3]=brch->kids[2];

    brch->kids[3]->keenfc[3]=4;
    brch->kids[4]->keenfc[3]=4;
    brch->kids[5]->keenfc[3]=4;

    brch->kids[3]->keenfcangle[3]=0;
    brch->kids[4]->keenfcangle[3]=0;
    brch->kids[5]->keenfcangle[3]=0;

    brch->kids[6]->keen[0]=brch->kids[2];
    brch->kids[6]->keen[1]=brch->kids[0];
    brch->kids[6]->keen[2]=brch->kids[1];
    brch->kids[6]->keen[4]=brch->kids[7];

    brch->kids[6]->keenfc[0]=0;
    brch->kids[6]->keenfc[1]=1;
    brch->kids[6]->keenfc[2]=2;
    brch->kids[6]->keenfc[4]=3;

    brch->kids[6]->keenfcangle[0]=3;
    brch->kids[6]->keenfcangle[1]=1;
    brch->kids[6]->keenfcangle[2]=1;
    brch->kids[6]->keenfcangle[4]=0;

    brch->kids[7]->keen[0]=brch->kids[5];
    brch->kids[7]->keen[1]=brch->kids[3];
    brch->kids[7]->keen[2]=brch->kids[4];
    brch->kids[7]->keen[3]=brch->kids[6];

    brch->kids[7]->keenfc[0]=0;
    brch->kids[7]->keenfc[1]=1;
    brch->kids[7]->keenfc[2]=2;
    brch->kids[7]->keenfc[3]=4;

    brch->kids[7]->keenfcangle[0]=3;
    brch->kids[7]->keenfcangle[1]=1;
    brch->kids[7]->keenfcangle[2]=1;
    brch->kids[7]->keenfcangle[3]=0;

  }

  if (brch->type==1&&brch->nsplit==8){
    for (ikid=0;ikid<brch->nkids;ikid++){
      for (ifc=0;ifc<brch->nlfc;ifc++){
        brch->kids[ikid]->keen[ifc]=brch;
        brch->kids[ikid]->keenfc[ifc]=ifc;
        brch->kids[ikid]->keenfcangle[ifc]=0;
      }
    }
    brch->kids[0]->keen[2]=brch->kids[5];
    brch->kids[1]->keen[3]=brch->kids[4];
    brch->kids[2]->keen[1]=brch->kids[6];
    brch->kids[3]->keen[0]=brch->kids[7];

    brch->kids[0]->keenfc[2]=2;
    brch->kids[1]->keenfc[3]=3;
    brch->kids[2]->keenfc[1]=1;
    brch->kids[3]->keenfc[0]=0;

    brch->kids[0]->keenfcangle[2]=1; //2
    brch->kids[1]->keenfcangle[3]=2; //1
    brch->kids[2]->keenfcangle[1]=0;
    brch->kids[3]->keenfcangle[0]=0;

    brch->kids[4]->keen[0]=brch;
    brch->kids[4]->keen[1]=brch->kids[6];
    brch->kids[4]->keen[2]=brch->kids[5];
    brch->kids[4]->keen[3]=brch->kids[1];

    brch->kids[4]->keenfc[0]=0;
    brch->kids[4]->keenfc[1]=3;
    brch->kids[4]->keenfc[2]=3;
    brch->kids[4]->keenfc[3]=3;

    brch->kids[4]->keenfcangle[0]=0;
    brch->kids[4]->keenfcangle[1]=1;
    brch->kids[4]->keenfcangle[2]=0; // to 5
    brch->kids[4]->keenfcangle[3]=2;


    brch->kids[5]->keen[0]=brch->kids[7];
    brch->kids[5]->keen[1]=brch;
    brch->kids[5]->keen[2]=brch->kids[0];
    brch->kids[5]->keen[3]=brch->kids[4];

    brch->kids[5]->keenfc[0]=2;
    brch->kids[5]->keenfc[1]=1;
    brch->kids[5]->keenfc[2]=2; //2
    brch->kids[5]->keenfc[3]=2;

    brch->kids[5]->keenfcangle[0]=2;
    brch->kids[5]->keenfcangle[1]=0;
    brch->kids[5]->keenfcangle[2]=1; // to 0
    brch->kids[5]->keenfcangle[3]=0;


    brch->kids[6]->keen[0]=brch->kids[7];
    brch->kids[6]->keen[1]=brch->kids[2];
    brch->kids[6]->keen[2]=brch;
    brch->kids[6]->keen[3]=brch->kids[4];

    brch->kids[6]->keenfc[0]=1;
    brch->kids[6]->keenfc[1]=1;
    brch->kids[6]->keenfc[2]=2;
    brch->kids[6]->keenfc[3]=1;

    brch->kids[6]->keenfcangle[0]=0;
    brch->kids[6]->keenfcangle[1]=0;
    brch->kids[6]->keenfcangle[2]=0;
    brch->kids[6]->keenfcangle[3]=1;

    brch->kids[7]->keen[0]=brch->kids[3];
    brch->kids[7]->keen[1]=brch->kids[6];
    brch->kids[7]->keen[2]=brch->kids[5];
    brch->kids[7]->keen[3]=brch;

    brch->kids[7]->keenfc[0]=0;
    brch->kids[7]->keenfc[1]=0;
    brch->kids[7]->keenfc[2]=0;
    brch->kids[7]->keenfc[3]=3;

    brch->kids[7]->keenfcangle[0]=0;
    brch->kids[7]->keenfcangle[1]=0;
    brch->kids[7]->keenfcangle[2]=2;
    brch->kids[7]->keenfcangle[3]=0;

  }

  if (brch->type==2&&brch->nsplit==4){
    for (ikid=0;ikid<brch->nkids;ikid++){
      for (ifc=0;ifc<brch->nlfc;ifc++){
        brch->kids[ikid]->keen[ifc]=brch;
        brch->kids[ikid]->keenfc[ifc]=ifc;
        brch->kids[ikid]->keenfcangle[ifc]=0;
      }
    }

    if (brch->fcqsplit==3){
      brch->kids[0]->keen[1]=brch->kids[3];
      brch->kids[2]->keen[2]=brch->kids[3];
      brch->kids[1]->keen[0]=brch->kids[3];

      brch->kids[0]->keenfc[1]=1;
      brch->kids[2]->keenfc[2]=2;
      brch->kids[1]->keenfc[0]=0;

      brch->kids[0]->keenfcangle[1]=1;
      brch->kids[2]->keenfcangle[2]=1;
      brch->kids[1]->keenfcangle[0]=3;

      brch->kids[3]->keen[0]=brch->kids[1];
      brch->kids[3]->keen[1]=brch->kids[0];
      brch->kids[3]->keen[2]=brch->kids[2];

      brch->kids[3]->keenfc[0]=0;
      brch->kids[3]->keenfc[1]=1;
      brch->kids[3]->keenfc[2]=2;

      brch->kids[3]->keenfcangle[0]=3;
      brch->kids[3]->keenfcangle[1]=1;
      brch->kids[3]->keenfcangle[2]=1;
    }


    if (brch->fcqsplit==4){
      brch->kids[0]->keen[1]=brch->kids[3];
      brch->kids[1]->keen[2]=brch->kids[3];
      brch->kids[2]->keen[0]=brch->kids[3];
      
      brch->kids[0]->keenfc[1]=1;
      brch->kids[1]->keenfc[2]=2;
      brch->kids[2]->keenfc[0]=0;
      
      brch->kids[0]->keenfcangle[1]=1;
      brch->kids[1]->keenfcangle[2]=1;
      brch->kids[2]->keenfcangle[0]=3;
      
      brch->kids[3]->keen[0]=brch->kids[2];
      brch->kids[3]->keen[1]=brch->kids[0];
      brch->kids[3]->keen[2]=brch->kids[1];

      brch->kids[3]->keenfc[0]=0;
      brch->kids[3]->keenfc[1]=1;
      brch->kids[3]->keenfc[2]=2;

      brch->kids[3]->keenfcangle[0]=3;
      brch->kids[3]->keenfcangle[1]=1;
      brch->kids[3]->keenfcangle[2]=1;
    }

  }

  if (brch->type==2&&brch->nsplit==2){
    for (ikid=0;ikid<brch->nkids;ikid++){
      for (ifc=0;ifc<brch->nlfc;ifc++){
        brch->kids[ikid]->keen[ifc]=brch;
        brch->kids[ikid]->keenfc[ifc]=ifc;
        brch->kids[ikid]->keenfcangle[ifc]=0;
      }
    }

    brch->kids[0]->keen[4]=brch->kids[1];

    brch->kids[0]->keenfc[4]=3;

    brch->kids[0]->keenfcangle[4]=0;

    brch->kids[1]->keen[3]=brch->kids[0];

    brch->kids[1]->keenfc[3]=4;

    brch->kids[1]->keenfcangle[3]=0;


  }

  if (brch->type==0&&brch->nsplit==2){
    for (ikid=0;ikid<brch->nkids;ikid++){
      for (ifc=0;ifc<brch->nlfc;ifc++){
        brch->kids[ikid]->keen[ifc]=brch;
        brch->kids[ikid]->keenfc[ifc]=ifc;
        brch->kids[ikid]->keenfcangle[ifc]=0;
      }
    }

    for (ikid=0;ikid<brch->nkids-1;ikid++){
      brch->kids[ikid]->keen[fcopp(brch->ctfc,brch->type)]=brch->kids[ikid+1];
      brch->kids[ikid]->keenfc[fcopp(brch->ctfc,brch->type)]=brch->ctfc;
    }
    for (ikid=1;ikid<brch->nkids;ikid++){
      brch->kids[ikid]->keen[brch->ctfc]=brch->kids[ikid-1];
      brch->kids[ikid]->keenfc[brch->ctfc]=fcopp(brch->ctfc,brch->type);
    }

  }

  for (ikid=0;ikid<brch->nkids;ikid++){
    for (ifc=0;ifc<brch->nlfc;ifc++){
      if (brch->kids[ikid]->keen[ifc]->nkids==0){brch->kids[ikid]->cl->fc[ifc].bc=0;}
    }
  }

}
