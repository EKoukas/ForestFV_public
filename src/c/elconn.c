#include "strdata.h"

void elconn (struct RUN * run, struct BRANCH * brch,int ifc) {

  struct BRANCH * curr; 
  struct BRANCH * neig; 
  int ifc2,iang,ilevel,ikid,ipnt,indnd,ifnd;

  curr=brch->keen[ifc];


  if (curr!=brch->prnt||curr->root==-1){
    brch->neigtr[ifc][0]=brch->keen[ifc];
 
    brch->lfc_neigs[ifc]=1;
    brch->neigfc[ifc]=brch->keenfc[ifc];
    brch->neigag[ifc]=brch->keenfcangle[ifc];

  }

  if ((curr==brch->prnt)&&(curr->root!=-1)){
    while(curr->level!=curr->keen[ifc]->level&&curr->keen[ifc]->root!=-1){
      curr=curr->keen[ifc];
    }

    if (curr->keen[ifc]->root!=-1){

      neig=curr->keen[ifc];
      ifc2=curr->keenfc[ifc];
      iang=curr->keenfcangle[ifc];
      ilevel=curr->level;
      while((neig->nkids!=0)&&(neig->level<brch->level)){
        ilevel=neig->level;

        if (neig->nsplit==8){
          neig=neig->kids[fcrefnd2elnd(ifc2,fcndrot(elnd2fcnd(ifc,brch->adrs[ilevel+1],brch->type),iang,brch->cl->fc[ifc].type),neig->type)];
        }
        if (brch->nsplit==4&&neig->nsplit==4){

          ipnt=fcndrot(elnd2fcnd(ifc,fcnd2elnd(brch->fcqsplit,brch->adrs[ilevel+1],brch->type),brch->type),iang,brch->cl->fc[ifc].type);
          indnd=fcrefnd2elnd(ifc2,ipnt,neig->type);
          ikid =-1;
          if (indnd==fcnd2elnd(neig->fcqsplit,0,neig->type)||indnd==fcopnd2elnd(neig->fcqsplit,0,neig->type)){ikid=0;}
          if (indnd==fcnd2elnd(neig->fcqsplit,1,neig->type)||indnd==fcopnd2elnd(neig->fcqsplit,1,neig->type)){ikid=1;}
          if (indnd==fcnd2elnd(neig->fcqsplit,2,neig->type)||indnd==fcopnd2elnd(neig->fcqsplit,2,neig->type)){ikid=2;}
          if (indnd==fcnd2elnd(neig->fcqsplit,3,neig->type)||indnd==fcopnd2elnd(neig->fcqsplit,3,neig->type)){ikid=3;}
          if (ikid==-1){printf("ikid -1\n");}
          neig=neig->kids[ikid];
        }
        if (brch->nsplit==2&&neig->nsplit==4){
          if (brch->type==0){
            if (brch->adrs[ilevel+1]==0){
              if (fcnd2elnd(0,ifc,brch->type)==fcnd2elnd(0,brch->ctfc,brch->type)||fcnd2elnd(0,ifc,brch->type)==fcnd2elnd(2,brch->ctfc,brch->type)){ifnd=0;}
              if (fcnd2elnd(1,ifc,brch->type)==fcnd2elnd(0,brch->ctfc,brch->type)||fcnd2elnd(1,ifc,brch->type)==fcnd2elnd(2,brch->ctfc,brch->type)){ifnd=1;}
              if (fcnd2elnd(2,ifc,brch->type)==fcnd2elnd(0,brch->ctfc,brch->type)||fcnd2elnd(2,ifc,brch->type)==fcnd2elnd(2,brch->ctfc,brch->type)){ifnd=2;}
              if (fcnd2elnd(3,ifc,brch->type)==fcnd2elnd(0,brch->ctfc,brch->type)||fcnd2elnd(3,ifc,brch->type)==fcnd2elnd(2,brch->ctfc,brch->type)){ifnd=3;}
            }
            if (brch->adrs[ilevel+1]==1){
              if (fcnd2elnd(0,ifc,brch->type)==fcnd2elnd(0,fcopp(brch->ctfc,brch->type),brch->type)||fcnd2elnd(0,ifc,brch->type)==fcnd2elnd(2,fcopp(brch->ctfc,brch->type),brch->type)){ifnd=0;}
              if (fcnd2elnd(1,ifc,brch->type)==fcnd2elnd(0,fcopp(brch->ctfc,brch->type),brch->type)||fcnd2elnd(1,ifc,brch->type)==fcnd2elnd(2,fcopp(brch->ctfc,brch->type),brch->type)){ifnd=1;}
              if (fcnd2elnd(2,ifc,brch->type)==fcnd2elnd(0,fcopp(brch->ctfc,brch->type),brch->type)||fcnd2elnd(2,ifc,brch->type)==fcnd2elnd(2,fcopp(brch->ctfc,brch->type),brch->type)){ifnd=2;}
              if (fcnd2elnd(3,ifc,brch->type)==fcnd2elnd(0,fcopp(brch->ctfc,brch->type),brch->type)||fcnd2elnd(3,ifc,brch->type)==fcnd2elnd(2,fcopp(brch->ctfc,brch->type),brch->type)){ifnd=3;}
            }

            ipnt=fcndrot(ifnd,iang,brch->cl->fc[ifc].type);
            indnd=fcrefnd2elnd(ifc2,ipnt,neig->type);
          }
          if (brch->type==2){
            if (brch->adrs[ilevel+1]==0){ifnd=0;}
            if (brch->adrs[ilevel+1]==1){ifnd=2;}

            ipnt=fcndrot(ifnd,iang,brch->cl->fc[ifc].type);
            indnd=fcrefnd2elnd(ifc2,ipnt,neig->type);
          }

          ikid =-1;
          if (indnd==fcnd2elnd(neig->fcqsplit,0,neig->type)||indnd==fcopnd2elnd(neig->fcqsplit,0,neig->type)){ikid=0;}
          if (indnd==fcnd2elnd(neig->fcqsplit,1,neig->type)||indnd==fcopnd2elnd(neig->fcqsplit,1,neig->type)){ikid=1;}
          if (indnd==fcnd2elnd(neig->fcqsplit,2,neig->type)||indnd==fcopnd2elnd(neig->fcqsplit,2,neig->type)){ikid=2;}
          if (indnd==fcnd2elnd(neig->fcqsplit,3,neig->type)||indnd==fcopnd2elnd(neig->fcqsplit,3,neig->type)){ikid=3;}
          if (ikid==-1){printf("ikid -1\n");}
          neig=neig->kids[ikid];

        }

        if (brch->nsplit==4&&neig->nsplit==2){
          ipnt=fcndrot(elnd2fcnd(ifc,fcnd2elnd(brch->fcqsplit,brch->adrs[ilevel+1],brch->type),brch->type),iang,brch->cl->fc[ifc].type);
          indnd=fcrefnd2elnd(ifc2,ipnt,neig->type);

          if (neig->type==0){
            if (indnd==fcnd2elnd(neig->ctfc,0,neig->type)||indnd==fcnd2elnd(neig->ctfc,1,neig->type)||indnd==fcnd2elnd(neig->ctfc,2,neig->type)||indnd==fcnd2elnd(neig->ctfc,3,neig->type)){ikid=0;}
            if (indnd==fcopnd2elnd(neig->ctfc,0,neig->type)||indnd==fcopnd2elnd(neig->ctfc,1,neig->type)||indnd==fcopnd2elnd(neig->ctfc,2,neig->type)||indnd==fcopnd2elnd(neig->ctfc,3,neig->type)){ikid=1;}
          }  
          if (neig->type==2){
            if (indnd==0||indnd==1||indnd==2){ikid=0;}
            if (indnd==3||indnd==4||indnd==5){ikid=1;}
          }  
          neig=neig->kids[ikid];
        }

        if (brch->nsplit==2&&neig->nsplit==2){
          if (brch->type==2){
            if (ifc<3){
              if(iang<=1){ikid=brch->adrs[ilevel+1];}
              if(iang>1){ikid=(1+brch->adrs[ilevel+1])%2;}
            }
            if (ifc2==3){ikid=0;}
            if (ifc2==4){ikid=1;}
            neig=neig->kids[ikid];
          }
          if (brch->type==0){
            if (ifc!=brch->ctfc||ifc!=fcopp(brch->ctfc,brch->type)){
              ikid=brch->adrs[ilevel+1];
            }
            if (ifc==brch->ctfc){ikid=1;}
            if (ifc==fcopp(brch->ctfc,brch->type)){ikid=0;}
            neig=neig->kids[ikid];
          }
        }
        if (brch->nsplit==1&&neig->nsplit==1){
          if (ifc<3){
            if(iang<=1){ikid=brch->adrs[ilevel+1];}
            if(iang>1){ikid=(1+brch->adrs[ilevel+1])%2;}
          }
          if (ifc2==3){ikid=0;}
          if (ifc2==4){ikid=1;}
          neig=neig->kids[ikid];
        }
      }
 
      brch->neigtr[ifc][0]=neig;
      brch->lfc_neigs[ifc]=1;
      brch->neigfc[ifc]=ifc2;
      brch->neigag[ifc]=iang;
    }
    if (curr->keen[ifc]->root==-1){
      brch->neigtr[ifc][0]=curr->keen[ifc];
      brch->lfc_neigs[ifc]=0;
      brch->neigfc[ifc]=curr->keenfc[ifc];
      brch->neigag[ifc]=curr->keenfcangle[ifc];
    }
  }

  brch->nsfc[ifc]=max(1,brch->lfc_neigs[ifc]);
}
