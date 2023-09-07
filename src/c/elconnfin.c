#include "strdata.h"

void elconnfin (BRANCH * brch,int ifc) {
  
  struct BRANCH * crnt;
  struct BRANCH * neig;
  int ifc2,iang,fctype,ifnd,ifnds,ineig,ipnt,indnd,ikid;

  neig=brch->neigtr[ifc][0];
  
  if (neig!=NULL&&neig->nkids!=0){

    ifc2=brch->neigfc[ifc];
    iang=brch->neigag[ifc];
    fctype=brch->cl->fc[ifc].type;

    if (neig->nsplit==8){
      
      brch->lfc_neigs[ifc]=4;
      brch->neigtr[ifc][0]=neig->kids[fcrefnd2elnd(ifc2,fcndrot(0,iang,fctype),neig->type)];
      brch->neigtr[ifc][1]=neig->kids[fcrefnd2elnd(ifc2,fcndrot(1,iang,fctype),neig->type)];
      brch->neigtr[ifc][2]=neig->kids[fcrefnd2elnd(ifc2,fcndrot(2,iang,fctype),neig->type)];
      brch->neigtr[ifc][3]=neig->kids[fcrefnd2elnd(ifc2,fcndrot(3,iang,fctype),neig->type)];
    
    } else if (brch->nsplit==4) {

      brch->lfc_neigs[ifc]=2;

      ineig=0;
      for (ifnds=0;ifnds<4;ifnds++) {

        if      (ifnds==0){ifnd=0;}
        else if (ifnds==1){ifnd=1;}
        else if (ifnds==2){ifnd=3;}
        else if (ifnds==3){ifnd=2;}

        if (elnd2fcnd(brch->fcqsplit,fcnd2elnd(ifc,ifnd,brch->type),brch->type)!=-1){

          ipnt=fcndrot(ifnd,iang,fctype);
          indnd=fcrefnd2elnd(ifc2,ipnt,neig->type);
          if (neig->nsplit==4){
            
            if      ( indnd==fcnd2elnd(neig->fcqsplit,0,neig->type) || indnd==fcopnd2elnd(neig->fcqsplit,0,neig->type) ) { ikid=0; }
            else if ( indnd==fcnd2elnd(neig->fcqsplit,1,neig->type) || indnd==fcopnd2elnd(neig->fcqsplit,1,neig->type) ) { ikid=1; }
            else if ( indnd==fcnd2elnd(neig->fcqsplit,2,neig->type) || indnd==fcopnd2elnd(neig->fcqsplit,2,neig->type) ) { ikid=2; }
            else if ( indnd==fcnd2elnd(neig->fcqsplit,3,neig->type) || indnd==fcopnd2elnd(neig->fcqsplit,3,neig->type) ) { ikid=3; }
          
          } else if (neig->nsplit==2){

            if      (indnd<3) {ikid=0;}
            else if (indnd>=3){ikid=1;}

          }

          brch->neigtr[ifc][ineig]=neig->kids[ikid];
          brch->neignd[ifc][ineig]=ifnd;
          ineig++;
        }
      }
    
    } else if (brch->nsplit==2) {
    
      ineig=0;
      fctype=brch->cl->fc[ifc].type;
      if (fctype==0){
        brch->lfc_neigs[ifc]=2;
        
        ineig=0;
        for (ifnds=0;ifnds<2;ifnds++){
          if      (ifnds==0){ifnd=0;}
          else if (ifnds==1){ifnd=2;}

          ipnt  = fcndrot(ifnd,iang,fctype);
          indnd = fcrefnd2elnd(ifc2,ipnt,neig->type);

          if (neig->nsplit==2){
            if(indnd<3) {ikid=0;}
            if(indnd>=3){ikid=1;}
          }

          if (neig->nsplit==4){
            if      (indnd==fcnd2elnd(neig->fcqsplit,0,neig->type)||indnd==fcopnd2elnd(neig->fcqsplit,0,neig->type)){ikid=0;}
            else if (indnd==fcnd2elnd(neig->fcqsplit,1,neig->type)||indnd==fcopnd2elnd(neig->fcqsplit,1,neig->type)){ikid=1;}
            else if (indnd==fcnd2elnd(neig->fcqsplit,2,neig->type)||indnd==fcopnd2elnd(neig->fcqsplit,2,neig->type)){ikid=2;}
            else if (indnd==fcnd2elnd(neig->fcqsplit,3,neig->type)||indnd==fcopnd2elnd(neig->fcqsplit,3,neig->type)){ikid=3;}
          }
          brch->neigtr[ifc][ineig]=neig->kids[ikid];
          ineig++;
        }
      
      } else if (fctype==1) {
        
        brch->lfc_neigs[ifc]=1;
        fctype=brch->cl->fc[ifc].type;
        if (ifc2==3){ikid=0;}
        if (ifc2==4){ikid=1;}
        brch->neigtr[ifc][0]=neig->kids[ikid];
      
      }

    }
  }
  brch->nsfc[ifc]=max(1,brch->lfc_neigs[ifc]);
}


void neigsfin (struct RUN * run) {
  
  int ifc;
  struct BRANCH * crnt;

  crnt=run->topo->locl; 
  while (crnt!=NULL){
    crnt->split=0;
    for (ifc=0;ifc<crnt->nlfc;ifc++){
      elconnfin(crnt,ifc);
	    crnt->nsfc[ifc]=max(1,crnt->lfc_neigs[ifc]);
    }
    crnt=crnt->lnxt;
  }

}