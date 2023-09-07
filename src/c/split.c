#include "strdata.h"

void split (int linked, RUN * run, BRANCH * brch) {
  
  int ikid,ifc,i,j,k,ilvl,onfc,onifc;
  int f,iv;
  int ind,ilnd,nlnd,imap;
  int ikidnd,indnd;
  int ikidnd2,indnd2;
  
  brch->nkids = brch->nsplit;

  brch->kids = malloc(8*sizeof(BRANCH *));

  //brch->split=0;
  for (ikid=0;ikid<brch->nsplit;ikid++){

    brch->kids[ikid]=createbranch(brch->type,ikid,0,0,brch->level+2);
    cellallocation(brch->kids[ikid]);
    createlfc(brch->kids[ikid]);
    brch->kids[ikid]->prnt=brch;
    brch->kids[ikid]->root=brch->root;
    brch->kids[ikid]->split=brch->split;

    brch->kids[ikid]->part=brch->part;
    brch->kids[ikid]->level=brch->level+1;
    brch->kids[ikid]->fcqsplit=brch->fcqsplit;
    brch->kids[ikid]->nsplit=brch->nsplit;

    for (ilvl=0;ilvl<brch->level+1;ilvl++){
     brch->kids[ikid]->adrs[ilvl]=brch->adrs[ilvl];
    }
    brch->kids[ikid]->adrs[brch->level+1]=ikid;
    brch->kids[ikid]->tag=tagaddress(brch->kids[ikid]->adrs,brch->kids[ikid]->level);

    if (linked==1){
      if (ikid==0) { // connect with previus leaf vis versa
      
        if (run->topo->locl==brch) {
          run->topo->locl=brch->kids[ikid];
          brch->kids[ikid]->lprv=NULL;
          brch->kids[ikid]->lnxt=NULL;
        } else {
          brch->kids[ikid]->lprv=brch->lprv;
          brch->kids[ikid]->lnxt=NULL;
          brch->lprv->lnxt=brch->kids[ikid];
	      }
      }
   
      if (ikid==brch->nsplit-1){  // connect with previus kids vis versa and update the lnxt leaf connection vis versa unless final
        brch->kids[ikid]->lprv=brch->kids[ikid-1];
        brch->kids[ikid-1]->lnxt=brch->kids[ikid];

        brch->kids[ikid]->lnxt=brch->lnxt;
        if(brch->lnxt!=NULL) brch->lnxt->lprv=brch->kids[ikid];
      }
   
      if ((ikid!=brch->nsplit-1)&&(ikid!=0)){ // connect with previus kids vis versa
        brch->kids[ikid]->lprv=brch->kids[ikid-1];
        brch->kids[ikid]->lnxt=NULL;
        brch->kids[ikid-1]->lnxt=brch->kids[ikid];
      }
  
    }
    if (run->topo->locl==brch){run->topo->locl-brch->kids[0];} 
  }

  for (ikid=0;ikid<brch->nsplit;ikid++) {

    if        (brch->type==0&&brch->nsplit==8) {
      
      for (ind=0;ind<8;ind++){
        onfc=0;
        for (ifc=0;ifc<brch->nlfc;ifc++){
          if ((elnd2fcnd (ifc,ind,0)!=-1)&&(elnd2fcnd(ifc,ikid,0)!=-1)){
            onfc++;onifc=ifc;
          }
        }
        if (onfc==3){
          brch->kids[ikid]->cl->nd[ind].x=brch->cl->nd[ind].x;
          brch->kids[ikid]->cl->nd[ind].y=brch->cl->nd[ind].y;
          brch->kids[ikid]->cl->nd[ind].z=brch->cl->nd[ind].z;
        } else if (onfc==2){
          brch->kids[ikid]->cl->nd[ind].x=0.5*(brch->cl->nd[ind].x+brch->cl->nd[ikid].x);
          brch->kids[ikid]->cl->nd[ind].y=0.5*(brch->cl->nd[ind].y+brch->cl->nd[ikid].y);
          brch->kids[ikid]->cl->nd[ind].z=0.5*(brch->cl->nd[ind].z+brch->cl->nd[ikid].z);
        } else if (onfc==1){
          brch->kids[ikid]->cl->nd[ind].x=0.25*(brch->cl->nd[fcnd2elnd(onifc,0,0)].x+brch->cl->nd[fcnd2elnd(onifc,1,0)].x+brch->cl->nd[fcnd2elnd(onifc,2,0)].x+brch->cl->nd[fcnd2elnd(onifc,3,0)].x);
          brch->kids[ikid]->cl->nd[ind].y=0.25*(brch->cl->nd[fcnd2elnd(onifc,0,0)].y+brch->cl->nd[fcnd2elnd(onifc,1,0)].y+brch->cl->nd[fcnd2elnd(onifc,2,0)].y+brch->cl->nd[fcnd2elnd(onifc,3,0)].y);
          brch->kids[ikid]->cl->nd[ind].z=0.25*(brch->cl->nd[fcnd2elnd(onifc,0,0)].z+brch->cl->nd[fcnd2elnd(onifc,1,0)].z+brch->cl->nd[fcnd2elnd(onifc,2,0)].z+brch->cl->nd[fcnd2elnd(onifc,3,0)].z);
        } else if (onfc==0){
          brch->kids[ikid]->cl->nd[ind].x=0.125*(brch->cl->nd[0].x+brch->cl->nd[1].x+brch->cl->nd[2].x+brch->cl->nd[3].x+brch->cl->nd[4].x+brch->cl->nd[5].x+brch->cl->nd[6].x+brch->cl->nd[7].x);
          brch->kids[ikid]->cl->nd[ind].y=0.125*(brch->cl->nd[0].y+brch->cl->nd[1].y+brch->cl->nd[2].y+brch->cl->nd[3].y+brch->cl->nd[4].y+brch->cl->nd[5].y+brch->cl->nd[6].y+brch->cl->nd[7].y);
          brch->kids[ikid]->cl->nd[ind].z=0.125*(brch->cl->nd[0].z+brch->cl->nd[1].z+brch->cl->nd[2].z+brch->cl->nd[3].z+brch->cl->nd[4].z+brch->cl->nd[5].z+brch->cl->nd[6].z+brch->cl->nd[7].z);
        }
      }
    
    } else if (brch->type==0&&brch->nsplit==4) {
      
      for (ind=0;ind<4;ind++){

        ikidnd=fcnd2elnd(brch->fcqsplit,ikid,brch->type);
        indnd= fcnd2elnd(brch->fcqsplit,ind,brch->type);

        if (((ikid-ind)%2!=0)||(ikid==ind)){
          
          brch->kids[ikid]->cl->nd[indnd].x = 0.5*(brch->cl->nd[ikidnd].x+brch->cl->nd[indnd].x);
          brch->kids[ikid]->cl->nd[indnd].y = 0.5*(brch->cl->nd[ikidnd].y+brch->cl->nd[indnd].y);
          brch->kids[ikid]->cl->nd[indnd].z = 0.5*(brch->cl->nd[ikidnd].z+brch->cl->nd[indnd].z);
        
        } else if (((ikid-ind)%2==0)&&(ikid!=ind)){

          ikidnd2 = fcnd2elnd(brch->fcqsplit,(ikid+1)%4,brch->type);
          indnd2  = fcnd2elnd(brch->fcqsplit,(ind+1)%4,brch->type);

          brch->kids[ikid]->cl->nd[indnd].x = 0.25*(brch->cl->nd[ikidnd].x+brch->cl->nd[indnd].x+brch->cl->nd[ikidnd2].x+brch->cl->nd[indnd2].x);
          brch->kids[ikid]->cl->nd[indnd].y = 0.25*(brch->cl->nd[ikidnd].y+brch->cl->nd[indnd].y+brch->cl->nd[ikidnd2].y+brch->cl->nd[indnd2].y);
          brch->kids[ikid]->cl->nd[indnd].z = 0.25*(brch->cl->nd[ikidnd].z+brch->cl->nd[indnd].z+brch->cl->nd[ikidnd2].z+brch->cl->nd[indnd2].z);
        }
      }

      for (ind=0;ind<4;ind++){

        ikidnd=fcopnd2elnd(brch->fcqsplit,ikid,brch->type);
        indnd= fcopnd2elnd(brch->fcqsplit,ind,brch->type);

        if (((ikid-ind)%2!=0)||(ikid==ind)){
          brch->kids[ikid]->cl->nd[indnd].x = 0.5*(brch->cl->nd[ikidnd].x+brch->cl->nd[indnd].x);
          brch->kids[ikid]->cl->nd[indnd].y = 0.5*(brch->cl->nd[ikidnd].y+brch->cl->nd[indnd].y);
          brch->kids[ikid]->cl->nd[indnd].z = 0.5*(brch->cl->nd[ikidnd].z+brch->cl->nd[indnd].z);
        }

        if (((ikid-ind)%2==0)&&(ikid!=ind)){
          ikidnd2=fcopnd2elnd(brch->fcqsplit,(ikid+1)%4,brch->type);
          indnd2= fcopnd2elnd(brch->fcqsplit,(ind+1)%4,brch->type);

          brch->kids[ikid]->cl->nd[indnd].x = 0.25*(brch->cl->nd[ikidnd].x+brch->cl->nd[indnd].x+brch->cl->nd[ikidnd2].x+brch->cl->nd[indnd2].x);
          brch->kids[ikid]->cl->nd[indnd].y = 0.25*(brch->cl->nd[ikidnd].y+brch->cl->nd[indnd].y+brch->cl->nd[ikidnd2].y+brch->cl->nd[indnd2].y);
          brch->kids[ikid]->cl->nd[indnd].z = 0.25*(brch->cl->nd[ikidnd].z+brch->cl->nd[indnd].z+brch->cl->nd[ikidnd2].z+brch->cl->nd[indnd2].z);
        }
      }

    } else if (brch->type==2&&brch->nsplit==8) {

      if (ikid<6){
        for (ind=0;ind<6;ind++){
          brch->kids[ikid]->cl->nd[ind].x=brch->cl->nd[ikid].x+0.5*(brch->cl->nd[ind].x-brch->cl->nd[ikid].x);
          brch->kids[ikid]->cl->nd[ind].y=brch->cl->nd[ikid].y+0.5*(brch->cl->nd[ind].y-brch->cl->nd[ikid].y);
          brch->kids[ikid]->cl->nd[ind].z=brch->cl->nd[ikid].z+0.5*(brch->cl->nd[ind].z-brch->cl->nd[ikid].z);
        }
      } else if (ikid>=6) {
        brch->kids[6]->cl->nd[0].x=brch->kids[2]->cl->nd[1].x;
        brch->kids[6]->cl->nd[0].y=brch->kids[2]->cl->nd[1].y;
        brch->kids[6]->cl->nd[0].z=brch->kids[2]->cl->nd[1].z;

        brch->kids[6]->cl->nd[1].x=brch->kids[0]->cl->nd[2].x;
        brch->kids[6]->cl->nd[1].y=brch->kids[0]->cl->nd[2].y;
        brch->kids[6]->cl->nd[1].z=brch->kids[0]->cl->nd[2].z;

        brch->kids[6]->cl->nd[2].x=brch->kids[0]->cl->nd[1].x;
        brch->kids[6]->cl->nd[2].y=brch->kids[0]->cl->nd[1].y;
        brch->kids[6]->cl->nd[2].z=brch->kids[0]->cl->nd[1].z;

        brch->kids[6]->cl->nd[3].x=brch->kids[2]->cl->nd[4].x;
        brch->kids[6]->cl->nd[3].y=brch->kids[2]->cl->nd[4].y;
        brch->kids[6]->cl->nd[3].z=brch->kids[2]->cl->nd[4].z;

        brch->kids[6]->cl->nd[4].x=brch->kids[0]->cl->nd[5].x;
        brch->kids[6]->cl->nd[4].y=brch->kids[0]->cl->nd[5].y;
        brch->kids[6]->cl->nd[4].z=brch->kids[0]->cl->nd[5].z;

        brch->kids[6]->cl->nd[5].x=brch->kids[0]->cl->nd[4].x;
        brch->kids[6]->cl->nd[5].y=brch->kids[0]->cl->nd[4].y;
        brch->kids[6]->cl->nd[5].z=brch->kids[0]->cl->nd[4].z;

        brch->kids[7]->cl->nd[0].x=brch->kids[5]->cl->nd[1].x;
        brch->kids[7]->cl->nd[0].y=brch->kids[5]->cl->nd[1].y;
        brch->kids[7]->cl->nd[0].z=brch->kids[5]->cl->nd[1].z;

        brch->kids[7]->cl->nd[1].x=brch->kids[3]->cl->nd[2].x;
        brch->kids[7]->cl->nd[1].y=brch->kids[3]->cl->nd[2].y;
        brch->kids[7]->cl->nd[1].z=brch->kids[3]->cl->nd[2].z;

        brch->kids[7]->cl->nd[2].x=brch->kids[3]->cl->nd[1].x;
        brch->kids[7]->cl->nd[2].y=brch->kids[3]->cl->nd[1].y;
        brch->kids[7]->cl->nd[2].z=brch->kids[3]->cl->nd[1].z;

        brch->kids[7]->cl->nd[3].x=brch->kids[5]->cl->nd[4].x;
        brch->kids[7]->cl->nd[3].y=brch->kids[5]->cl->nd[4].y;
        brch->kids[7]->cl->nd[3].z=brch->kids[5]->cl->nd[4].z;

        brch->kids[7]->cl->nd[4].x=brch->kids[3]->cl->nd[5].x;
        brch->kids[7]->cl->nd[4].y=brch->kids[3]->cl->nd[5].y;
        brch->kids[7]->cl->nd[4].z=brch->kids[3]->cl->nd[5].z;

        brch->kids[7]->cl->nd[5].x=brch->kids[3]->cl->nd[4].x;
        brch->kids[7]->cl->nd[5].y=brch->kids[3]->cl->nd[4].y;
        brch->kids[7]->cl->nd[5].z=brch->kids[3]->cl->nd[4].z;
      }
    
    } else if (brch->type==2&&brch->nsplit==4) {

      if (ikid<3){
        for (ind=0;ind<3;ind++){
          ikidnd=fcnd2elnd(brch->fcqsplit,ikid,brch->type);
          indnd= fcnd2elnd(brch->fcqsplit,ind,brch->type);
            
          brch->kids[ikid]->cl->nd[indnd].x=brch->cl->nd[ikidnd].x+0.5*(brch->cl->nd[indnd].x-brch->cl->nd[ikidnd].x);
          brch->kids[ikid]->cl->nd[indnd].y=brch->cl->nd[ikidnd].y+0.5*(brch->cl->nd[indnd].y-brch->cl->nd[ikidnd].y);
          brch->kids[ikid]->cl->nd[indnd].z=brch->cl->nd[ikidnd].z+0.5*(brch->cl->nd[indnd].z-brch->cl->nd[ikidnd].z);
        }

        for (ind=0;ind<3;ind++) {   
          ikidnd=fcopnd2elnd(brch->fcqsplit,ikid,brch->type);
          indnd= fcopnd2elnd(brch->fcqsplit,ind,brch->type);

          brch->kids[ikid]->cl->nd[indnd].x=brch->cl->nd[ikidnd].x+0.5*(brch->cl->nd[indnd].x-brch->cl->nd[ikidnd].x);
          brch->kids[ikid]->cl->nd[indnd].y=brch->cl->nd[ikidnd].y+0.5*(brch->cl->nd[indnd].y-brch->cl->nd[ikidnd].y);
          brch->kids[ikid]->cl->nd[indnd].z=brch->cl->nd[ikidnd].z+0.5*(brch->cl->nd[indnd].z-brch->cl->nd[ikidnd].z);
        }
      }

      if (ikid==3){
        if (brch->fcqsplit==3){
          brch->kids[3]->cl->nd[0].x=brch->kids[1]->cl->nd[1].x;
          brch->kids[3]->cl->nd[0].y=brch->kids[1]->cl->nd[1].y;
          brch->kids[3]->cl->nd[0].z=brch->kids[1]->cl->nd[1].z;

          brch->kids[3]->cl->nd[1].x=brch->kids[0]->cl->nd[2].x;
          brch->kids[3]->cl->nd[1].y=brch->kids[0]->cl->nd[2].y;
          brch->kids[3]->cl->nd[1].z=brch->kids[0]->cl->nd[2].z;

          brch->kids[3]->cl->nd[2].x=brch->kids[0]->cl->nd[1].x;
          brch->kids[3]->cl->nd[2].y=brch->kids[0]->cl->nd[1].y;
          brch->kids[3]->cl->nd[2].z=brch->kids[0]->cl->nd[1].z;

          brch->kids[3]->cl->nd[3].x=brch->kids[1]->cl->nd[4].x;
          brch->kids[3]->cl->nd[3].y=brch->kids[1]->cl->nd[4].y;
          brch->kids[3]->cl->nd[3].z=brch->kids[1]->cl->nd[4].z;

          brch->kids[3]->cl->nd[4].x=brch->kids[0]->cl->nd[5].x;
          brch->kids[3]->cl->nd[4].y=brch->kids[0]->cl->nd[5].y;
          brch->kids[3]->cl->nd[4].z=brch->kids[0]->cl->nd[5].z;

          brch->kids[3]->cl->nd[5].x=brch->kids[0]->cl->nd[4].x;
          brch->kids[3]->cl->nd[5].y=brch->kids[0]->cl->nd[4].y;
          brch->kids[3]->cl->nd[5].z=brch->kids[0]->cl->nd[4].z;
        
        } else if (brch->fcqsplit==4){
          
          brch->kids[3]->cl->nd[0].x=brch->kids[2]->cl->nd[1].x;
          brch->kids[3]->cl->nd[0].y=brch->kids[2]->cl->nd[1].y;
          brch->kids[3]->cl->nd[0].z=brch->kids[2]->cl->nd[1].z;

          brch->kids[3]->cl->nd[1].x=brch->kids[0]->cl->nd[2].x;
          brch->kids[3]->cl->nd[1].y=brch->kids[0]->cl->nd[2].y;
          brch->kids[3]->cl->nd[1].z=brch->kids[0]->cl->nd[2].z;

          brch->kids[3]->cl->nd[2].x=brch->kids[0]->cl->nd[1].x;
          brch->kids[3]->cl->nd[2].y=brch->kids[0]->cl->nd[1].y;
          brch->kids[3]->cl->nd[2].z=brch->kids[0]->cl->nd[1].z;

          brch->kids[3]->cl->nd[3].x=brch->kids[2]->cl->nd[4].x;
          brch->kids[3]->cl->nd[3].y=brch->kids[2]->cl->nd[4].y;
          brch->kids[3]->cl->nd[3].z=brch->kids[2]->cl->nd[4].z;

          brch->kids[3]->cl->nd[4].x=brch->kids[0]->cl->nd[5].x;
          brch->kids[3]->cl->nd[4].y=brch->kids[0]->cl->nd[5].y;
          brch->kids[3]->cl->nd[4].z=brch->kids[0]->cl->nd[5].z;

          brch->kids[3]->cl->nd[5].x=brch->kids[0]->cl->nd[4].x;
          brch->kids[3]->cl->nd[5].y=brch->kids[0]->cl->nd[4].y;
          brch->kids[3]->cl->nd[5].z=brch->kids[0]->cl->nd[4].z;
        }
      }
    
    } else if (brch->type==1&&brch->nsplit==8) {
      
      if (ikid<4){
        for (ind=0;ind<4;ind++){
          brch->kids[ikid]->cl->nd[ind].x=brch->cl->nd[ikid].x+0.5*(brch->cl->nd[ind].x-brch->cl->nd[ikid].x);
          brch->kids[ikid]->cl->nd[ind].y=brch->cl->nd[ikid].y+0.5*(brch->cl->nd[ind].y-brch->cl->nd[ikid].y);
          brch->kids[ikid]->cl->nd[ind].z=brch->cl->nd[ikid].z+0.5*(brch->cl->nd[ind].z-brch->cl->nd[ikid].z);
        }
    
      } else if (ikid>=4){

        brch->kids[4]->cl->nd[0].x=brch->kids[1]->cl->nd[2].x;
        brch->kids[4]->cl->nd[0].y=brch->kids[1]->cl->nd[2].y;
        brch->kids[4]->cl->nd[0].z=brch->kids[1]->cl->nd[2].z;

        brch->kids[4]->cl->nd[1].x=brch->kids[2]->cl->nd[0].x;
        brch->kids[4]->cl->nd[1].y=brch->kids[2]->cl->nd[0].y;
        brch->kids[4]->cl->nd[1].z=brch->kids[2]->cl->nd[0].z;

        brch->kids[4]->cl->nd[2].x=brch->kids[1]->cl->nd[0].x;
        brch->kids[4]->cl->nd[2].y=brch->kids[1]->cl->nd[0].y;
        brch->kids[4]->cl->nd[2].z=brch->kids[1]->cl->nd[0].z;

        brch->kids[4]->cl->nd[3].x=brch->kids[1]->cl->nd[3].x;
        brch->kids[4]->cl->nd[3].y=brch->kids[1]->cl->nd[3].y;
        brch->kids[4]->cl->nd[3].z=brch->kids[1]->cl->nd[3].z;

        brch->kids[5]->cl->nd[0].x=brch->kids[1]->cl->nd[3].x;
        brch->kids[5]->cl->nd[0].y=brch->kids[1]->cl->nd[3].y;
        brch->kids[5]->cl->nd[0].z=brch->kids[1]->cl->nd[3].z;

        brch->kids[5]->cl->nd[1].x=brch->kids[0]->cl->nd[3].x;
        brch->kids[5]->cl->nd[1].y=brch->kids[0]->cl->nd[3].y;
        brch->kids[5]->cl->nd[1].z=brch->kids[0]->cl->nd[3].z;

        brch->kids[5]->cl->nd[2].x=brch->kids[0]->cl->nd[2].x;
        brch->kids[5]->cl->nd[2].y=brch->kids[0]->cl->nd[2].y;
        brch->kids[5]->cl->nd[2].z=brch->kids[0]->cl->nd[2].z;

        brch->kids[5]->cl->nd[3].x=brch->kids[0]->cl->nd[1].x;
        brch->kids[5]->cl->nd[3].y=brch->kids[0]->cl->nd[1].y;
        brch->kids[5]->cl->nd[3].z=brch->kids[0]->cl->nd[1].z;

        brch->kids[6]->cl->nd[0].x=brch->kids[2]->cl->nd[0].x;
        brch->kids[6]->cl->nd[0].y=brch->kids[2]->cl->nd[0].y;
        brch->kids[6]->cl->nd[0].z=brch->kids[2]->cl->nd[0].z;

        brch->kids[6]->cl->nd[1].x=brch->kids[2]->cl->nd[3].x;
        brch->kids[6]->cl->nd[1].y=brch->kids[2]->cl->nd[3].y;
        brch->kids[6]->cl->nd[1].z=brch->kids[2]->cl->nd[3].z;

        brch->kids[6]->cl->nd[2].x=brch->kids[1]->cl->nd[3].x;
        brch->kids[6]->cl->nd[2].y=brch->kids[1]->cl->nd[3].y;
        brch->kids[6]->cl->nd[2].z=brch->kids[1]->cl->nd[3].z;

        brch->kids[6]->cl->nd[3].x=brch->kids[1]->cl->nd[2].x;
        brch->kids[6]->cl->nd[3].y=brch->kids[1]->cl->nd[2].y;
        brch->kids[6]->cl->nd[3].z=brch->kids[1]->cl->nd[2].z;

        brch->kids[7]->cl->nd[0].x=brch->kids[3]->cl->nd[2].x;
        brch->kids[7]->cl->nd[0].y=brch->kids[3]->cl->nd[2].y;
        brch->kids[7]->cl->nd[0].z=brch->kids[3]->cl->nd[2].z;

        brch->kids[7]->cl->nd[1].x=brch->kids[3]->cl->nd[1].x;
        brch->kids[7]->cl->nd[1].y=brch->kids[3]->cl->nd[1].y;
        brch->kids[7]->cl->nd[1].z=brch->kids[3]->cl->nd[1].z;

        brch->kids[7]->cl->nd[2].x=brch->kids[3]->cl->nd[0].x;
        brch->kids[7]->cl->nd[2].y=brch->kids[3]->cl->nd[0].y;
        brch->kids[7]->cl->nd[2].z=brch->kids[3]->cl->nd[0].z;

        brch->kids[7]->cl->nd[3].x=brch->kids[0]->cl->nd[2].x;
        brch->kids[7]->cl->nd[3].y=brch->kids[0]->cl->nd[2].y;
        brch->kids[7]->cl->nd[3].z=brch->kids[0]->cl->nd[2].z;

      }
      
    } else if (brch->type==2&&brch->nsplit==2) {
          
      brch->kids[0]->cl->nd[0].x=brch->cl->nd[0].x;
      brch->kids[0]->cl->nd[0].y=brch->cl->nd[0].y;
      brch->kids[0]->cl->nd[0].z=brch->cl->nd[0].z;

      brch->kids[0]->cl->nd[1].x=brch->cl->nd[1].x;
      brch->kids[0]->cl->nd[1].y=brch->cl->nd[1].y;
      brch->kids[0]->cl->nd[1].z=brch->cl->nd[1].z;

      brch->kids[0]->cl->nd[2].x=brch->cl->nd[2].x;
      brch->kids[0]->cl->nd[2].y=brch->cl->nd[2].y;
      brch->kids[0]->cl->nd[2].z=brch->cl->nd[2].z;

      brch->kids[0]->cl->nd[3].x=0.5*(brch->cl->nd[0].x+brch->cl->nd[3].x);
      brch->kids[0]->cl->nd[3].y=0.5*(brch->cl->nd[0].y+brch->cl->nd[3].y);
      brch->kids[0]->cl->nd[3].z=0.5*(brch->cl->nd[0].z+brch->cl->nd[3].z);

      brch->kids[0]->cl->nd[4].x=0.5*(brch->cl->nd[1].x+brch->cl->nd[4].x);
      brch->kids[0]->cl->nd[4].y=0.5*(brch->cl->nd[1].y+brch->cl->nd[4].y);
      brch->kids[0]->cl->nd[4].z=0.5*(brch->cl->nd[1].z+brch->cl->nd[4].z);

      brch->kids[0]->cl->nd[5].x=0.5*(brch->cl->nd[2].x+brch->cl->nd[5].x);
      brch->kids[0]->cl->nd[5].y=0.5*(brch->cl->nd[2].y+brch->cl->nd[5].y);
      brch->kids[0]->cl->nd[5].z=0.5*(brch->cl->nd[2].z+brch->cl->nd[5].z);

      brch->kids[1]->cl->nd[3].x=brch->cl->nd[3].x;
      brch->kids[1]->cl->nd[3].y=brch->cl->nd[3].y;
      brch->kids[1]->cl->nd[3].z=brch->cl->nd[3].z;

      brch->kids[1]->cl->nd[4].x=brch->cl->nd[4].x;
      brch->kids[1]->cl->nd[4].y=brch->cl->nd[4].y;
      brch->kids[1]->cl->nd[4].z=brch->cl->nd[4].z;

      brch->kids[1]->cl->nd[5].x=brch->cl->nd[5].x;
      brch->kids[1]->cl->nd[5].y=brch->cl->nd[5].y;
      brch->kids[1]->cl->nd[5].z=brch->cl->nd[5].z;

      brch->kids[1]->cl->nd[0].x=0.5*(brch->cl->nd[0].x+brch->cl->nd[3].x);
      brch->kids[1]->cl->nd[0].y=0.5*(brch->cl->nd[0].y+brch->cl->nd[3].y);
      brch->kids[1]->cl->nd[0].z=0.5*(brch->cl->nd[0].z+brch->cl->nd[3].z);

      brch->kids[1]->cl->nd[1].x=0.5*(brch->cl->nd[1].x+brch->cl->nd[4].x);
      brch->kids[1]->cl->nd[1].y=0.5*(brch->cl->nd[1].y+brch->cl->nd[4].y);
      brch->kids[1]->cl->nd[1].z=0.5*(brch->cl->nd[1].z+brch->cl->nd[4].z);

      brch->kids[1]->cl->nd[2].x=0.5*(brch->cl->nd[2].x+brch->cl->nd[5].x);
      brch->kids[1]->cl->nd[2].y=0.5*(brch->cl->nd[2].y+brch->cl->nd[5].y);
      brch->kids[1]->cl->nd[2].z=0.5*(brch->cl->nd[2].z+brch->cl->nd[5].z);

	  }
	   
    for (ifc=0;ifc<brch->nlfc;ifc++){
      brch->kids[ikid]->cl->fc[ifc].bc = brch->cl->fc[ifc].bc;
    }

    brch->kids[ikid]->part=brch->part;
    leafallocation(run,brch->kids[ikid]);
    brch->kids[ikid]->el->crith=brch->el->crith;
    brch->kids[ikid]->split=brch->split;
     
    run->topo->nleavestmp++;
    run->topo->lleaves[brch->part]++;
    normalvectorcalc(run,brch->kids[ikid]);
    tagvectorcalc(run,brch->kids[ikid]);
    volmproperties(brch->kids[ikid],run);
    
    //dx=brch->kids[ikid]->cl->xc-brch->cl->xc;
    //dy=brch->kids[ikid]->cl->yc-brch->cl->yc;
    //dz=brch->kids[ikid]->cl->zc-brch->cl->zc;

    for(iv=0;iv<NEQ;iv++){ 
      brch->kids[ikid]->el->S[iv]     = brch->el->S[iv];//+dx*brch->el->SG[iv][0]+dy*brch->el->SG[iv][1]+dz*brch->el->SG[iv][2];
      //brch->kids[ikid]->el->SG[iv][0] = brch->el->SG[iv][0];
      //brch->kids[ikid]->el->SG[iv][1] = brch->el->SG[iv][1];
      //brch->kids[ikid]->el->SG[iv][2] = brch->el->SG[iv][2];
    }

    for(iv=0;iv<2;iv++) {
      brch->kids[ikid]->el->SI[0]=brch->el->SI[0];
    }
  }  
  
  leafdeallocation(run,brch);
  brchconn(brch);

}