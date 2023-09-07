#include "strdata.h"

void least_squares_grad (struct RUN* run) {

  int iv,ing,bc_temp,f;
  double temp_WG[3];
  struct BRANCH * crnt;

  // Check mesh change
  if (g_mesh_change_grad_scheme==1) {

    crnt=run->topo->locl; 
    while (crnt!=NULL){

      // d matrix with distances

      
      // w matrix with weights


      // Compute G and G-1 matrix


      //


      crnt=crnt->lnxt;
    }
    g_mesh_change_grad_scheme=0;

  } 


  // Matrix multiplication (G-1)(dT)(WT)(W)(Tn-Tp)
  crnt=run->topo->locl; 
  while (crnt!=NULL){

    m0_facevalues(run,crnt,ing,0,0,0);                // crnt values

    for(f=0;f<crnt->nlfc;f++){                        // Loop faces of crnt cell

      if (crnt->cl->fc[f].bc==0) {                    // Check if boudary face

        for (ing=0;ing<crnt->lfc_neigs[f];ing++) {  // Loop neigs of face (1/2/4)
          
          bc_temp = crnt->cl->fc[f].bc;
          
          m0_facevalues(run,crnt,ing,f,bc_temp,1);



          
          for(iv=0;iv<NEQ;iv++){
            crnt->el->WG[iv][0] = temp_WG[0];
            crnt->el->WG[iv][1] = temp_WG[1];
            crnt->el->WG[iv][2] = temp_WG[2];
          }

        }
      }


    }


    crnt=crnt->lnxt;
  }



}