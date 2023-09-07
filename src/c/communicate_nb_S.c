#include "strdata.h"

void communicate_nb_S (struct RUN * run,int ifield,int i_action) {

  struct BRANCH * crnt;
  int ipt,ipr,iv,itag,iel,ilvl,nlvl,ifc,idr,ing,ib;
  int ibufsize;
  int tag;
  int * adr;

  MPI_Status status;
  MPI_Request request;
  adr=malloc(10*sizeof(int));

	if (i_action==0) {	// Sends data

		for (ipr=0;ipr<run->con->size;ipr++){ // loop proxies to send

			ib=0;
			if 			(ifield==0){ ibufsize = run->topo->nprox[ipr]*(NEQ);       }  // S
			else if (ifield==1){ ibufsize = run->topo->nprox[ipr]*6*(NEQ);     }  // SF
			
			if (ibufsize!=0) {
	
				for (itag=0;itag<run->topo->nprox[ipr];itag++){
					iel=run->topo->proxdrys[ipr][itag];
					tag=run->topo->proxtags[ipr][itag];
					tag2adr(adr,tag,&nlvl);
					crnt=run->topo->drys[iel]->brch;
					ilvl=0;
					while (crnt->nkids!=0){
						crnt=crnt->kids[adr[ilvl]-1];
						ilvl++;
					}
					
					if (ifield==0) {
						
						for (iv=0;iv<NEQ;iv++){
							g_sendbuff_S[ipr][ib]=crnt->el->S[iv];
							ib++;
						}
					
					}	else if (ifield==1) {
						
						for (ifc=0;ifc<6;ifc++){    
							for (iv=0;iv<NEQ;iv++){								
								g_sendbuff_SF[ipr][ib]=crnt->el->SF[ifc][iv];
								ib++;
							}
						}

					}
										
				}

				if (ifield==0){ 
					MPI_Isend(g_sendbuff_S[ipr], ibufsize,MPI_DOUBLE,ipr,(run->con->rank+(2*run->con->size*ifield)),MPI_COMM_WORLD,&request);
					MPI_Request_free(&request);
				} else if (ifield==1){ 
					MPI_Isend(g_sendbuff_SF[ipr],ibufsize,MPI_DOUBLE,ipr,(run->con->rank+(2*run->con->size*ifield)),MPI_COMM_WORLD,&request);
					MPI_Request_free(&request);
				} 

			}
				
		}

	} else if (i_action==1) {

		for (ipt=0;ipt<run->con->size;ipt++) { // loop of proxies to receive

			if ((ipt!=run->con->rank)&&(run->topo->nbuff[ipt]!=0)) {

				if 			(ifield==0){ ibufsize = run->topo->nbuff[ipt]*  (NEQ); }
				else if (ifield==1){ ibufsize = run->topo->nbuff[ipt]*6*(NEQ); }
				
				if (ifield==0){ 
					MPI_Irecv(g_recvbuff_S[ipt],ibufsize,MPI_DOUBLE,ipt,(ipt+(2*run->con->size*ifield)),MPI_COMM_WORLD,&request);
					MPI_Wait(&request, &status);
					MPI_Request_free(&request);
				} else if (ifield==1){ 
					MPI_Irecv(g_recvbuff_SF[ipt],ibufsize,MPI_DOUBLE,ipt,(ipt+(2*run->con->size*ifield)),MPI_COMM_WORLD,&request);
					MPI_Wait(&request, &status);
					MPI_Request_free(&request);
				} 

				ib=0;
				for (itag=0;itag<run->topo->nbuff[ipt];itag++) {
					iel=run->topo->buffdrys[ipt][itag];
					tag=run->topo->bufftags[ipt][itag];
					tag2adr(adr,tag,&nlvl);
					crnt=run->topo->drys[iel]->brch;
					ilvl=0;
					
					while (crnt->nkids!=0){
						crnt=crnt->kids[adr[ilvl]-1];
						ilvl++;
					}       
					
					if (ifield==0){
						
						for (iv=0;iv<(NEQ);iv++){
							crnt->el->S[iv] = g_recvbuff_S[ipt][ib];
							ib++;
						}

					} else if (ifield==1) {
						
						for (ifc=0;ifc<6;ifc++){    
							for (iv=0;iv<NEQ;iv++){
								crnt->el->SF[ifc][iv] = g_recvbuff_SF[ipt][ib];
								ib++;
							}
						}

					}

				}					
				
			}

		}

	} else {
		printf("Wrong send/recv option in communicate_nb_S: %d, should be 0 or 1 for send or receive:EXITING \n",i_action);
		exit(0);
	}

  free(adr);
}