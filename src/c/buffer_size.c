#include "strdata.h"

void buffer_size_malloc(struct RUN * run) {

  int ibufsize,ipr,ipt;

  g_sendbuff_S  = malloc(run->con->size*sizeof(double *));
  g_sendbuff_SF = malloc(run->con->size*sizeof(double *));
  
  for (ipr=0;ipr<run->con->size;ipr++){
    ibufsize = run->topo->nprox[ipr];   
    
    g_sendbuff_S[ipr]  = malloc(max(1,ibufsize*(NEQ))   *sizeof(double));
    g_sendbuff_SF[ipr] = malloc(max(1,ibufsize*(NEQ*18)) *sizeof(double));
    
  }

  g_recvbuff_S  = malloc(run->con->size*sizeof(double *));
  g_recvbuff_SF = malloc(run->con->size*sizeof(double *));

  for (ipt=0;ipt<run->con->size;ipt++){
    ibufsize = run->topo->nbuff[ipt];   
    
    g_recvbuff_S[ipt]  = malloc(max(1,ibufsize*(NEQ))   *sizeof(double));
    g_recvbuff_SF[ipt] = malloc(max(1,ibufsize*(NEQ*18)) *sizeof(double));
    
  }

}


void buffer_size_free(struct RUN * run) {

  int ibufsize,ipr,ipt;

  for (ipr=0;ipr<run->con->size;ipr++) {
    ibufsize = run->topo->nprox[ipr]; 
    
    free(g_sendbuff_S[ipr]);
    free(g_sendbuff_SF[ipr]);
    
  }
  free(g_sendbuff_S);
  free(g_sendbuff_SF);


  for (ipt=0;ipt<run->con->size;ipt++){
    ibufsize = run->topo->nbuff[ipt]; 
    
    free(g_recvbuff_S[ipt]);
    free(g_recvbuff_SF[ipt]);
    
  }
  free(g_recvbuff_S);
  free(g_recvbuff_SF);
  
}