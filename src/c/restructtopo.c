#include "strdata.h"

void restructtopo (RUN * run) {
    
  //repartition(run);  // rearrange branches on processors, rearrange keen buffers
  repartition_parmetis(run);  // rearrange branches on processors, rearrange keen buffers

  rebuildlist(run);  // restitch local liked list
  
  rebuildprox(run);  // rearrange neig buffers
  
  rebuildconn(run);  // populate buffer
  
}