#include "strdata.h"

void init_domain (struct RUN * run) {

  switch(MODEL){
    case 0: m0_init_domain(run); break;
    case 1: m1_init_domain(run); break;
    case 2: m2_init_domain(run); break;
    case 3: m3_init_domain(run); break;
    case 4: m4_init_domain(run); break;
    case 5: m5_init_domain(run); break;
  }

}