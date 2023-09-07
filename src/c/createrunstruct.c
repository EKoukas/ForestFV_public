#include "strdata.h"

void  createrunstruct (struct RUN * run) {

  run->con   = malloc(sizeof(CON));
  run->topo  = malloc(sizeof(FOREST));
  run->sol_L = malloc(sizeof(SOL));
  run->sol_R = malloc(sizeof(SOL));

}