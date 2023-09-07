#include "strdata.h"

void check_int(int i_temp, int i_bound_low, int i_bound_high, int line) {

  if ((i_temp<i_bound_low)||(i_temp>i_bound_high)) {
    printf("Wrong input file: variable out of bounds \n"); 
    printf("line %d \n ",line); 
    printf("Exiting \n "); 
    exit(0);
  } 

}

void check_db(double temp, double bound_low, double bound_high, int line) {

  if ((temp<bound_low)||(temp>bound_high)) {
    printf("Wrong input file: variable out of bounds \n"); 
    printf("line %d \n",line); 
    printf("Exiting \n "); 
    exit(0);
  } 

}