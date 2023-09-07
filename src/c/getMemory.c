#include "strdata.h"

void getMemory(struct RUN * run,
    int* currRealMem, int* peakRealMem,
    int* currVirtMem, int* peakVirtMem) {
    
    int ierr;
    
    // stores each word in status file
    char buffer[1024] = "";

    // linux file contaiNS this-process info
    FILE* file = fopen("/proc/self/status", "r");

    // read the entire file
    while (fscanf(file, " %1023s", buffer) == 1) {

        if (strcmp(buffer, "VmRSS:") == 0) {
            ierr=fscanf(file, " %d", currRealMem);
        }
        if (strcmp(buffer, "VmHWM:") == 0) {
            ierr=fscanf(file, " %d", peakRealMem);
        }
        if (strcmp(buffer, "VmSize:") == 0) {
            ierr=fscanf(file, " %d", currVirtMem);
        }
        if (strcmp(buffer, "VmPeak:") == 0) {
            ierr=fscanf(file, " %d", peakVirtMem);
        }
    }
    fclose(file);

    double memory;
    memory = ((double) *currVirtMem)*0.001;
 
    printf("Memory usage of proc: %d size: %d :  %f [MB] \n",run->con->rank,run->topo->pleaves,memory);
    MPI_Barrier(MPI_COMM_WORLD);
    if(run->con->rank==0){ printf (" \n"); }

}
