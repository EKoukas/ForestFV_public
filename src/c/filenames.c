#include"strdata.h"

void filenames(struct RUN * run) {

	run->con->cnamelength=strlen(run->con->casename);

	run->con->filend  = (char *)malloc(sizeof(run->con->cnamelength)+4);
	run->con->fileel  = (char *)malloc(sizeof(run->con->cnamelength)+4);
	run->con->filecon = (char *)malloc(sizeof(run->con->cnamelength)+4);
	run->con->filebr1 = (char *)malloc(sizeof(run->con->cnamelength)+4);
	run->con->filebr2 = (char *)malloc(sizeof(run->con->cnamelength)+4);
	run->con->filebr3 = (char *)malloc(sizeof(run->con->cnamelength)+4);
	run->con->filebr4 = (char *)malloc(sizeof(run->con->cnamelength)+4);

	strcpy(run->con->filend, run->con->casename); // Copy casename to filend
	strcpy(run->con->fileel, run->con->casename);
	strcpy(run->con->filecon,run->con->casename);
	strcpy(run->con->filebr1,run->con->casename);
	strcpy(run->con->filebr2,run->con->casename);
	strcpy(run->con->filebr3,run->con->casename);
	strcpy(run->con->filebr4,run->con->casename);

	strcat(run->con->filend, ".nd");	// Add extention to filend
	strcat(run->con->fileel, ".el");	// Add extention to fileel
	strcat(run->con->filecon,".ct");
	strcat(run->con->filebr1,".ba");
	strcat(run->con->filebr2,".bb");
	strcat(run->con->filebr3,".bc");
	strcat(run->con->filebr4,".bd");

}