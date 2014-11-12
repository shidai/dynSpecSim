#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
#include <fftw3.h>
//#include "fitsio.h"

typedef struct acfStruct {
	int ns; // sampling number of spatial scale 
	int nf; // sampling number of frequency scale
	double *s; // spatial scale
	double *f; // bw scale
	double *acf; 
	double *acf2d; 
} acfStruct;

int calACF (acfStruct *acfStructure);
