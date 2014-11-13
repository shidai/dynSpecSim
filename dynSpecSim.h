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
	double *acf2d;  // ACF 
	double *psrt;  // power spactrum
	fftw_complex *eField; // complex electric field
	fftw_complex *intensity;  // intensity 
	double **dynSpec; // dynamic spectrum 
} acfStruct;

int idft2d (acfStruct *acfStructure);
int dft2d (acfStruct *acfStructure, fftw_complex *out);
int calACF (acfStruct *acfStructure, double step);
int power (acfStruct *acfStructure);
void deallocateMemory (acfStruct *acfStructure);
void allocateMemory (acfStruct *acfStructure);
int simDynSpec (acfStruct *acfStructure);
