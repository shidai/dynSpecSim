#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
#include <fftw3.h>
//#include "fitsio.h"

typedef struct acfStruct {
	double phaseGradient;
	double bw; // observing bandwidth
	double f0; // scintillation bandwidth
	double tint; // integration time
	double t0; // scintillation time-scale
	int nchn;
	int nsubint;
	int ns; // sampling number of spatial scale 
	int nf; // sampling number of frequency scale
	double size[2]; // sampling boundary
	double steps;
	double stepf;
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
int calACF (acfStruct *acfStructure);
int power (acfStruct *acfStructure);
void deallocateMemory (acfStruct *acfStructure);
void allocateMemory (acfStruct *acfStructure);
int simDynSpec (acfStruct *acfStructure);
int windowSize (acfStruct *acfStructure, double *size);
int calSize (acfStruct *acfStructure, double *size, double *ratio);
double find_peak_value (int n, double *s);
void preAllocateMemory (acfStruct *acfStructure);

void palett(int TYPE, float CONTRA, float BRIGHT);
void heatMap (float *tab, acfStruct *acfStructure);
