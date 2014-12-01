// functions used to simulate dynamic spectrum  
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
#include "dynSpecSim.h"
#include "T2toolkit.h"

int calACF (acfStruct *acfStructure, double step)
{
	int i,j;
	int ns = acfStructure->ns;
	int nf = acfStructure->nf;
	double acf[ns*nf];

	for (i = 0; i < ns; i++)
	{
		//acfStructure->s[i] = i*step;
		acfStructure->s[i] = (i-ns/2)*step;
	}

  for (i = 0; i < nf; i++)
	{
		//acfStructure->f[i] = i*step;
		acfStructure->f[i] = (i-nf/2)*step;
	}

	long seed;
	double rand;
	seed = TKsetSeed();
	rand = TKgaussDev(&seed);
	printf ("%lf\n",rand);
	int n = 0;
	for (i = 0; i < nf; i++)
	{
		for (j = 0; j < ns; j++)
		{
			//acf[n] = acfStructure->s[j]+acfStructure->f[i];
			//acf[n] = exp(-pow((pow(acfStructure->s[j],2.5)+pow(acfStructure->f[i],1.5)),2.0/3.0));
			//acf[n] = exp(-pow((pow(fabs(acfStructure->s[j]),2.5)+pow(fabs(acfStructure->f[i]),1.5)),2.0/3.0));
			acf[n] = exp(-pow((pow(fabs(acfStructure->s[j]+2.0*rand*0.4*acfStructure->f[i]),2.5)+pow(fabs(acfStructure->f[i]),1.5)),2.0/3.0));
			n++;
		}
	}

	//n = 0;
	//for (i = 0; i < nf; i++)
	//{
	//	for (j = 0; j < ns; j++)
	//	{
	//		printf ("%lf  ", acf[n]);
	//		n++;
	//	}
	//	printf ("\n");
	//}

	//n = 0;
	//for (i = 0; i < (2*nf-2); i++)
	//{
	//	for (j = 0; j < (2*ns-2); j++)
	//	{
	//		if (i >= nf-1 && j >= ns-1)
	//		{
	//			acfStructure->acf2d[n] = acf[ns*(2*nf-2-i)+(2*ns-2-j)];
	//			//acf2d[n] = acf[ns*(i-nf+1)+(j-ns+1)];
	//		}
	//		else if  (i >= nf-1 && j < ns-1)
	//		{
	//			acfStructure->acf2d[n] = acf[ns*(2*nf-2-i)+j];
	//			//acf2d[n] = acf[ns*(i-nf+1)+(ns-j-1)];
	//		}
	//		else if  (i < nf-1 && j < ns-1)
	//		{
	//			acfStructure->acf2d[n] = acf[ns*i+j];
	//			//acf2d[n] = acf[ns*(nf-i-1)+(ns-j-1)];
	//		}
	//		else if  (i < nf-1 && j >= ns-1)
	//		{
	//			acfStructure->acf2d[n] = acf[ns*i+(2*ns-2-j)];
	//			//acf2d[n] = acf[ns*(nf-i-1)+(j-ns+1)];
	//		}
	//		n++;
	//	}
	//}

	n = 0;
	for (i = 0; i < nf; i++)
	{
		for (j = 0; j < ns; j++)
		{
			if (i >= nf/2-1 && j >= ns/2-1)
			{
				acfStructure->acf2d[i*ns+j] = acf[ns*(i-nf/2+1)+(j-ns/2+1)];
			}
			else if  (i >= nf/2-1 && j < ns/2-1)
			{
				acfStructure->acf2d[i*ns+j] = acf[ns*(i-nf/2+1)+(j+ns/2-1)];
			}
			else if  (i < nf/2-1 && j < ns/2-1)
			{
				acfStructure->acf2d[i*ns+j] = acf[ns*(i+nf/2-1)+(j+ns/2-1)];
			}
			else if  (i < nf/2-1 && j >= ns/2-1)
			{
				acfStructure->acf2d[i*ns+j] = acf[ns*(i+nf/2-1)+(j-ns/2+1)];
			}
			n++;
		}
	}

	//n = 0;
	//for (i = 0; i < nf; i++)
	//{
	//	for (j = 0; j < ns; j++)
	//	{
	//		printf ("%.0lf  ", acfStructure->acf2d[n]);
	//		n++;
	//	}
	//	printf ("\n");
	//}

	//n = 0;
	//for (i = 0; i < (2*nf-2); i++)
	//{
	//	for (j = 0; j < (2*ns-2); j++)
	//	{
	//		printf ("%.0lf  ", acfStructure->acf2d[n]);
	//		n++;
	//	}
	//	printf ("\n");
	//}

	return 0;
}


int dft2d (acfStruct *acfStructure, fftw_complex *out)
{
	int n0 = acfStructure->nf;
	int n1 = acfStructure->ns;
	//int n0 = 2*acfStructure->nf-2;
	//int n1 = 2*acfStructure->ns-2;
	double *in;

	in = acfStructure->acf2d;

	fftw_plan p;
	
	p = fftw_plan_dft_r2c_2d (n0, n1, in, out, FFTW_MEASURE);

	fftw_execute(p);

	fftw_destroy_plan(p);
  
	return 0;
}

int idft2d (acfStruct *acfStructure)
{
	int n0 = acfStructure->nf;
	int n1 = acfStructure->ns;
	//int n0 = 2*acfStructure->nf-2;
	//int n1 = 2*acfStructure->ns-2;

	fftw_plan p;
	
	p = fftw_plan_dft_2d (n0, n1, acfStructure->eField, acfStructure->intensity, FFTW_BACKWARD, FFTW_MEASURE);

	fftw_execute(p);

	fftw_destroy_plan(p);
  
	return 0;
}

int power (acfStruct *acfStructure)
{
	int i;
	int nf = acfStructure->nf;
	int ns = acfStructure->ns;
	
	/////////////////////////////////////////////////////////////////////////////////
	//printf ("intializing fft...\n");
	acfStruct acfTest; // initialize the system, don't know why....
	acfTest.ns = acfStructure->ns;
	acfTest.nf = acfStructure->nf;

	allocateMemory (&acfTest);
	for (i = 0; i < nf*ns; i++)
	//for (i = 0; i < (2*nf-2)*(2*ns-2); i++)
	{
		//printf ("acf2d %lf\n",acfStructure->acf2d[i]);
		acfTest.acf2d[i] = acfStructure->acf2d[i];
	}

	fftw_complex *out_t;
	out_t = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nf*(ns/2+1));
	//out_t = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (2*nf-2)*ns);
	dft2d (&acfTest, out_t);

	deallocateMemory (&acfTest);
	//////////////////////////////////////////////////////////////////////////////

  fftw_complex *out;
	
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nf*(ns/2+1));
	//out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (2*nf-2)*ns);
	
	dft2d (acfStructure, out);

	//for (i = 0; i < (2*nf-2)*ns; i++)
	//{
	//	acfStructure->psrt[i] = sqrt(fabs(out[i][0]));
	//	//printf ("%lf\n", acfStructure->psrt[i]);
	//}

	int n, j;
	n = 0;
	for (i = 0; i < nf; i++)
	//for (i = 0; i < 2*nf-2; i++)
	{
		for (j = 0; j < ns; j++)
		//for (j = 0; j < 2*ns-2; j++)
		{
			if (j < ns/2+1)
			{
				acfStructure->psrt[n] = sqrt(fabs(out[i*(ns/2+1)+j][0]));
			}
			else
			{
				acfStructure->psrt[n] = sqrt(fabs(out[i*(ns/2+1)+ns-j][0]));
			}
			//printf ("%lf ", acfStructure->psrt[n]);
			n++;
		}
		//printf ("\n");
	}

	fftw_free(out); 
	fftw_free(out_t); 

	return 0;
}

void allocateMemory (acfStruct *acfStructure)
{
	int ns = acfStructure->ns;
	int nf = acfStructure->nf;
	
	acfStructure->s = (double *)malloc(sizeof(double)*ns);
	acfStructure->f = (double *)malloc(sizeof(double)*nf);
	acfStructure->acf2d = (double *)malloc(sizeof(double)*ns*nf);
	acfStructure->psrt = (double *)malloc(sizeof(double)*ns*nf);
	acfStructure->eField = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nf*ns);
	acfStructure->intensity = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nf*ns);
	//acfStructure->acf2d = (double *)malloc(sizeof(double)*(2*ns-2)*(2*nf-2));
	//acfStructure->psrt = (double *)malloc(sizeof(double)*(2*ns-2)*(2*nf-2));
	//acfStructure->eField = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (2*nf-2)*(2*ns-2));
	//acfStructure->intensity = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (2*nf-2)*(2*ns-2));
	
	acfStructure->dynSpec = (double **)fftw_malloc(sizeof(double *)*nf);
	//acfStructure->dynSpec = (double **)fftw_malloc(sizeof(double *) * (2*nf-2));

	int i;
	for (i = 0; i < nf; i++)
	{
		acfStructure->dynSpec[i] = (double *)fftw_malloc(sizeof(double)*ns);
	}
}

void deallocateMemory (acfStruct *acfStructure)
{
	int ns = acfStructure->ns;
	
	fftw_free(acfStructure->s); 
	fftw_free(acfStructure->f); 
	fftw_free(acfStructure->eField); 
	fftw_free(acfStructure->intensity); 
	free(acfStructure->acf2d);
	free(acfStructure->psrt);

	//fftw_free(acfStructure->dynSpec); 
	int i;
	for (i = 0; i < ns; i++)
	//for (i = 0; i < 2*ns-2; i++)
	{
		free(acfStructure->dynSpec[i]);
	}

	free(acfStructure->dynSpec);
}

int simDynSpec (acfStruct *acfStructure)
{
	int nf = acfStructure->nf;
	int ns = acfStructure->ns;
	long seed;
	
	int i;
	int j;
	seed = TKsetSeed();

	//for (i = 0; i < (2*nf-2)*(2*ns-2); i++)
	for (i = 0; i < nf*ns; i++)
	{
		//acfStructure->eField[i][0] = acfStructure->psrt[i];
		//acfStructure->eField[i][1] = acfStructure->psrt[i];
		acfStructure->eField[i][0] = acfStructure->psrt[i]*TKgaussDev(&seed);
		acfStructure->eField[i][1] = acfStructure->psrt[i]*TKgaussDev(&seed);
		//printf ("%lf\n",TKgaussDev(&seed));
	}

	/////////////////////////////////////////////////////////////////////////////////
	//printf ("intializing ifft...\n");
	acfStruct acfTest; // initialize the system, don't know why....
	acfTest.ns = acfStructure->ns;
	acfTest.nf = acfStructure->nf;

	allocateMemory (&acfTest);
	//for (i = 0; i < (2*nf-2)*(2*ns-2); i++)
	for (i = 0; i < nf*ns; i++)
	{
		acfTest.eField[i][0] = acfStructure->eField[i][0];
		acfTest.eField[i][1] = acfStructure->eField[i][0];
	}

	idft2d (&acfTest);

	deallocateMemory (&acfTest);
	//////////////////////////////////////////////////////////////////////////////

	 FILE *fp;
	 if ((fp = fopen("test.dat", "w+")) == NULL)
	 {
		 fprintf (stdout, "Can't open file\n");
		 exit(1);
	 }

	// ifft
	idft2d (acfStructure);

	int n;
	n = 0;
	for (i = 0; i < nf; i++)
	//for (i = 0; i < (2*nf-2); i++)
	{
		for (j = 0; j < ns; j++)
		//for (j = 0; j < (2*ns-2); j++)
		{
			acfStructure->dynSpec[i][j] = pow(acfStructure->intensity[n][1]/(nf*ns),2.0)+pow(acfStructure->intensity[n][0]/(nf*ns),2.0);
			//acfStructure->dynSpec[i][j] = pow(acfStructure->intensity[n][1]/((2*nf-2)*(2*ns-2)),2.0)+pow(acfStructure->intensity[n][0]/((2*nf-2)*(2*ns-2)),2.0);
			//printf ("%lf  ", acfStructure->intensity[n][0]/((2*nf-2)*(2*ns-2)));
			//printf ("%lf  ", acfStructure->dynSpec[i][j]);
			fprintf (fp, "%lf  ", acfStructure->dynSpec[i][j]);
			n++;
		}
		//printf ("\n");
		fprintf (fp, "\n");
	}

	if (fclose (fp) != 0)
		fprintf (stderr, "Error closing\n");

	return 0;
}
