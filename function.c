// functions used to simulate dynamic spectrum  
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
#include "dynSpecSim.h"

int calACF (acfStruct *acfStructure)
{
	int i,j;
	int ns = acfStructure->ns;
	int nf = acfStructure->nf;
	double acf[ns*nf];
	double acf2d[(2*ns-2)*(2*nf-2)];

	int n = 0;
	for (i = 0; i < nf; i++)
	{
		for (j = 0; j < ns; j++)
		{
			acf[n] = acfStructure->s[j]+acfStructure->f[i];
			//acf[n] = exp(-pow((pow(acfStructure->s[j],2.5)+pow(acfStructure->f[i],1.5)),2.0/3.0));
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

	n = 0;
	for (i = 0; i < (2*nf-2); i++)
	{
		for (j = 0; j < (2*ns-2); j++)
		{
			if (i >= nf-1 && j >= ns-1)
			{
				acf2d[n] = acf[ns*(i-nf+1)+(j-ns+1)];
			}
			else if  (i >= nf-1 && j < ns-1)
			{
				acf2d[n] = acf[ns*(i-nf+1)+(ns-j-1)];
			}
			else if  (i < nf-1 && j < ns-1)
			{
				acf2d[n] = acf[ns*(nf-i-1)+(ns-j-1)];
			}
			else if  (i < nf-1 && j >= ns-1)
			{
				acf2d[n] = acf[ns*(nf-i-1)+(j-ns+1)];
			}
			n++;
		}
	}

	n = 0;
	for (i = 0; i < (2*nf-2); i++)
	{
		for (j = 0; j < (2*ns-2); j++)
		{
			printf ("%.0lf  ", acf2d[n]);
			n++;
		}
		printf ("\n");
	}

	acfStructure->acf = acf;
	return 0;
}


/*
int dft_profiles (int N, double *in, fftw_complex *out)
// dft of profiles
{
	//  dft of profiles 
	///////////////////////////////////////////////////////////////////////
	
	//printf ("%lf\n", in[0]);
	//double *in;
	//fftw_complex *out;
	fftw_plan p;
	
	//in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	//out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_MEASURE);

	fftw_execute(p);

	fftw_destroy_plan(p);
	//fftw_free(in); 
	//fftw_free(out);
  
	return 0;
}

int preA7 (double *s, double *p, int nphase, int nchn, params *param)
//int preA7 (int *k, double amp_s[][NP], double amp_p[][NP], double phi_s[][NP], double phi_p[][NP], double *s, double *p, int nphase, int nchn)
// preparation for calculating A7 of Talyor 1992  
{
	// nphase is the dimention of one profile, nchn is number of profiles
	// k is the dimention of amp of one profile 
	param->nchn = nchn;
	int i,j;
	
	/////////////////////////////////////////////////////////////////////////////////
	double test[nphase];  // initialize the system, don't know why....

	for (i=0;i<nphase;i++)
	{
		test[i]=s[i];
	}
	fftw_complex *out_t;
	out_t = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nphase);
	dft_profiles(nphase,test,out_t);
	//////////////////////////////////////////////////////////////////////////////

    fftw_complex *out_s;
	fftw_complex *out_p;
	
	out_s = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nphase);
	out_p = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nphase);
	
	double s_temp[nphase];  // store one template and profile
	double p_temp[nphase];  

	int n;
	double r_s[nphase/2],im_s[nphase/2];
	double r_p[nphase/2],im_p[nphase/2];
	for (i = 0; i < nchn; i++)
	{
	    for (j=0;j<nphase;j++)
	    {
		    s_temp[j]=s[i*nphase + j];
		    p_temp[j]=p[i*nphase + j];
	    }

	    dft_profiles(nphase,s_temp,out_s);
	    //printf ("%lf %lf\n", out_s[1][0], out_s[1][1]);

	    dft_profiles(nphase,p_temp,out_p);

	    //double amp_s[N/2],phi_s[N/2];
	    //double amp_p[N/2],phi_p[N/2];

		n = 0;
	    for (j = 0; j <= nphase/2-1; j++)
	    {
		    r_s[j]=out_s[j+1][0];
		    im_s[j]=out_s[j+1][1];
		    r_p[j]=out_p[j+1][0];
		    im_p[j]=out_p[j+1][1];
		    //printf ("%lf %lf\n", r_p[i], im_p[i]);
		    //printf ("%lf %lf\n", out_s[i][0], out_s[i][1]);
		    n++;
	    }
	    //printf ("%d\n", n);
	    //printf ("%d %d\n", nphase, nchn);

	    for (j = 0; j < n; j++)
	    {
		    param->a_s[i][j]=sqrt(r_s[j]*r_s[j]+im_s[j]*im_s[j]);
		    param->a_p[i][j]=sqrt(r_p[j]*r_p[j]+im_p[j]*im_p[j]);
		    param->p_s[i][j]=atan2(im_s[j],r_s[j]);
		    param->p_p[i][j]=atan2(im_p[j],r_p[j]);
		    //printf ("%lf %lf %lf\n", r_s[i], im_s[i], amp_s[i]);
		    //printf ("%lf %lf %lf\n", r_p[i], im_p[i], amp_p[i]);
		    //printf ("%lf\n", amp_s[i]);
		    //printf ("%lf\n", amp_p[i]);
	    }
	}
	//(*k)=n;
	param->num = n;

	fftw_free(out_s); 
	fftw_free(out_p); 
	fftw_free(out_t); 

	return 0;
}
*/

