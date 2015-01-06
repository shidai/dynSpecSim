//dynamic spectrum simulator
//first acf defined by spatial scale s0 and bandwidth f0
//the time scale t0 is given by s0 = V*t0 where V is the scintillation
//velocity. This is the velocity that the line of sight moves through the
//scattering region.
//
//
//the spatial behavior is acf(s) = exp(-(s/s0)^(5/3)) and the frequency
//behavior is acf(f) = exp(-f/f0) we have modeled the combined behavior as
//acf(s,f) = exp(-((s/s0)^(5/2) + (f/f0)^(3/2))^(2/3)) which gives a
//reasonably good approximation to an exact calculation for isotropic
//scattering.
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "dynSpecSim.h"
#include "T2toolkit.h"

int main (int argc, char *argv[])
{
	int i;
	double bw, f0;
	double tint, t0;
	for (i = 0; i < argc; i++)
	{
		if (strcmp(argv[i],"-f") == 0)
		{
			bw = atof(argv[++i]);
			f0 = atof(argv[++i]);
			//printf ("Observing bandwidth: %lf (MHz)\n", bw);
			//printf ("Scintillation bandwidth: %lf (MHz)\n", f0);
			//strcpy(fname,argv[++i]);
		}
		else if (strcmp(argv[i],"-t")==0)
		{
			tint = atof(argv[++i]);
			t0 = atof(argv[++i]);
			//printf ("Integration time: %lf (s)\n", tint);
			//printf ("Scintillation time-scale: %lf (s)\n", t0);
		}
	}

	acfStruct acfStructure;

	acfStructure.bw = bw;
	acfStructure.f0 = f0;
	acfStructure.tint = tint;
	acfStructure.t0 = t0;

	acfStructure.nchn = 1024;
	acfStructure.nsubint = 64;

	allocateMemory (&acfStructure);

	calACF (&acfStructure);
	power (&acfStructure);
	simDynSpec (&acfStructure);

	deallocateMemory (&acfStructure);

	return 0;
}
