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
	/*
	int index, n;
	for (i=0;i<argc;i++)
	{
		if (strcmp(argv[i],"-f") == 0)
		{
            index = i + 1;
			n = 0;
			while ( (index + n) < argc && strcmp(argv[index+n],"-std") != 0 && strcmp(argv[index+n],"-pt") != 0 && strcmp(argv[index+n],"-o") != 0 && strcmp(argv[index+n],"-single") != 0)
			{
				n++;
		    }
			//strcpy(fname,argv[++i]);
		}
		else if (strcmp(argv[i],"-std")==0)
		{
			strcpy(tname,argv[++i]);
			mode = 0; // standard template format
			printf ("standard template format\n");
			//sscanf(argv[++i],"%d",&nbin);
		}
		else if (strcmp(argv[i],"-pt")==0)
		{
			strcpy(tname,argv[++i]);
			mode = 1; // ptime template
			printf ("ptime template format\n");
		}
		else if (strcmp(argv[i],"-o")==0)
		{
			strcpy(oname,argv[++i]);
		}
		else if (strcmp(argv[i],"-single")==0)
		{
			tmode = 0; // do freq-dependent matching, and get one TOA
		}
		else if (strcmp(argv[i],"-multi")==0)
		{
			tmode = 1; // do freq-dependent matching, and get TOA for each channel
		}
		else if (strcmp(argv[i],"-fitDM")==0)
		{
			fitDM = 1; // do freq-dependent matching, and get TOA for each channel
		}
	}
	*/

	double step = atof(argv[1]);  // sampling
	int ns = 2*(int)((3+step)/step);
	int nf = 2*(int)((6+step)/step);
	acfStruct acfStructure;

	acfStructure.ns = ns;
	acfStructure.nf = nf;
	allocateMemory (&acfStructure);

	calACF (&acfStructure, step);
	power (&acfStructure);
	simDynSpec (&acfStructure);

	deallocateMemory (&acfStructure);

	return 0;
}
