#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cpgplot.h"
#include "dynSpecSim.h"

//#define MAX(a,b) ( (a)>(b) ? (a) : (b))
//#define MIN(a,b) ( (a)<(b) ? (a) : (b))


void heatMap (float *tab, int dimx, int dimy)
{
  //int i,j;                     
  //int dimx = acfStructure.ns;
	//int dimy = acfStructure.nf; // dimensions 
  //float tab[dimx*dimy];       // value
  float xmin,xmax;            /* intervalle. en x */
  float ymin,ymax;            /* intervalle. en y */
  float zmin,zmax;            /* min et max des valeurs de la fonction */
  float tr[6];                /* matrice utilisee par pgimag */

  xmin=-dimx/2; xmax=dimx/2;
  ymin=-dimy/2; ymax=dimy/2;
  zmin=0; zmax=1;

  /* Matrice de passage pixels -> coordonnees user
     --------------------------------------------- */
  //tr[0]=0;
  //tr[1]=1;
  //tr[2]=0;
  //tr[3]=0;
  //tr[4]=0;
  //tr[5]=1;
  tr[0]=xmin;
  tr[1]=(xmax-xmin)/dimx;
  tr[2]=0;
  tr[3]=ymin;
  tr[4]=0;
  tr[5]=(ymax-ymin)/dimy;
 
	// plot 
  cpgbeg(0,"/xs",1,1);
  //cpgsch(1.5);
  //cpgenv(0,dimy,0,dimy,1,0);        
  cpgenv(xmin,xmax,ymin,ymax,1,0);        
	cpgsvp(0.1,0.9,0.1,0.9);
	//cpgswin(0.05,0.95,0.05,0.95);
  cpglab("Subintegration","Frequency (MHz)","Dynamic spectrum");
	palett(3, -1.0, 0.5);
  cpgimag(tab,dimx,dimy,1,dimx,1,dimy,zmin,zmax,tr);

  cpgend();
} 

void palett(int TYPE, float CONTRA, float BRIGHT)
{
//-----------------------------------------------------------------------
// Set a "palette" of colors in the range of color indices used by
// PGIMAG.
//-----------------------------------------------------------------------
	float GL[] = {0.0, 1.0};
	float GR[] = {0.0, 1.0};
	float GG[] = {0.0, 1.0};
	float GB[] = {0.0, 1.0};
	float RL[] = {-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7};
	float RR[] = { 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0};
	float RG[] = { 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0};
	float RB[] = { 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0};
	float HL[] = {0.0, 0.2, 0.4, 0.6, 1.0};
	float HR[] = {0.0, 0.5, 1.0, 1.0, 1.0};
	float HG[] = {0.0, 0.0, 0.5, 1.0, 1.0};
	float HB[] = {0.0, 0.0, 0.0, 0.3, 1.0};
	float WL[] = {0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0};
	float WR[] = {0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0};
	float WG[] = {0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0};
	float WB[] = {0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0};
	float AL[] = {0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5, 0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0};
	float AR[] = {0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	float AG[] = {0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8, 0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0};
	float AB[] = {0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
      
	if (TYPE == 1)
	{   
		//-- gray scale
		cpgctab(GL, GR, GG, GB, 2, CONTRA, BRIGHT);
	}
	else if (TYPE == 2) 
	{
		//-- rainbow
		cpgctab(RL, RR, RG, RB, 9, CONTRA, BRIGHT);
	}
	else if (TYPE == 3) 
	{
		//-- heat
		cpgctab(HL, HR, HG, HB, 5, CONTRA, BRIGHT);
	}
	else if (TYPE == 4) 
	{
		//-- weird IRAF
		cpgctab(WL, WR, WG, WB, 10, CONTRA, BRIGHT);
	}
	else if (TYPE == 5) 
	{
		//-- AIPS
		cpgctab(AL, AR, AG, AB, 20, CONTRA, BRIGHT);
	}
}

