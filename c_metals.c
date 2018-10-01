#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <mpi.h>

#include "allvars.h"
#include "proto.h"

/*! \file c_metals.c 
 *  
 *  This file contains some of the routines for the
 *  chemical model, feedback and decoupling. 
 */

#ifdef SFR_METALS

int *Nlines, cont_sm;
double **logLambda_i;
double **logT_i;
double FeHgas;

#define  Nab 8
#define  delta 0.05
double metal;
double XH, yhelium;

#define Nmasas 10		/* Number of mass intervals on the tables */
#define Nelements 13		/* Number of chemical elements on the tables */

float **nsimfww;
double **yield1, **yield2, **yield3, **yield4, **yield5;

/* This functions reads the cooling tables from 
      Sutherland & Dopita 1993 */

void read_coolrate_table(void)
{
  int i;


  FILE *fpcool;

  Nlines = calloc(Nab, sizeof(int));
  logT_i = calloc(Nab, sizeof(double *));
  logLambda_i = calloc(Nab, sizeof(double *));

  fpcool = fopen("tablas_cooling/mzero.cie", "r");	/*  primordial abundance */
  if(fpcool == NULL)
    {
      printf("can't open Cooling Table for primordial abundance\n");
      return;
    }
  fscanf(fpcool, "%d\n", Nlines + 0);

  *(logT_i + 0) = mymalloc(*(Nlines + 0) * sizeof(double));
  *(logLambda_i + 0) = mymalloc(*(Nlines + 0) * sizeof(double));
  for(i = 0; i < (*(Nlines + 0)); i++)
    fscanf(fpcool, "%lf %*f %*f %*f %*f %lf %*f %*f %*f %*f %*f %*f\n", *(logT_i + 0) + i,
	   *(logLambda_i + 0) + i);


  fpcool = fopen("tablas_cooling/m-30.cie", "r");
  if(fpcool == NULL)
    {
      printf("can't open Cooling Table -3.0 \n");
      return;
    }
  fscanf(fpcool, "%d\n", Nlines + 1);
  *(logT_i + 1) = mymalloc(*(Nlines + 1) * sizeof(double));
  *(logLambda_i + 1) = mymalloc(*(Nlines + 1) * sizeof(double));
  for(i = 0; i < (*(Nlines + 1)); i++)
    fscanf(fpcool, "%lf %*f %*f %*f %*f %lf %*f %*f %*f %*f %*f %*f\n", *(logT_i + 1) + i,
	   *(logLambda_i + 1) + i);

  fpcool = fopen("tablas_cooling/m-20.cie", "r");
  if(fpcool == NULL)
    {
      printf("can't open Cooling Table -2.0 \n");
      return;
    }
  fscanf(fpcool, "%d\n", Nlines + 2);
  *(logT_i + 2) = mymalloc(*(Nlines + 2) * sizeof(double));
  *(logLambda_i + 2) = mymalloc(*(Nlines + 2) * sizeof(double));
  for(i = 0; i < (*(Nlines + 2)); i++)
    fscanf(fpcool, "%lf %*f %*f %*f %*f %lf %*f %*f %*f %*f %*f %*f\n", *(logT_i + 2) + i,
	   *(logLambda_i + 2) + i);


  fpcool = fopen("tablas_cooling/m-15.cie", "r");
  if(fpcool == NULL)
    {
      printf("can't open Cooling Table -1.5\n");
      return;
    }
  fscanf(fpcool, "%d\n", Nlines + 3);
  *(logT_i + 3) = mymalloc(*(Nlines + 3) * sizeof(double));
  *(logLambda_i + 3) = mymalloc(*(Nlines + 3) * sizeof(double));
  for(i = 0; i < (*(Nlines + 3)); i++)
    fscanf(fpcool, "%lf %*f %*f %*f %*f %lf %*f %*f %*f %*f %*f %*f\n", *(logT_i + 3) + i,
	   *(logLambda_i + 3) + i);


  fpcool = fopen("tablas_cooling/m-10.cie", "r");
  if(fpcool == NULL)
    {
      printf("can't open Cooling Table -1.0\n");
      return;
    }
  fscanf(fpcool, "%d\n", Nlines + 4);
  *(logT_i + 4) = mymalloc(*(Nlines + 4) * sizeof(double));
  *(logLambda_i + 4) = mymalloc(*(Nlines + 4) * sizeof(double));
  for(i = 0; i < (*(Nlines + 4)); i++)
    fscanf(fpcool, "%lf %*f %*f %*f %*f %lf %*f %*f %*f %*f %*f %*f\n", *(logT_i + 4) + i,
	   *(logLambda_i + 4) + i);

  fpcool = fopen("tablas_cooling/m-05.cie", "r");
  if(fpcool == NULL)
    {
      printf("can't open Cooling Table -0.5\n");
      return;
    }
  fscanf(fpcool, "%d\n", Nlines + 5);
  *(logT_i + 5) = mymalloc(*(Nlines + 5) * sizeof(double));
  *(logLambda_i + 5) = mymalloc(*(Nlines + 5) * sizeof(double));
  for(i = 0; i < (*(Nlines + 5)); i++)
    fscanf(fpcool, "%lf %*f %*f %*f %*f %lf %*f %*f %*f %*f %*f %*f\n", *(logT_i + 5) + i,
	   *(logLambda_i + 5) + i);

  fpcool = fopen("tablas_cooling/m-00.cie", "r");
  if(fpcool == NULL)
    {
      printf("can't open Cooling Table 0.0\n");
      return;
    }
  fscanf(fpcool, "%d\n", Nlines + 6);
  *(logT_i + 6) = mymalloc(*(Nlines + 6) * sizeof(double));
  *(logLambda_i + 6) = mymalloc(*(Nlines + 6) * sizeof(double));
  for(i = 0; i < (*(Nlines + 6)); i++)
    fscanf(fpcool, "%lf %*f %*f %*f %*f %lf %*f %*f %*f %*f %*f %*f\n", *(logT_i + 6) + i,
	   *(logLambda_i + 6) + i);

  fpcool = fopen("tablas_cooling/m+05.cie", "r");
  if(fpcool == NULL)
    {
      printf("can't open Cooling Table 1.5\n");
      return;
    }
  fscanf(fpcool, "%d\n", Nlines + 7);
  *(logT_i + 7) = mymalloc(*(Nlines + 7) * sizeof(double));
  *(logLambda_i + 7) = mymalloc(*(Nlines + 7) * sizeof(double));
  for(i = 0; i < (*(Nlines + 7)); i++)
    fscanf(fpcool, "%lf %*f %*f %*f %*f %lf %*f %*f %*f %*f %*f %*f\n", *(logT_i + 7) + i,
	   *(logLambda_i + 7) + i);

}



/* This function selects the correct cooling rate for a given temperature
and abundance */

double get_Lambda_SD(double logT, double abund)
{
  double fhi = 0., flow = 0., t = 0.;
  double Lambda = 0., logLambda = 0.;
  int itab = 0, j = 0;


  if(abund <= -3.5)
    itab = 0;
  if(abund > -3.5 && abund <= -2.5)
    itab = 1;
  if(abund > -2.5 && abund <= -1.75)
    itab = 2;
  if(abund > -1.75 && abund <= -1.25)
    itab = 3;
  if(abund > -1.25 && abund <= -0.75)
    itab = 4;
  if(abund > -0.75 && abund <= -0.25)
    itab = 5;
  if(abund > -0.25 && abund <= 0.25)
    itab = 6;
  if(abund > 0.25)
    itab = 7;

  if(logT < *(*(logT_i + itab) + 0))
    Lambda = 0;
  else
    {
      if(logT > *(*(logT_i + itab) + *(Nlines + itab) - 1))
	Lambda = pow(10.0, *(*(logLambda_i + itab) + *(Nlines + itab) - 1));
      else
	{
/* Interpolation to get the correct Cooling Rate */

	  t = (logT - *(*(logT_i + itab) + 0)) / delta;
	  j = (int) t;

	  fhi = t - j;
	  flow = 1 - fhi;

	  logLambda = flow * (*(*(logLambda_i + itab) + j)) + fhi * (*(*(logLambda_i + itab) + j + 1));
	  Lambda = pow(10.0, logLambda);
	}

    }

  return Lambda;
}


/* This function reads the SNII yields tables from Woosley & Weaver 1995 
   Chemical elements: 3He, 12C, 24Mg, 16O, 56Fe, 28Si, H, 14N, 20Ne, 32S, 40Ca, 62Zn, 56N
   As 56N decays into 56Fe in a very short time-scale, we add its contribution to
   56Fe. Finally we take into account 12 elements */

#ifdef SFR_SNII
void read_yield_table(void)
{
  int i, j;

  FILE *fpyield;

  yield1 = calloc(Nelements, sizeof(double *));
  for(i = 0; i < Nelements; i++)
    *(yield1 + i) = calloc(Nmasas, sizeof(double));

  fpyield = fopen("yieldsww/yield1.dat-ww", "r");	/* primordial abundance */
  if(fpyield == NULL)
    {
      printf("can't open Yield Table for primordial abundance \n");
      return;
    }
  for(i = 0; i < Nelements; i++)
    for(j = 0; j < Nmasas; j++)
      {
	fscanf(fpyield, "   %lf\n", (*(yield1 + i) + j));
	if(i == 12)
	  *(*(yield1 + 4) + j) += *(*(yield1 + 12) + j);

	*(*(yield1 + 4) + j) /= 2;	/* Correction of Fe yield */
      }

  yield2 = calloc(Nelements, sizeof(double *));
  for(i = 0; i < Nelements; i++)
    *(yield2 + i) = calloc(Nmasas, sizeof(double));

  fpyield = fopen("yieldsww/yield2.dat-ww", "r");
  if(fpyield == NULL)
    {
      printf("can't open Yield Table 2 \n");
      return;
    }
  for(i = 0; i < Nelements; i++)
    for(j = 0; j < Nmasas; j++)
      {
	fscanf(fpyield, "   %lf\n", (*(yield2 + i) + j));
	if(i == 12)
	  *(*(yield2 + 4) + j) += *(*(yield2 + 12) + j);

	*(*(yield2 + 4) + j) /= 2;
      }


  yield3 = calloc(Nelements, sizeof(double *));
  for(i = 0; i < Nelements; i++)
    *(yield3 + i) = calloc(Nmasas, sizeof(double));

  fpyield = fopen("yieldsww/yield3.dat-ww", "r");
  if(fpyield == NULL)
    {
      printf("can't open Yield Table 3 \n");
      return;
    }
  for(i = 0; i < Nelements; i++)
    for(j = 0; j < Nmasas; j++)
      {
	fscanf(fpyield, "   %lf\n", (*(yield3 + i) + j));
	if(i == 12)
	  *(*(yield3 + 4) + j) += *(*(yield3 + 12) + j);

	*(*(yield3 + 4) + j) /= 2;
      }

  yield4 = calloc(Nelements, sizeof(double *));
  for(i = 0; i < Nelements; i++)
    *(yield4 + i) = calloc(Nmasas, sizeof(double));

  fpyield = fopen("yieldsww/yield4.dat-ww", "r");
  if(fpyield == NULL)
    {
      printf("can't open Yield Table 4 \n");
      return;
    }
  for(i = 0; i < Nelements; i++)
    for(j = 0; j < Nmasas; j++)
      {
	fscanf(fpyield, "   %lf\n", (*(yield4 + i) + j));
	if(i == 12)
	  *(*(yield4 + 4) + j) += *(*(yield4 + 12) + j);

	*(*(yield4 + 4) + j) /= 2;
      }

  yield5 = calloc(Nelements, sizeof(double *));
  for(i = 0; i < Nelements; i++)
    *(yield5 + i) = calloc(Nmasas, sizeof(double));

  fpyield = fopen("yieldsww/yield5.dat-ww", "r");
  if(fpyield == NULL)
    {
      printf("can't open Yield Table 5 \n");
      return;
    }
  for(i = 0; i < Nelements; i++)
    for(j = 0; j < Nmasas; j++)
      {
	fscanf(fpyield, "   %lf\n", (*(yield5 + i) + j));
	if(i == 12)
	  *(*(yield5 + 4) + j) += *(*(yield5 + 12) + j);

	*(*(yield5 + 4) + j) /= 2;
      }

  fclose(fpyield);
}



/* This function selects the correct  SNII yields for a given metallicity and
   stellar mass. */
double SNII_yields(int indice)
{
  int i, j, ik;
  int kmetal = 0;

  double yield[12][10];
  double delta_metals = 0;
  double metals_in_element[12];
  double check = 0.;
  double metal = 0;

#define  Zsol 0.02

# ifdef SFR_FEEDBACK
  double energy_sn = 0;
# endif

  for(i = 0; i < Nelements - 1; i++)	/* Now we have 12 elements */
    for(j = 0; j < Nmasas; j++)
      yield[i][j] = 0;

  metal = 0;
  for(ik = 1; ik < Nelements - 1; ik++)	/* all chemical elements but H & He */
    if(ik != 6)
      metal += P[indice].Zm[ik];

  metal /= (P[indice].Mass * Zsol);	/* metallicity in solar units */

  if(metal < 0)
    {
      printf("ERROR: Negative metallicity %g\n", metal);
      endrun(330);
    }

  for(i = 0; i < Nelements - 1; i++)
    check += P[indice].ZmReservoir[i];
  if(check > 0)
#ifdef SFR_FEEDBACK
    if(Flag_phase == 1)
#endif
      {
	printf("ERROR Part=%d enters SNII_yields with Reservoir=%g\n", P[indice].ID, check);
	endrun(331);
      }

  if(metal < 0.0001)
    {
      kmetal = 0;
      for(i = 0; i < Nelements - 1; i++)
	for(j = 0; j < Nmasas; j++)
	  yield[i][j] = *(*(yield1 + i) + j);
    }

  if(metal >= 0.0001 && metal < 0.01)
    {
      kmetal = 1;
      for(i = 0; i < Nelements - 1; i++)
	for(j = 0; j < Nmasas; j++)
	  yield[i][j] = *(*(yield2 + i) + j);
    }

  if(metal >= 0.01 && metal < 0.1)
    {
      kmetal = 2;
      for(i = 0; i < Nelements - 1; i++)
	for(j = 0; j < Nmasas; j++)
	  yield[i][j] = *(*(yield3 + i) + j);
    }

  if(metal >= 0.1 && metal < 1)
    {
      kmetal = 3;
      for(i = 0; i < Nelements - 1; i++)
	for(j = 0; j < Nmasas; j++)
	  yield[i][j] = *(*(yield4 + i) + j);
    }

  if(metal >= 1)
    {
      kmetal = 4;
      for(i = 0; i < Nelements - 1; i++)
	for(j = 0; j < Nmasas; j++)
	  yield[i][j] = *(*(yield5 + i) + j);
    }

  /* Calculation of the production of metals by the stellar mass considered. */

  for(i = 0; i < Nelements - 1; i++)	/* Note: 56N has already decayed into 56Fe */
    {
      metals_in_element[i] = 0;

      for(j = 0; j < Nmasas; j++)
	{
	  metals_in_element[i] += P[indice].Mass * (*(*(nsimfww + kmetal) + j)) * (yield[i][j]);
#ifdef SFR_FEEDBACK
	  if(i == 0)
	    energy_sn += (*(*(nsimfww + kmetal) + j));
#endif
	}
      P[indice].ZmReservoir[i] = metals_in_element[i];

      delta_metals += metals_in_element[i];	/* total mass in metals to be ejected */
    }

  for(i = 0; i < Nelements - 1; i++)
    P[indice].Zm[i] *= (1 - delta_metals / P[indice].Mass);

  P[indice].Mass -= delta_metals;


#ifdef SFR_FEEDBACK

  if(P[indice].EnergySN != 0)
    {
      printf("ERROR - Star comes with EnergySN = %g\n", P[indice].EnergySN);
      endrun(23487686);
    }

  /* if the stars was transformed from a gas particle that already had a reservoir
   *  of SN energy , then the star keeps the energy and adds the new energy */
  P[indice].EnergySN += energy_sn * P[indice].Mass * ESN / SOLAR_MASS * All.UnitMass_in_g / All.HubbleParam;	/*conversion to internal units. IMF comes in units of 1M_sun */


  /*Energy Test */
  DEnergy_feedback += P[indice].EnergySN;
#endif

  return delta_metals;
}



/* This function construct the Initial Mass Function compatible with yield tables.  */
void imf(void)
{
  int i;

  /* Salpeter IMF */
#define pend 1.35		/* Slope  */
#define cnorm 0.1706		/* Normalization */

  int *Nfila;			/* Number of files of each table */

  Nfila = calloc(5, sizeof(int));
  *(Nfila + 0) = 10;
  *(Nfila + 1) = 10;
  *(Nfila + 2) = 10;
  *(Nfila + 3) = 10;
  *(Nfila + 4) = 12;

  nsimfww = calloc(5, sizeof(float *));
  for(i = 0; i < 5; i++)
    *(nsimfww + i) = calloc(*(Nfila + i), sizeof(float));

  /*  Z=0  in solar units */
  *(*(nsimfww + 0) + 0) = cnorm / (1 - pend) * (pow(11, (1 - pend)) - pow(10, (1 - pend))) / 10.5;
  *(*(nsimfww + 0) + 1) = cnorm / (1 - pend) * (pow(13, (1 - pend)) - pow(11, (1 - pend))) / 12.0;
  *(*(nsimfww + 0) + 2) = cnorm / (1 - pend) * (pow(14, (1 - pend)) - pow(13, (1 - pend))) / 13.5;
  *(*(nsimfww + 0) + 3) = cnorm / (1 - pend) * (pow(16, (1 - pend)) - pow(14, (1 - pend))) / 15.0;
  *(*(nsimfww + 0) + 4) = cnorm / (1 - pend) * (pow(18, (1 - pend)) - pow(16, (1 - pend))) / 17.0;
  *(*(nsimfww + 0) + 5) = cnorm / (1 - pend) * (pow(23, (1 - pend)) - pow(18, (1 - pend))) / 20.5;
  *(*(nsimfww + 0) + 6) = cnorm / (1 - pend) * (pow(24, (1 - pend)) - pow(23, (1 - pend))) / 23.5;
  *(*(nsimfww + 0) + 7) = cnorm / (1 - pend) * (pow(28, (1 - pend)) - pow(24, (1 - pend))) / 26.0;
  *(*(nsimfww + 0) + 8) = cnorm / (1 - pend) * (pow(38, (1 - pend)) - pow(28, (1 - pend))) / 32.0;
  *(*(nsimfww + 0) + 9) = cnorm / (1 - pend) * (pow(40, (1 - pend)) - pow(38, (1 - pend))) / 39.0;


  /*  Z=0.0001 in solar units  */
  *(*(nsimfww + 1) + 0) = cnorm / (1 - pend) * (pow(11, (1 - pend)) - pow(10, (1 - pend))) / 10.5;
  *(*(nsimfww + 1) + 1) = cnorm / (1 - pend) * (pow(13, (1 - pend)) - pow(11, (1 - pend))) / 12.0;
  *(*(nsimfww + 1) + 2) = cnorm / (1 - pend) * (pow(16, (1 - pend)) - pow(13, (1 - pend))) / 14.5;
  *(*(nsimfww + 1) + 3) = cnorm / (1 - pend) * (pow(18, (1 - pend)) - pow(16, (1 - pend))) / 17.0;
  *(*(nsimfww + 1) + 4) = cnorm / (1 - pend) * (pow(20, (1 - pend)) - pow(18, (1 - pend))) / 19.0;
  *(*(nsimfww + 1) + 5) = cnorm / (1 - pend) * (pow(23, (1 - pend)) - pow(20, (1 - pend))) / 21.5;
  *(*(nsimfww + 1) + 6) = cnorm / (1 - pend) * (pow(26, (1 - pend)) - pow(23, (1 - pend))) / 24.5;
  *(*(nsimfww + 1) + 7) = cnorm / (1 - pend) * (pow(36, (1 - pend)) - pow(26, (1 - pend))) / 31.0;
  *(*(nsimfww + 1) + 8) = cnorm / (1 - pend) * (pow(38, (1 - pend)) - pow(36, (1 - pend))) / 37.0;
  *(*(nsimfww + 1) + 9) = cnorm / (1 - pend) * (pow(40, (1 - pend)) - pow(38, (1 - pend))) / 39.0;


  /*  Z=0.01 in solar units */
  *(*(nsimfww + 2) + 0) = cnorm / (1 - pend) * (pow(11, (1 - pend)) - pow(10, (1 - pend))) / 10.5;
  *(*(nsimfww + 2) + 1) = cnorm / (1 - pend) * (pow(13, (1 - pend)) - pow(11, (1 - pend))) / 12.0;
  *(*(nsimfww + 2) + 2) = cnorm / (1 - pend) * (pow(16, (1 - pend)) - pow(13, (1 - pend))) / 14.5;
  *(*(nsimfww + 2) + 3) = cnorm / (1 - pend) * (pow(18, (1 - pend)) - pow(16, (1 - pend))) / 17.0;
  *(*(nsimfww + 2) + 4) = cnorm / (1 - pend) * (pow(20, (1 - pend)) - pow(18, (1 - pend))) / 19.0;
  *(*(nsimfww + 2) + 5) = cnorm / (1 - pend) * (pow(23, (1 - pend)) - pow(20, (1 - pend))) / 21.5;
  *(*(nsimfww + 2) + 6) = cnorm / (1 - pend) * (pow(31, (1 - pend)) - pow(23, (1 - pend))) / 27.0;
  *(*(nsimfww + 2) + 7) = cnorm / (1 - pend) * (pow(35, (1 - pend)) - pow(31, (1 - pend))) / 33.0;
  *(*(nsimfww + 2) + 8) = cnorm / (1 - pend) * (pow(38, (1 - pend)) - pow(35, (1 - pend))) / 36.5;
  *(*(nsimfww + 2) + 9) = cnorm / (1 - pend) * (pow(40, (1 - pend)) - pow(38, (1 - pend))) / 39.0;


  /*  Z=0.1 in solar units */
  *(*(nsimfww + 3) + 0) = cnorm / (1 - pend) * (pow(11, (1 - pend)) - pow(10, (1 - pend))) / 10.5;
  *(*(nsimfww + 3) + 1) = cnorm / (1 - pend) * (pow(13, (1 - pend)) - pow(11, (1 - pend))) / 12.0;
  *(*(nsimfww + 3) + 2) = cnorm / (1 - pend) * (pow(16, (1 - pend)) - pow(13, (1 - pend))) / 14.5;
  *(*(nsimfww + 3) + 3) = cnorm / (1 - pend) * (pow(18, (1 - pend)) - pow(16, (1 - pend))) / 17.0;
  *(*(nsimfww + 3) + 4) = cnorm / (1 - pend) * (pow(20, (1 - pend)) - pow(18, (1 - pend))) / 19.0;
  *(*(nsimfww + 3) + 5) = cnorm / (1 - pend) * (pow(23, (1 - pend)) - pow(20, (1 - pend))) / 21.5;
  *(*(nsimfww + 3) + 6) = cnorm / (1 - pend) * (pow(31, (1 - pend)) - pow(23, (1 - pend))) / 27.0;
  *(*(nsimfww + 3) + 7) = cnorm / (1 - pend) * (pow(35, (1 - pend)) - pow(31, (1 - pend))) / 33.0;
  *(*(nsimfww + 3) + 8) = cnorm / (1 - pend) * (pow(38, (1 - pend)) - pow(35, (1 - pend))) / 36.5;
  *(*(nsimfww + 3) + 9) = cnorm / (1 - pend) * (pow(40, (1 - pend)) - pow(38, (1 - pend))) / 39.0;


  /*  Z=1 in solar units */
  *(*(nsimfww + 4) + 0) = cnorm / (1 - pend) * (pow(10, (1 - pend)) - pow(9, (1 - pend))) / 9.5;
  *(*(nsimfww + 4) + 1) = cnorm / (1 - pend) * (pow(11, (1 - pend)) - pow(10, (1 - pend))) / 10.5;
  *(*(nsimfww + 4) + 2) = cnorm / (1 - pend) * (pow(13, (1 - pend)) - pow(11, (1 - pend))) / 12.0;
  *(*(nsimfww + 4) + 3) = cnorm / (1 - pend) * (pow(16, (1 - pend)) - pow(13, (1 - pend))) / 14.5;
  *(*(nsimfww + 4) + 4) = cnorm / (1 - pend) * (pow(17, (1 - pend)) - pow(16, (1 - pend))) / 16.5;
  *(*(nsimfww + 4) + 5) = cnorm / (1 - pend) * (pow(18, (1 - pend)) - pow(17, (1 - pend))) / 17.5;
  *(*(nsimfww + 4) + 6) = cnorm / (1 - pend) * (pow(20, (1 - pend)) - pow(18, (1 - pend))) / 19.0;
  *(*(nsimfww + 4) + 7) = cnorm / (1 - pend) * (pow(23, (1 - pend)) - pow(20, (1 - pend))) / 21.5;
  *(*(nsimfww + 4) + 8) = cnorm / (1 - pend) * (pow(27, (1 - pend)) - pow(23, (1 - pend))) / 25.0;
  *(*(nsimfww + 4) + 9) = cnorm / (1 - pend) * (pow(32, (1 - pend)) - pow(27, (1 - pend))) / 29.5;
  *(*(nsimfww + 4) + 10) = cnorm / (1 - pend) * (pow(36, (1 - pend)) - pow(32, (1 - pend))) / 34.0;
  *(*(nsimfww + 4) + 11) = cnorm / (1 - pend) * (pow(40, (1 - pend)) - pow(36, (1 - pend))) / 38.0;

}
#endif /* end of  SFR_SNII */



#endif /* end of SFR_METALS */
