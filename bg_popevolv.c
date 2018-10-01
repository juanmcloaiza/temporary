#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "allvars.h"
#include "proto.h"

#ifdef BG_SFR


/*
 * This routine does the stellar evolution, and yields the metal release 
 */
void bg_stellar_evolution(double age_of_star, double dtime, double mass, FLOAT * metals_star,
			  FLOAT * metals_released, FLOAT * metals_change)
{
  int k;
  double tau = 500.0;		/* release the metals over a fiducial timescale */
  double mz, frac;

  /* for the moment, we have a dummy function here */

  /* convert the times to Myrs */

  age_of_star *= All.UnitTime_in_s / All.HubbleParam;	/* now in seconds */
  age_of_star /= SEC_PER_MEGAYEAR;	/* now in megayears */

  dtime *= All.UnitTime_in_s / All.HubbleParam;	/* now in seconds */
  dtime /= SEC_PER_MEGAYEAR;	/* now in megayears */

  /* compute the amount of metals released by the stellar population */

  mz = METAL_YIELD * mass * (exp(-age_of_star / tau) - exp(-(age_of_star + dtime) / tau));
  mz /= (1 - METAL_YIELD * (1 - exp(-age_of_star / tau)));

  for(k = 0; k < BG_NELEMENTS; k++)
    metals_released[k] = 0;

  metals_released[2] = mz;

  frac = mz / mass;

  for(k = 0; k < BG_NELEMENTS; k++)
    metals_change[k] = -frac * metals_star[k];
}



/* returns a global metallicity for the particle */
double bg_get_metallicity(int j)
{
  double met;
  int k;

  for(k = 2, met = 0; k < BG_NELEMENTS; k++)
    met += P[j].Metals[k];

  return met / P[j].Mass;
}

/* This function returns the elapsed time between a0,a1 */

double bg_get_elapsed_time(double a0, double a1)
{
#define WORKSIZE 100000
  gsl_function F;
  gsl_integration_workspace *workspace;
  double result, abserr, age;

  if(All.ComovingIntegrationOn)
    {
      F.function = &bg_time_integ;
      workspace = gsl_integration_workspace_alloc(WORKSIZE);
      gsl_integration_qag(&F, a0, a1, 1 / All.Hubble,	/* note: absolute error just a dummy */
			  1.0e-7, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
      gsl_integration_workspace_free(workspace);


      age = result * All.UnitTime_in_s / All.HubbleParam;	/* now in seconds */
      age /= SEC_PER_MEGAYEAR;	/* now in megayears */
    }
  else
    {
      age = (a1 - a0) * All.UnitTime_in_s / All.HubbleParam;	/* now in seconds */
      age /= SEC_PER_MEGAYEAR;	/* now in megayears */
    }

  return age;
}


double bg_time_integ(double a, void *param)
{
  double hubble_a;

  hubble_a = hubble_function(All.Time);

  return 1 / (a * hubble_a);
}

#endif
