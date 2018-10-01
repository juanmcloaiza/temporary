#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"

#ifdef BG_SFR

#include "bg_cooling.h"

/*
 * These routines does metal-enrichment and wind formation for
 * the STELLA project on BlueGene.
 */




static FLOAT MetalsReleased[BG_NELEMENTS];
static FLOAT StarMass;


/* Main driver routine that deals with neighbour finding and parallelization 
   (export to other CPUs if needed) */
void bg_enrich(void)
{
  int *noffset, *nbuffer, *nsend, *nsend_local;
  int i, j, k, n, ndone, ndonetot, ntot, ntotleft, npleft;
  int maxfill, numstars, level, ngrp, sendTask, recvTask, nexport;
  double hubble_a, time_hubble_a, a3inv, mass;
  double time_begstep, dt, dtime, age_of_star;
  FLOAT MetalsChangeInStar[BG_NELEMENTS];
  MPI_Status status;


  if(All.ComovingIntegrationOn)	/* Factors for comoving integration of hydro */
    {
      a3inv = 1 / (All.Time * All.Time * All.Time);
      hubble_a = All.Hubble * sqrt(All.Omega0 / (All.Time * All.Time * All.Time)
				   + (1 - All.Omega0 - All.OmegaLambda) / (All.Time * All.Time) +
				   All.OmegaLambda);
      time_hubble_a = All.Time * hubble_a;
    }
  else
    a3inv = time_hubble_a = 1;



  noffset = mymalloc(sizeof(int) * NTask);	/* offsets of bunches in common list */
  nbuffer = mymalloc(sizeof(int) * NTask);
  nsend_local = mymalloc(sizeof(int) * NTask);
  nsend = mymalloc(sizeof(int) * NTask * NTask);


  for(n = 0, numstars = 0; n < NumPart; n++)
    {
      if(P[n].Ti_endstep == All.Ti_Current)
        if(P[n].Type == 4)	/* count stars  */
          numstars++;
    }

  MPI_Allreduce(&numstars, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  i = 0;			/* begin with this index */
  npleft = numstars;		/* particles left for this task */
  ntotleft = ntot;		/* particles left for all tasks together */


  while(ntotleft > 0)
    {
      for(j = 0; j < NTask; j++)
	nsend_local[j] = 0;

      /* do local particles and prepare export list */

      for(nexport = 0, ndone = 0; i < NumPart && nexport < All.BunchSizeMetal - NTask; i++)
	if(P[i].Type == 4)
	  {
	    if(P[i].Ti_endstep == All.Ti_Current)
	      {
		ndone++;

		if(All.ComovingIntegrationOn)
		  time_begstep = All.TimeBegin * exp(P[i].Ti_begstep * All.Timebase_interval);
		else
		  time_begstep = All.TimeBegin + P[i].Ti_begstep * All.Timebase_interval;

		/*  the actual time-step */
		dt = (P[i].Ti_endstep - P[i].Ti_begstep) * All.Timebase_interval;
		if(All.ComovingIntegrationOn)
		  dtime = All.Time * dt / time_hubble_a;
		else
		  dtime = dt;

		age_of_star = bg_get_elapsed_time(P[i].StellarAge, time_begstep);


		/* this function encapsulates the stellar evolution. It returns how much
		   mass in each element is returned during stellar evolution of the stellar
		   population between [age_of_star, age_of_star+ dtime] */

		bg_stellar_evolution(age_of_star, dtime, P[i].Mass, P[i].Metals, MetalsReleased,
				     MetalsChangeInStar);

		/* Let's change the mass of the star accordingly */
		for(j = 0; j < BG_NELEMENTS; j++)
		  {
		    P[i].Metals[j] += MetalsChangeInStar[j];
		    P[i].Mass += MetalsChangeInStar[j];
		  }

		/* now we check whether the star has just been created in the last step,
		   in which case it is eligible to spawn a wind */

		if(fabs(P[i].StellarAge - time_begstep) < 1.0e-6 * time_begstep)	/* note: P[i].StellarAge may be single precision, that's why we may miss identity due to round off */
		  StarMass = P[i].Mass;	/* the star is in its first step after creation */
		else
		  StarMass = 0;	/* this signals that the star has been around for a while -> no wind creation */


		for(j = 0; j < NTask; j++)
		  Exportflag[j] = 0;

		bg_enrich_evaluate(i, 0);

		for(j = 0; j < NTask; j++)
		  {
		    if(Exportflag[j])
		      {
			MetalDataIn[nexport].Pos[0] = P[i].Pos[0];
			MetalDataIn[nexport].Pos[1] = P[i].Pos[1];
			MetalDataIn[nexport].Pos[2] = P[i].Pos[2];

			MetalDataIn[nexport].Hsml = PPP[i].Hsml;
			MetalDataIn[nexport].SolidAngleWeightSum = P[i].SolidAngleWeightSum;
			MetalDataIn[nexport].StarMass = StarMass;

			for(k = 0; k < BG_NELEMENTS; k++)
			  MetalDataIn[nexport].Metals[k] = MetalsReleased[k];

			MetalDataIn[nexport].ID = P[i].ID;
			MetalDataIn[nexport].Index = i;
			MetalDataIn[nexport].Task = j;
			nexport++;
			nsend_local[j]++;
		      }
		  }
	      }
	  }

      qsort(MetalDataIn, nexport, sizeof(struct metaldata_in), bg_metals_compare_key);

      for(j = 1, noffset[0] = 0; j < NTask; j++)
	noffset[j] = noffset[j - 1] + nsend_local[j - 1];

      MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);


      /* now do the particles that need to be exported */

      for(level = 1; level < (1 << PTask); level++)
	{
	  for(j = 0; j < NTask; j++)
	    nbuffer[j] = 0;
	  for(ngrp = level; ngrp < (1 << PTask); ngrp++)
	    {
	      maxfill = 0;
	      for(j = 0; j < NTask; j++)
		{
		  if((j ^ ngrp) < NTask)
		    if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
		      maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		}
	      if(maxfill >= All.BunchSizeMetal)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      /* get the particles */
		      MPI_Sendrecv(&MetalDataIn[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct metaldata_in), MPI_BYTE,
				   recvTask, TAG_ENRICH_A,
				   &MetalDataGet[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct metaldata_in),
				   MPI_BYTE, recvTask, TAG_ENRICH_A, MPI_COMM_WORLD, &status);
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }

	  for(j = 0; j < nbuffer[ThisTask]; j++)
	    bg_enrich_evaluate(j, 1);

	  level = ngrp - 1;
	}

      MPI_Allreduce(&ndone, &ndonetot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      npleft -= ndone;
      ntotleft -= ndonetot;
    }

  myfree(nsend);
  myfree(nsend_local);
  myfree(nbuffer);
  myfree(noffset);



  /* As a final act, we need to recompute the total mass of the gas particles, as they may have been enriched */
  /* CDV
  for(n = 0; n < N_gas; n++)
    {
      if(P[n].Type == 0)
	{
	  for(k = 0, mass = 0; k < BG_NELEMENTS; k++)
	    mass += P[n].Metals[k];

	  P[n].Mass = mass;
	}
    }
  */

  for (n = 0; n < N_gas; n++)
    {
      if (P[n].Type == 0)
	{
	  for (k = 0, mass = 0; k < BG_NELEMENTS; k++)
	    {
	      mass += P[n].Metals[k];
	      if (strcmp(SPH_Element_Name[k], "Silicon") == 0)
		{
		  mass += SULPSILI * P[n].Metals[k];
		  mass += CALCSILI * P[n].Metals[k];
		}
	    }

	  P[n].Mass = mass;
	}
    }

}



/*! This function represents the core of the enrichment calculation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
void bg_enrich_evaluate(int target, int mode)
{
  int j, k, n;
  unsigned int id;
  int startnode, numngb_inbox;
  double h, h2, omega, dm, solidAngleWeightSum, mass_star;
  double dx, dy, dz, r, r2;

#ifdef BG_WINDS
  double prob, dir[3], theta, phi, norm, ascale;
#endif
  FLOAT *pos, *metals;

  if(mode == 0)
    {
      pos = P[target].Pos;
      h = PPP[target].Hsml;
      solidAngleWeightSum = P[target].SolidAngleWeightSum;
      metals = MetalsReleased;
      mass_star = StarMass;
      id = P[target].ID;
    }
  else
    {
      pos = MetalDataGet[target].Pos;
      h = MetalDataGet[target].Hsml;
      solidAngleWeightSum = MetalDataGet[target].SolidAngleWeightSum;
      metals = MetalDataGet[target].Metals;
      mass_star = MetalDataGet[target].StarMass;
      id = MetalDataGet[target].ID;
    }

  h2 = h * h;

  startnode = All.MaxPart;
  do
    {
      numngb_inbox = ngb_treefind_variable(&pos[0], h, &startnode);

      for(n = 0; n < numngb_inbox; n++)
	{
	  j = Ngblist[n];

	  dx = pos[0] - P[j].Pos[0];
	  dy = pos[1] - P[j].Pos[1];
	  dz = pos[2] - P[j].Pos[2];

#ifdef PERIODIC			/*  now find the closest image in the given box size  */
	  if(dx > boxHalf_X)
	    dx -= boxSize_X;
	  if(dx < -boxHalf_X)
	    dx += boxSize_X;
	  if(dy > boxHalf_Y)
	    dy -= boxSize_Y;
	  if(dy < -boxHalf_Y)
	    dy += boxSize_Y;
	  if(dz > boxHalf_Z)
	    dz -= boxSize_Z;
	  if(dz < -boxHalf_Z)
	    dz += boxSize_Z;
#endif
	  r2 = dx * dx + dy * dy + dz * dz;

	  if(r2 < h2)
	    {
	      r = sqrt(r2);

	      omega = 2 * M_PI * (1 - r / sqrt(r2 + PPP[j].Hsml * PPP[j].Hsml));

	      for(k = 0; k < BG_NELEMENTS; k++)
		{
		  dm = metals[k] * omega * P[j].Mass / solidAngleWeightSum;
		  P[j].Metals[k] += dm;
		}

#ifdef BG_WINDS
	      prob = omega * All.WindMassLoading * mass_star / solidAngleWeightSum;

	      if(get_random_number(P[j].ID + id) < prob)	/* put this gas particle into the wind */
		{
		  for(k = 0, norm = 0; k < 3; k++)
		    {
		      dir[k] = P[j].Pos[k] - pos[k];
		      norm += dir[k] * dir[k];
		    }

		  if(norm > 0)
		    {
		      norm = sqrt(norm);
		      for(k = 0; k < 3; k++)
			dir[k] /= norm;
		    }
		  else
		    {
		      theta = acos(2 * get_random_number(P[j].ID + 1 + id) - 1);
		      phi = 2 * M_PI * get_random_number(P[j].ID + 2 + id);

		      dir[0] = sin(theta) * cos(phi);
		      dir[1] = sin(theta) * sin(phi);
		      dir[2] = cos(theta);
		    }

		  if(All.ComovingIntegrationOn)
		    ascale = All.Time;
		  else
		    ascale = 1;

		  for(k = 0; k < 3; k++)
		    {
		      P[j].Vel[k] += All.WindSpeed * ascale * dir[k];
		      SphP[j].VelPred[k] += All.WindSpeed * ascale * dir[k];
		    }
		}
#endif
	    }
	}
    }
  while(startnode >= 0);

}


/*! This routine is a comparison kernel used in a sort routine to group
 *  particles that are exported to the same processor.
 */
int bg_metals_compare_key(const void *a, const void *b)
{
  if(((struct metaldata_in *) a)->Task < (((struct metaldata_in *) b)->Task))
    return -1;

  if(((struct metaldata_in *) a)->Task > (((struct metaldata_in *) b)->Task))
    return +1;

  return 0;
}

#endif
