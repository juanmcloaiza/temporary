#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"

/*! \file blackhole.c
 *  \brief routines for gas accretion onto black holes, and black hole mergers
 */

#ifdef BLACK_HOLES

static double hubble_a, ascale;


static int N_gas_swallowed, N_BH_swallowed;

void blackhole_accretion(void)
{
  int i, j, k, n, ngrp, maxfill, source;
  int ntot, ntotleft, ndone, ndonetot;
  int *nbuffer, *noffset, *nsend_local, *nsend;
  int level, sendTask, recvTask;
  int nexport, place, num_blackholes;
  int Ntot_gas_swallowed, Ntot_BH_swallowed;
  int total_num_blackholes;
  double mdot, rho, bhvel, soundspeed, meddington, dt, mdot_in_msun_per_year;
  double mass_real, total_mass_real;
  double mass_holes, total_mass_holes, total_mdot;
  double fac;
  double mdoteddington, total_mdoteddington;
  MPI_Status status;

  if(ThisTask == 0)
    {
      printf("Beginning black-hole accretion\n");
      fflush(stdout);
    }

  if(All.ComovingIntegrationOn)
    {
      ascale = All.Time;
      hubble_a = hubble_function(All.Time);
    }
  else
    hubble_a = ascale = 1;


  /* Let's first compute the Mdot values */

  for(n = 0; n < NumPart; n++)
    if(P[n].Type == 5)
      if(P[n].Ti_endstep == All.Ti_Current)
	{
	  mdot = 0;		/* if no accretion model is enabled, we have mdot=0 */

	  rho = P[n].b1.BH_Density;

	  bhvel = sqrt(pow(P[n].Vel[0] - P[n].b3.BH_SurroundingGasVel[0], 2) +
		       pow(P[n].Vel[1] - P[n].b3.BH_SurroundingGasVel[1], 2) +
		       pow(P[n].Vel[2] - P[n].b3.BH_SurroundingGasVel[2], 2));

	  if(All.ComovingIntegrationOn)
	    {
	      bhvel /= All.Time;
	      rho /= pow(All.Time, 3);
	    }

	  soundspeed = sqrt(GAMMA * P[n].b2.BH_Entropy * pow(rho, GAMMA_MINUS1));

	  /* Note: we take here a radiative efficiency of 0.1 for Eddington accretion */

	  meddington = (4 * M_PI * GRAVITY * C * PROTONMASS / (0.1 * C * C * THOMPSON)) * P[n].BH_Mass
	    * All.UnitTime_in_s;

#ifdef BONDI
	  mdot = 4. * M_PI * All.BlackHoleAccretionFactor * All.G * All.G *
	    P[n].BH_Mass * P[n].BH_Mass * rho / pow((pow(soundspeed, 2) + pow(bhvel, 2)), 1.5);
#endif


#ifdef ENFORCE_EDDINGTON_LIMIT
	  if(mdot > All.BlackHoleEddingtonFactor * meddington)
	    mdot = All.BlackHoleEddingtonFactor * meddington;
#endif
	  P[n].BH_Mdot = mdot;

	  if(meddington > 0)
	    P[n].BH_MdotEddington = mdot / meddington;	/* this stores the accretion rate in units of the Eddington rate */
	  else
	    P[n].BH_MdotEddington = 0;

	  if(P[n].BH_Mass > 0)
	    {
#ifndef LONGIDS
	      fprintf(FdBlackHolesDetails, "BH=%u %g %g %g %g %g\n",
		      P[n].ID, All.Time, P[n].BH_Mass, mdot, rho, soundspeed);
#else
	      fprintf(FdBlackHolesDetails, "BH=%llu %g %g %g %g %g\n",
		      P[n].ID, All.Time, P[n].BH_Mass, mdot, rho, soundspeed);
#endif
	    }

	  dt = (P[n].Ti_endstep - P[n].Ti_begstep) * All.Timebase_interval / hubble_a;

#ifdef BH_DRAG
	  /* add a drag force for the black-holes,
	     accounting for the accretion */

	  if(P[n].BH_Mass > 0)
	    {
	      /*
		 fac = P[n].BH_Mdot * dt / P[n].BH_Mass;
	       */
	      fac = meddington * dt / P[n].BH_Mass;

	      if(fac > 1)
		fac = 1;

	      if(dt > 0)
		for(k = 0; k < 3; k++)
		  P[n].g.GravAccel[k] +=
		    -ascale * ascale * fac / dt * (P[n].Vel[k] - P[n].b3.BH_SurroundingGasVel[k]) / ascale;
	    }
#endif

	  P[n].BH_Mass += P[n].BH_Mdot * dt;


#ifdef BH_KINETICFEEDBACK
	  if(mdot >= 0.99 * meddington)
	    {
	      P[n].ActiveTime += dt;
	      P[n].ActiveEnergy += All.BlackHoleFeedbackFactor * 0.1 * mdot * dt *
		pow(C / All.UnitVelocity_in_cm_per_s, 2);
	    }
#endif
	}


  /* Now let's invoke the functions that stochasticall swallow gas
   * and deal with black hole mergers.
   */

  if(ThisTask == 0)
    {
      printf("Start swallowing of gas particles and black holes\n");
      fflush(stdout);
    }


  N_gas_swallowed = N_BH_swallowed = 0;

  for(n = 0, num_blackholes = 0; n < NumPart; n++)
    {
      if(P[n].Type == 5)
	if(P[n].Ti_endstep == All.Ti_Current)
	  num_blackholes++;
    }

  MPI_Allreduce(&num_blackholes, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  noffset = (int *) mymalloc(sizeof(int) * NTask);	/* offsets of bunches in common list */
  nbuffer = (int *) mymalloc(sizeof(int) * NTask);
  nsend_local = (int *) mymalloc(sizeof(int) * NTask);
  nsend = (int *) mymalloc(sizeof(int) * NTask * NTask);








  /** Let's first spread the feedback energy, and determine which particles may be swalled by whom */


  i = 0;			/* first particle for this task */
  ntotleft = ntot;		/* particles left for all tasks together */

  while(ntotleft > 0)
    {
      for(j = 0; j < NTask; j++)
	nsend_local[j] = 0;

      /* do local particles and prepare export list */

      for(nexport = 0, ndone = 0; i < NumPart && nexport < All.BunchSizeBlackhole - NTask; i++)
	if(P[i].Type == 5)
	  if(P[i].Ti_endstep == All.Ti_Current)
	    {
	      ndone++;

	      for(j = 0; j < NTask; j++)
		Exportflag[j] = 0;

	      blackhole_evaluate(i, 0);

	      for(j = 0; j < NTask; j++)
		{
		  if(Exportflag[j])
		    {
		      for(k = 0; k < 3; k++)
			{
			  BlackholeDataIn[nexport].Pos[k] = P[i].Pos[k];
			  BlackholeDataIn[nexport].Vel[k] = P[i].Vel[k];
			}

		      BlackholeDataIn[nexport].Hsml = PPP[i].Hsml;
		      BlackholeDataIn[nexport].Mass = P[i].Mass;
		      BlackholeDataIn[nexport].BH_Mass = P[i].BH_Mass;
#ifdef BH_KINETICFEEDBACK
		      BlackholeDataIn[nexport].ActiveTime = P[i].ActiveTime;
		      BlackholeDataIn[nexport].ActiveEnergy = P[i].ActiveEnergy;
#endif
		      BlackholeDataIn[nexport].Density = P[i].b1.BH_Density;
		      BlackholeDataIn[nexport].Mdot = P[i].BH_Mdot;
		      BlackholeDataIn[nexport].Csnd =
			sqrt(GAMMA * P[i].b2.BH_Entropy *
			     pow(P[i].b1.BH_Density / (ascale * ascale * ascale), GAMMA_MINUS1));
		      BlackholeDataIn[nexport].Dt =
			(P[i].Ti_endstep - P[i].Ti_begstep) * All.Timebase_interval / hubble_a;
		      BlackholeDataIn[nexport].ID = P[i].ID;
		      BlackholeDataIn[nexport].Index = i;
		      BlackholeDataIn[nexport].Task = j;
		      nexport++;
		      nsend_local[j]++;
		    }
		}
	    }

      qsort(BlackholeDataIn, nexport, sizeof(struct blackholedata_in), blackhole_compare_key);

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
	      if(maxfill >= All.BunchSizeBlackhole)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      /* get the particles */
		      MPI_Sendrecv(&BlackholeDataIn[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct blackholedata_in), MPI_BYTE,
				   recvTask, TAG_BH_A,
				   &BlackholeDataGet[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct blackholedata_in),
				   MPI_BYTE, recvTask, TAG_BH_A, MPI_COMM_WORLD, &status);
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }

	  /* now do the imported particles */
	  for(j = 0; j < nbuffer[ThisTask]; j++)
	    blackhole_evaluate(j, 1);

	  /* get the result */

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
	      if(maxfill >= All.BunchSizeBlackhole)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      /* send the results */
		      MPI_Sendrecv(&BlackholeDataResult[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct blackholedata_out),
				   MPI_BYTE, recvTask, TAG_BH_B,
				   &BlackholeDataPartialResult[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct blackholedata_out),
				   MPI_BYTE, recvTask, TAG_BH_B, MPI_COMM_WORLD, &status);

		      /* add the result to the particles */
		      for(j = 0; j < nsend_local[recvTask]; j++)
			{
			  source = j + noffset[recvTask];
			  place = BlackholeDataIn[source].Index;
#ifdef REPOSITION_ON_POTMIN
			  if(P[place].BH_MinPot > BlackholeDataPartialResult[source].BH_MinPot)
			    {
			      P[place].BH_MinPot = BlackholeDataPartialResult[source].BH_MinPot;
			      for(k = 0; k < 3; k++)
				P[place].BH_MinPotPos[k] = BlackholeDataPartialResult[source].BH_MinPotPos[k];
			    }
#endif
			}
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }

	  level = ngrp - 1;
	}

      MPI_Allreduce(&ndone, &ndonetot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      ntotleft -= ndonetot;
    }










  /* Now do the swallowing of particles */



  i = 0;			/* first particle for this task */
  ntotleft = ntot;		/* particles left for all tasks together */

  while(ntotleft > 0)
    {
      for(j = 0; j < NTask; j++)
	nsend_local[j] = 0;

      /* do local particles and prepare export list */

      for(nexport = 0, ndone = 0; i < NumPart && nexport < All.BunchSizeBlackhole - NTask; i++)
	if(P[i].Type == 5)
	  if(P[i].Ti_endstep == All.Ti_Current)
	    {
	      ndone++;

	      if(P[i].SwallowID == 0)
		{
		  for(j = 0; j < NTask; j++)
		    Exportflag[j] = 0;

		  blackhole_evaluate_swallow(i, 0);

		  for(j = 0; j < NTask; j++)
		    {
		      if(Exportflag[j])
			{
			  for(k = 0; k < 3; k++)
			    BlackholeDataIn[nexport].Pos[k] = P[i].Pos[k];

			  BlackholeDataIn[nexport].Hsml = PPP[i].Hsml;
			  BlackholeDataIn[nexport].BH_Mass = P[i].BH_Mass;
			  BlackholeDataIn[nexport].ID = P[i].ID;
			  BlackholeDataIn[nexport].Index = i;
			  BlackholeDataIn[nexport].Task = j;
			  nexport++;
			  nsend_local[j]++;
			}
		    }
		}
	    }

      qsort(BlackholeDataIn, nexport, sizeof(struct blackholedata_in), blackhole_compare_key);

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
	      if(maxfill >= All.BunchSizeBlackhole)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      /* get the particles */
		      MPI_Sendrecv(&BlackholeDataIn[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct blackholedata_in), MPI_BYTE,
				   recvTask, TAG_BH_A,
				   &BlackholeDataGet[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct blackholedata_in),
				   MPI_BYTE, recvTask, TAG_BH_A, MPI_COMM_WORLD, &status);
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }

	  /* now do the imported particles */
	  for(j = 0; j < nbuffer[ThisTask]; j++)
	    blackhole_evaluate_swallow(j, 1);

	  /* get the result */

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
	      if(maxfill >= All.BunchSizeBlackhole)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      /* send the results */
		      MPI_Sendrecv(&BlackholeDataResult[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct blackholedata_out),
				   MPI_BYTE, recvTask, TAG_BH_B,
				   &BlackholeDataPartialResult[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct blackholedata_out),
				   MPI_BYTE, recvTask, TAG_BH_B, MPI_COMM_WORLD, &status);

		      /* add the result to the particles */
		      for(j = 0; j < nsend_local[recvTask]; j++)
			{
			  source = j + noffset[recvTask];
			  place = BlackholeDataIn[source].Index;

			  P[place].b4.dBH_accreted_Mass += BlackholeDataPartialResult[source].Mass;
			  P[place].b5.dBH_accreted_BHMass += BlackholeDataPartialResult[source].BH_Mass;
			  for(k = 0; k < 3; k++)
			    P[place].b6.dBH_accreted_momentum[k] +=
			      BlackholeDataPartialResult[source].AccretedMomentum[k];

			}
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }

	  level = ngrp - 1;
	}

      MPI_Allreduce(&ndone, &ndonetot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      ntotleft -= ndonetot;
    }


  myfree(nsend);
  myfree(nsend_local);
  myfree(nbuffer);
  myfree(noffset);

  MPI_Reduce(&N_gas_swallowed, &Ntot_gas_swallowed, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&N_BH_swallowed, &Ntot_BH_swallowed, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      printf("Accretion done: %d gas particles swallowed, %d BH particles swallowed\n",
	     Ntot_gas_swallowed, Ntot_BH_swallowed);
      fflush(stdout);
    }



#ifdef REPOSITION_ON_POTMIN
  for(n = 0; n < NumPart; n++)
    if(P[n].Ti_endstep == All.Ti_Current)
      if(P[n].Type == 5)
	for(k = 0; k < 3; k++)
	  P[n].Pos[k] = P[n].BH_MinPotPos[k];
#endif



  for(n = 0, num_blackholes = 0, mdot = 0, mass_holes = 0, mass_real = 0, mdoteddington = 0; n < NumPart; n++)
    {
#ifdef FLTROUNDOFFREDUCTION
      if(P[n].Type == 0)
	SphP[n].i.Injected_BH_Energy = FLT(SphP[n].i.dInjected_BH_Energy);
#endif

      if(P[n].Type == 5)
	{
	  if(P[n].Ti_endstep == All.Ti_Current)
	    {
#ifdef FLTROUNDOFFREDUCTION
	      P[n].b4.BH_accreted_Mass = FLT(P[n].b4.dBH_accreted_Mass);
	      P[n].b5.BH_accreted_BHMass = FLT(P[n].b5.dBH_accreted_BHMass);
	      for(k = 0; k < 3; k++)
		P[n].b6.BH_accreted_momentum[k] = FLT(P[n].b6.dBH_accreted_momentum[k]);
#endif
	      if(P[n].b4.BH_accreted_Mass > 0)
		{
		  for(k = 0; k < 3; k++)
		    P[n].Vel[k] =
		      (P[n].Vel[k] * P[n].Mass + P[n].b6.BH_accreted_momentum[k]) /
		      (P[n].Mass + P[n].b4.BH_accreted_Mass);

		  P[n].Mass += P[n].b4.BH_accreted_Mass;
		  P[n].BH_Mass += P[n].b5.BH_accreted_BHMass;
		  P[n].b4.BH_accreted_Mass = 0;
		}
	    }


#ifdef BH_KINETICFEEDBACK
	  if(P[n].ActiveTime > All.BlackHoleActiveTime)
	    {
	      P[n].ActiveTime = 0;
	      P[n].ActiveEnergy = 0;
	    }
#endif

	  mass_holes += P[n].BH_Mass;
	  mass_real += P[n].Mass;
	  mdot += P[n].BH_Mdot;
	  mdoteddington += P[n].BH_MdotEddington;
	  num_blackholes++;

	}
    }

  MPI_Reduce(&mass_holes, &total_mass_holes, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&mass_real, &total_mass_real, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&mdot, &total_mdot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&mdoteddington, &total_mdoteddington, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&num_blackholes, &total_num_blackholes, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      /* convert to solar masses per yr */
      mdot_in_msun_per_year =
	total_mdot * (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

      fprintf(FdBlackHoles, "%g %d %g %g %g %g %g\n",
	      All.Time, total_num_blackholes, total_mass_holes, total_mdot, mdot_in_msun_per_year,
	      total_mass_real, total_mdoteddington);
      fflush(FdBlackHoles);
    }


  fflush(FdBlackHolesDetails);
}






void blackhole_evaluate(int target, int mode)
{
  int startnode, numngb, j, k, n, index, id;
  FLOAT *pos, *velocity, h_i, dt, mdot, rho, mass, bh_mass, csnd;
  double w, dx, dy, dz, h_i2, r2, r, u, hinv, hinv3, wk;
  double p, vrel;

#ifdef BH_KINETICFEEDBACK
  /*  double deltavel; */
  double activetime, activeenergy;
#endif
#ifdef BH_THERMALFEEDBACK
  double energy;
#endif
#ifdef REPOSITION_ON_POTMIN
  FLOAT minpotpos[3] = { 0, 0, 0 }, minpot = 1.0e30;
#endif

  if(mode == 0)
    {
      pos = P[target].Pos;
      rho = P[target].b1.BH_Density;
      mdot = P[target].BH_Mdot;
      dt = (P[target].Ti_endstep - P[target].Ti_begstep) * All.Timebase_interval / hubble_a;
      h_i = PPP[target].Hsml;
      mass = P[target].Mass;
      bh_mass = P[target].BH_Mass;
      velocity = P[target].Vel;
      csnd =
	sqrt(GAMMA * P[target].b2.BH_Entropy *
	     pow(P[target].b1.BH_Density / (ascale * ascale * ascale), GAMMA_MINUS1));
      index = target;
      id = P[target].ID;
#ifdef BH_KINETICFEEDBACK
      activetime = P[target].ActiveTime;
      activeenergy = P[target].ActiveEnergy;
#endif

    }
  else
    {
      pos = BlackholeDataGet[target].Pos;
      rho = BlackholeDataGet[target].Density;
      mdot = BlackholeDataGet[target].Mdot;
      dt = BlackholeDataGet[target].Dt;
      h_i = BlackholeDataGet[target].Hsml;
      mass = BlackholeDataGet[target].Mass;
      bh_mass = BlackholeDataGet[target].BH_Mass;
      velocity = BlackholeDataGet[target].Vel;
      csnd = BlackholeDataGet[target].Csnd;
      index = BlackholeDataGet[target].Index;
      id = BlackholeDataGet[target].ID;
#ifdef BH_KINETICFEEDBACK
      activetime = BlackholeDataGet[target].ActiveTime;
      activeenergy = BlackholeDataGet[target].ActiveEnergy;
#endif
    }

  /* initialize variables before SPH loop is started */
  h_i2 = h_i * h_i;

  /* Now start the actual SPH computation for this particle */
  startnode = All.MaxPart;
  do
    {
      numngb = ngb_treefind_blackhole(&pos[0], h_i, &startnode);

      for(n = 0; n < numngb; n++)
	{
	  j = Ngblist[n];

	  if(P[j].Mass > 0)
	    {
	      if(mass > 0)
		{
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

		  if(r2 < h_i2)
		    {
#ifdef REPOSITION_ON_POTMIN
		      /* if this option is switched on, we may also encounter dark matter particles or stars */
		      if(P[j].p.Potential < minpot)
			{
			  /* compute relative velocities */

			  for(k = 0, vrel = 0; k < 3; k++)
			    vrel += (P[j].Vel[k] - velocity[k]) * (P[j].Vel[k] - velocity[k]);

			  vrel = sqrt(vrel) / ascale;

			  if(vrel <= 0.25 * csnd)
			    {
			      minpot = P[j].p.Potential;
			      for(k = 0; k < 3; k++)
				minpotpos[k] = P[j].Pos[k];
			    }
			}
#endif
		      if(P[j].Type == 5)	/* we have a black hole merger */
			{
			  if(r2 > 0)
			    {
			      /* compute relative velocity of BHs */

			      for(k = 0, vrel = 0; k < 3; k++)
				vrel += (P[j].Vel[k] - velocity[k]) * (P[j].Vel[k] - velocity[k]);

			      vrel = sqrt(vrel) / ascale;

			      if(vrel > 0.5 * csnd)
				{
#ifndef LONGIDS
				  fprintf(FdBlackHolesDetails,
					  "ThisTask=%d, time=%g: id=%u would like to swallow %u, but vrel=%g csnd=%g\n",
					  ThisTask, All.Time, id, P[j].ID, vrel, csnd);
#else
				  fprintf(FdBlackHolesDetails,
					  "ThisTask=%d, time=%g: id=%llu would like to swallow %llu, but vrel=%g csnd=%g\n",
					  ThisTask, All.Time, id, P[j].ID, vrel, csnd);
#endif
				}
			      else
				{
				  if(P[j].SwallowID < id)
				    P[j].SwallowID = id;
				}
			    }
			}
		      if(P[j].Type == 0)
			{
			  /* here we have a gas particle */

			  r = sqrt(r2);
			  hinv = 1 / h_i;
#ifndef  TWODIMS
			  hinv3 = hinv * hinv * hinv;
#else
			  hinv3 = hinv * hinv / boxSize_Z;
#endif

			  u = r * hinv;

			  if(u < 0.5)
			    wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
			  else
			    wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

#ifdef SWALLOWGAS
			  /* compute accretion probability */

			  if((bh_mass - mass) > 0)
			    p = (bh_mass - mass) * wk / rho;
			  else
			    p = 0;

			  /* compute random number, uniform in [0,1] */
			  w = get_random_number(P[j].ID);
			  if(w < p)
			    {
			      if(P[j].SwallowID < id)
				P[j].SwallowID = id;
			    }
#endif

			  if(P[j].Mass > 0)
			    {
#ifdef BH_THERMALFEEDBACK
			      energy = All.BlackHoleFeedbackFactor * 0.1 * mdot * dt *
				pow(C / All.UnitVelocity_in_cm_per_s, 2);

			      SphP[j].i.dInjected_BH_Energy += FLT(energy * P[j].Mass * wk / rho);
#endif

#ifdef BH_KINETICFEEDBACK
			      if(activetime > All.BlackHoleActiveTime)
				{
				  SphP[j].i.dInjected_BH_Energy += FLT(activeenergy * P[j].Mass * wk / rho);
				}
#endif
			    }

			}
		    }
		}
	    }
	}
    }
  while(startnode >= 0);

  /* Now collect the result at the right place */
  if(mode == 0)
    {
#ifdef REPOSITION_ON_POTMIN
      for(k = 0; k < 3; k++)
	P[target].BH_MinPotPos[k] = minpotpos[k];
      P[target].BH_MinPot = minpot;
#endif
    }
  else
    {
#ifdef REPOSITION_ON_POTMIN
      for(k = 0; k < 3; k++)
	BlackholeDataResult[target].BH_MinPotPos[k] = minpotpos[k];
      BlackholeDataResult[target].BH_MinPot = minpot;
#endif
    }
}


void blackhole_evaluate_swallow(int target, int mode)
{
  int startnode, numngb, j, k, n, id;
  DOUBLE accreted_mass, accreted_BH_mass, accreted_momentum[3];
  FLOAT *pos, h_i, bh_mass;


  if(mode == 0)
    {
      pos = P[target].Pos;
      h_i = PPP[target].Hsml;
      id = P[target].ID;
      bh_mass = P[target].BH_Mass;
    }
  else
    {
      pos = BlackholeDataGet[target].Pos;
      h_i = BlackholeDataGet[target].Hsml;
      id = BlackholeDataGet[target].ID;
      bh_mass = BlackholeDataGet[target].BH_Mass;
    }

  accreted_mass = 0;
  accreted_BH_mass = 0;
  accreted_momentum[0] = accreted_momentum[1] = accreted_momentum[2] = 0;

  startnode = All.MaxPart;
  do
    {
      numngb = ngb_treefind_blackhole(&pos[0], h_i, &startnode);

      for(n = 0; n < numngb; n++)
	{
	  j = Ngblist[n];

	  if(P[j].SwallowID == id)
	    {
	      if(P[j].Type == 5)	/* we have a black hole merger */
		{
#ifndef LONGIDS
		  fprintf(FdBlackHolesDetails,
			  "ThisTask=%d, time=%g: id=%u swallows %u (%g %g)\n",
			  ThisTask, All.Time, id, P[j].ID, bh_mass, P[j].BH_Mass);
#else
		  fprintf(FdBlackHolesDetails,
			  "ThisTask=%d, time=%g: id=%llu swallows %llu (%g %g)\n",
			  ThisTask, All.Time, id, P[j].ID, bh_mass, P[j].BH_Mass);
#endif

		  accreted_mass += FLT(P[j].Mass);
		  accreted_BH_mass += FLT(P[j].BH_Mass);

		  for(k = 0; k < 3; k++)
		    accreted_momentum[k] += FLT(P[j].Mass * P[j].Vel[k]);

		  P[j].Mass = 0;
		  P[j].BH_Mass = 0;

		  N_BH_swallowed++;
		}
	    }

	  if(P[j].Type == 0)
	    {
	      if(P[j].SwallowID == id)
		{
		  accreted_mass += FLT(P[j].Mass);

		  for(k = 0; k < 3; k++)
		    accreted_momentum[k] += FLT(P[j].Mass * P[j].Vel[k]);

		  P[j].Mass = 0;

		  N_gas_swallowed++;
		}
	    }
	}
    }
  while(startnode >= 0);

  /* Now collect the result at the right place */
  if(mode == 0)
    {
      P[target].b4.dBH_accreted_Mass = accreted_mass;
      P[target].b5.dBH_accreted_BHMass = accreted_BH_mass;
      for(k = 0; k < 3; k++)
	P[target].b6.dBH_accreted_momentum[k] = accreted_momentum[k];
    }
  else
    {
      BlackholeDataResult[target].Mass = accreted_mass;
      BlackholeDataResult[target].BH_Mass = accreted_BH_mass;
      for(k = 0; k < 3; k++)
	BlackholeDataResult[target].AccretedMomentum[k] = accreted_momentum[k];
    }
}





/*! This is a comparison kernel for a sort routine, which is used to group
 *  particles that are going to be exported to the same CPU.
 */
int blackhole_compare_key(const void *a, const void *b)
{
  if(((struct blackholedata_in *) a)->Task < (((struct blackholedata_in *) b)->Task))
    return -1;
  if(((struct blackholedata_in *) a)->Task > (((struct blackholedata_in *) b)->Task))
    return +1;
  return 0;
}

#endif
