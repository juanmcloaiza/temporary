#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"

#ifdef SFR_METALS
#include "c_metals.h"
#endif


/*! \file accel.c
 *  \brief driver routines to carry out force computation
 */
#ifdef REIONIZATION
void heating(void);
#endif


/*! This routine computes the accelerations for all active particles.  First, the gravitational forces are
 * computed. This also reconstructs the tree, if needed, otherwise the drift/kick operations have updated the
 * tree to make it fullu usable at the current time.
 *
 * If gas particles are presented, the `interior' of the local domain is determined. This region is guaranteed
 * to contain only particles local to the processor. This information will be used to reduce communication in
 * the hydro part.  The density for active SPH particles is computed next. If the number of neighbours should
 * be outside the allowed bounds, it will be readjusted by the function ensure_neighbours(), and for those
 * particle, the densities are recomputed accordingly. Finally, the hydrodynamical forces are added.
 */
void compute_accelerations(int mode)
{
  double tstart, tend;

#if defined(BUBBLES) || defined(MULTI_BUBBLES)
  double hubble_a;
#endif

#ifdef LT_STELLAREVOLUTION
  double tstartSn, tendSn;
  double CPU_sev, sumCPU_sev;
  int SnEv;
#endif

#ifdef SFR_METALS
  int i;
  double u_i, a3inv;

  if(All.ComovingIntegrationOn)	/* Factors for comoving integration of hydro */
    a3inv = 1 / (All.Time * All.Time * All.Time);
  else
    a3inv = 1;

  TotalEnergy = 0;
  DEnergy_spawned = 0;
  DEnergy_converted = 0;
  DEnergy_radiation = 0;
  DEnergy_promotion = 0;
  DEnergy_feedback = 0;


/*Not only for active particles! */
  for(i = 0; i < NumPart; i++)
    if(P[i].Type == 0)
      {
	u_i = SphP[i].Entropy / GAMMA_MINUS1 * pow(SphP[i].Density * a3inv, GAMMA_MINUS1);
	TotalEnergy += u_i * P[i].Mass;
      }

#endif


#ifdef REIONIZATION
  heating();
#endif


  if(ThisTask == 0)
    {
      printf("Start force computation...\n");
      fflush(stdout);
    }

#ifdef PMGRID
  if(All.PM_Ti_endstep == All.Ti_Current)
    {
      tstart = second();
      long_range_force();
      tend = second();
      All.CPU_PM += timediff(tstart, tend);
      CPU_Step[CPU_MESH] += timediff(tstart, tend);
    }
#endif

#ifndef ONLY_PM

  tstart = second();		/* measure the time for the full force computation */

  gravity_tree();		/* computes gravity accel. */

  if(All.TypeOfOpeningCriterion == 1 && All.Ti_Current == 0)
    gravity_tree();		/* For the first timestep, we redo it
				 * to allow usage of relative opening
				 * criterion for consistent accuracy.
				 */
  tend = second();
  All.CPU_Gravity += timediff(tstart, tend);

#endif


#ifdef FORCETEST
  gravity_forcetest();
#endif


#ifdef HPM
  hpm_find_eqs();
#else

  if(All.TotN_gas > 0)
    {
      if(ThisTask == 0)
	{
	  printf("Start density computation...\n");
	  fflush(stdout);
	}

      tstart = second();

#ifdef SFR_METALS
#ifdef SFR_ENRICH		/* equivalent to SNI OR SNII */
      flag_SN_starparticles();	/* mark SNI star particles */
#endif
#endif

#if defined(SFR_METALS) && defined(SFR_DECOUPLING)
      copy_densities();
      find_low_density_tail();
#endif
      density();		/* computes density, and pressure */


#if defined(MAGNETIC) && defined(BSMOOTH)
      bsmooth();
#endif


#if defined(SFR_METALS) && defined(SFR_DECOUPLING) && defined(SFR_PROMOTION)
      find_hot_neighbours();
      promote_particles();
      copy_densities();		/* optional ? */
      density();
#endif

#if ( defined(CONDUCTION) || defined(CR_DIFFUSION) || defined(SMOOTH_PHI) || defined(VOLUME_CORRECTION))
      compute_smoothed_values();

#ifdef STOP_AFTER_VOL_CORRECTION
   savepositions(All.SnapshotFileCount++);
   MPI_Barrier(MPI_COMM_WORLD);
   if(ThisTask ==0) endrun(9965);
#endif
#endif

      tend = second();
      All.CPU_Hydro += timediff(tstart, tend);

      tstart = second();

      force_update_hmax();

      tend = second();
      All.CPU_Predict += timediff(tstart, tend);

      if(ThisTask == 0)
	{
	  printf("Start hydro-force computation...\n");
	  fflush(stdout);
	}

      tstart = second();

      hydro_force();		/* adds hydrodynamical accelerations 
				   and computes du/dt  */
      tend = second();
      All.CPU_Hydro += timediff(tstart, tend);

#ifdef MHM
      kinetic_feedback_mhm();
#endif


#ifdef BLACK_HOLES
      blackhole_accretion();

#ifdef FOF
      /* this will find new black hole seed halos */
      if(All.Time >= All.TimeNextBlackHoleCheck)
	{
	  fof_fof(-1);

	  if(All.ComovingIntegrationOn)
	    All.TimeNextBlackHoleCheck *= All.TimeBetBlackHoleSearch;
	  else
	    All.TimeNextBlackHoleCheck += All.TimeBetBlackHoleSearch;
	}
#endif
#endif


#ifdef BG_STELLAR_EVOLUTION
      bg_enrich();		/* do chemical enrichment and wind generation in blue-gene model */
#endif


#ifdef COOLING
      tstart = second();

      cooling_and_starformation();	/* do radiative cooling and star formation */

#ifdef LT_STELLAREVOLUTION
      if(ThisTask == 0)
	{
	  printf("Start supernovae computation...\n");
	  fflush(stdout);
	}

      tstartSn = second();

      SnEv = evolve_SN(EVOLVE_SN);

      tendSn = second();
      CPU_sev = timediff(tstartSn, tendSn);

      MPI_Reduce(&CPU_sev, &sumCPU_sev, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      if(ThisTask == 0)
	All.CPU_SEv += sumCPU_sev / NTask;
#endif

#ifdef SFR_METALS
#ifdef SFR_ENRICH
#ifndef SFR_FEEDBACK
      if(ThisTask == 0)
	{
	  printf("... updating weights ...\n");
	  fflush(stdout);
	}
      update_weights();

      if(ThisTask == 0)
	{
	  printf("Start metal enrichment computation...\n");
	  fflush(stdout);
	}
      enrichment();		/* distribute metals within neighbourhood */

#else
/* HOT phase Flag_phase =1*/
      Flag_phase = 1;
      update_weights();
      if(ThisTask == 0)
	{
	  printf("Start metal enrichment computation phase I...\n");
	  fflush(stdout);
	}

      enrichment();		/* distribute metals within neighbourhood */

/* COLD phase Flag_phase =2*/
      Flag_phase = 2;
      update_weights();
      if(ThisTask == 0)
	{
	  printf("Start metal enrichment computation phase II..Cleaning...\n");
	  fflush(stdout);
	}

      enrichment();		/* distribute metals within neighbourhood */

      Flag_phase = 0;

#endif
#endif

#ifdef SFR_FEEDBACK
      /*      phase_mass(); */
      energy_test();
#endif

#endif

#ifdef CR_DIFFUSION_GREEN
      greenf_diffusion();
#endif


      tend = second();
      All.CPU_Hydro += timediff(tstart, tend);
      All.CPU_SfrCool += timediff(tstart, tend);
      CPU_Step[CPU_COOLINGSFR] += timediff(tstart, tend);
#endif



#ifdef BUBBLES
      if(All.Time >= All.TimeOfNextBubble)
	{
#ifdef FOF
	  fof_fof(-1);
	  bubble();
#else
	  bubble();
#endif
	  if(All.ComovingIntegrationOn)
	    {
              hubble_a = hubble_function(All.Time);
	      All.TimeOfNextBubble *= (1.0 + All.BubbleTimeInterval * hubble_a);
	    }
	  else
	    All.TimeOfNextBubble += All.BubbleTimeInterval / All.UnitTime_in_Megayears;
	  if(ThisTask == 0)
	    printf("Time of the bubble generation: %g\n", 1. / All.TimeOfNextBubble - 1.);
	}
#endif


#if defined(MULTI_BUBBLES) && defined(FOF)
      if(All.Time >= All.TimeOfNextBubble)
	{
	  fof_fof(-1);

	  if(All.ComovingIntegrationOn)
	    {
              hubble_a = Hubble_func(All.Time);
	      All.TimeOfNextBubble *= (1.0 + All.BubbleTimeInterval * hubble_a);
	    }
	  else
	    All.TimeOfNextBubble += All.BubbleTimeInterval / All.UnitTime_in_Megayears;
	  if(ThisTask == 0)
	    printf("Time of the bubble generation: %g\n", 1. / All.TimeOfNextBubble - 1.);
	}
#endif



    }

#endif /* end of HPM */

  if(ThisTask == 0)
    {
      printf("force computation done.\n");
      fflush(stdout);
    }
}

#ifdef REIONIZATION
void heating(void)
{
  int i;
  double u, temp, meanweight, a3inv;

  /* reionization: Tmin = 10^4 Kelvin @ z = 10: */

  if(Flag_FullStep)
    {
      if(All.not_yet_reionized)
	{
	  if(1 / All.Time - 1 < 10.0)
	    {
	      All.not_yet_reionized = 0;

	      meanweight = 4.0 / (8 - 5 * (1 - HYDROGEN_MASSFRAC)) * PROTONMASS;	/* fully reionized */

	      a3inv = 1 / (All.Time * All.Time * All.Time);

	      for(i = 0; i < N_gas; i++)
		{
		  u = SphP[i].Entropy / GAMMA_MINUS1 * pow(SphP[i].a2.Density * a3inv, GAMMA_MINUS1);

		  temp =
		    meanweight / BOLTZMANN * GAMMA_MINUS1 * u * All.UnitEnergy_in_cgs / All.UnitMass_in_g;

		  if(temp < 1.0e4)
		    temp = 1.0e4;

		  u =
		    temp / (meanweight / BOLTZMANN * GAMMA_MINUS1 * All.UnitEnergy_in_cgs /
			    All.UnitMass_in_g);

		  SphP[i].Entropy = u * GAMMA_MINUS1 / pow(SphP[i].a2.Density * a3inv, GAMMA_MINUS1);
		}
	    }
	}
    }
}
#endif
