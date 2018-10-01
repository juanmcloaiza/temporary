#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "allvars.h"
#include "proto.h"
#ifdef COSMIC_RAYS
#include "cosmic_rays.h"
#endif

/*! \file timestep.c 
 *  \brief routines for 'kicking' particles in
 *  momentum space and assigning new timesteps
 */

static double fac1, fac2, fac3, hubble_a, atime, a3inv;
static double dt_displacement = 0;

/*! This function advances the system in momentum space, i.e. it does apply the 'kick' operation after the
 *  forces have been computed. Additionally, it assigns new timesteps to particles. At start-up, a
 *  half-timestep is carried out, as well as at the end of the simulation. In between, the half-step kick that
 *  ends the previous timestep and the half-step kick for the new timestep are combined into one operation.
 */
void advance_and_find_timesteps(void)
{
  int i, j, ti_step, ti_min, tend, tstart;
  double dt_entr, dt_entr2, dt_gravkick, dt_hydrokick, dt_gravkick2, dt_hydrokick2, t0, t1;
  double minentropy, aphys;
  FLOAT dv[3];

#ifdef FLEXSTEPS
  int ti_grp;
#endif
#if defined(PSEUDOSYMMETRIC) && !defined(FLEXSTEPS)
  double apred, prob;
  int ti_step2;
#endif
#ifdef PMGRID
  double dt_gravkickA, dt_gravkickB;
#endif
#ifdef MAKEGLASS
  double disp, dispmax, globmax, dmean, fac, disp2sum, globdisp2sum;
#endif
#if defined(TIME_DEP_ART_VISC) || defined(TIME_DEP_MAGN_DISP)
  double dmin1, dmin2;
#endif

#ifdef CHEMISTRY
  int ifunc, mode;
  double a_start, a_end;
#endif

  t0 = second();

  if(All.ComovingIntegrationOn)
    {
      fac1 = 1 / (All.Time * All.Time);
      fac2 = 1 / pow(All.Time, 3 * GAMMA - 2);
      fac3 = pow(All.Time, 3 * (1 - GAMMA) / 2.0);
      hubble_a = hubble_function(All.Time);

      a3inv = 1 / (All.Time * All.Time * All.Time);
      atime = All.Time;
    }
  else
    fac1 = fac2 = fac3 = hubble_a = a3inv = atime = 1;

  if(Flag_FullStep || dt_displacement == 0)
    find_dt_displacement_constraint(hubble_a * atime * atime);

#ifdef PMGRID
  if(All.ComovingIntegrationOn)
    dt_gravkickB = get_gravkick_factor(All.PM_Ti_begstep, All.Ti_Current) -
      get_gravkick_factor(All.PM_Ti_begstep, (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2);
  else
    dt_gravkickB = (All.Ti_Current - (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2) * All.Timebase_interval;

  if(All.PM_Ti_endstep == All.Ti_Current)	/* need to do long-range kick */
    {
      /* make sure that we reconstruct the domain/tree next time because we don't kick the tree nodes in this case */
      All.NumForcesSinceLastDomainDecomp = 1 + All.TotNumPart * All.TreeDomainUpdateFrequency;
    }
#endif


#ifdef MAKEGLASS
  for(i = 0, dispmax = 0, disp2sum = 0; i < NumPart; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  P[i].g.GravAccel[j] *= -1;
#ifdef PMGRID
	  P[i].GravPM[j] *= -1;
	  P[i].g.GravAccel[j] += P[i].GravPM[j];
	  P[i].GravPM[j] = 0;
#endif
	}

      disp = sqrt(P[i].g.GravAccel[0] * P[i].g.GravAccel[0] +
		  P[i].g.GravAccel[1] * P[i].g.GravAccel[1] + P[i].g.GravAccel[2] * P[i].g.GravAccel[2]);

      disp *= 2.0 / (3 * All.Hubble * All.Hubble);

      disp2sum += disp * disp;

      if(disp > dispmax)
	dispmax = disp;
    }

  MPI_Allreduce(&dispmax, &globmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&disp2sum, &globdisp2sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  dmean = pow(P[0].Mass / (All.Omega0 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)), 1.0 / 3);

  if(globmax > dmean)
    fac = dmean / globmax;
  else
    fac = 1.0;

  if(ThisTask == 0)
    {
      printf("\nglass-making:  dmean= %g  global disp-maximum= %g  rms= %g\n\n",
	     dmean, globmax, sqrt(globdisp2sum / All.TotNumPart));
      fflush(stdout);
    }

  for(i = 0, dispmax = 0; i < NumPart; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  P[i].Vel[j] = 0;
	  P[i].Pos[j] += fac * P[i].g.GravAccel[j] * 2.0 / (3 * All.Hubble * All.Hubble);
	  P[i].g.GravAccel[j] = 0;
	}
    }
#endif




  /* Now assign new timesteps and kick */

#ifdef FLEXSTEPS
  if((All.Ti_Current % (4 * All.PresentMinStep)) == 0)
    if(All.PresentMinStep < TIMEBASE)
      All.PresentMinStep *= 2;

  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Ti_endstep == All.Ti_Current)
	{
	  ti_step = get_timestep(i, &aphys, 0);

	  /* make it a power 2 subdivision */
	  ti_min = TIMEBASE;
	  while(ti_min > ti_step)
	    ti_min >>= 1;
	  ti_step = ti_min;

	  if(ti_step < All.PresentMinStep)
	    All.PresentMinStep = ti_step;
	}
    }

  ti_step = All.PresentMinStep;
  MPI_Allreduce(&ti_step, &All.PresentMinStep, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  ti_step = (int) (All.MaxSizeTimestep / All.Timebase_interval);
  /* make it a power 2 subdivision */
  ti_min = TIMEBASE;
  while(ti_min > ti_step)
    ti_min >>= 1;
  All.PresentMaxStep = ti_min;


  if(ThisTask == 0)
    printf("Syn Range = %g  PresentMinStep = %d  PresentMinStep = %d \n",
	   (double) All.PresentMaxStep / All.PresentMinStep, All.PresentMinStep, All.PresentMaxStep);

#endif


  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Ti_endstep == All.Ti_Current)
	{
	  ti_step = get_timestep(i, &aphys, 0);

	  /* make it a power 2 subdivision */
	  ti_min = TIMEBASE;
	  while(ti_min > ti_step)
	    ti_min >>= 1;
	  ti_step = ti_min;

#ifdef FLEXSTEPS
	  ti_grp = P[i].FlexStepGrp % All.PresentMaxStep;
	  ti_grp = (ti_grp / All.PresentMinStep) * All.PresentMinStep;
	  ti_step = ((P[i].Ti_endstep + ti_grp + ti_step) / ti_step) * ti_step - (P[i].Ti_endstep + ti_grp);
#else

#ifdef PSEUDOSYMMETRIC
	  if(P[i].Type != 0)
	    {
	      if(P[i].Ti_endstep > P[i].Ti_begstep)
		{
		  apred = aphys + ((aphys - P[i].AphysOld) / (P[i].Ti_endstep - P[i].Ti_begstep)) * ti_step;
		  if(fabs(apred - aphys) < 0.5 * aphys)
		    {
		      ti_step2 = get_timestep(i, &apred, -1);
		      ti_min = TIMEBASE;
		      while(ti_min > ti_step2)
			ti_min >>= 1;
		      ti_step2 = ti_min;

		      if(ti_step2 < ti_step)
			{
			  get_timestep(i, &apred, ti_step);
			  prob =
			    ((apred - aphys) / (aphys - P[i].AphysOld) * (P[i].Ti_endstep -
									  P[i].Ti_begstep)) / ti_step;
			  if(prob < get_random_number(P[i].ID))
			    ti_step /= 2;
			}
		      else if(ti_step2 > ti_step)
			{
			  get_timestep(i, &apred, 2 * ti_step);
			  prob =
			    ((apred - aphys) / (aphys - P[i].AphysOld) * (P[i].Ti_endstep -
									  P[i].Ti_begstep)) / ti_step;
			  if(prob < get_random_number(P[i].ID + 1))
			    ti_step *= 2;
			}
		    }
		}
	      P[i].AphysOld = aphys;
	    }
#endif

#ifdef SYNCHRONIZATION
	  if(ti_step > (P[i].Ti_endstep - P[i].Ti_begstep))	/* timestep wants to increase */
	    {
	      if(((TIMEBASE - P[i].Ti_endstep) % ti_step) > 0)
		ti_step = P[i].Ti_endstep - P[i].Ti_begstep;	/* leave at old step */
	    }
#endif
#endif /* end of FLEXSTEPS */


	  if(All.Ti_Current == TIMEBASE)	/* we here finish the last timestep. */
	    ti_step = 0;

	  if((TIMEBASE - All.Ti_Current) < ti_step)	/* check that we don't run beyond the end */
	    ti_step = TIMEBASE - All.Ti_Current;



	  tstart = (P[i].Ti_begstep + P[i].Ti_endstep) / 2;	/* midpoint of old step */
	  tend = P[i].Ti_endstep + ti_step / 2;	/* midpoint of new step */

	  if(All.ComovingIntegrationOn)
	    {
	      dt_entr = (tend - tstart) * All.Timebase_interval;
	      dt_entr2 = (tend - P[i].Ti_endstep) * All.Timebase_interval;
	      dt_gravkick = get_gravkick_factor(tstart, tend);
	      dt_hydrokick = get_hydrokick_factor(tstart, tend);
	      dt_gravkick2 = get_gravkick_factor(P[i].Ti_endstep, tend);
	      dt_hydrokick2 = get_hydrokick_factor(P[i].Ti_endstep, tend);
	    }
	  else
	    {
	      dt_entr = dt_gravkick = dt_hydrokick = (tend - tstart) * All.Timebase_interval;
	      dt_gravkick2 = dt_hydrokick2 = dt_entr2 = (tend - P[i].Ti_endstep) * All.Timebase_interval;
	    }

	  P[i].Ti_begstep = P[i].Ti_endstep;
	  P[i].Ti_endstep = P[i].Ti_begstep + ti_step;



	  /* do the kick */

	  for(j = 0; j < 3; j++)
	    {
	      dv[j] = P[i].g.GravAccel[j] * dt_gravkick;
	      P[i].Vel[j] += dv[j];
	    }

	  if(P[i].Type == 0)	/* SPH stuff */
	    {
	      for(j = 0; j < 3; j++)
		{
		  dv[j] += SphP[i].a.HydroAccel[j] * dt_hydrokick;
		  P[i].Vel[j] += SphP[i].a.HydroAccel[j] * dt_hydrokick;

		  SphP[i].VelPred[j] =
		    P[i].Vel[j] - dt_gravkick2 * P[i].g.GravAccel[j] -
		    dt_hydrokick2 * SphP[i].a.HydroAccel[j];
#ifdef PMGRID
		  SphP[i].VelPred[j] += P[i].GravPM[j] * dt_gravkickB;
#endif
#ifdef MAGNETIC
		  SphP[i].B[j] += SphP[i].DtB[j] * dt_entr;
		  SphP[i].BPred[j] = SphP[i].B[j] - SphP[i].DtB[j] * dt_entr2;
#endif
		}
#if defined(MAGNETIC) && defined(DIVBCLEANING_DEDNER)
	      SphP[i].Phi += SphP[i].DtPhi * dt_entr;
	      SphP[i].PhiPred = SphP[i].Phi - SphP[i].DtPhi * dt_entr2;
#endif
#ifdef TIME_DEP_ART_VISC
	      SphP[i].alpha += SphP[i].Dtalpha * dt_entr;
	      SphP[i].alpha = DMIN(SphP[i].alpha, All.ArtBulkViscConst);
	      if(SphP[i].alpha < All.AlphaMin)
		SphP[i].alpha = All.AlphaMin;
#endif
#ifdef TIME_DEP_MAGN_DISP
	      SphP[i].Balpha += SphP[i].DtBalpha * dt_entr;
	      SphP[i].Balpha = DMIN(SphP[i].Balpha, All.ArtMagDispConst);
	      if(SphP[i].Balpha < All.ArtMagDispMin)
		SphP[i].Balpha = All.ArtMagDispMin;
#endif
	      /* In case of cooling, we prevent that the entropy (and
	         hence temperature decreases by more than a factor 0.5 */

	      if(SphP[i].e.DtEntropy * dt_entr > -0.5 * SphP[i].Entropy)
		SphP[i].Entropy += SphP[i].e.DtEntropy * dt_entr;
	      else
		SphP[i].Entropy *= 0.5;




#ifdef CHEMISTRY
	      /* update the chemical abundances for the new density and temperature */
	      a_start = All.TimeBegin * exp(P[i].Ti_begstep * All.Timebase_interval);
	      a_end = All.TimeBegin * exp(P[i].Ti_endstep * All.Timebase_interval);

	      /* time in cosmic expansion parameter */
	      ifunc = compute_abundances(mode = 1, i, a_start, a_end);
#endif

	      if(All.MinEgySpec)
		{
		  minentropy = All.MinEgySpec * GAMMA_MINUS1 / pow(SphP[i].a2.Density * a3inv, GAMMA_MINUS1);
		  if(SphP[i].Entropy < minentropy)
		    {
		      SphP[i].Entropy = minentropy;
		      SphP[i].e.DtEntropy = 0;
		    }
		}

	      /* In case the timestep increases in the new step, we
	         make sure that we do not 'overcool' when deriving
	         predicted temperatures. The maximum timespan over
	         which prediction can occur is ti_step/2, i.e. from
	         the middle to the end of the current step */

	      dt_entr = ti_step / 2 * All.Timebase_interval;
	      if(SphP[i].Entropy + SphP[i].e.DtEntropy * dt_entr < 0.5 * SphP[i].Entropy)
		SphP[i].e.DtEntropy = -0.5 * SphP[i].Entropy / dt_entr;
	    }
	}
    }



#ifdef PMGRID
  if(All.PM_Ti_endstep == All.Ti_Current)	/* need to do long-range kick */
    {
      ti_step = TIMEBASE;
      while(ti_step > (dt_displacement / All.Timebase_interval))
	ti_step >>= 1;

      if(ti_step > (All.PM_Ti_endstep - All.PM_Ti_begstep))	/* PM-timestep wants to increase */
	{
	  /* we only increase if an integer number of steps will bring us to the end */
	  if(((TIMEBASE - All.PM_Ti_endstep) % ti_step) > 0)
	    ti_step = All.PM_Ti_endstep - All.PM_Ti_begstep;	/* leave at old step */
	}

      if(All.Ti_Current == TIMEBASE)	/* we here finish the last timestep. */
	ti_step = 0;

      tstart = (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2;
      tend = All.PM_Ti_endstep + ti_step / 2;

      if(All.ComovingIntegrationOn)
	dt_gravkick = get_gravkick_factor(tstart, tend);
      else
	dt_gravkick = (tend - tstart) * All.Timebase_interval;

      All.PM_Ti_begstep = All.PM_Ti_endstep;
      All.PM_Ti_endstep = All.PM_Ti_begstep + ti_step;

      if(All.ComovingIntegrationOn)
	dt_gravkickB = -get_gravkick_factor(All.PM_Ti_begstep, (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2);
      else
	dt_gravkickB =
	  -((All.PM_Ti_begstep + All.PM_Ti_endstep) / 2 - All.PM_Ti_begstep) * All.Timebase_interval;

      for(i = 0; i < NumPart; i++)
	{
	  for(j = 0; j < 3; j++)	/* do the kick */
	    P[i].Vel[j] += P[i].GravPM[j] * dt_gravkick;

	  if(P[i].Type == 0)
	    {
	      if(All.ComovingIntegrationOn)
		{
		  dt_gravkickA = get_gravkick_factor(P[i].Ti_begstep, All.Ti_Current) -
		    get_gravkick_factor(P[i].Ti_begstep, (P[i].Ti_begstep + P[i].Ti_endstep) / 2);
		  dt_hydrokick = get_hydrokick_factor(P[i].Ti_begstep, All.Ti_Current) -
		    get_hydrokick_factor(P[i].Ti_begstep, (P[i].Ti_begstep + P[i].Ti_endstep) / 2);
		}
	      else
		dt_gravkickA = dt_hydrokick =
		  (All.Ti_Current - (P[i].Ti_begstep + P[i].Ti_endstep) / 2) * All.Timebase_interval;

	      for(j = 0; j < 3; j++)
		SphP[i].VelPred[j] = P[i].Vel[j]
		  + P[i].g.GravAccel[j] * dt_gravkickA
		  + SphP[i].a.HydroAccel[j] * dt_hydrokick + P[i].GravPM[j] * dt_gravkickB;
	    }
	}
    }
#endif

  t1 = second();
  All.CPU_TimeLine += timediff(t0, t1);
  CPU_Step[CPU_TIMELINE] += timediff(t0, t1);
}




/*! This function normally (for flag==0) returns the maximum allowed timestep of a particle, expressed in
 *  terms of the integer mapping that is used to represent the total simulated timespan. The physical
 *  acceleration is returned in aphys. The latter is used in conjunction with the PSEUDOSYMMETRIC integration
 *  option, which also makes of the second function of get_timestep. When it is called with a finite timestep
 *  for flag, it returns the physical acceleration that would lead to this timestep, assuming timestep
 *  criterion 0.
 */
int get_timestep(int p,		/*!< particle index */
		 double *aphys,	/*!< acceleration (physical units) */
		 int flag	/*!< either 0 for normal operation, or finite timestep to get corresponding
				   aphys */ )
{
  double ax, ay, az, ac;
  double csnd = 0, dt = 0, dt_courant = 0;
  int ti_step;

#ifdef BLACK_HOLES
  double dt_accr;
#endif

#ifdef CONDUCTION
  double dt_cond;
#endif

#ifdef NS_TIMESTEP
  double dt_NS = 0;
#endif

#ifdef CR_DIFFUSION
  double dt_CR;
#endif

#ifdef NONEQUILIBRIUM
  double dt_cool, dt_elec;
#endif

  if(flag == 0)
    {
      ax = fac1 * P[p].g.GravAccel[0];
      ay = fac1 * P[p].g.GravAccel[1];
      az = fac1 * P[p].g.GravAccel[2];

#ifdef PMGRID
      ax += fac1 * P[p].GravPM[0];
      ay += fac1 * P[p].GravPM[1];
      az += fac1 * P[p].GravPM[2];
#endif

      if(P[p].Type == 0)
	{
	  ax += fac2 * SphP[p].a.HydroAccel[0];
	  ay += fac2 * SphP[p].a.HydroAccel[1];
	  az += fac2 * SphP[p].a.HydroAccel[2];
	}

      ac = sqrt(ax * ax + ay * ay + az * az);	/* this is now the physical acceleration */
      *aphys = ac;
    }
  else
    ac = *aphys;

  if(ac == 0)
    ac = 1.0e-30;


  switch (All.TypeOfTimestepCriterion)
    {
    case 0:
      if(flag > 0)
	{
	  dt = flag * All.Timebase_interval;

	  dt /= hubble_a;	/* convert dloga to physical timestep  */

	  ac = 2 * All.ErrTolIntAccuracy * atime * All.SofteningTable[P[p].Type] / (dt * dt);
	  *aphys = ac;
	  return flag;
	}
      dt = sqrt(2 * All.ErrTolIntAccuracy * atime * All.SofteningTable[P[p].Type] / ac);
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
#ifdef ADAPTIVE_GRAVSOFT_FORGAS_HSML
      if(P[p].Type == 0)
	dt = sqrt(2 * All.ErrTolIntAccuracy * atime * PPP[p].Hsml * ADAPTIVE_GRAVSOFT_FORGAS_HSML / 2.8 / ac);
#else
      if(P[p].Type == 0)
	dt =
	  sqrt(2 * All.ErrTolIntAccuracy * atime * All.SofteningTable[P[p].Type] *
	       pow(P[p].Mass / All.ReferenceGasMass, 1.0 / 3) / ac);
#endif
#endif
      break;

    default:
      endrun(888);
      break;
    }


  if(P[p].Type == 0)
    {
      csnd = sqrt(GAMMA * SphP[p].Pressure / SphP[p].a2.Density);

      if(All.ComovingIntegrationOn)
	dt_courant = 2 * All.CourantFac * All.Time * PPP[p].Hsml / (fac3 * SphP[p].MaxSignalVel);
      else
	dt_courant = 2 * All.CourantFac * PPP[p].Hsml / SphP[p].MaxSignalVel;

      if(dt_courant < dt)
	dt = dt_courant;

#ifdef MYFALSE
      dt_viscous = All.CourantFac * SphP[p].MaxViscStep / hubble_a;	/* to convert dloga to physical dt */

      if(dt_viscous < dt)
	dt = dt_viscous;
#endif

#ifdef NS_TIMESTEP
      if(fabs(SphP[p].ViscEntropyChange))
	{
	  dt_NS = VISC_TIMESTEP_PARAMETER * SphP[p].Entropy / 
	    SphP[p].ViscEntropyChange / hubble_a;
 	  
	  if(dt_NS < dt)
	    dt = dt_NS;
	}
#endif

#ifdef CONDUCTION
      if(fabs(SphP[p].CondEnergyChange))
	{
	  dt_cond =
	    COND_TIMESTEP_PARAMETER * SphP[p].Entropy / GAMMA_MINUS1 * pow(SphP[p].a2.Density * a3inv,
									   GAMMA_MINUS1) /
	    fabs(SphP[p].CondEnergyChange);

	  /* convert from dloga back to dt */
	  dt_cond /= hubble_a;

	  if(dt_cond < dt)
	    dt = dt_cond;
	}
#endif

#if (defined(COSMIC_RAYS) && defined(COOLING) && defined(CR_TIMESTEP))
      if(SphP[p].CR_C0 != 0.0)
	{
	  CR_Particle_GetThermalizationRate(SphP + p, &dE_CR, &dN_CR);

	  dt_CR = 0.1 * (SphP[p].CR_E0 / dE_CR);


	  if((dt_CR < dt) && (dt_CR > 0.0))
	    {
	      dt = dt_CR;
	    }

	  dt_CR = 0.1 * (SphP[p].CR_n0 / dN_CR);

	  if((dt_CR < dt) && (dt_CR > 0.0))
	    {
	      dt = dt_CR;
	    }
	}
#endif


#ifdef CR_DIFFUSION
      /* For the diffusion timestep criterion, we look at the green's
       * function for diffusive effects
       * We find that in the exp, we have something like
       * r^2 / 4 alpha t
       * when alpha is the diffusivity.
       * Things should be safe if this factor is smaller than one.
       *
       * The approx. value of r is of the order of (m/rho)^1/3
       *
       */
      if(All.CR_DiffusionCoeff > 0.0)
	{
	  dt_CR = All.CR_DiffusionTimeScale * pow(P[p].Mass / (SphP[p].a2.Density * a3inv), 2.0 / 3.0) / All.CR_DiffusionCoeff;

	  if(All.CR_DiffusionDensScaling != 0.0)
	    {
	      dt_CR *=
		pow(SphP[p].a2.Density * a3inv / All.CR_DiffusionDensZero, -All.CR_DiffusionDensScaling);
	    }

	  if(All.CR_DiffusionEntropyScaling != 0.0)
	    {
	      dt_CR *= pow(SphP[p].Entropy / All.CR_DiffusionEntropyZero, -All.CR_DiffusionEntropyScaling);
	    }


	  if(dt_CR < dt)
	    {
	      dt = dt_CR;

#ifdef CR_DIFFUSION_VERBOSE
	      printf("CR_DIFF timestep had to be applied for particle %i.\n"
		     "dt = %g -- rho_phys = %g -- entropy = %g -- a = %g\n",
		     dt_CR, SphP[p].a2.Density * a3inv, SphP[p].Entropy, atime);
#endif
	    }

	}
#endif




    }

#ifdef BLACK_HOLES
  if(P[p].Type == 5)
    {
      if(P[p].BH_Mdot > 0 && P[p].BH_Mass > 0)
	{
	  dt_accr = 0.25 * P[p].BH_Mass / P[p].BH_Mdot;
	  if(dt_accr < dt)
	    dt = dt_accr;
	}
    }
#endif

#ifdef NONEQUILIBRIUM
  /* another criterion given by the local cooling time */

  if(P[p].Type == 0)
    {
      dt_cool = fabs(SphP[p].t_cool);	/* still in yrs */
      dt_cool *= SEC_PER_YEAR;	/* in seconds */
      dt_cool /= All.UnitTime_in_s;
      dt_cool *= All.HubbleParam;	/* internal units */

      dt_cool = All.Epsilon * dt_cool;


      if(dt_cool > 0 && dt_cool < dt)
	dt = dt_cool;


      /* yet another criterion given by the electron number density change */

      dt_elec = fabs(SphP[p].t_elec);	/* still in yrs */
      dt_elec *= SEC_PER_YEAR;	/* in seconds */
      dt_elec /= All.UnitTime_in_s;
      dt_elec *= All.HubbleParam;	/* internal units */

      dt_elec = All.Epsilon * dt_elec;

      if(dt_elec > 0 && dt_elec < dt)
	dt = dt_elec;
    }
#endif



  /* convert the physical timestep to dloga if needed. Note: If comoving integration has not been selected,
     hubble_a=1.
   */
  dt *= hubble_a;

#ifdef ONLY_PM
  dt = All.MaxSizeTimestep;
#endif



  if(dt >= All.MaxSizeTimestep)
    dt = All.MaxSizeTimestep;


  if(dt >= dt_displacement)
    dt = dt_displacement;

  if(dt < All.MinSizeTimestep)
    {
#ifndef NOSTOP_WHEN_BELOW_MINTIMESTEP
      printf("warning: Timestep wants to be below the limit `MinSizeTimestep'\n");

      if(P[p].Type == 0)
	{
	  printf
	    ("Part-ID=%d  dt=%g dtc=%g ac=%g xyz=(%g|%g|%g)  hsml=%g  maxcsnd=%g dt0=%g eps=%g\n",
	     (int) P[p].ID, dt, dt_courant * hubble_a, ac, P[p].Pos[0], P[p].Pos[1], P[p].Pos[2], PPP[p].Hsml,
	     csnd,
	     sqrt(2 * All.ErrTolIntAccuracy * atime * All.SofteningTable[P[p].Type] / ac) * hubble_a,
	     All.SofteningTable[P[p].Type]);
#ifdef CONDUCTION
	  printf
	    ("Part-ID=%d  dt_cond=%g  A=%g  rho=%g  condench=%g  dtold=%g\n",
	     (int) P[p].ID, dt_cond * hubble_a, SphP[p].Entropy, SphP[p].a2.Density,
	     SphP[p].CondEnergyChange, (P[p].Ti_endstep - P[p].Ti_begstep) * All.Timebase_interval);
#endif

#ifdef NS_TIMESTEP
	  printf
	    ("Part-ID=%d  dt_NS=%g  A=%g  rho=%g  dotAvisc=%g  dtold=%g, meanpath=%g \n",
	     (int) P[p].ID, dt_NS * hubble_a, SphP[p].Entropy, SphP[p].a2.Density,
	     SphP[p].ViscEntropyChange, (P[p].Ti_endstep - P[p].Ti_begstep) * All.Timebase_interval, All.IonMeanFreePath * pow((SphP[p].Entropy * pow(SphP[p].a2.Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1), 2.0) / SphP[p].a2.Density);

	  printf("Stressd=(%g|%g|%g) \n", SphP[p].u.s.StressDiag[0],SphP[p].u.s.StressDiag[1],SphP[p].u.s.StressDiag[2]);
	  printf("Stressoffd=(%g|%g|%g) \n", SphP[p].u.s.StressOffDiag[0],SphP[p].u.s.StressOffDiag[1],SphP[p].u.s.StressOffDiag[2]);
#endif


	}
      else
	{
	  printf("Part-ID=%d  dt=%g ac=%g xyz=(%g|%g|%g)\n", (int) P[p].ID, dt, ac, P[p].Pos[0], P[p].Pos[1],
		 P[p].Pos[2]);
	}
      fflush(stdout);
      endrun(888);
#endif
      dt = All.MinSizeTimestep;
    }

  ti_step = (int) (dt / All.Timebase_interval);

#ifdef CHEMISTRY
  if(ti_step == 0)
    {
      printf("\nError: A timestep of size zero was assigned on the integer timeline!\n"
	     "We better stop.\n"
	     "Task=%d Part-ID=%d dt=%g dt_elec=%g dt_cool=%g tibase=%g ti_step=%d ac=%g xyz=(%g|%g|%g)\n\n",
	     ThisTask, P[p].ID, dt, SphP[p].t_elec, SphP[p].t_cool, All.Timebase_interval, ti_step, ac,
	     P[p].Pos[0], P[p].Pos[1], P[p].Pos[2]);
      fflush(stdout);
      endrun(818);
    }
#endif


  if(!(ti_step > 0 && ti_step < TIMEBASE))
    {
      printf("\nError: A timestep of size zero was assigned on the integer timeline!\n"
	     "We better stop.\n"
	     "Task=%d Part-ID=%d dt=%g dtc=%g dtv=%g dtdis=%g tibase=%g ti_step=%d ac=%g xyz=(%g|%g|%g) tree=(%g|%g|%g)\n\n",
	     ThisTask, (int) P[p].ID, dt, dt_courant, dt, dt_displacement,
	     All.Timebase_interval, ti_step, ac,
	     P[p].Pos[0], P[p].Pos[1], P[p].Pos[2], P[p].g.GravAccel[0], P[p].g.GravAccel[1],
	     P[p].g.GravAccel[2]);
#ifdef PMGRID
      printf("pm_force=(%g|%g|%g)\n", P[p].GravPM[0], P[p].GravPM[1], P[p].GravPM[2]);
#endif
      if(P[p].Type == 0)
	printf("hydro-frc=(%g|%g|%g)\n", SphP[p].a.HydroAccel[0], SphP[p].a.HydroAccel[1],
	       SphP[p].a.HydroAccel[2]);

#ifdef COSMIC_RAYS
      if(P[p].Type == 0)
	printf("Cosmic Ray Properties: C0: %g -- q0  : %g -- P  : %g\n"
	       "                       Rho: %g\n",
	       SphP[p].CR_C0, SphP[p].CR_q0, CR_Particle_Pressure(SphP + p), SphP[p].a2.Density);
#endif

      fflush(stdout);
      endrun(818);
    }

  return ti_step;
}


/*! This function computes an upper limit ('dt_displacement') to the global timestep of the system based on
 *  the rms velocities of particles. For cosmological simulations, the criterion used is that the rms
 *  displacement should be at most a fraction MaxRMSDisplacementFac of the mean particle separation. Note that
 *  the latter is estimated using the assigned particle masses, separately for each particle type. If comoving
 *  integration is not used, the function imposes no constraint on the timestep.
 */
void find_dt_displacement_constraint(double hfac /*!<  should be  a^2*H(a)  */ )
{
  int i, j, type, *temp;
  int count[6];
  long long count_sum[6];
  double v[6], v_sum[6], mim[6], min_mass[6];
  double dt, dmean, asmth = 0;

  dt_displacement = All.MaxSizeTimestep;

  if(All.ComovingIntegrationOn)
    {
      for(type = 0; type < 6; type++)
	{
	  count[type] = 0;
	  v[type] = 0;
	  mim[type] = 1.0e30;
	}

      for(i = 0; i < NumPart; i++)
	{
	  v[P[i].Type] += P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2];
	  if(mim[P[i].Type] > P[i].Mass)
	    mim[P[i].Type] = P[i].Mass;
	  count[P[i].Type]++;
	}

      MPI_Allreduce(v, v_sum, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(mim, min_mass, 6, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

      temp = (int *) mymalloc(NTask * 6 * sizeof(int));
      MPI_Allgather(count, 6, MPI_INT, temp, 6, MPI_INT, MPI_COMM_WORLD);
      for(i = 0; i < 6; i++)
	{
	  count_sum[i] = 0;
	  for(j = 0; j < NTask; j++)
	    count_sum[i] += temp[j * 6 + i];
	}
      myfree(temp);

#ifdef SFR
      /* add star and gas particles together to treat them on equal footing, using the original gas particle
         spacing. */
      v_sum[0] += v_sum[4];
      count_sum[0] += count_sum[4];
#ifdef BLACK_HOLES
      v_sum[0] += v_sum[5];
      count_sum[0] += count_sum[5];
      v_sum[5] = v_sum[0];
      count_sum[5] = count_sum[0];
      min_mass[5] = All.OrigGasMass;
#endif
      v_sum[4] = v_sum[0];
      count_sum[4] = count_sum[0];
      min_mass[4] = All.OrigGasMass;
      min_mass[0] = All.OrigGasMass;
#endif

      for(type = 0; type < 6; type++)
	{
	  if(count_sum[type] > 0)
	    {
	      if(type == 0 || (type == 4 && All.StarformationOn))
		dmean =
		  pow(min_mass[type] / (All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)),
		      1.0 / 3);
	      else
		dmean =
		  pow(min_mass[type] /
		      ((All.Omega0 - All.OmegaBaryon) * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)),
		      1.0 / 3);

	      dt = All.MaxRMSDisplacementFac * hfac * dmean / sqrt(v_sum[type] / count_sum[type]);

#ifdef PMGRID
	      asmth = All.Asmth[0];
#ifdef PLACEHIGHRESREGION
	      if(((1 << type) & (PLACEHIGHRESREGION)))
		asmth = All.Asmth[1];
#endif
	      if(asmth < dmean)
		dt = All.MaxRMSDisplacementFac * hfac * asmth / sqrt(v_sum[type] / count_sum[type]);
#endif

	      if(ThisTask == 0)
		printf("type=%d  dmean=%g asmth=%g minmass=%g a=%g  sqrt(<p^2>)=%g  dlogmax=%g\n",
		       type, dmean, asmth, min_mass[type], All.Time, sqrt(v_sum[type] / count_sum[type]), dt);

	      if(dt < dt_displacement)
		dt_displacement = dt;
	    }
	}

      if(ThisTask == 0)
	printf("displacement time constraint: %g  (%g)\n", dt_displacement, All.MaxSizeTimestep);
    }
}
