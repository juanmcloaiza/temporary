#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "allvars.h"
#include "proto.h"
#include "forcetree.h"
#ifdef COSMIC_RAYS
#include "cosmic_rays.h"
#endif
#include "bg_cooling.h"

#ifdef BG_SFR

/*
 * This routine does cooling and star formation for
 * the STELLA project on BlueGene.
 */


void cooling_and_starformation(void)
/* cooling routine when star formation is enabled */
{
  int i, j, k, flag;
  int stars_spawned, tot_spawned, stars_converted, tot_converted, number_of_stars_generated;
  unsigned int bits;
  double dt, dtime, ascale = 1, hubble_a = 0, a3inv;
  double time_hubble_a, unew, mass_of_star, frac;
  double sum_sm, total_sm, rate, sum_mass_stars, total_sum_mass_stars;
  double p, prob, egycurrent, rate_in_msunperyear;
  double sfrrate, totsfrrate, dmax1, dmax2, u_to_temp_fac;
  double rho, redshift, pressure, element_metallicity[BG_NELEMENTS];

  static int z_index_old = -3;
  float d_z;

  u_to_temp_fac = (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS / BOLTZMANN * GAMMA_MINUS1
    * All.UnitEnergy_in_cgs / All.UnitMass_in_g;


  if(All.ComovingIntegrationOn)
    {
      /* Factors for comoving integration of hydro */
      a3inv = 1 / (All.Time * All.Time * All.Time);
      hubble_a = hubble_function(All.Time);

      time_hubble_a = All.Time * hubble_a;
      ascale = All.Time;
      redshift = 1 / ascale - 1;
    }
  else
    {
      a3inv = ascale = time_hubble_a = 1;
      redshift = 0;
    }

  stars_spawned = stars_converted = 0;
  sum_sm = sum_mass_stars = 0;

  for(bits = 0; GENERATIONS > (1 << bits); bits++);

  for(i = 0; i < N_gas; i++)
    {
      if(P[i].Type == 0)
	if(P[i].Ti_endstep == All.Ti_Current)
	  {
	    dt = (P[i].Ti_endstep - P[i].Ti_begstep) * All.Timebase_interval;
	    /*  the actual time-step */

	    if(All.ComovingIntegrationOn)
	      dtime = All.Time * dt / time_hubble_a;
	    else
	      dtime = dt;

	    /* check whether conditions for star formation are fulfilled.
	     *  
	     * f=1  normal cooling
	     * f=0  star formation
	     */
	    flag = 1;		/* default is normal cooling */

	    if(SphP[i].a2.Density * a3inv >= All.PhysDensThresh)
	      flag = 0;

	    if(All.ComovingIntegrationOn)
	      if(SphP[i].a2.Density < All.OverDensThresh)
		flag = 1;


	    if(flag == 1)	/* normal implicit isochoric cooling */
	      {
		SphP[i].Sfr = 0;

#ifndef BG_COOLING

		ne = SphP[i].Ne;	/* electron abundance (gives ionization state and mean molecular weight) */

		unew = DMAX(All.MinEgySpec,
			    (SphP[i].Entropy + SphP[i].e.DtEntropy * dt) /
			    GAMMA_MINUS1 * pow(SphP[i].a2.Density * a3inv, GAMMA_MINUS1));

		unew = DoCooling(unew, SphP[i].a2.Density * a3inv, dtime, &ne);
		SphP[i].Ne = ne;

#else /* BG_COOLING */

		/* compute element metallicity as element mass fractions */
		for(j = 0; j < BG_NELEMENTS; j++)
		  element_metallicity[j] = P[i].Metals[j] / P[i].Mass;

		get_redshift_index(redshift, z_index_old, &d_z);

		if(cooling_redshift_index < 0 && z_index_old == -3)
		  {
		    if(ThisTask == 0)
		      GetCollisTable();

		    BroadcastCoolingTables();
		  }
		else if(cooling_redshift_index != z_index_old)
		  {
		    if(ThisTask == 0)
		      GetCoolingTables(z_index_old);

		    BroadcastCoolingTables();
		  }

		z_index_old = cooling_redshift_index;

		/* convert to CGS units for the cooling routine */
		unew = DMAX(All.MinEgySpec,
			    (SphP[i].Entropy + SphP[i].e.DtEntropy * dt) /
			    GAMMA_MINUS1 * pow(SphP[i].a2.Density * a3inv, GAMMA_MINUS1));
		unew *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;

		rho = SphP[i].a2.Density * a3inv;
		rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;

		dtime *= All.UnitTime_in_s / All.HubbleParam;

		/* do the cooling */
		unew = DoCooling(unew, rho, dtime, redshift, d_z, element_metallicity);

		/* convert energy back to code units */
		unew *= All.UnitDensity_in_cgs / All.UnitPressure_in_cgs;

#endif /* BG_COOLING */

		if(P[i].Ti_endstep > P[i].Ti_begstep)	/* upon start-up, we need to protect against dt==0 */
		  {
		    /* note: the adiabatic rate has been already added in ! */

		    if(dt > 0)
		      {
			SphP[i].e.DtEntropy = (unew * GAMMA_MINUS1 /
					       pow(SphP[i].a2.Density * a3inv,
						   GAMMA_MINUS1) - SphP[i].Entropy) / dt;

			if(SphP[i].e.DtEntropy < -0.5 * SphP[i].Entropy / dt)
			  SphP[i].e.DtEntropy = -0.5 * SphP[i].Entropy / dt;
		      }
		  }
	      }

	    if(flag == 0)	/* active star formation  - we are governed by an effective equation of state */
	      {
		egycurrent = TEMP_THRESH / u_to_temp_fac;

		pressure = GAMMA_MINUS1 * All.PhysDensThresh * egycurrent *
		  pow(SphP[i].a2.Density * a3inv / All.PhysDensThresh, All.GammaEffective);

		/* convert pressure to comoving pressure */
		pressure /= pow(a3inv, GAMMA);

		SphP[i].Entropy = pressure / pow(SphP[i].a2.Density, GAMMA);
		SphP[i].e.DtEntropy = 0;


		/* convert back to physical pressure */
		pressure *= pow(a3inv, GAMMA);
		/* convert to cgs units */
		pressure *= All.UnitPressure_in_cgs;

		p = KENNICUTT_NORM * KENNICUTT_COEFF *
		  pow(pressure * GAMMA / GRAVITY, (KENNICUTT_EXP - 1) / 2) * dtime * All.UnitTime_in_s;

		if(p > 1.0)
		  printf("SFR probability greater than 1!!! p = %g", p);

		sum_sm += P[i].Mass * (1 - exp(-p));

		/* the upper bits of the gas particle ID store how man stars this gas
		   particle gas already generated */

		number_of_stars_generated = (P[i].ID >> (32 - bits));

		mass_of_star = P[i].Mass / (GENERATIONS - number_of_stars_generated);

		SphP[i].Sfr = KENNICUTT_NORM * KENNICUTT_COEFF * P[i].Mass * All.UnitMass_in_g *
		  pow(pressure * GAMMA / GRAVITY, (KENNICUTT_EXP - 1) / 2) * (SEC_PER_YEAR / SOLAR_MASS);


		prob = P[i].Mass / mass_of_star * (1 - exp(-p));


		if(get_random_number(P[i].ID + 1) < prob)	/* ok, make a star */
		  {
		    if(number_of_stars_generated == (GENERATIONS - 1))
		      {
			/* here we turn the gas particle itself into a star */
			Stars_converted++;
			stars_converted++;

			sum_mass_stars += P[i].Mass;

			P[i].Type = 4;

			P[i].StellarAge = All.Time;
		      }
		    else
		      {
			/* here we spawn a new star particle */

			if(NumPart + stars_spawned >= All.MaxPart)
			  {
			    printf
			      ("On Task=%d with NumPart=%d we try to spawn %d particles. Sorry, no space left...(All.MaxPart=%d)\n",
			       ThisTask, NumPart, stars_spawned, All.MaxPart);
			    fflush(stdout);
			    endrun(8888);
			  }

			P[NumPart + stars_spawned] = P[i];
			P[NumPart + stars_spawned].Type = 4;

			P[i].ID += (1 << (32 - bits));

			P[NumPart + stars_spawned].Mass = mass_of_star;

			frac = mass_of_star / P[i].Mass;

			for(k = 0; k < BG_NELEMENTS; k++)
			  {
			    P[NumPart + stars_spawned].Metals[k] = frac * P[i].Metals[k];
			    P[i].Metals[k] -= P[NumPart + stars_spawned].Metals[k];
			  }

			P[i].Mass -= P[NumPart + stars_spawned].Mass;
			sum_mass_stars += P[NumPart + stars_spawned].Mass;

			P[NumPart + stars_spawned].StellarAge = All.Time;

			force_add_star_to_tree(i, NumPart + stars_spawned);

			stars_spawned++;
		      }
		  }

	      }
	  }

    }				/* end of main loop over active particles */


  MPI_Allreduce(&stars_spawned, &tot_spawned, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&stars_converted, &tot_converted, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(tot_spawned > 0 || tot_converted > 0)
    {
      if(ThisTask == 0)
	{
	  printf("\n----> spawned %d stars, converted %d gas particles into stars\n\n",
		 tot_spawned, tot_converted);
	  fflush(stdout);
	}


      All.TotNumPart += tot_spawned;
      All.TotN_gas -= tot_converted;
      NumPart += stars_spawned;
      NumForceUpdate += stars_spawned;

      /* Note: N_gas is only reduced once rearrange_particle_sequence is called */

      /* Note: New tree construction can be avoided because of  `force_add_star_to_tree()' */
    }

  for(i = 0, sfrrate = 0; i < N_gas; i++)
    if(P[i].Type == 0)
      sfrrate += SphP[i].Sfr;

  MPI_Allreduce(&sfrrate, &totsfrrate, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  MPI_Reduce(&sum_sm, &total_sm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sum_mass_stars, &total_sum_mass_stars, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(ThisTask == 0)
    {
      if(All.TimeStep > 0)
	rate = total_sm / (All.TimeStep / time_hubble_a);
      else
	rate = 0;

      /* convert to solar masses per yr */

      rate_in_msunperyear = rate * (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

      fprintf(FdSfr, "%g %g %g %g %g\n", All.Time, total_sm, totsfrrate, rate_in_msunperyear,
	      total_sum_mass_stars);
      fflush(FdSfr);
    }
}


void set_units_sfr(void)
{
  All.OverDensThresh =
    All.CritOverDensity * All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

  All.PhysDensThresh = All.CritPhysDensity * PROTONMASS / HYDROGEN_MASSFRAC / All.UnitDensity_in_cgs;
}


double get_starformation_rate(int i)
{
  int flag;
  double rateOfSF, a3inv, pressure, egycurrent, u_to_temp_fac;

  u_to_temp_fac = (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS / BOLTZMANN * GAMMA_MINUS1
    * All.UnitEnergy_in_cgs / All.UnitMass_in_g;

  if(All.ComovingIntegrationOn)
    a3inv = 1 / (All.Time * All.Time * All.Time);
  else
    a3inv = 1;

  flag = 1;			/* default is normal cooling */

  if(SphP[i].a2.Density * a3inv >= All.PhysDensThresh)
    flag = 0;

  if(All.ComovingIntegrationOn)
    if(SphP[i].a2.Density < All.OverDensThresh)
      flag = 1;

  if(flag == 1)
    return 0;

  egycurrent = TEMP_THRESH / u_to_temp_fac;

  /* compute to physical pressure */
  pressure = GAMMA_MINUS1 * All.PhysDensThresh * egycurrent *
    pow(SphP[i].a2.Density * a3inv / All.PhysDensThresh, All.GammaEffective);
  /* convert to cgs units */
  pressure *= All.UnitPressure_in_cgs;

  rateOfSF = KENNICUTT_NORM * KENNICUTT_COEFF * P[i].Mass * All.UnitMass_in_g *
    pow(pressure * GAMMA / GRAVITY, (KENNICUTT_EXP - 1) / 2) * (SEC_PER_YEAR / SOLAR_MASS);

  return rateOfSF;
}

#endif /* closing of BG_SFR-conditional */
