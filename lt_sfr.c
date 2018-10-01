#ifdef LT_STELLAREVOLUTION

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "allvars.h"
#include "proto.h"
#include "forcetree.h"
#ifdef COSMIC_RAYS
#include "cosmic_rays.h"
#endif

#ifdef LT_SEv_INFO
/* temporary */
/* static int SEvInfo_grain; */
/* */
double mZ[LT_NMet], tot_mZ[LT_NMet];

#define STAT_INUM  5
#define ERatio     0
#define EFRatio    1
#define XRatio     2
#define OEgy       3
#define NoREFratio 4

#define REScount   5
#define NoREScount 6

double StatSum[STAT_INUM + 2], tot_StatSum[STAT_INUM + 2];
double StatMin[STAT_INUM], tot_StatMin[STAT_INUM];
double StatMax[STAT_INUM], tot_StatMax[STAT_INUM];

#if defined(LT_EXTEGY_INFO) && defined(LT_MOD_EFFM)
static int EEInfo_grain;
#endif

#endif

#ifdef LT_CharT_INFO
static int CharTInfo_grain;
int ncool, nsf;
int tot_ncool, tot_nsf;
double tcool_tcross[22], tcool_tff[22];
double *tot_tcool_tcross, *tot_tcool_tff, tot_tc_tcr[22], tot_tc_tff[22], tot_div;
#endif

double CPU_eff_iter, CPU_ee_info, CPU_sn_info;
double sumCPU_eff_iter, sumCPU_ee_info, sumCPU_sn_info;

static double *sum_sm, *sum_mass_stars, *total_sm, *total_sum_mass_stars;

#ifdef LT_MOD_EFFM

#define PERTURB_MAX_ITER 100
#define PERTURB_REL_PREC 0.00005

struct
{
  double rho, PhysDensTh, FEVP, fEVP, EpsilonExt, tsfr, deltat, FactorSN, ehot_old, ne, Z;
} perturbp;

static const gsl_root_fsolver_type *T;
static gsl_root_fsolver *s;
static gsl_function F;


int INLINE_FUNC perturb(double *, double *);


double xfx(double x, void *params)
{
  double factorEVP, egyI, egyII, egyhot, tcool, y, x_;

  if(All.Mod_SEff)
    {
      factorEVP = pow(perturbp.rho / perturbp.PhysDensTh, -0.8) * perturbp.FEVP *
	(1 + perturbp.EpsilonExt * x / perturbp.EgySpecSN * perturbp.tsfr / perturbp.deltat);
    }
  perturbp.fEVP = factorEVP;

  if(x > 0)
    {
      egyI = perturbp.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;
      if(All.Mod_SEff)
	egyII = perturbp.EpsilonExt * (1 - x) /
	  (perturbp.FactorSN * (1 + factorEVP)) * (1 - x) / x * (perturbp.tsfr / perturbp.deltat);
      else
	egyII = perturbp.EpsilonExt * (1 - x) /
	  (perturbp.FactorSN * (1 + factorEVP)) / x * (perturbp.tsfr / perturbp.deltat);
    }
  else
    {
      egyI = (perturbp.EgySpecSN) / (1 + factorEVP) + All.EgySpecCold + perturbp.EpsilonExt;
      egyII = 0;
    }

  egyhot = (egyI + egyII);
  perturbp.ehot_old = egyhot;

#ifdef LT_METAL_COOLING
  if((tcool = GetCoolingTime(egyhot, perturbp.rho, &(perturbp.ne), perturbp.Z)) == 0)
    tcool = 1e13;
#else
  if((tcool = GetCoolingTime(egyhot, perturbp.rho, &(perturbp.ne))) == 0)
    tcool = 1e13;
#endif

  if(x > 0)
    {
      if(All.Mod_SEff)
	y = perturbp.tsfr * egyhot / tcool /
	  ((perturbp.FactorSN * perturbp.EgySpecSN - (1 - perturbp.FactorSN) * All.EgySpecCold) +
	   perturbp.EpsilonExt * perturbp.tsfr / perturbp.deltat * ((1 - x) / x * (1 - x)));
      else
	y = perturbp.tsfr * egyhot / tcool /
	  ((perturbp.FactorSN * All.EgySpecSN - (1 - perturbp.FactorSN) * All.EgySpecCold) +
	   perturbp.EpsilonExt * perturbp.tsfr / perturbp.deltat * ((1 - x) / x));
    }
  else
    y = perturbp.tsfr * egyhot / tcool /
      (perturbp.FactorSN * perturbp.EgySpecSN - (1 - perturbp.FactorSN) * All.EgySpecCold);

  if(y < 1e-2)
    x_ = y;
  else
    x_ = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

  return x_ - x;

}
#endif /* refers to LT_MOD_EFFM */


void cooling_and_starformation(void)	/* cooling routine when star formation is enabled */
{
  int i, flag, stars_spawned, tot_spawned, stars_converted, tot_converted;
  int number_of_stars_generated, bits, Yset, YZbin;
  double dt, dtime, ascale = 1, hubble_a = 0, a3inv, ne = 1;
  double time_hubble_a, unew, mass_of_star;
  double Sum_sm, Sum_mass_stars, Total_sm, Total_sum_mass_stars;
  double sm, rate;
  double p, prob;
  double cloudmass;
  double factorEVP;
  double tsfr, trelax;
  double egyhot, egyeff, egycurrent, tcool, x, y, rate_in_msunperyear;
  double sfrrate, totsfrrate, dmax1, dmax2, dmin1, dmin2;
  int j, IMFi;

#ifdef WINDS
  double v;
  double norm, dir[3];

#ifdef ISOTROPICWINDS
  double theta, phi;
#endif

#ifdef LT_HOT_WINDS
  FLOAT u_hotwinds;
#endif

#ifdef LT_SEv_INFO
  /* temporary */
/*   double all_mass, all_num, tot_allmass, tot_allnum; */
  /* */
  int windn, tot_windn;
  double windv_min, windv_max, windv, tot_windv_min, tot_windv_max, tot_windv;
#endif
#endif

  double myFactorEVP, myPhysDensThresh, myFactorSN, myEgySpecSN;

#if defined(LT_CharT) || defined(LT_CharT_INFO)
  double sound_crossing_time, freefall_time;
#endif
  double SNEgy;
  double mz, mstar, factor;

#ifdef LT_MOD_EFFM
#ifdef LT_EXTEGY_INFO
  double egy_new;
  double xrel, xold;
#endif
  double x_lo, x_hi, xoldold;
  int niter
#endif
  int GENERATIONS;
  double tstart, tend;
  double Z, Zsol, exp_minus_p;
  int ti_step_Ia, ti_step_II;

#ifdef LT_SNII
  double tot_mass;
#endif

#if defined(LT_EXTEGY_INFO) && defined(LT_MOD_EFFM)
  double r, egyhot_old;
#endif

#ifdef COSMIC_RAYS
  double rCR_dE, rCR_dN;
#endif

#if defined(QUICK_LYALPHA) || defined(BH_THERMALFEEDBACK) || defined (BH_KINETICFEEDBACK)
  double temp, u_to_temp_fac;

  u_to_temp_fac = (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS / BOLTZMANN * GAMMA_MINUS1
    * All.UnitEnergy_in_cgs / All.UnitMass_in_g;
#endif

  if(sum_sm == 0x0)
    {
      sum_sm = (double *) malloc(IMFs_dim * sizeof(double));
      total_sm = (double *) malloc(IMFs_dim * sizeof(double));
      sum_mass_stars = (double *) malloc(IMFs_dim * sizeof(double));
      total_sum_mass_stars = (double *) malloc(IMFs_dim * sizeof(double));
    }

  if(All.ComovingIntegrationOn)
    {
      /* Factors for comoving integration of hydro */
      a3inv = 1 / (All.Time * All.Time * All.Time);
      hubble_a = hubble_function(All.Time);

      time_hubble_a = All.Time * hubble_a;
      ascale = All.Time;
    }
  else
    a3inv = ascale = time_hubble_a = 1;



  for(i = 0; i < IMFs_dim; i++)
    {
      sum_sm[i] = sum_mass_stars[i] = 0;
      total_sm[i] = total_sum_mass_stars[i] = 0;
      Sum_sm = Sum_mass_stars = Total_sm = Total_sum_mass_stars = 0;
    }
  stars_spawned = stars_converted = 0;

#ifdef LT_CharT_INFO
  if(++CharTInfo_grain > CharT_Info_GRAIN)
    {
      tot_tcool_tcross = (double *) calloc(NTask * 22, sizeof(double));
      tot_tcool_tff = (double *) calloc(NTask * 22, sizeof(double));

      for(i = 0; i < 22; i++)
	tcool_tcross[i] = tcool_tff[i] = 0;
      ncool = nsf = tot_ncool = tot_nsf = 0;
    }
#endif

#ifdef LT_SEv_INFO

#ifdef LT_SNII
  for(i = 0; i < LT_NMet; i++)
    mZ[i] = tot_mZ[i] = 0;
#endif

#ifdef WINDS
  windv_min = 1e6;
  windv_max = windv = windn = 0;
#endif

#endif

#if defined(LT_EXTEGY_INFO) && defined(LT_MOD_EFFM)
  EEInfo_grain++;
#endif

  /* temporary */
/* #ifdef LT_SEv_INFO */
/*   SEvInfo_grain++; */
/*   if(SEvInfo_grain > SEvInfo_GRAIN) */
/*     {   */
/*       all_mass = all_num = 0; */
/*       for(j = 0; j < Trun_size; j++) */
/* 	Tpop_num[j] = Tpop_mass[j] = 0; */
/*       for(j = 0; j < TCrun_size; j++) */
/* 	TCpop_num[j] = TCpop_mass[j] = 0; */
/*     } */
/* #endif */
  /* */

  for(i = 0; i < N_gas; i++)
    if((P[i].Type == 0) && (P[i].Ti_endstep == All.Ti_Current))
      {
	dt = (P[i].Ti_endstep - P[i].Ti_begstep) * All.Timebase_interval;
	/*  the actual time-step */

	if(All.ComovingIntegrationOn)
	  dtime = All.Time * dt / time_hubble_a;
	else
	  dtime = dt;

#ifndef LT_LOCAL_IRA
	SphP[i].mstar = 0;
#endif

#ifndef LT_SMOOTH_Z
	Z = get_metallicity(i, 2);
#else
	Z = SphP[i].Zsmooth * 0.0885;
#endif

/*         Z = get_metallicity(i, 2); */
/* #ifdef LT_SMOOTH_Z */
/*         if(Z > 0) */
/*           Z = (Z +  SphP[i].Zsmooth * 0.0885) / 2; */
/* #endif */

	Zsol = get_metallicity_solarunits(Z);

	IMFi = get_IMF_index(i);

	GENERATIONS = IMFs[IMFi].Generations;
	for(bits = 0; GENERATIONS > (1 << bits); bits++);

	Yset = IMFs[IMFi].YSet;
	IMFp = &IMFs[IMFi];
	myFactorSN = IMFs[IMFi].FactorSN;
	myEgySpecSN = IMFs[IMFi].EgySpecSN;
	for(YZbin = IIZbins_dim[Yset] - 1; Z < IIZbins[Yset][YZbin] && YZbin > 0; YZbin--)
	  ;

	if(IMFs[IMFi].SFTh_Zdep)
	  {
	    /*
	       if the effective model is allowed to depend on metallicity,
	       thresholds and evaporation factors differ for different Z.
	     */
	    getindex(&Zvalue[0], 0, ZBins - 1, &Zsol, &flag);

	    if(flag == 0 || flag == ZBins - 1)
	      {
		myFactorEVP = IMFs[IMFi].FEVP[flag];
		myPhysDensThresh = IMFs[IMFi].PhysDensThresh[flag];
	      }
	    else
	      {
		p = (Zvalue[flag + 1] - Z) / (Zvalue[flag + 1] - Zvalue[flag]);
		/* interpolate */
		myFactorEVP = IMFs[IMFi].FEVP[flag] * p + IMFs[IMFi].FEVP[flag + 1] * (1 - p);
		myPhysDensThresh =
		  IMFs[IMFi].PhysDensThresh[flag] * p + IMFs[IMFi].PhysDensThresh[flag + 1] * (1 - p);
	      }
	  }
	else
	  {
	    myFactorEVP = IMFs[IMFi].FEVP[0];
	    myPhysDensThresh = IMFs[IMFi].PhysDensThresh[0];
	  }

	/* check whether conditions for star formation are fulfilled.
	 *
	 * f=1  normal cooling
	 * f=0  star formation
	 */
	flag = 1;		/* default is normal cooling */

	if(SphP[i].a2.Density * a3inv >= myPhysDensThresh)
	  flag = 0;

	if(All.ComovingIntegrationOn)
	  if(SphP[i].a2.Density < All.OverDensThresh)
	    flag = 1;


	/* temporary */
/* #ifdef LT_SEv_INFO */
/* 	if(SEvInfo_grain > SEvInfo_GRAIN) */
/* 	  { */
/* 	    egyhot = DMAX(All.MinEgySpec, (SphP[i].Entropy + SphP[i].e.DtEntropy * dt) / */
/*                           GAMMA_MINUS1 * pow(SphP[i].a2.Density * a3inv, GAMMA_MINUS1)); */

/* #ifdef LT_METAL_COOLING */
/* 	    if((tcool = GetCoolingTime(egyhot, SphP[i].a2.Density * a3inv, &ne, Z)) == 0) */
/* 	      tcool = 1e13; */
/* #else */
/* 	    if((tcool = GetCoolingTime(egyhot, SphP[i].a2.Density * a3inv, &ne)) == 0) */
/* 	      tcool = 1e13; */
/* #endif */

/* 	    j = DMAX(DMIN((log10(tcool) + 5.0) / delta_TCrun, TCrun_size-1), 0); */
/* 	    TCpop_num[j]++; */
/* 	    TCpop_mass[j] += P[i].Mass; */

/* 	    egyhot *=  All.UnitEnergy_in_cgs / All.UnitMass_in_g; */
/* 	    j = DMAX(DMIN((LogTemp(egyhot, SphP[i].Ne) - 3.0) / delta_Trun, Trun_size-1), 0); */
/* 	    Tpop_num[j]++; */
/* 	    Tpop_mass[j] += P[i].Mass; */

/* 	    all_mass += P[i].Mass; */
/* 	    all_num++; */
/* 	  } */

/* #endif */
	/* */


#if defined(LT_CharT) || defined(LT_CharT_INFO)

	egyhot = DMAX(All.MinEgySpec, (SphP[i].Entropy + SphP[i].e.DtEntropy * dt) /
		      GAMMA_MINUS1 * pow(SphP[i].a2.Density * a3inv, GAMMA_MINUS1));
#ifdef LT_METAL_COOLING
	if((tcool = GetCoolingTime(egyhot, SphP[i].a2.Density * a3inv, &ne, Z)) == 0)
	  tcool = 1e13;
#else
	if((tcool = GetCoolingTime(egyhot, SphP[i].a2.Density * a3inv, &ne)) == 0)
	  tcool = 1e13;
#endif

	sound_crossing_time = 1.0 / sqrt(GAMMA * SphP[i].Pressure / SphP[i].a2.Density) *
	  pow(P[i].Mass / SphP[i].a2.Density, 1.0 / 3.0);
	freefall_time = 1.0 / sqrt(All.G * SphP[i].a2.Density);

#ifdef LT_CharT_INFO
	if(CharTInfo_grain > CharT_Info_GRAIN)
	  {
	    if((p = log10(SphP[i].a2.Density / 1e-8)) < 0)
	      p = 0;
	    if((p /= (6.0 / 20.0)) > 19)
	      p = 19;
	    tcool_tcross[(int) p] += tcool / sound_crossing_time;
	    tcool_tff[(int) p] += tcool / freefall_time;
	  }
#endif
#ifdef LT_CharT

	if(sound_crossing_time < tcool)
	  flag = 1;
	if(freefall_time > sound_crossing_time)
	  flag = 1;
#endif

#endif


#ifdef BLACK_HOLES
	if(P[i].Mass == 0)
	  flag = 1;
#endif

#ifdef WINDS
	if(SphP[i].DelayTime > 0)
	  flag = 1;		/* only normal cooling for particles in the wind */

	if(SphP[i].DelayTime > 0)
	  SphP[i].DelayTime -= dtime;

	if(SphP[i].DelayTime > 0)
	  if(SphP[i].a2.Density * a3inv < All.WindFreeTravelDensFac * All.PhysDensThresh)
	    SphP[i].DelayTime = 0;

	if(SphP[i].DelayTime < 0)
	  SphP[i].DelayTime = 0;

	/* exclude winds from cooling ?? */
	if(SphP[i].DelayTime > 0)
	  continue;
#endif


#ifdef QUICK_LYALPHA
	temp = u_to_temp_fac * (SphP[i].Entropy + SphP[i].e.DtEntropy * dt) /
	  GAMMA_MINUS1 * pow(SphP[i].a2.Density * a3inv, GAMMA_MINUS1);

	if(SphP[i].a2.Density > All.OverDensThresh && temp < 1.0e5)
	  flag = 0;
	else
	  flag = 1;
#endif


#if !defined(NOISMPRESSURE) && !defined(QUICK_LYALPHA)
	if(flag == 1)		/* normal implicit isochoric cooling */
#endif
	  {
	    SphP[i].Sfr = 0;

	    ne = SphP[i].Ne;	/* electron abundance (gives ionization state and mean molecular weight) */

#ifdef LT_CharT_INFO
	    if(CharTInfo_grain > CharT_Info_GRAIN)
	      {
		ncool++;
		tcool_tcross[20] += tcool / sound_crossing_time;
		tcool_tff[20] += tcool / freefall_time;
	      }
#endif

	    unew = DMAX(All.MinEgySpec,
			(SphP[i].Entropy + SphP[i].e.DtEntropy * dt) /
			GAMMA_MINUS1 * pow(SphP[i].a2.Density * a3inv, GAMMA_MINUS1));

#if defined(BH_THERMALFEEDBACK) || defined(BH_KINETICFEEDBACK)
	    if(SphP[i].i.Injected_BH_Energy > 0)
	      {
		if(P[i].Mass == 0)
		  SphP[i].i.Injected_BH_Energy = 0;
		else
		  unew += SphP[i].i.Injected_BH_Energy / P[i].Mass;

		temp = u_to_temp_fac * unew;


		if(temp > 5.0e9)
		  unew = 5.0e9 / u_to_temp_fac;

		SphP[i].i.Injected_BH_Energy = 0;

	      }
#endif

	    SNEgy = SphP[i].EgyRes / P[i].Mass;
#ifdef LT_METAL_COOLING
	    unew = DoCooling(unew, SphP[i].a2.Density * a3inv, dtime, &ne, Z);
#else
	    unew = DoCooling(unew, SphP[i].a2.Density * a3inv, dtime, &ne);
#endif
	    unew += SNEgy;;

	    SphP[i].Ne = ne;


	    if(P[i].Ti_endstep > P[i].Ti_begstep)	/* upon start-up, we need to protect against dt==0 */
	      {
		/* note: the adiabatic rate has been already added in ! */

		if(dt > 0)
		  {
#ifdef COSMIC_RAYS
		    unew += CR_Particle_ThermalizeAndDissipate(SphP + i, dtime);
#endif /* COSMIC_RAYS */


		    SphP[i].e.DtEntropy = (unew * GAMMA_MINUS1 /
					   pow(SphP[i].a2.Density * a3inv,
					       GAMMA_MINUS1) - SphP[i].Entropy) / dt;

		    if(SphP[i].e.DtEntropy < -0.5 * SphP[i].Entropy / dt)
		      SphP[i].e.DtEntropy = -0.5 * SphP[i].Entropy / dt;

		    SphP[i].EgyRes = 0;
		  }
	      }

#ifdef LT_WINDS_EXT_NOSF_BRANCH
	    /* Here comes the wind model also for non-starforming particles */

	    if(SNEgy > 0.5 * IMFs[IMFi].WindEnergy)
	      {
		p = IMFs[IMFi].WindEfficiency * (SNEgy / myEgySpecSN);
		/*p =  All.WindEnergyFraction * SNEgy / uold; */

		prob = 1 - exp(-p);

		if(drand48() < prob)	/* ok, make the particle go into the wind */
		  {
#ifndef LT_WIND_VELOCITY
		    v = sqrt(2 * SNEgy / IMFs[IMFi].WindEfficiency);
#else
		    v = LT_WIND_VELOCITY;
#endif

#ifdef LT_SEv_INFO

		    if(v < windv_min)
		      windv_min = v;
		    if(v > windv_max)
		      windv_max = v;
		    windv += v;
		    windn++;

#endif

#ifdef ISOTROPICWINDS
		    theta = acos(2 * drand48() - 1);
		    phi = 2 * M_PI * drand48();

		    dir[0] = sin(theta) * cos(phi);
		    dir[1] = sin(theta) * sin(phi);
		    dir[2] = cos(theta);
#else
		    dir[0] = P[i].g.GravAccel[1] * P[i].Vel[2] - P[i].g.GravAccel[2] * P[i].Vel[1];
		    dir[1] = P[i].g.GravAccel[2] * P[i].Vel[0] - P[i].g.GravAccel[0] * P[i].Vel[2];
		    dir[2] = P[i].g.GravAccel[0] * P[i].Vel[1] - P[i].g.GravAccel[1] * P[i].Vel[0];
#endif

		    for(j = 0, norm = 0; j < 3; j++)
		      norm += dir[j] * dir[j];

		    norm = sqrt(norm);
		    if(drand48() < 0.5)
		      norm = -norm;

		    if(norm != 0)
		      {
			for(j = 0; j < 3; j++)
			  dir[j] /= norm;

			for(j = 0; j < 3; j++)
			  {
			    P[i].Vel[j] += v * ascale * dir[j];
			    SphP[i].VelPred[j] += v * ascale * dir[j];
			  }

			SphP[i].DelayTime = DMAX(2*PPP[i].hsml,All.WindFreeTravelLength) / v;
		      }
		  }
	      }
#endif
	  }

	if(flag == 0)		/* active star formation */
	  {
#if !defined(QUICK_LYALPHA)

	    tsfr = sqrt(myPhysDensThresh / (SphP[i].a2.Density * a3inv)) * All.MaxSfrTimescale;
	    factorEVP = pow(SphP[i].a2.Density * a3inv / myPhysDensThresh, -0.8) * myFactorEVP;
	    SNEgy = SphP[i].EgyRes / P[i].Mass;

	    egyhot = myEgySpecSN / (1 + factorEVP) + All.EgySpecCold;

#if defined(LT_EXTEGY_INFO) && defined(LT_MOD_EFFM)
	    egyhot_old = egyhot;
#endif
	    ne = SphP[i].Ne;
#ifdef LT_METAL_COOLING
	    /*Z = get_metallicity(i); already done */
	    if((tcool = GetCoolingTime(egyhot, SphP[i].a2.Density * a3inv, &ne, Z)) == 0)
	      tcool = 1e13;
#else
	    tcool = GetCoolingTime(egyhot, SphP[i].a2.Density * a3inv, &ne);
#endif
	    SphP[i].Ne = ne;

	    y = tsfr / tcool * egyhot / (myFactorSN * myEgySpecSN - (1 - myFactorSN) * All.EgySpecCold);
	    if(y < 1e-4)
	      x = y;
	    else
	      x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

#ifdef LT_CharT_INFO
	    if(CharTInfo_grain > CharT_Info_GRAIN)
	      {
		nsf++;
		tcool_tcross[21] += tcool / sound_crossing_time;
		tcool_tff[21] += tcool / freefall_time;
	      }
#endif

#ifdef LT_MOD_EFFM
	    tstart = second();
	    SNEgy = 0;
	    if(SphP[i].EgyRes > 0 && dtime > 0 && x < 1)
	      {
		xold = x;
		perturbp.rho = SphP[i].a2.Density * a3inv;
		perturbp.PhysDensTh = myPhysDensThresh;
		perturbp.FEVP = myFactorEVP;
		perturbp.EpsilonExt = SNEgy;
		perturbp.tsfr = tsfr;
		perturbp.deltat = dtime;
		perturbp.FactorSN = myFactorSN;
		perturbp.ehot_old = egyhot;
		perturbp.ne = ne;
		perturbp.Z = Z;

		F.function = &xfx;
		F.params = NULL;
		my_gslstatus = 0;

		if(s == NULL)
		  {
		    T = gsl_root_fsolver_brent;
		    s = gsl_root_fsolver_alloc(T);
		  }

		x_lo = 0.0001;
		x_hi = 0.9999;
		gsl_root_fsolver_set(s, &F, x_lo, x_hi);
		if(!my_gslstatus)
		  {
		    niter = 0;
		    xoldold = xold = x;
		    do
		      {
			xold = x;
			niter++;

			my_gslstatus = gsl_root_fsolver_iterate(s);
			x = gsl_root_fsolver_root(s);
			x_lo = gsl_root_fsolver_x_lower(s);
			x_hi = gsl_root_fsolver_x_upper(s);
			my_gslstatus = gsl_root_test_interval(x_lo, x_hi, 0, PERTURB_REL_PREC);

			xrel = fabs(x - xold) / x;
		      }
		    while((my_gslstatus == GSL_CONTINUE) && (niter < PERTURB_MAX_ITER));
		    xold = xoldold;
		  }
		else
		  {
		    perturbp.rho = SphP[i].a2.Density * a3inv;
		    perturbp.PhysDensTh = myPhysDensThresh;
		    perturbp.FEVP = myFactorEVP;
		    perturbp.EpsilonExt = SNEgy;
		    perturbp.tsfr = tsfr;
		    perturbp.deltat = dtime;
		    perturbp.FactorSN = myFactorSN;
		    perturbp.ehot_old = egyhot;
		    perturbp.ne = ne;
		    perturbp.Z = Z;
		    niter = perturb(&x, &xrel);

		    egyhot = perturbp.ehot_old;
		    factorEVP = perturbp.fEVP;
		  }

		if((niter == PERTURB_MAX_ITER) && (xrel > PERTURB_REL_PREC))
		  {
		    if(xrel > 0.01)
		      {
			printf("  >>1 %8.6g %d %d %8.6g %8.6g %8.6g %8.6g %8.6g %8.6g %8.6g %8.6g %i\n",
			       All.Time, P[i].ID, ThisTask,
			       SphP[i].a2.Density * a3inv * All.UnitDensity_in_cgs,
			       convert_u_to_temp(SphP[i].Entropy *
						 pow(SphP[i].a2.Density * a3inv,
						     GAMMA_MINUS1) / GAMMA_MINUS1 *
						 All.UnitPressure_in_cgs / All.UnitDensity_in_cgs,
						 SphP[i].a2.Density * All.UnitDensity_in_cgs, &ne), Z,
			       dtime, SNEgy * All.UnitEnergy_in_cgs / All.UnitMass_in_g, xrel, x, xold,
			       niter);
			fflush(stdout);
		      }
		    /*x = (x + xold) / 2; */
		  }

	      }
	    tend = second();
	    CPU_eff_iter = timediff(tstart, tend);

#endif /* closes LT_MOD_EFFM */


	    egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

	    if(dt > 0)
	      {
		if(P[i].Ti_endstep > P[i].Ti_begstep)	/* upon start-up, we need to protect against dt==0 */
		  {
		    trelax = tsfr * (1 - x) / x / (myFactorSN * (1 + factorEVP));
		    egycurrent =
		      SphP[i].Entropy * pow(SphP[i].a2.Density * a3inv, GAMMA_MINUS1) / GAMMA_MINUS1;


#ifdef COSMIC_RAYS
		    egycurrent += CR_Particle_ThermalizeAndDissipate(SphP + i, dtime);
#endif /* COSMIC_RAYS */

#if !defined(LT_MOD_EFFM) && !defined(LT_SNegy_IN_HOTPHASE)
		    egycurrent += SphP[i].EgyRes / P[i].Mass;
#endif

#if defined(BH_THERMALFEEDBACK) || defined(BH_KINETICFEEDBACK)
		    if(SphP[i].i.Injected_BH_Energy > 0)
		      {
			egycurrent += SphP[i].i.Injected_BH_Energy / P[i].Mass;

			temp = u_to_temp_fac * egycurrent;

			if(temp > 5.0e9)
			  egycurrent = 5.0e9 / u_to_temp_fac;

			if(egycurrent > egyeff)
			  {
#ifndef LT_METAL_COOLING
			    tcool = GetCoolingTime(egycurrent, SphP[i].a2.Density * a3inv, &ne);
#else
			    if((tcool = GetCoolingTime(egyhot, SphP[i].a2.Density * a3inv, &ne, Z)) == 0)
			      tcool = 1e13;
#endif
			    if(tcool < trelax && tcool > 0)
			      trelax = tcool;
			  }

			SphP[i].i.Injected_BH_Energy = 0;
		      }
#endif



#if !defined(NOISMPRESSURE)
		    SphP[i].Entropy =
		      (egyeff +
		       (egycurrent -
			egyeff) * exp(-dtime / trelax)) * GAMMA_MINUS1 /
		      pow(SphP[i].a2.Density * a3inv, GAMMA_MINUS1);

		    SphP[i].e.DtEntropy = 0;
#endif
#if defined(LT_EXTEGY_INFO) && defined(LT_MOD_EFFM)
		    if(EEInfo_grain > SEvInfo_GRAIN)
		      {
			tstart = second();
			if(SphP[i].EgyRes > 0)
			  {
			    StatSum[REScount]++;
			    r = x / xold;
			    if(xold > 0)
			      {
				if(r < StatMin[XRatio] && r > 0)
				  StatMin[XRatio] = r;
				if(r > StatMax[XRatio])
				  StatMax[XRatio] = r;
				StatSum[XRatio] += r;
			      }
			    r = egyhot_old / egyhot;
			    if(r < StatMin[ERatio] && r > 0)
			      StatMin[ERatio] = r;
			    if(r > StatMax[ERatio])
			      StatMax[ERatio] = r;
			    StatSum[ERatio] += r;

			    r = (SphP[i].EgyRes / P[i].Mass) / All.IRA_erg_per_g;
			    if(r < StatMin[OEgy] && r > 0)
			      StatMin[OEgy] = r;
			    if(r > StatMax[OEgy])
			      StatMax[OEgy] = r;
			    StatSum[OEgy] += r;
			  }

			r = egycurrent / egy_new;
			if(SphP[i].EgyRes > 0)
			  {
			    if(r < StatMin[EFRatio] && r > 0)
			      StatMin[EFRatio] = r;
			    if(r > StatMax[EFRatio])
			      StatMax[EFRatio] = r;
			    StatSum[EFRatio] += r;
			  }
			else
			  {
			    StatSum[NoREScount]++;
			    if(r < StatMin[NoREFratio] && r > 0)
			      StatMin[NoREFratio] = r;
			    if(r > StatMax[NoREFratio])
			      StatMax[NoREFratio] = r;
			    StatSum[NoREFratio] += r;
			  }

			tend = second();
			CPU_ee_info = timediff(tstart, tend);
		      }
#endif

		    SphP[i].EgyRes = 0;
		  }
	      }

	    cloudmass = x * P[i].Mass;
#ifdef LT_EJECTA_IN_HOTPHASE
	    SphP[i].x = x;
#endif

	    if(tsfr < dtime)
	      tsfr = dtime;

	    sm = dtime / tsfr * cloudmass;

	    p = sm / P[i].Mass;

	    sum_sm[IMFi] += P[i].Mass * (1 - exp(-p));

	    SphP[i].Sfr = (1 - myFactorSN) * cloudmass / tsfr *
	      (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);


	    exp_minus_p = exp(-p);
#ifndef LT_STOCHASTIC_IRA
	    mstar = sm * exp_minus_p;
#else
	    mstar = P[i].Mass * (1 - exp_minus_p);
#endif


#if defined(LT_SNII) && defined(LT_LOCAL_IRA)
	    if(IMFs[IMFi].nonZeroIRA)
	      {
		for(j = tot_mass = 0; j < LT_NMetP; j++)
		  {
		    /* add the newly formed elements */

		    mz = mstar * SnII_ShortLiv_Yields[Yset][j][YZbin];

		    /* subtract the fraction locked in stars at birth
		     *   m_star * (All.MassFrac_inIRA - myFactorSN) * SphP[i].Metals[j] / P[i].Mass
		     */
		    /*
		       mz -= (1 - exp_minus_p) * SphP[i].Metals[j] * (All.MassFrac_inIRA - myFactorSN);
		     */

#ifdef LT_ACCOUNT_NONPROC_METALS
		    /* account for non-processed elements */
		    mz += sm * exp_minus_p * SphP[i].Metals[j] / P[i].Mass *
		      SnII_ShortLiv_Yields[Yset][Hyd][YZbin];
#endif
		    /*
		     *                     MassFrac_inIRA = mass fraction of stars in IRA Sn
		     *                     FactorSN = restored fraction by IRA Sn
		     *                     MassFra_inIRA - FactorSN = locked fraction
		     */

		    tot_mass += mz;

		    if(SphP[i].Metals[j] + mz < 0)
		      {
			printf(" \n\n !!! GULP !!! \n %i %i %i :: %g %g %g %g\n", ThisTask, P[i].ID, j,
			       SphP[i].Metals[j], mz, sm, p);
			endrun(101);
		      }
		    else
		      {
			SphP[i].Metals[j] += mz;
#ifdef LT_SEv_INFO
			mZ[j] += mz;
#endif
		      }
		  }

		if(tot_mass > P[i].Mass)
		  {
		    printf(" \n\n !!! GULP 2 !!! \n %i %i %i :: %g %g %g %g\n", ThisTask, P[i].ID, j,
			   tot_mass, P[i].Mass, sm, p);
		    endrun(102);
		  }
		tot_mass += mstar * SnII_ShortLiv_Yields[Yset][LT_NMetP][YZbin];
	      }
#endif

#ifndef LT_LOCAL_IRA
	    if(IMFs[IMFi].nonZeroIRA)
	      SphP[i].mstar = mstar;
#endif


	    /* the upper bits of the gas particle ID store how man stars this gas
	       particle gas already generated */

	    number_of_stars_generated = (P[i].ID >> (32 - bits));

	    mass_of_star = P[i].Mass / (GENERATIONS - number_of_stars_generated);

	    prob = P[i].Mass * (1 - exp(-p)) / mass_of_star;

#else /* belongs to ifndef(QUICK_LYALPHA) */

	    prob = 2.0;		/* this will always cause a star creation event */

#endif /* ends to QUICK_LYALPHA */

	    if(get_random_number(P[i].ID + 1) < prob)	/* ok, make a star */
	      {

#if defined(LT_SNIa) || defined(LT_AGB)
		ti_step_Ia = get_chemstep(1, 0, All.Time, All.Ti_Current, All.Ti_Current, IMFi);
#endif
#ifdef LT_SNII
		ti_step_II = get_chemstep(2, 0, All.Time, All.Ti_Current, All.Ti_Current, IMFi);
#endif

		if(number_of_stars_generated == (GENERATIONS - 1))
		  {
		    /* here we turn the gas particle itself into a star */
		    Stars_converted++;
		    stars_converted++;

		    sum_mass_stars[IMFi] += P[i].Mass;

		    P[i].Type = 4;

		    P[i].StellarAge = All.Time;

#ifdef LT_SEvDbg
		    if(ThisTask == 0 && FirstID == 0)
		      if(drand48() > 0.5)
			FirstID = P[i].ID;
#endif

		    if(N_star + 1 >= All.MaxPartMet)
		      {
			printf
			  ("On Task=%d with NumPart=%d we try to convert %d particles. Sorry, no space left...(All.MaxPartMet=%d)\n",
			   ThisTask, NumPart, stars_converted, All.MaxPartMet);
			fflush(stdout);
			endrun(8889);
		      }

		    P[i].MetID = N_star;
		    MetP[N_star].PID = i;
#if defined(LT_SNIa) || defined(LT_AGB)
		    if(ti_step_Ia > 0)
		      {
			if((MetP[N_star].NextChemStepIa = All.Ti_Current + ti_step_Ia) > TIMEBASE)
			  MetP[N_star].NextChemStepIa = TIMEBASE - 1;
		      }
		    else
		      MetP[N_star].NextChemStepIa = TIMEBASE + 2;
		    MetP[N_star].LastIaTime = All.mean_lifetime;

#endif
#ifdef LT_SNII
		    if(ti_step_II > 0)
		      {
			if((MetP[N_star].NextChemStepII = All.Ti_Current + ti_step_II) > TIMEBASE)
			  MetP[N_star].NextChemStepII = TIMEBASE - 1;
		      }
		    else
		      MetP[N_star].NextChemStepII = TIMEBASE + 2;
		    MetP[N_star].LastIITime = IMFs[IMFi].ShortLiv_TimeTh;

		    P[i].Mass *= (1 - IMFs[IMFi].FactorSN);
#endif
		    MetP[N_star].iMass = P[i].Mass;
		    MetP[N_star].SLguess = SphP[i].Hsml * pow(All.SpreadNeighCoeff, 0.3333333);

		    for(j = 0; j < LT_NMetP; j++)
		      MetP[N_star].Metals[j] = SphP[i].Metals[j];

#ifdef LT_SMOOTH_Z
		    MetP[N_star].Zsmooth = SphP[i].Zsmooth;
#endif
#ifdef LT_TRACK_CONTRIBUTES
		    MetP[N_star].contrib = SphP[i].contrib;
#endif

		    N_star++;
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

		    if(N_star + 1 >= All.MaxPartMet)
		      {
			printf
			  ("On Task=%d with NumPart=%d we try to spawn %d particles. Sorry, no space left...(All.MaxPartMet=%d)\n",
			   ThisTask, NumPart, stars_spawned, All.MaxPartMet);
			fflush(stdout);
			endrun(8889);
		      }

		    P[NumPart + stars_spawned] = P[i];
		    P[NumPart + stars_spawned].Type = 4;
		    P[i].ID += (1 << (32 - bits));
#ifdef LT_SEvDbg
		    if(ThisTask == 0 && FirstID == 0)
		      if(ThisTask == 0 && FirstID == 0)
			if(drand48() > 0.5)
			  FirstID = P[NumPart + stars_spawned].ID;
#endif

		    P[NumPart + stars_spawned].Mass = mass_of_star;

		    MetP[N_star].iMass = P[NumPart + stars_spawned].Mass;
		    factor = P[NumPart + stars_spawned].Mass / P[i].Mass;

		    P[i].Mass -= P[NumPart + stars_spawned].Mass;
		    sum_mass_stars[IMFi] += P[NumPart + stars_spawned].Mass;

		    P[NumPart + stars_spawned].StellarAge = All.Time;

		    P[NumPart + stars_spawned].MetID = N_star;
		    MetP[N_star].PID = NumPart + stars_spawned;

#ifdef LT_SNIa
		    if(ti_step_Ia > 0)
		      {
			if((MetP[N_star].NextChemStepIa = All.Ti_Current + ti_step_Ia) > TIMEBASE)
			  MetP[N_star].NextChemStepIa = TIMEBASE - 1;
		      }
		    else
		      MetP[N_star].NextChemStepIa = TIMEBASE + 2;
		    MetP[N_star].LastIaTime = All.mean_lifetime;
#endif
#ifdef LT_SNII
		    if(ti_step_II > 0)
		      {
			if((MetP[N_star].NextChemStepII = All.Ti_Current + ti_step_II) > TIMEBASE)
			  MetP[N_star].NextChemStepII = TIMEBASE - 1;
		      }
		    else
		      MetP[N_star].NextChemStepII = TIMEBASE + 2;
		    MetP[N_star].LastIITime = IMFs[IMFi].ShortLiv_TimeTh;
#endif


		    for(j = 0; j < LT_NMetP; j++)
		      {
#ifdef LT_LOCAL_IRA
			SphP[i].Metals[j] = DMAX(SphP[i].Metals[j] - mstar * SnII_ShortLiv_Yields[Yset][j][YZbin], 0);	/* to avoid round-off errors */
#endif
			MetP[N_star].Metals[j] = SphP[i].Metals[j] * factor;
			SphP[i].Metals[j] *= (1 - factor);
#ifdef LT_LOCAL_IRA
			SphP[i].Metals[j] += mstar * SnII_ShortLiv_Yields[Yset][j][YZbin];
#endif
		      }

		    P[NumPart + stars_spawned].Mass *= (1 - IMFs[IMFi].FactorSN);

		    MetP[N_star].SLguess = SphP[i].Hsml * pow(All.SpreadNeighCoeff, 0.3333333);

#ifdef LT_SMOOTH_Z
		    MetP[N_star].Zsmooth = SphP[i].Zsmooth;
#endif
#ifdef LT_TRACK_CONTRIBUTES
		    MetP[N_star].contrib = SphP[i].contrib;
#endif

		    N_star++;

		    force_add_star_to_tree(i, NumPart + stars_spawned);

		    stars_spawned++;
		  }
	      }

#ifdef COSMIC_RAYS
	    CR_Particle_SupernovaFeedback(&SphP[i], p * All.FeedbackEnergy * All.CR_SNEff, dtime);
#endif

#ifdef WINDS
	    /* Here comes the wind model */

	    if(P[i].Type == 0)	/* to protect using a particle that has been turned into a star */
	      {
		p = IMFs[IMFi].WindEfficiency * sm / P[i].Mass;

		prob = 1 - exp(-p);

		if(get_random_number(P[i].ID + 2) < prob)	/* ok, make the particle go into the wind */
		  {
#ifndef LT_WIND_VELOCITY
		    v =
		      sqrt(2 * IMFs[IMFi].WindEnergyFraction / IMFs[IMFi].WindEfficiency *
			   (IMFs[IMFi].totFactorSN / (1 - IMFs[IMFi].totFactorSN) * myEgySpecSN + SNEgy));

#ifdef LT_SEv_INFO

		    if(v < windv_min)
		      windv_min = v;
		    if(v > windv_max)
		      windv_max = v;
		    windv += v;
		    windn++;

#endif
#else
		    v = LT_WIND_VELOCITY;
#ifdef LT_SEv_INFO
		    IMFs[IMFi].WindEnergyFraction = v * v * IMFs[IMFi].WindEfficiency /
		      (2 * (IMFs[IMFi].totFactorSN / (1 - IMFs[IMFi].totFactorSN) * myEgySpecSN + SNEgy));

		    if(IMFs[IMFi].WindEnergyFraction < windv_min)
		      windv_min = IMFs[IMFi].WindEnergyFraction;
		    if(IMFs[IMFi].WindEnergyFraction > windv_max)
		      windv_max = IMFs[IMFi].WindEnergyFraction;
		    windv += IMFs[IMFi].WindEnergyFraction;
		    windn++;
#endif
#endif

#ifdef ISOTROPICWINDS
		    theta = acos(2 * get_random_number(P[i].ID + 3) - 1);
		    phi = 2 * M_PI * get_random_number(P[i].ID + 4);

		    dir[0] = sin(theta) * cos(phi);
		    dir[1] = sin(theta) * sin(phi);
		    dir[2] = cos(theta);
#else
		    dir[0] = P[i].g.GravAccel[1] * P[i].Vel[2] - P[i].g.GravAccel[2] * P[i].Vel[1];
		    dir[1] = P[i].g.GravAccel[2] * P[i].Vel[0] - P[i].g.GravAccel[0] * P[i].Vel[2];
		    dir[2] = P[i].g.GravAccel[0] * P[i].Vel[1] - P[i].g.GravAccel[1] * P[i].Vel[0];
#endif

		    for(j = 0, norm = 0; j < 3; j++)
		      norm += dir[j] * dir[j];

		    norm = sqrt(norm);
		    if(get_random_number(P[i].ID + 5) < 0.5)
		      norm = -norm;

		    if(norm != 0)
		      {
			for(j = 0; j < 3; j++)
			  dir[j] /= norm;

			for(j = 0; j < 3; j++)
			  {
			    P[i].Vel[j] += v * ascale * dir[j];
			    SphP[i].VelPred[j] += v * ascale * dir[j];
			  }
			/* Try to avoid problems for early haloes forming stars via H2 */
			SphP[i].DelayTime = DMAX(2*PPP[i].hsml,All.WindFreeTravelLength) / v;
#ifdef LT_HOT_WINDS
			u_hotwinds = (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * 1.16e7;
			u_hotwinds *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;
			u_hotwinds /= 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));

			SphP[i].Entropy =
			  GAMMA_MINUS1 * u_hotwinds / pow(SphP[i].a2.Density * a3inv, GAMMA_MINUS1);
#endif
		      }
		  }
	      }
#endif
	  }
      }				/* end of main loop over active particles */

#ifdef LT_SEvDbg
  if(checkFirstID == 0)
    {
      MPI_Bcast(&FirstID, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
      if(FirstID > 0)
	checkFirstID = 1;
    }
#endif

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
      All.TotN_star += tot_spawned + tot_converted;

      /* Note: N_gas is only reduced once rearrange_particle_sequence is called */

      /* Note: New tree construction can be avoided because of  `force_add_star_to_tree()' */
    }

  for(i = 0, sfrrate = 0; i < N_gas; i++)
    if(P[i].Type == 0)
      sfrrate += SphP[i].Sfr;

  MPI_Allreduce(&sfrrate, &totsfrrate, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for(i = 0; i < IMFs_dim; i++)
    {
      Sum_sm += sum_sm[i];
      Sum_mass_stars += sum_mass_stars[i];
    }
#ifdef LT_SEv_INFO
  MPI_Allreduce(&Sum_sm, &Total_sm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  MPI_Reduce(&Sum_sm, &Total_sm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
  MPI_Reduce(&Sum_mass_stars, &Total_sum_mass_stars, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sum_sm[0], &total_sm[0], IMFs_dim, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sum_mass_stars[0], &total_sum_mass_stars[0], IMFs_dim, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      if(All.TimeStep > 0)
	rate = Total_sm / (All.TimeStep / time_hubble_a);
      else
	rate = 0;

      /* convert to solar masses per yr */

      rate_in_msunperyear = rate * (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

      fprintf(FdSfr, "%g %g %g %g %g ", All.Time, Total_sm, totsfrrate, rate_in_msunperyear,
	      Total_sum_mass_stars);
      for(i = 0; i < IMFs_dim; i++)
	fprintf(FdSfr, "%g ", total_sm[i] / (All.TimeStep / time_hubble_a) *
		(All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR));
      for(i = 0; i < IMFs_dim; i++)
	fprintf(FdSfr, "%g ", total_sum_mass_stars[i]);
      fprintf(FdSfr, "\n");

      fflush(FdSfr);
    }
#ifdef LT_CharT_INFO
  if(CharTInfo_grain > CharT_Info_GRAIN)
    {
      for(i = 0; i < 20; i++)
	{
	  tcool_tcross[i] /= N_gas;
	  tcool_tff[i] /= N_gas;
	}
      if(ncool > 0)
	{
	  tcool_tcross[20] /= ncool;
	  tcool_tff[20] /= ncool;
	}
      if(nsf > 0)
	{
	  tcool_tcross[21] /= nsf;
	  tcool_tff[21] /= nsf;
	}

      MPI_Gather(&tcool_tcross[0], 22, MPI_DOUBLE, &tot_tcool_tcross[0], 22, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Gather(&tcool_tff[0], 22, MPI_DOUBLE, &tot_tcool_tff[0], 22, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      if(ThisTask == 0)
	{
	  for(j = 0; j < 22; j++)
	    {
	      tot_tc_tcr[j] = 0;
	      tot_div = 0;
	      for(i = 0; i < NTask; i++)
		if(tot_tcool_tcross[i * 22 + j] != 0)
		  {
		    tot_tc_tcr[j] += tot_tcool_tcross[i * 22 + j];
		    tot_div++;
		  }
	      if(tot_div > 0)
		tot_tc_tcr[j] /= tot_div;
	    }
	  for(j = 0; j < 22; j++)
	    {
	      tot_tc_tff[j] = 0;
	      tot_div = 0;
	      for(i = 0; i < NTask; i++)
		if(tot_tcool_tff[i * 22 + j] != 0)
		  {
		    tot_tc_tff[j] += tot_tcool_tff[i * 22 + j];
		    tot_div++;
		  }
	      if(tot_div > 0)
		tot_tc_tff[j] /= tot_div;
	    }

	  fprintf(FdCharT, "%g ", All.Time);
	  for(i = 0; i < 22; i++)
	    fprintf(FdCharT, "%g ", tot_tc_tcr[i]);
	  for(i = 0; i < 22; i++)
	    fprintf(FdCharT, "%g ", tot_tc_tff[i]);
	  fprintf(FdCharT, "\n");
	}
      CharTInfo_grain = 0;
    }
#endif

/* temporary */

/* #ifdef LT_SEv_INFO */
/*   if(SEvInfo_grain > SEvInfo_GRAIN) */
/*     { */
/*       SEvInfo_grain = 0; */
/*       MPI_Reduce(&Tpop_num[0], &tot_Tpop_num[0], Trun_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); */
/*       MPI_Reduce(&Tpop_mass[0], &tot_Tpop_mass[0], Trun_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); */

/*       MPI_Reduce(&TCpop_num[0], &tot_TCpop_num[0], TCrun_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); */
/*       MPI_Reduce(&TCpop_mass[0], &tot_TCpop_mass[0], TCrun_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); */

/*       MPI_Reduce(&all_mass, &tot_allmass, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); */
/*       MPI_Reduce(&all_num, &tot_allnum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); */

/*       if(ThisTask == 0) */
/* 	{ */
/* 	  fprintf(FdTrun, "> %g ", All.Time); */
/* 	  for(j = 0; j < Trun_size; j++) */
/* 	    fprintf(FdTrun, "%g ", tot_Tpop_num[j] / tot_allnum); */
/* 	  fprintf(FdTrun, "\n"); */
/* 	  fprintf(FdTrun, "< %g ", All.Time); */
/* 	  for(j = 0; j < Trun_size; j++) */
/* 	    fprintf(FdTrun, "%g ", tot_Tpop_mass[j] / tot_allmass); */
/* 	  fprintf(FdTrun, "\n"); */
/* 	  fflush(FdTrun); */

/* 	  fprintf(FdTCrun, "> %g ", All.Time); */
/* 	  for(j = 0; j < TCrun_size; j++) */
/* 	    fprintf(FdTCrun, "%g ", tot_TCpop_num[j] / tot_allnum); */
/* 	  fprintf(FdTCrun, "\n"); */
/* 	  fprintf(FdTCrun, "< %g ", All.Time); */
/* 	  for(j = 0; j < TCrun_size; j++) */
/* 	    fprintf(FdTCrun, "%g ", tot_TCpop_mass[j] / tot_allmass); */
/* 	  fprintf(FdTCrun, "\n"); */
/* 	  fflush(FdTCrun); */
/* 	} */

/*       MPI_Barrier(MPI_COMM_WORLD); */
/*     } */
/* #endif */

  /* */



#if defined(LT_EXTEGY_INFO) && defined(LT_MOD_EFFM)

  /* collect some infos about external energy in the reservoir SphP[].EgyRes.
     could be used for any egy source (AGN etc).
   */
  if(EEInfo_grain > SEvInfo_GRAIN)
    {
      tstart = second();
      MPI_Reduce(&StatSum[0], &tot_StatSum[0], STAT_INUM + 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&StatMin[0], &tot_StatMin[0], STAT_INUM, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
      MPI_Reduce(&StatMax[0], &tot_StatMax[0], STAT_INUM, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

      if(ThisTask == 0)
	{
	  fprintf(FdExtEgy, "%8.6lg ", All.Time);

	  if(tot_StatSum[REScount])
	    {
	      if(tot_StatSum[REScount] > 2)
		for(i = 0; i < STAT_INUM - 1; i++)
		  tot_StatSum[i] = (tot_StatSum[i] -=
				    tot_StatMin[i] + tot_StatMax[i]) / (tot_StatSum[REScount] - 2);
	      else
		for(i = 0; i < STAT_INUM - 1; i++)
		  tot_StatSum[i] /= tot_StatSum[REScount];

	      fprintf(FdExtEgy, "%g "
		      "%8.6lg %8.6lg %8.6lg "
		      "%8.6lg %8.6lg %8.6lg " "%8.6lg %8.6lg %8.6lg " "%8.6lg %8.6lg %8.6lg ",
		      /* (2) # of particles with nonzero reservoir */
		      tot_StatSum[REScount],
		      /* (3-5) {min,max,mean} ratio between egyhot before and after modified eff model */
		      tot_StatMin[ERatio], tot_StatMax[ERatio], tot_StatSum[ERatio],
		      /* (6-8) {min,max,mean} ratio between egy before and after effective model + modif.eff.model */
		      tot_StatMin[EFRatio], tot_StatMax[EFRatio], tot_StatSum[EFRatio],
		      /* (9-11) {min,max,mean} ratio between cold fraction before and after modified eff model */
		      tot_StatMin[XRatio], tot_StatMax[XRatio], tot_StatSum[XRatio],
		      /* (12-14) {min,max,mean} ratio between egy of reservoir and egy due to IRA SnII */
		      tot_StatMin[OEgy], tot_StatMax[OEgy], tot_StatSum[OEgy]);
	    }
	  if(tot_StatSum[NoREScount])
	    {
	      if(tot_StatSum[NoREScount] > 2)
		tot_StatSum[NoREFratio] = (tot_StatSum[NoREFratio] -=
					   tot_StatMin[NoREFratio] +
					   tot_StatMax[NoREFratio]) / (tot_StatSum[NoREScount] - 2);
	      else
		tot_StatSum[NoREFratio] /= tot_StatSum[NoREScount];

	      fprintf(FdExtEgy, "%g %8.6lg %8.6lg %8.6lg",
		      /* # of particles with reservoir = 0 */
		      /* {min, max,mean} ratio between egy before and after effective model */
		      tot_StatSum[NoREScount], tot_StatMin[NoREFratio], tot_StatMax[NoREFratio],
		      tot_StatSum[NoREFratio]);
	    }

	  fprintf(FdExtEgy, "\n");
	  fflush(FdExtEgy);
	}
      tend = second();
      CPU_ee_info += timediff(tstart, tend);

      MPI_Reduce(&CPU_ee_info, &sumCPU_ee_info, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      if(ThisTask == 0)
	All.CPU_EE_info += sumCPU_ee_info / NTask;

      EEInfo_grain = 0;
    }
#endif /* closes LT_EXTEGY_INFO */

#ifdef LT_SEv_INFO

  /* collect some informations */

  tstart = second();
#ifdef LT_SNII

  /* write metals produced by SnII in IRA */
#ifdef LT_LOCAL_IRA
  if(Total_sm > 0)
    {
      MPI_Reduce(&mZ[0], &tot_mZ[0], LT_NMet, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      if(ThisTask == 0)
	{
	  for(j = 0; j < LT_NIMFs; j++)
	    {
	      fprintf(FdSn, "IRA %1d %.8lg 0 0 %g", j, All.Time,
		      total_sm[j] * All.UnitMass_in_g / SOLAR_MASS * 1);
	      for(i = 0; i < LT_NMet; i++)
		fprintf(FdSn, "%.8lg ", tot_mZ[i] * All.UnitMass_in_g / SOLAR_MASS);
	      fprintf(FdSn, "\n");
	    }
	  fflush(FdSn);
	}
    }
#endif

#endif
#ifdef WINDS
  /* write the minimum, maximum and mean wind velocity (to check the external egy contrib) */
  MPI_Allreduce(&windn, &tot_windn, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(tot_windn > 0)
    {
      if(windv_min == 1e6)
	windv_min = 0;
      MPI_Reduce(&windv, &tot_windv, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&windv_min, &tot_windv_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
      MPI_Reduce(&windv_max, &tot_windv_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      if(ThisTask == 0)
	{
	  if(tot_windn > 1)
	    tot_windv = (tot_windv - tot_windv_max) / (tot_windn - 1);

	  fprintf(FdWinds, "%g %i %g %g %g\n", All.Time, tot_windn, tot_windv_min, tot_windv_max, tot_windv);
	  fflush(FdWinds);
	}
    }
#endif

  tend = second();
  CPU_sn_info = timediff(tstart, tend);

  MPI_Reduce(&CPU_sn_info, &sumCPU_sn_info, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(ThisTask == 0)
    All.CPU_SN_info += sumCPU_sn_info / NTask;

#endif /* closes LT_SEv_INFO  */

#ifdef LT_MOD_EFFM
  gsl_set_error_handler(old_error_handler);
#endif

  MPI_Reduce(&CPU_eff_iter, &sumCPU_eff_iter, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(ThisTask == 0)
    All.CPU_Eff_Iter += sumCPU_eff_iter / NTask;

}


double get_starformation_rate(int i)
{
  double rateOfSF;
  double a3inv;
  int IMFi, flag, Yset;
  double tsfr;
  double factorEVP, egyhot, ne, tcool, y, x, cloudmass;

  double Z, Zsol;

  double myFactorEVP, myPhysDensThresh, myEgySpecSN, myFactorSN;

#ifdef LT_MOD_EFFM
  double SNEgy, xrel, xold;
  int niter;
  double dt, dtime, hubble_a, time_hubble_a;
#endif

  xclouds = 0;

  IMFi = get_IMF_index(i);
  Yset = IMFs[IMFi].YSet;
  IMFp = &IMFs[IMFi];
  myFactorSN = IMFs[IMFi].FactorSN;
  myEgySpecSN = IMFs[IMFi].EgySpecSN;

  if(IMFs[IMFi].SFTh_Zdep)
    {
#ifndef LT_SMOOTH_Z
      Z = get_metallicity(i, 2);
#else
      Z = SphP[i].Zsmooth * 0.0885;
#endif
/*       Z = get_metallicity(i, 2); */
/* #ifdef LT_SMOOTH_Z */
/*       Z = (Z + SphP[i].Zsmooth * 0.0885) / 2; */
/* #endif */

      Zsol = get_metallicity_solarunits(Z);
      getindex(&Zvalue[0], 0, ZBins - 1, &Zsol, &flag);
      if(flag == 0 || flag == ZBins - 1)
	{
	  myFactorEVP = IMFs[IMFi].FEVP[flag];
	  myPhysDensThresh = IMFs[IMFi].PhysDensThresh[flag];
	}
      else
	{
	  x = (Zvalue[flag + 1] - Z) / (Zvalue[flag + 1] - Zvalue[flag]);
	  myFactorEVP = IMFs[IMFi].FEVP[flag] * x + IMFs[IMFi].FEVP[flag + 1] * (1 - x);
	  myPhysDensThresh =
	    IMFs[IMFi].PhysDensThresh[flag] * x + IMFs[IMFi].PhysDensThresh[flag + 1] * (1 - x);
	}
    }
  else
    {
      myFactorEVP = IMFs[IMFi].FEVP[0];
      myPhysDensThresh = IMFs[IMFi].PhysDensThresh[0];
      Z = 0;
    }

  if(All.ComovingIntegrationOn)
    a3inv = 1 / (All.Time * All.Time * All.Time);
  else
    a3inv = 1;

  flag = 1;			/* default is normal cooling */

  if(SphP[i].a2.Density * a3inv >= myPhysDensThresh)

    flag = 0;

  if(All.ComovingIntegrationOn)
    if(SphP[i].a2.Density < All.OverDensThresh)
      flag = 1;

  if(flag == 1)
    return 0;


  tsfr = sqrt(myPhysDensThresh / (SphP[i].a2.Density * a3inv)) * All.MaxSfrTimescale;

  factorEVP = pow(SphP[i].a2.Density * a3inv / myPhysDensThresh, -0.8) * myFactorEVP;

  egyhot = myEgySpecSN / (1 + factorEVP) + All.EgySpecCold;

  ne = SphP[i].Ne;
#ifdef LT_METAL_COOLING
  if((tcool = GetCoolingTime(egyhot, SphP[i].a2.Density * a3inv, &ne, Z)) == 0)
    tcool = 1e13;
#else
  tcool = GetCoolingTime(egyhot, SphP[i].a2.Density * a3inv, &ne);
#endif

  y = tsfr / tcool * egyhot / (myFactorSN * myEgySpecSN - (1 - myFactorSN) * All.EgySpecCold);
  if(y < 1e-4)
    x = y;
  else
    x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

#ifdef LT_MOD_EFFM

  dt = (P[i].Ti_endstep - P[i].Ti_begstep) * All.Timebase_interval;

  if(All.ComovingIntegrationOn)
    {
      hubble_a = All.Hubble * sqrt(All.Omega0 / (All.Time * All.Time * All.Time)
				   + (1 - All.Omega0 - All.OmegaLambda) / (All.Time * All.Time) +
				   All.OmegaLambda);

      time_hubble_a = All.Time * hubble_a;
      dtime = All.Time * dt / time_hubble_a;
    }
  else
    dtime = dt;

  SNEgy = SphP[i].EgyRes / P[i].Mass;
  if(SNEgy > 0 && x < 1 && dtime > 0)
    {
      niter = 0;
      xold = x;
      perturbp.rho = SphP[i].a2.Density * a3inv;
      perturbp.PhysDensTh = myPhysDensThresh;
      perturbp.FEVP = myFactorEVP;
      perturbp.EpsilonExt = SNEgy;
      perturbp.tsfr = tsfr;
      perturbp.deltat = dtime;
      perturbp.FactorSN = myFactorSN;
      perturbp.ehot_old = egyhot;
      perturbp.ne = ne;
#ifdef LT_METAL_COOLING
      perturbp.Z = Z;
#endif

      niter = perturb(&x, &xrel);

      if(niter > PERTURB_MAX_ITER && xrel > PERTURB_REL_PREC)
	{
	  if(xrel > 0.1 && x > 1e-2)
	    {
	      printf("  >>2 %8.6g %d %d %8.6g %8.6g %8.6g %8.6g %8.6g %8.6g %i\n",
		     All.Time, P[i].ID, ThisTask,
		     SphP[i].a2.Density * All.UnitDensity_in_cgs,
		     convert_u_to_temp(SphP[i].Entropy * pow(SphP[i].a2.Density * a3inv, GAMMA_MINUS1) /
				       GAMMA_MINUS1 * All.UnitPressure_in_cgs / All.UnitDensity_in_cgs,
				       SphP[i].a2.Density * All.UnitDensity_in_cgs, &ne), dtime,
		     SNEgy * All.UnitEnergy_in_cgs / All.UnitMass_in_g, xrel, x, niter);
	      fflush(stdout);
	    }
	  /*x = (x + xold) / 2; */
	}

    }

#endif

  cloudmass = x * P[i].Mass;

  xclouds = x;


  rateOfSF = (1 - myFactorSN) * cloudmass / tsfr;

  /* convert to solar masses per yr */

  rateOfSF *= (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);

  return rateOfSF;
}


void rearrange_particle_sequence(void)
{
  int i, j;
  struct particle_data psave;

#ifdef BLACK_HOLES
  int count_elim, count_gaselim, tot_elim, tot_gaselim;
#endif

  if(Stars_converted)
    {
      N_gas -= Stars_converted;
      Stars_converted = 0;

      for(i = 0; i < N_gas; i++)
	if(P[i].Type != 0)
	  {
	    for(j = N_gas; j < NumPart; j++)
	      if(P[j].Type == 0)
		break;

	    if(j >= NumPart)
	      endrun(181170);

	    if(P[i].Type == 4)
	      MetP[P[i].MetID].PID = j;

	    psave = P[i];
	    P[i] = P[j];
	    SphP[i] = SphP[j];
	    P[j] = psave;
	  }
    }

#ifdef BLACK_HOLES
  count_elim = 0;
  count_gaselim = 0;

  for(i = 0; i < NumPart; i++)
    if(P[i].Mass == 0)
      {
	if(P[i].Type == 0)
	  {
	    P[i] = P[N_gas - 1];
	    SphP[i] = SphP[N_gas - 1];

	    P[N_gas - 1] = P[NumPart - 1];

	    N_gas--;

	    count_gaselim++;
	  }
	else
	  {
	    P[i] = P[NumPart - 1];
	  }

	NumPart--;
	i--;

	count_elim++;
      }

  MPI_Allreduce(&count_elim, &tot_elim, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&count_gaselim, &tot_gaselim, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      printf("Blackholes: Eliminated %d gas particles and merged away %d black holes.\n",
	     tot_gaselim, tot_elim - tot_gaselim);
      fflush(stdout);
    }

  All.TotNumPart -= tot_elim;
  All.TotN_gas -= tot_gaselim;
#endif
}

int write_eff_model(int num, int IMFi)
{
  char string[200];
  FILE *file;
  int i;

  if(ThisTask == 0)
    {
      if(num < 0)
	sprintf(string, "%s.%02d", "eff_model.dat", IMFi);
      else
	sprintf(string, "eff_model.dat.%02d.%03d", IMFi, num);


      /* write data on a file to save time in the next run */
      if((file = fopen(string, "w")) == NULL)
	{
	  printf("it's impossible to write in <%s> \n", string);
	  fclose(file);
	  return -1;
	}
      else
	{
	  for(i = 0; i < ZBins; i++)
	    fprintf(file, "%lg %lg %lg\n", Zvalue[i], PhysDensTh[i], FEVP[i]);
	  fclose(file);
	}
    }
  return 0;
}

int read_eff_model(int num, int IMFi)
{
  int i, j, k, Zbins;
  char string[200];
  FILE *file;

  k = 0;

  if(num < 0)
    sprintf(string, "%s.%02d", "eff_model.dat", IMFi);
  else
    sprintf(string, "eff_model.dat.%02d.%03d", IMFi, num);

  if((file = fopen(string, "r")) == NULL)
    k = 1;
  else
    {
      fscanf(file, "%d", &Zbins);
      for(i = j = 0; i < Zbins; i++)
	{
	  j += fscanf(file, "%*g %lg %lg", &IMFs[IMFi].PhysDensThresh[i], &IMFs[IMFi].FEVP[i]);
	  printf("[Fe/H]: %-4.3lg - RhoTh: %-9.7lg  - fEVP: %-9.7lg\n",
		 Zvalue[i], IMFs[IMFi].PhysDensThresh[i], IMFs[IMFi].FEVP[i]);
	}
      if(j != Zbins * 2)
	k = 1;
    }
  return k;
}

void init_clouds_cm(double *PhDTh, double *fEVP, double EgySN, double FSN, int Zbins, double *ZArray)
{
  int i;
  int *offset, *counts;

  if(ThisTask == 0)
    printf("\n\ninitialize effective model.. \n");

  /* no file or file corrupted. need to recalculate */

  if(ThisTask == 0)
    printf("it is needed to recalculate metallicity dependence for effective model.. it will take a while\n");

  /* distribute work over processor */
  offset = (int *) malloc(sizeof(int) * NTask);
  counts = (int *) malloc(sizeof(int) * NTask);

  offset[0] = 0;
  counts[0] = Zbins / NTask + (0 < (int) fmod(Zbins, NTask));
  for(i = 1; i < NTask; i++)
    {
      if((counts[i] = Zbins / NTask + (i < (int) fmod(Zbins, NTask))))
	offset[i] = offset[i - 1] + counts[i - 1];
      else
	offset[i] = 0;
    }

  for(i = 0; i < Zbins; i++)
    PhDTh[i] = fEVP[i] = 0;

  /* ThisTask will calculate params for its range of Z */
  for(i = 0; i < counts[ThisTask]; i++)
    init_clouds(1, EgySN, FSN, ZArray[i], &PhDTh[i], &fEVP[i]);

  /* collect informations */
  MPI_Allgatherv((void *) &PhDTh[offset[ThisTask]], counts[ThisTask], MPI_DOUBLE,
		 (void *) &PhDTh[0], &counts[0], &offset[0], MPI_DOUBLE, MPI_COMM_WORLD);

  MPI_Allgatherv((void *) &fEVP[offset[ThisTask]], counts[ThisTask], MPI_DOUBLE,
		 (void *) &fEVP[0], &counts[0], &offset[0], MPI_DOUBLE, MPI_COMM_WORLD);

  return;
}

void init_clouds(int mode, double egySN, double fSN, double Z, double *PhDTh, double *fEVP)
{
  int Zbin;
  double Zsol, A0, dens, tcool, ne, coolrate, egyhot, x, u4, meanweight;
  double tsfr, y, peff, fac, neff, egyeff, factorEVP, sigma, thresholdStarburst;

  /*
   * calculating effective model parameters you have different choices
   * in the case you have settled on LT_METAL_COOLING:
   *   (a) make A0 parameter (onset of thermal instability) dependent on
   *       metal cooling
   *   (b) leave A0 parameter independent of metal cooling
   * the following code trigger this choice as stated by compiler directives
   */

  Zsol = get_metallicity_solarunits(Z);
  for(Zbin = ZBins - 1; Zsol < Zvalue[Zbin] && Zbin > 0; Zbin--)
    ;

  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* note: assuming FULL ionization */
  u4 = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * 1.0e4;
  u4 *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;


  if(mode)
    {
      if(All.SFTh_Zdep == 1)
	*fEVP = All.TempSupernova / ThInst_onset[Zbin];
      A0 = *fEVP;

      egyhot = egySN / A0;

      if(All.ComovingIntegrationOn)
	dens = 1.0e6 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
      else
	dens = 1.0e6 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

      if(All.ComovingIntegrationOn)
	{
	  All.Time = 1.0;	/* to be guaranteed to get z=0 rate */
	  IonizeParams();
	}

      ne = 1.0;
      SetZeroIonization();
#ifdef LT_METAL_COOLING
      if((tcool = GetCoolingTime(egyhot, dens, &ne, Z)) == 0)
	tcool = 1e13;
#else
      tcool = GetCoolingTime(egyhot, dens, &ne);
#endif

      coolrate = egyhot / tcool / dens;

      x = (egyhot - u4) / (egyhot - All.EgySpecCold);

      *PhDTh = x / pow(1 - x, 2) *
	(fSN * egySN - (1 - fSN) * All.EgySpecCold) / (All.MaxSfrTimescale * coolrate);

      if(mode == 2)
	return;
    }

  dens = *PhDTh * 10;

  do
    {
      tsfr = sqrt(*PhDTh / (dens)) * All.MaxSfrTimescale;
      factorEVP = pow(dens / *PhDTh, -0.8) * *fEVP;
      egyhot = egySN / (1 + factorEVP) + All.EgySpecCold;

      ne = 0.5;

#ifdef LT_METAL_COOLING
      if((tcool = GetCoolingTime(egyhot, dens, &ne, Z)) == 0)
	tcool = 1e13;
#else
      tcool = GetCoolingTime(egyhot, dens, &ne);
#endif

      y = tsfr / tcool * egyhot / (fSN * egySN - (1 - fSN) * All.EgySpecCold);
      if(y < 1e-4)
	x = y;
      else
	x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

      egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

      peff = GAMMA_MINUS1 * dens * egyeff;

      fac = 1 / (log(dens * 1.025) - log(dens));
      dens *= 1.025;

      neff = -log(peff) * fac;


      tsfr = sqrt(*PhDTh / (dens)) * All.MaxSfrTimescale;
      factorEVP = pow(dens / *PhDTh, -0.8) * *fEVP;
      egyhot = egySN / (1 + factorEVP) + All.EgySpecCold;

      ne = 0.5;

#ifdef LT_METAL_COOLING
      if((tcool = GetCoolingTime(egyhot, dens, &ne, Z)) == 0)
	tcool = 1e13;
#else
      tcool = GetCoolingTime(egyhot, dens, &ne);
#endif

      y = tsfr / tcool * egyhot / (fSN * egySN - (1 - fSN) * All.EgySpecCold);
      if(y < 1e-4)
	x = y;
      else
	x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

      egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

      peff = GAMMA_MINUS1 * dens * egyeff;

      neff += log(peff) * fac;
    }
  while(neff > 4.0 / 3);

  thresholdStarburst = dens;

  integrate_sfr(*PhDTh, *fEVP, egySN, fSN, Z);

  sigma = 10.0 / All.Hubble * 1.0e-10 / pow(1.0e-3, 2);

  printf("\n[%i] .: [Fe] = %.3g :. \n"
	 "\nA0= %g  \n"
	 "Computed: PhysDensThresh= %g  (int units)         %g h^2 cm^-3\n"
	 "EXPECTED FRACTION OF COLD GAS AT THRESHOLD = %g\n\n"
	 "tcool=%g dens=%g egyhot=%g\n"
	 "Run-away sets in for dens=%g\n"
	 "Dynamic range for quiescent star formation= %g\n\n"
	 "Isotherm sheet central density: %g   z0=%g\n",
	 ThisTask,
	 Zvalue[Zbin],
	 A0,
	 *PhDTh,
	 *PhDTh / (PROTONMASS / HYDROGEN_MASSFRAC / All.UnitDensity_in_cgs), x, tcool, dens,
	 egyhot, thresholdStarburst, thresholdStarburst / *PhDTh,
	 M_PI * All.G * sigma * sigma / (2 * GAMMA_MINUS1) / u4,
	 GAMMA_MINUS1 * u4 / (2 * M_PI * All.G * sigma));
  fflush(stdout);

  if(All.ComovingIntegrationOn)
    {
      All.Time = All.TimeBegin;
      IonizeParams();
    }

  return;
}

void integrate_sfr(double PhDTh, double fEVP, double egySN, double fSN, double Z)
{
  double rho0, rho, rho2, q, dz, gam, sigma = 0, sigma_u4, sigmasfr = 0, ne, P1;
  double x = 0, y, P, P2, x2, y2, tsfr2, factorEVP2, egyhot2, tcool2, drho, dq;
  double meanweight, u4, z, tsfr, tcool, egyhot, factorEVP, egyeff, egyeff2;
  char buff[30];
  FILE *fd;

  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* note: assuming FULL ionization */
  u4 = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * 1.0e4;
  u4 *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

  if(All.ComovingIntegrationOn)
    {
      All.Time = 1.0;		/* to be guaranteed to get z=0 rate */
      IonizeParams();
    }

  if(ThisTask == 0)
    fd = fopen("eos.txt", "w");
  else
    fd = 0;

  for(rho = PhDTh; rho <= 1000 * PhDTh; rho *= 1.1)
    {
      tsfr = sqrt(PhDTh / rho) * All.MaxSfrTimescale;

      factorEVP = pow(rho / PhDTh, -0.8) * fEVP;

      egyhot = egySN / (1 + factorEVP) + All.EgySpecCold;

      ne = 1.0;
#ifdef LT_METAL_COOLING
      if((tcool = GetCoolingTime(egyhot, rho, &ne, Z)) == 0)
	tcool = 1e13;
#else
      tcool = GetCoolingTime(egyhot, rho, &ne);
#endif
      y = tsfr / tcool * egyhot / (fSN * egySN - (1 - fSN) * All.EgySpecCold);
      if(y < 1e-4)
	x = y;
      else
	x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

      egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

      P = GAMMA_MINUS1 * rho * egyeff;

      if(ThisTask == 0)
	{
	  fprintf(fd, "%g %g\n", rho, P);
	}
    }

  if(ThisTask == 0)
    fclose(fd);


  sprintf(buff, "sfrrate.%03d.txt", sfrrate_filenum);
  if(ThisTask == 0)
    fd = fopen(buff, "w");
  else
    fd = 0;

  for(rho0 = PhDTh; rho0 <= 10000 * PhDTh; rho0 *= 1.02)
    {
      z = 0;
      rho = rho0;
      q = 0;
      dz = 0.001;

      sigma = sigmasfr = sigma_u4 = 0;

      while(rho > 0.0001 * rho0)
	{
	  if(rho > PhDTh)
	    {
	      tsfr = sqrt(PhDTh / rho) * All.MaxSfrTimescale;

	      factorEVP = pow(rho / PhDTh, -0.8) * fEVP;

	      egyhot = egySN / (1 + factorEVP) + All.EgySpecCold;

	      ne = 1.0;
#ifdef LT_METAL_COOLING
	      if((tcool = GetCoolingTime(egyhot, rho, &ne, Z)) == 0)
		tcool = 1e13;
#else
	      tcool = GetCoolingTime(egyhot, rho, &ne);
#endif

	      y = tsfr / tcool * egyhot / (fSN * egySN - (1 - fSN) * All.EgySpecCold);
	      if(y < 1e-4)
		x = y;
	      else
		x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

	      egyeff = egyhot * (1 - x) + All.EgySpecCold * x;

	      P = P1 = GAMMA_MINUS1 * rho * egyeff;

	      rho2 = 1.1 * rho;
	      tsfr2 = sqrt(PhDTh / rho2) * All.MaxSfrTimescale;
	      factorEVP2 = pow(rho2 / PhDTh, -0.8) * fEVP;
	      egyhot2 = egySN / (1 + factorEVP) + All.EgySpecCold;
#ifdef LT_METAL_COOLING
	      if((tcool2 = GetCoolingTime(egyhot2, rho2, &ne, Z)) == 0)
		tcool2 = 1e13;
#else
	      tcool2 = GetCoolingTime(egyhot2, rho2, &ne);
#endif
	      y2 = tsfr2 / tcool2 * egyhot2 / (fSN * egySN - (1 - fSN) * All.EgySpecCold);
	      if(y2 < 1e-4)
		x2 = y2;
	      else
		x2 = 1 + 1 / (2 * y2) - sqrt(1 / y2 + 1 / (4 * y2 * y2));
	      egyeff2 = egyhot2 * (1 - x2) + All.EgySpecCold * x2;
	      P2 = GAMMA_MINUS1 * rho2 * egyeff2;

	      gam = log(P2 / P1) / log(rho2 / rho);
	    }
	  else
	    {
	      tsfr = 0;

	      P = GAMMA_MINUS1 * rho * u4;
	      gam = 1.0;


	      sigma_u4 += rho * dz;
	    }



	  drho = q;
	  dq = -(gam - 2) / rho * q * q - 4 * M_PI * All.G / (gam * P) * rho * rho * rho;

	  sigma += rho * dz;
	  if(tsfr > 0)
	    {
	      sigmasfr += (1 - IMFs[sfrrate_filenum].FactorSN) * rho * x / tsfr * dz;
	    }

	  rho += drho * dz;
	  q += dq * dz;
	}


      sigma *= 2;		/* to include the other side */
      sigmasfr *= 2;
      sigma_u4 *= 2;


      if(ThisTask == 0)
	{
	  fprintf(fd, "%g %g %g %g\n", rho0, sigma, sigmasfr, sigma_u4);
	}
    }


  if(All.ComovingIntegrationOn)
    {
      All.Time = All.TimeBegin;
      IonizeParams();
    }

  if(ThisTask == 0)
    fclose(fd);
}



#ifdef LT_MOD_EFFM
int INLINE_FUNC perturb(double *x, double *xrel)
     /*
        use iterative method to converge to the new value of x.
        use Wegstein correction to be faster and more robust.
      */
{
#define old 0
#define now 1
  double factorEVP, egyI, egyII, egyhot, tcool, y, x_temp;
  double ox;
  double x_value[2], y_value[2];
  int iiter, oiter;

  oiter = 0;
  do
    {

      iiter = 0;
      ox = *x;
      y_value[old] = *x;
      y_value[now] = *x;

      x_value[old] = *x;
      x_value[now] = *x;
      do
	{
	  if(All.Mod_SEff)
	    factorEVP = pow(perturbp.rho / perturbp.PhysDensTh, -0.8) * perturbp.FEVP *
	      (1 + perturbp.EpsilonExt * *x / perturbp.EgySpecSN * perturbp.tsfr / perturbp.deltat);
	  perturbp.fEVP = factorEVP;

	  if(*x > 0)
	    {
	      egyI = perturbp.EgySpecSN / (1 + factorEVP) + All.EgySpecCold;
	      if(All.Mod_SEff)
		egyII = perturbp.EpsilonExt * (1 - *x) /
		  (perturbp.FactorSN * (1 + factorEVP)) * (1 - *x) / *x * (perturbp.tsfr / perturbp.deltat);
	      else
		egyII = perturbp.EpsilonExt * (1 - *x) /
		  (perturbp.FactorSN * (1 + factorEVP)) / *x * (perturbp.tsfr / perturbp.deltat);
	    }
	  else
	    {
	      egyI = (perturbp.EgySpecSN + perturbp.EpsilonExt) / (1 + factorEVP) + All.EgySpecCold;
	      egyII = 0;
	    }

	  egyhot = (egyI + egyII);
	  perturbp.ehot_old = egyhot;

#ifdef LT_METAL_COOLING
	  if((tcool = GetCoolingTime(egyhot, perturbp.rho, &(perturbp.ne), perturbp.Z)) == 0)
	    tcool = 1e13;
#else
	  tcool = GetCoolingTime(egyhot, perturbp.rho, &(perturbp.ne));
#endif

	  if(*x > 0)
	    {
	      if(All.Mod_SEff)
		y = perturbp.tsfr * egyhot / tcool /
		  ((perturbp.FactorSN * perturbp.EgySpecSN - (1 - perturbp.FactorSN) * All.EgySpecCold) +
		   perturbp.EpsilonExt * perturbp.tsfr / perturbp.deltat * ((1 - *x) / *x * (1 - *x)));
	      else
		y = perturbp.tsfr * egyhot / tcool /
		  ((perturbp.FactorSN * perturbp.EgySpecSN - (1 - perturbp.FactorSN) * All.EgySpecCold) +
		   perturbp.EpsilonExt * perturbp.tsfr / perturbp.deltat * ((1 - *x) / *x));
	    }
	  else
	    y = perturbp.tsfr * egyhot / tcool /
	      (perturbp.FactorSN * perturbp.EgySpecSN - (1 - perturbp.FactorSN) * All.EgySpecCold);

	  if(y < 5e-3)
	    *x = y;
	  else
	    *x = 1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y));

	  if(iiter == 0)
	    {
	      y_value[now] = *x;
	      x_value[now] = *x;
	    }
	  else
	    {
	      y_value[old] = y_value[now];
	      y_value[now] = *x;

	      x_temp = (y_value[old] * x_value[now] - y_value[now] * x_value[old]) /
		(x_value[now] - x_value[old] - y_value[now] + y_value[old]);
	      if(x_temp < 0 || x_temp > 1)
		*x = (y_value[now] + x_value[now]) / 2;
	      else
		*x = x_temp;

	      x_value[old] = x_value[now];
	      x_value[now] = *x;
	      y_value[old] = y_value[now];

	    }

	  *xrel = fabs(x_value[old] - x_value[now]) / x_value[old];
	  /*
	     xrel_II = fabs(x_value[old] - y_value[now]) / x_value[old];
	     if(xrel_II < *xrel)
	     {
	     x_value[now] = y_value[now];
	     *xrel = xrel_II;
	     }
	   */
	  iiter++;
	}
      while(*xrel > PERTURB_REL_PREC && iiter < PERTURB_MAX_ITER);
      oiter++;
    }
  while(((*x > 1) || (*x < 0)) && oiter < 2);

  if(iiter == PERTURB_MAX_ITER)
    {
      if(fabs(x_value[old] - ox) < fabs(x_value[now] - ox))
	*x = x_value[old];
      else
	*x = x_value[now];
    }

  if(oiter == 2)
    {
      if(*x < (1 + 0.01))
	{
	  printf("acc: %d %d %g\n", ThisTask, iiter, *x);
	  *x = 1.0;
	  return iiter;
	}
      if(*x > -0.01)
	{
	  printf("acc: %d %d %g\n", ThisTask, iiter, *x);
	  *x = 0;
	  return iiter;
	}

      printf("ACC! :: Task : %d has not got viable solution in perturb()\n"
	     "%g %d %d\n", ThisTask, *x, iiter, oiter);
      fflush(stdout);
      endrun(901901);
    }


#undef old
#undef now
  return iiter;
}



#endif /* refers to LT_MOD_EFFM */

#endif
