#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_sf_gamma.h>

#include "allvars.h"
#include "proto.h"
#ifdef COSMIC_RAYS
#include "cosmic_rays.h"
#endif

#ifdef SFR_METALS
#include "c_metals.h"
#endif

#ifdef MACHNUM
#ifdef COSMIC_RAYS
#define h  All.HubbleParam
#define cm (h/All.UnitLength_in_cm)
#define s  (h/All.UnitTime_in_s)
#define LightSpeed (2.9979e10*cm/s)
#define c2   ( LightSpeed * LightSpeed )
#endif
#endif



/*! \file init.c
 *  \brief code for initialisation of a simulation from initial conditions
 */


/*! This function reads the initial conditions, and allocates storage for the
 *  tree(s). Various variables of the particle data are initialised and An
 *  intial domain decomposition is performed. If SPH particles are present,
 *  the inial SPH smoothing lengths are determined.
 */
void init(void)
{
  int i, j;
  double a3, atime;

#if defined(COSMIC_RAYS) && defined(MACHNUM)
  double Pth1, PCR1, rBeta, C_phys, q_phys;
#endif
#ifdef CR_INITPRESSURE
  double cr_pressure, q_phys, C_phys;
#endif
#ifdef CHEMISTRY
  int ifunc;
  double min_t_cool, max_t_cool;
  double min_t_elec, max_t_elec;
  double a_start, a_end;
#endif

  All.Time = All.TimeBegin;

  switch (All.ICFormat)
    {
    case 1:
    case 2:
    case 3:
      read_ic(All.InitCondFile);
      break;
    case 4:
      if(RestartFlag == 2)
	read_ic(All.InitCondFile);
      else
	read_ic_cluster(All.InitCondFile);
      break;
    case 5:			/* this is a special format for bepis' cluster ICs */
      if(RestartFlag == 2)
	read_ic(All.InitCondFile);
      else
	read_ic_cluster_gas(All.InitCondFile);
      break;
    default:
      if(ThisTask == 0)
	printf("ICFormat=%d not supported.\n", All.ICFormat);
      endrun(0);
    }

  All.Time = All.TimeBegin;

#ifdef CR_DIFFUSION_GREEN
  All.TimeOfLastDiffusion = All.Time;
#endif

#ifdef COOLING
#ifndef BG_COOLING
#ifdef SFR_METALS
  if(RestartFlag == 0)
    XH = HYDROGEN_MASSFRAC;
#endif
  IonizeParams();
#endif
#endif

#ifdef CHEMISTRY
  InitChem();
#endif

  if(All.ComovingIntegrationOn)
    {
      All.Timebase_interval = (log(All.TimeMax) - log(All.TimeBegin)) / TIMEBASE;
      All.Ti_Current = 0;
      a3 = All.Time * All.Time * All.Time;
      atime = All.Time;
    }
  else
    {
      All.Timebase_interval = (All.TimeMax - All.TimeBegin) / TIMEBASE;
      All.Ti_Current = 0;
      a3 = 1;
      atime = 1;
    }

  set_softenings();

  All.NumCurrentTiStep = 0;	/* setup some counters */
  All.SnapshotFileCount = 0;
  if(RestartFlag == 2)
    All.SnapshotFileCount = atoi(All.InitCondFile + strlen(All.InitCondFile) - 3) + 1;

#ifdef OUTPUTLINEOFSIGHT
  All.Ti_nextlineofsight = (int) (log(All.TimeFirstLineOfSight / All.TimeBegin) / All.Timebase_interval);
  if(RestartFlag == 2)
    endrun(78787);
#endif

  All.TotNumOfForces = 0;
  All.NumForcesSinceLastDomainDecomp = 0;
#if defined(MAGNETIC) && defined(BSMOOTH)
  All.MainTimestepCounts = 0;
#endif

  All.Cadj_Cost = 1.0e-30;
  All.Cadj_Cpu = 1.0e-3;

  if(All.ComovingIntegrationOn)
    if(All.PeriodicBoundariesOn == 1)
      check_omega();

  All.TimeLastStatistics = All.TimeBegin - All.TimeBetStatistics;
#ifdef BLACK_HOLES
  All.TimeNextBlackHoleCheck = All.TimeBegin;
#endif



#ifdef BUBBLES
  if(All.ComovingIntegrationOn)
    All.TimeOfNextBubble = 1. / (1. + All.FirstBubbleRedshift);
  else
    All.TimeOfNextBubble = All.TimeBegin + All.BubbleTimeInterval / All.UnitTime_in_Megayears;
  
  if(ThisTask == 0)
    printf("Initial time: %g and first bubble time %g \n",
	   All.TimeBegin, All.TimeOfNextBubble);
  
  if(RestartFlag == 2 && All.TimeBegin > All.TimeOfNextBubble)
    {
      printf("Restarting from the snapshot file with the wrong FirstBubbleRedshift! \n");
      endrun(0);
    }
#endif

#ifdef MULTI_BUBBLES
  if(All.ComovingIntegrationOn)
    All.TimeOfNextBubble = 1. / (1. + All.FirstBubbleRedshift);
  else
    All.TimeOfNextBubble = All.TimeBegin + All.BubbleTimeInterval / All.UnitTime_in_Megayears;

  if(ThisTask == 0)
    printf("Initial time: %g and time of the first bubbles %g \n", All.TimeBegin, All.TimeOfNextBubble);

  if(RestartFlag == 2 && All.TimeBegin > All.TimeOfNextBubble)
    {
      printf("Restarting from the snapshot file with the wrong FirstBubbleRedshift! \n");
      endrun(0);
    }
#endif

  if(All.ComovingIntegrationOn)	/*  change to new velocity variable */
    {
      for(i = 0; i < NumPart; i++)
	for(j = 0; j < 3; j++)
	  P[i].Vel[j] *= sqrt(All.Time) * All.Time;
    }

  for(i = 0; i < NumPart; i++)	/*  start-up initialization */
    {
      for(j = 0; j < 3; j++)
	P[i].g.GravAccel[j] = 0;
#ifdef PMGRID
      for(j = 0; j < 3; j++)
	P[i].GravPM[j] = 0;
#endif
      P[i].Ti_endstep = 0;
      P[i].Ti_begstep = 0;

      P[i].OldAcc = 0;
      P[i].GravCost = 1;
      P[i].p.Potential = 0;
#ifdef STELLARAGE
      if(RestartFlag == 0)
	P[i].StellarAge = 0;
#endif


#ifdef BG_SFR
      /*
       * This should be read from a file containing the mass fractions of the elements we are
       * tracking: Hydrogen, Helium, Carbon, Nitrogen, Oxygen, Neon, Magnesium, Silicon and Iron.
       * Sulphur and Calcium mass fractions are proportional to Silicon mass fraction.
       */

      /* solar mass fractions */
      // const FLOAT mass_fractions[] = {0.70681476, 0.2806, 2.0626E-03, 8.3549E-04, 5.4951E-03, 1.4145E-03, 5.9115E-04, 6.8310E-04, 1.1039E-03};

      /* primordial mass fractions */
      const FLOAT mass_fractions[] = { 0.24, 0, 0, 0, 0, 0, 0, 0 };

      for(j = 0; j < BG_NELEMENTS; j++)
	P[i].Metals[j] = mass_fractions[j] * P[i].Mass;
#endif


#ifdef METALS
#ifndef SFR_METALS
      if(RestartFlag == 0)
	P[i].Metallicity = 0;
#else
      if(RestartFlag == 0)
	{
	  for(j = 0; j < 12; j++)
	    {
	      P[i].Zm[j] = 0;
	      P[i].ZmReservoir[j] = 0;
	    }

	  P[i].Zm[6] = HYDROGEN_MASSFRAC * (P[i].Mass);
	  P[i].Zm[0] = (1 - HYDROGEN_MASSFRAC) * (P[i].Mass);
	  /*     Order of chemical elements:   He, Carbon,Mg,O,Fe,Si,H,N,Ne,S,Ca,Zn */

#ifndef READ_HSML
	  PPP[i].Hsml = 0;
#endif
#ifdef SFR_FEEDBACK
	  PPP[i].EnergySN = 0;
	  PPP[i].EnergySNCold = 0;
#endif
	}
#endif

#endif

#ifdef BLACK_HOLES
      if(RestartFlag == 0 && P[i].Type == 5)
	P[i].BH_Mass = All.SeedBlackHoleMass;
#ifdef BH_KINETICFEEDBACK
      if(P[i].Type == 5)
	{
	  P[i].ActiveTime = 0;
	  P[i].ActiveEnergy = 0;
	}
#endif
#endif
    }

#ifdef PMGRID
  All.PM_Ti_endstep = All.PM_Ti_begstep = 0;
#endif


#ifdef FLEXSTEPS
  All.PresentMinStep = TIMEBASE;
  for(i = 0; i < NumPart; i++)	/*  start-up initialization */
    {
      P[i].FlexStepGrp = (int) (TIMEBASE * get_random_number(P[i].ID));
    }
#endif


  for(i = 0; i < N_gas; i++)	/* initialize sph_properties */
    {
      for(j = 0; j < 3; j++)
	{
	  SphP[i].VelPred[j] = P[i].Vel[j];
	  SphP[i].a.HydroAccel[j] = 0;
	}

      SphP[i].e.DtEntropy = 0;

#ifdef CHEMISTRY
      SphP[i].Gamma = GAMMA;	/* set universal value */
      SphP[i].t_cool = 0;
      SphP[i].t_elec = 0;
#endif

#ifdef SFR_DECOUPLING
      SphP[i].DensityOld = 0;
#endif
#ifdef SFR_PROMOTION
      SphP[i].DensityAvg = 0;
      SphP[i].EntropyAvg = 0;
      SphP[i].DensPromotion = 0;
      SphP[i].TempPromotion = 0;

#endif

      if(RestartFlag == 0)
	{
#ifndef READ_HSML
	  PPP[i].Hsml = 0;
#endif
	  SphP[i].a2.Density = -1;
#ifdef COOLING
	  SphP[i].Ne = 1.0;
#endif
	}
#ifdef WINDS
      SphP[i].DelayTime = 0;
#endif
#ifdef SFR
      SphP[i].Sfr = 0;
#endif
#ifdef MAGNETIC
      for(j = 0; j < 3; j++)
	{
	  SphP[i].DtB[j] = 0;
	  SphP[i].BPred[j] = SphP[i].B[j];
	}
#ifdef BINISET
      SphP[i].B[0] = All.BiniX;
      SphP[i].B[1] = All.BiniY;
      SphP[i].B[2] = All.BiniZ;
      SphP[i].BPred[0] = All.BiniX;
      SphP[i].BPred[1] = All.BiniY;
      SphP[i].BPred[2] = All.BiniZ;
#endif
#ifdef TIME_DEP_MAGN_DISP
#ifdef HIGH_MAGN_DISP_START
      SphP[i].Balpha = All.ArtMagDispConst;
#else
      SphP[i].Balpha = All.ArtMagDispMin;
#endif
      SphP[i].DtBalpha = 0.0;
#endif
#ifdef DIVBCLEANING_DEDNER
      SphP[i].Phi = SphP[i].PhiPred = SphP[i].DtPhi = 0;
#endif
#endif

#ifdef TIME_DEP_ART_VISC
#ifdef HIGH_ART_VISC_START
      if(HIGH_ART_VISC_START == 0)
	SphP[i].alpha = All.ArtBulkViscConst;
      if(HIGH_ART_VISC_START > 0)
	if(P[i].Pos[0] > HIGH_ART_VISC_START)
	  SphP[i].alpha = All.ArtBulkViscConst;
	else
	  SphP[i].alpha = All.AlphaMin;
      if(HIGH_ART_VISC_START < 0)
	if(P[i].Pos[0] < -HIGH_ART_VISC_START)
	  SphP[i].alpha = All.ArtBulkViscConst;
	else
	  SphP[i].alpha = All.AlphaMin;
#else
      SphP[i].alpha = All.AlphaMin;
#endif
      SphP[i].Dtalpha = 0.0;
#endif

#if defined(BH_THERMALFEEDBACK) || defined(BH_KINETICFEEDBACK)
      SphP[i].i.Injected_BH_Energy = 0;
#endif
    }

#ifdef HPM
  All.HPM_entr0 = All.HPM_entr1 = All.HPM_ne0 = All.HPM_ne1 = 0;
#endif

  ngb_treeallocate(MAX_NGB);

  force_treeallocate((int) (All.TreeAllocFactor * All.MaxPart), All.MaxPart);

  All.NumForcesSinceLastDomainDecomp = (long long) (1 + All.TotNumPart * All.TreeDomainUpdateFrequency);

  Flag_FullStep = 1;		/* to ensure that Peano-Hilber order is done */

  DomainDecomposition();	/* do initial domain decomposition (gives equal numbers of particles) */

  set_softenings();

  /* will build tree */
  ngb_treebuild();

  All.Ti_Current = 0;

  setup_smoothinglengths();

  TreeReconstructFlag = 1;

  /* at this point, the entropy variable actually contains the 
   * internal energy, read in from the initial conditions file. 
   * Once the density has been computed, we can convert to entropy.
   */

  for(i = 0; i < N_gas; i++)	/* initialize sph_properties */
    {
#ifdef ISOTHERM_EQS
	  if(ThisTask == 0 && i == 0)
	    printf("Using ISOTHERMAL EOS, GAMMA = %e !\n",GAMMA);
#else
      if(header.flag_entropy_instead_u == 0)
	{
	  if(ThisTask == 0 && i == 0)
	    printf("Converting u -> entropy, GAMMA = %e !\n", GAMMA);
	  SphP[i].Entropy = GAMMA_MINUS1 * SphP[i].Entropy / pow(SphP[i].a2.Density / a3, GAMMA_MINUS1);
	}
#endif
      SphP[i].e.DtEntropy = 0;

#ifdef MACHNUM
      SphP[i].Shock_MachNumber = 1.0;
#ifdef COSMIC_RAYS
      Pth1 = SphP[i].Entropy * pow(SphP[i].a2.Density / a3, GAMMA);

#ifdef CR_IC_PHYSICAL
      C_phys = SphP[i].CR_C0;
      q_phys = SphP[i].CR_q0;
#else
      C_phys = SphP[i].CR_C0 * pow(SphP[i].a2.Density, (All.CR_Alpha - 1.0) / 3.0);
      q_phys = SphP[i].CR_q0 * pow(SphP[i].a2.Density, 1.0 / 3.0);
#endif

      rBeta = gsl_sf_beta((All.CR_Alpha - 2.0) * 0.5, (3.0 - All.CR_Alpha) * 0.5) *
	gsl_sf_beta_inc((All.CR_Alpha - 2.0) * 0.5, (3.0 - All.CR_Alpha) * 0.5,
			1.0 / (1.0 + q_phys * q_phys));

      PCR1 = C_phys * c2 * SphP[i].a2.Density * rBeta / 6.0;
      PCR1 *= pow(atime, -3.0 * GAMMA);
      SphP[i].PreShock_XCR = PCR1 / Pth1;

      SphP[i].PreShock_PhysicalDensity = SphP[i].a2.Density / a3;
      SphP[i].PreShock_PhysicalEnergy =
	SphP[i].Entropy / GAMMA_MINUS1 * pow(SphP[i].a2.Density / a3, GAMMA_MINUS1);

      SphP[i].Shock_DensityJump = 1.0001;
      SphP[i].Shock_EnergyJump = 1.0;
#endif /* COSMIC_RAYS */
#endif /* MACHNUM */

#ifdef REIONIZATION
      All.not_yet_reionized = 1;
#endif

#ifdef CR_IC_PHYSICAL
      /* Scale CR variables so that values from IC file are now the
       * physical values, not the adiabatic invariants
       */

      SphP[i].CR_C0 *= pow(SphP[i].a2.Density, (1.0 - All.CR_Alpha) / 3.0);
      SphP[i].CR_q0 *= pow(SphP[i].a2.Density, -1.0 / 3.0);
#endif

#ifdef CR_INITPRESSURE

      cr_pressure = CR_INITPRESSURE * SphP[i].Entropy * pow(SphP[i].a2.Density / a3, GAMMA);
      SphP[i].Entropy *= (1 - CR_INITPRESSURE);
      q_phys = 1.685;
      C_phys =
	cr_pressure / (SphP[i].a2.Density / a3 * CR_Tab_Beta(q_phys) * (C / All.UnitVelocity_in_cm_per_s) *
		       (C / All.UnitVelocity_in_cm_per_s) / 6.0);

      SphP[i].CR_C0 = C_phys * pow(SphP[i].a2.Density, (1.0 - All.CR_Alpha) / 3.0);
      SphP[i].CR_q0 = q_phys * pow(SphP[i].a2.Density, -1.0 / 3.0);
#endif
    }


#ifdef CHEMISTRY

  if(ThisTask == 0)
    {
      printf("Initial abundances: \n");
      printf("HI=%g, HII=%g, HeI=%g, HeII=%g, HeIII=%g \n",
	     SphP[1].HI, SphP[1].HII, SphP[1].HeI, SphP[1].HeII, SphP[1].HeIII);

      printf("HM=%g, H2I=%g, H2II=%g, elec=%g, %d\n",
	     SphP[1].HM, SphP[1].H2I, SphP[1].H2II, SphP[1].elec, P[1].ID);

      printf("x=%g, y=%g, z=%g, vx=%g, vy=%g, vz=%g, density=%g, entropy=%g\n",
	     P[N_gas - 1].Pos[0], P[N_gas - 1].Pos[1], P[N_gas - 1].Pos[2], P[N_gas - 1].Vel[0],
	     P[N_gas - 1].Vel[1], P[N_gas - 1].Vel[2], SphP[N_gas - 1].Density, SphP[N_gas - 1].Entropy);
    }

  /* need predict the cooling time and elec_dot here */
  min_t_cool = min_t_elec = 1.0e30;
  max_t_cool = max_t_elec = -1.0e30;

  for(i = 0; i < N_gas; i++)
    {
      a_start = All.Time;
      a_end = All.Time + 0.001;	/* 0.001 as an arbitrary value */

      ifunc = compute_abundances(0, i, a_start, a_end);


      if(fabs(SphP[i].t_cool) < min_t_cool)
	min_t_cool = fabs(SphP[i].t_cool);
      if(fabs(SphP[i].t_cool) > max_t_cool)
	max_t_cool = fabs(SphP[i].t_cool);

      if(fabs(SphP[i].t_elec) < min_t_elec)
	min_t_elec = fabs(SphP[i].t_elec);
      if(fabs(SphP[i].t_elec) > max_t_elec)
	max_t_elec = fabs(SphP[i].t_elec);

    }

  fprintf(stdout, "PE %d t_cool min= %g, max= %g in yrs \n", ThisTask, min_t_cool, max_t_cool);
  fflush(stdout);
  fprintf(stdout, "PE %d t_elec min= %g, max= %g in yrs \n", ThisTask, min_t_elec, max_t_elec);
  fflush(stdout);

#endif



}


/*! This routine computes the mass content of the box and compares it to the
 * specified value of Omega-matter.  If discrepant, the run is terminated.
 */
void check_omega(void)
{
  double mass = 0, masstot, omega;
  int i;

  for(i = 0; i < NumPart; i++)
    mass += P[i].Mass;

  MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  omega =
    masstot / (All.BoxSize * All.BoxSize * All.BoxSize) / (3 * All.Hubble * All.Hubble / (8 * M_PI * All.G));

  if(fabs(omega - All.Omega0) > 1.0e-3)
    {
      if(ThisTask == 0)
	{
	  printf("\n\nI've found something odd!\n");
	  printf
	    ("The mass content accounts only for Omega=%g,\nbut you specified Omega=%g in the parameterfile.\n",
	     omega, All.Omega0);
	  printf("\nI better stop.\n");

	  fflush(stdout);
	}
      endrun(1);
    }
}



/*! This function is used to find an initial smoothing length for each SPH
 *  particle. It guarantees that the number of neighbours will be between
 *  desired_ngb-MAXDEV and desired_ngb+MAXDEV. For simplicity, a first guess
 *  of the smoothing length is provided to the function density(), which will
 *  then iterate if needed to find the right smoothing length.
 */
void setup_smoothinglengths(void)
{
  int i, no, p;

  if(RestartFlag == 0)
    {
      for(i = 0; i < N_gas; i++)
	{
	  no = Father[i];

	  while(10 * All.DesNumNgb * P[i].Mass > Nodes[no].u.d.mass)
	    {
	      p = Nodes[no].u.d.father;

	      if(p < 0)
		break;

	      no = p;
	    }

#ifndef READ_HSML
#ifndef TWODIMS
	  PPP[i].Hsml =
	    pow(3.0 / (4 * M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 3) * Nodes[no].len;
#else
	  PPP[i].Hsml =
	    pow(1.0 / (M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 2) * Nodes[no].len;
#endif
#endif
	}
    }

#ifdef BLACK_HOLES
  if(RestartFlag == 0 || RestartFlag == 2)
    {
      for(i = 0; i < NumPart; i++)
	if(P[i].Type == 5)
	  PPP[i].Hsml = All.SofteningTable[5];
    }
#endif

#ifdef BG_SFR
  if(RestartFlag == 0)
    {
      for(i = 0; i < NumPart; i++)
	if(P[i].Type == 4)
	  {
	    no = Father[i];

	    while(10 * All.DesNumNgb * P[i].Mass > Nodes[no].u.d.mass)
	      {
		p = Nodes[no].u.d.father;

		if(p < 0)
		  break;

		no = p;
	      }
	    PPP[i].Hsml =
	      pow(3.0 / (4 * M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 3) * Nodes[no].len;
	  }
    }
#endif

  density();

#if defined(MAGNETIC) && defined(BFROMROTA)
  if(RestartFlag == 0)
    {
      if(ThisTask == 0)
	printf("Converting: Vector Potential -> Bfield\n");
      rot_a();
    }
#endif

#if defined(MAGNETIC) && defined(BSMOOTH)
  if(RestartFlag == 0)
    bsmooth();
#endif
}
