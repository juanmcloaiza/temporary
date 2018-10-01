#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <sys/types.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>

#include "allvars.h"
#include "proto.h"
#ifdef COSMIC_RAYS
#include "cosmic_rays.h"
#endif
#ifdef BG_COOLING
#include "bg_cooling.h"
#endif

#if defined(BG_COOLING) || defined (BG_STELLAR_EVOLUTION)
#include "bg_proto.h"
#include "bg_vars.h"
#endif

/*! \file begrun.c
 *  \brief initial set-up of a simulation run
 *
 *  This file contains various functions to initialize a simulation run. In
 *  particular, the parameterfile is read in and parsed, the initial
 *  conditions or restart files are read, and global variables are initialized
 *  to their proper values.
 */



/*! This function performs the initial set-up of the simulation. First, the
 *  parameterfile is set, then routines for setting units, reading
 *  ICs/restart-files are called, auxialiary memory is allocated, etc.
 */
void begrun(void)
{
  struct global_data_all_processes all;

#ifdef LT_METAL_EFF_MODEL
  int i;
#endif

  if(ThisTask == 0)
    {
      //      printf("\nThis is P-Gadget, version `%s', svn-revision `%s'.\n", GADGETVERSION, svn_version());
      printf("\nThis is P-Gadget, version `%s'.\n", GADGETVERSION);
      printf("\nRunning on %d processors.\n", NTask);
    }

#if defined(X86FIX) && defined(SOFTDOUBLEDOUBLE)
  x86_fix();			/* disable 80bit treatment of internal FPU registers in favour of proper IEEE 64bit double precision arithmetic */
#endif

#ifdef PEDANTIC_MEMORY_HANDLER
  mymalloc_init(PEDANTIC_MEMORY_CEILING * 1024.0 * 1024.0);
#endif

  read_parameter_file(ParameterFile);	/* ... read in parameters for this run */

#ifdef DEBUG
  write_pid_file();
  enable_core_dumps_and_fpu_exceptions();
#endif

  allocate_commbuffers();	/* ... allocate buffer-memory for particle 
				   exchange during force computation */
#ifdef VOLUME_CORRECTION
  vol_weights_init();
#endif

  set_units();

#if defined(BG_COOLING) || defined(BG_STELLAR_EVOLUTION)
  InitChemistry();

  if(ThisTask == 0)
    {
      printf("Chemistry initialized.\n");
      fflush(stdout);
    }

#endif

#ifdef BG_STELLAR_EVOLUTION
  init_imf();

  if(ThisTask == 0)
    {
      printf("IMF initialized.\n");
      fflush(stdout);
    }

  init_yields();

  if(ThisTask == 0)
    {
      printf("Yields initialized.\n");
      fflush(stdout);
    }
#endif

#ifdef COOLING
  All.Time = All.TimeBegin;
  InitCool();
#endif

#ifdef CHEMISTRY
  InitChem();
#endif

#if defined(SFR) && !defined(LT_STELLAREVOLUTION)
#ifndef SFR_METALS
#ifndef BG_SFR
  init_clouds();
#endif
#else
  Flag_phase = 0;
  Flag_promotion = 0;
#endif
  Stars_converted = 0;
#endif

#ifdef PERIODIC
  ewald_init();
#ifdef LONG_X
  boxSize_X = All.BoxSize * LONG_X;
  boxHalf_X = 0.5 * All.BoxSize * LONG_X;
#else
  boxSize_X = All.BoxSize;
  boxHalf_X = 0.5 * All.BoxSize;
#endif
#ifdef LONG_Y
  boxSize_Y = All.BoxSize * LONG_Y;
  boxHalf_Y = 0.5 * All.BoxSize * LONG_Y;
#else
  boxSize_Y = All.BoxSize;
  boxHalf_Y = 0.5 * All.BoxSize;
#endif
#ifdef LONG_Z
  boxSize_Z = All.BoxSize * LONG_Z;
  boxHalf_Z = 0.5 * All.BoxSize * LONG_Z;
#else
  boxSize_Z = All.BoxSize;
  boxHalf_Z = 0.5 * All.BoxSize;
#endif
#endif

#ifdef DARKENERGY
#ifdef TIMEDEPDE
  fwa_init();
#endif
#endif

#ifdef TIME_DEP_ART_VISC
#ifdef ISOTHERM_EQS
  All.ViscSource = All.ViscSource0;
  All.DecayTime = 0.1;
  if(ThisTask == 0){
    printf("TIME_DEP_ART_VISC + ISOTHERM_EQS, ignoring ViscosityDecayLength in parameterfile\n");
    printf("All.ViscSource = %e, All.DecayTime = %e\n", All.ViscSource, All.DecayTime);
  }

#else
  All.ViscSource = All.ViscSource0 / log((GAMMA + 1) / (GAMMA - 1));
  All.DecayTime = 1.0 / All.DecayLength * sqrt((GAMMA - 1) / (2 * GAMMA));
  if(ThisTask == 0){
    printf("TIME_DEP_ART_VISC \n");
    printf("All.ViscSource = %e, All.DecayTime = %e\n", All.ViscSource, All.DecayTime);
  }
#endif
#endif

  open_outputfiles();

#ifdef BG_STELLAR_EVOLUTION
  if(ThisTask == 0)
    write_outputfiles_header();
#endif

  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

  gsl_rng_set(random_generator, 42);	/* start-up seed */

#ifdef PMGRID
  long_range_init();
#endif

  All.TimeLastRestartFile = CPUThisRun;

  if(RestartFlag == 0 || RestartFlag == 2)
    {
      set_random_numbers();

      init();			/* ... read in initial model */
    }
  else
    {
      all = All;		/* save global variables. (will be read from restart file) */

      restart(RestartFlag);	/* ... read restart file. Note: This also resets 
				   all variables in the struct `All'. 
				   However, during the run, some variables in the parameter
				   file are allowed to be changed, if desired. These need to 
				   copied in the way below.
				   Note:  All.PartAllocFactor is treated in restart() separately.  
				 */

      All.MinSizeTimestep = all.MinSizeTimestep;
      All.MaxSizeTimestep = all.MaxSizeTimestep;
      All.BufferSize = all.BufferSize;
      All.BunchSizeForce = all.BunchSizeForce;
      All.BunchSizeDensity = all.BunchSizeDensity;
      All.BunchSizeHydro = all.BunchSizeHydro;
      All.BunchSizeDomain = all.BunchSizeDomain;
#ifdef SFR_METALS
      All.BunchSizeMetal = all.BunchSizeMetal;
#ifdef SFR_FEEDBACK
      All.BunchSizeHotNgbs = all.BunchSizeHotNgbs;
#endif
#endif
#ifdef BLACK_HOLES
      All.BunchSizeBlackhole = all.BunchSizeBlackhole;
#endif
#ifdef FOF
      All.BunchSizeFoF = all.BunchSizeFoF;
#endif
#ifdef MHM
      All.BunchSizeKinetic = all.BunchSizeKinetic;
#endif
      All.TimeLimitCPU = all.TimeLimitCPU;
      All.ResubmitOn = all.ResubmitOn;
      All.TimeBetSnapshot = all.TimeBetSnapshot;
      All.TimeBetStatistics = all.TimeBetStatistics;
      All.CpuTimeBetRestartFile = all.CpuTimeBetRestartFile;
      All.ErrTolIntAccuracy = all.ErrTolIntAccuracy;
      All.MaxRMSDisplacementFac = all.MaxRMSDisplacementFac;

      All.ErrTolForceAcc = all.ErrTolForceAcc;
      All.TypeOfTimestepCriterion = all.TypeOfTimestepCriterion;
      All.TypeOfOpeningCriterion = all.TypeOfOpeningCriterion;
      All.NumFilesWrittenInParallel = all.NumFilesWrittenInParallel;
      All.TreeDomainUpdateFrequency = all.TreeDomainUpdateFrequency;

      All.OutputListOn = all.OutputListOn;
      All.CourantFac = all.CourantFac;

      All.OutputListLength = all.OutputListLength;
      memcpy(All.OutputListTimes, all.OutputListTimes, sizeof(double) * All.OutputListLength);

#ifdef TIME_DEP_ART_VISC
      All.ViscSource = all.ViscSource;
      All.ViscSource0 = all.ViscSource0;
      All.DecayTime = all.DecayTime;
      All.DecayLength = all.DecayLength;
      All.AlphaMin = all.AlphaMin;
#endif

#ifdef MAGNETIC_DISSIPATION
      All.ArtMagDispConst = all.ArtMagDispConst;
#ifdef TIME_DEP_MAGN_DISP
      All.ArtMagDispMin = all.ArtMagDispMin;
      All.ArtMagDispSource = all.ArtMagDispSource;
      All.ArtMagDispTime = all.ArtMagDispTime;
#endif
#endif

#ifdef DIVBCLEANING_DEDNER
      All.DivBcleanParabolicSigma = all.DivBcleanParabolicSigma;
      All.DivBcleanHyperbolicSigma = all.DivBcleanHyperbolicSigma;
#endif

#ifdef DARKENERGY
      All.DarkEnergyParam = all.DarkEnergyParam;
#endif

      strcpy(All.ResubmitCommand, all.ResubmitCommand);
      strcpy(All.OutputListFilename, all.OutputListFilename);
      strcpy(All.OutputDir, all.OutputDir);
      strcpy(All.RestartFile, all.RestartFile);
      strcpy(All.EnergyFile, all.EnergyFile);
      strcpy(All.InfoFile, all.InfoFile);
      strcpy(All.CpuFile, all.CpuFile);
      strcpy(All.TimingsFile, all.TimingsFile);
      strcpy(All.SnapshotFileBase, all.SnapshotFileBase);

      if(All.TimeMax != all.TimeMax)
	readjust_timebase(All.TimeMax, all.TimeMax);

#ifdef NO_TREEDATA_IN_RESTART
      /* if this is not activated, the tree was stored in the restart-files,
         which also allocated the storage for it already */

      ngb_treeallocate(MAX_NGB);
      force_treeallocate(All.TreeAllocFactor * All.MaxPart, All.MaxPart);

      /* ensures that domain reconstruction will be done and new tree will be constructed */
      All.NumForcesSinceLastDomainDecomp = 1 + All.TotNumPart * All.TreeDomainUpdateFrequency;
#endif
    }

#ifdef PMGRID
  long_range_init_regionsize();
#endif


#ifdef COSMIC_RAYS
  CR_initialize_beta_tabs(All.CR_Alpha);
  CR_Tab_Initialize();
#ifdef COSMIC_RAY_TEST
  CR_test_routine();
#endif

#endif

  if(All.ComovingIntegrationOn)
    init_drift_table();

  if(RestartFlag == 2)
    All.Ti_nextoutput = find_next_outputtime(All.Ti_Current + 1);
  else
    All.Ti_nextoutput = find_next_outputtime(All.Ti_Current);


  All.TimeLastRestartFile = CPUThisRun;
}




/*! Computes conversion factors between internal code units and the
 *  cgs-system.
 */
void set_units(void)
{
  double meanweight;

#ifdef CONDUCTION
#ifndef CONDUCTION_CONSTANT
  double coulomb_log;
#endif
#endif
#ifdef STATICNFW
  double Mtot;
#endif

  All.UnitTime_in_s = All.UnitLength_in_cm / All.UnitVelocity_in_cm_per_s;
  All.UnitTime_in_Megayears = All.UnitTime_in_s / SEC_PER_MEGAYEAR;

  if(All.GravityConstantInternal == 0)
    All.G = GRAVITY / pow(All.UnitLength_in_cm, 3) * All.UnitMass_in_g * pow(All.UnitTime_in_s, 2);
  else
    All.G = All.GravityConstantInternal;
#ifdef TIMEDEPGRAV
  All.Gini=All.G;
  All.G = All.Gini * dGfak(All.TimeBegin);
#endif

  All.UnitDensity_in_cgs = All.UnitMass_in_g / pow(All.UnitLength_in_cm, 3);
  All.UnitPressure_in_cgs = All.UnitMass_in_g / All.UnitLength_in_cm / pow(All.UnitTime_in_s, 2);
  All.UnitCoolingRate_in_cgs = All.UnitPressure_in_cgs / All.UnitTime_in_s;
  All.UnitEnergy_in_cgs = All.UnitMass_in_g * pow(All.UnitLength_in_cm, 2) / pow(All.UnitTime_in_s, 2);

  /* convert some physical input parameters to internal units */

  All.Hubble = HUBBLE * All.UnitTime_in_s;

  if(ThisTask == 0)
    {
      printf("\nHubble (internal units) = %g\n", All.Hubble);
      printf("G (internal units) = %g\n", All.G);
      printf("UnitMass_in_g = %g \n", All.UnitMass_in_g);
      printf("UnitTime_in_s = %g \n", All.UnitTime_in_s);
      printf("UnitVelocity_in_cm_per_s = %g \n", All.UnitVelocity_in_cm_per_s);
      printf("UnitDensity_in_cgs = %g \n", All.UnitDensity_in_cgs);
      printf("UnitEnergy_in_cgs = %g \n", All.UnitEnergy_in_cgs);
#ifdef MYSWITCH
      printf("Black Hole Mass (internal units) = %g \n", my_Mbh);
#endif
      printf("\n");
    }

  meanweight = 4.0 / (1 + 3 * HYDROGEN_MASSFRAC);	/* note: assuming NEUTRAL GAS */

#ifdef ISOTHERM_EQS
  All.MinEgySpec = 0;
#else
  All.MinEgySpec = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.MinGasTemp;
  All.MinEgySpec *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;
#endif

#ifdef LT_STELLAREVOLUTION
  All.MinChemTimeStep /= 86400 * 365 * 1e9;
  /* 
     Add here Spreading Length in the case of fixed distance
   */
  All.SnIaEgy /= All.UnitEnergy_in_cgs;
  All.SnIIEgy /= All.UnitEnergy_in_cgs;

#ifdef LT_HOT_EJECTA
  /* All.EgySpecEjecta is supposed to be given in paramfile as Km/sec        */
  /* then we have to transform in erg / g through the 1e10 conversion factor */
  All.EgySpecEjecta = pow(All.EgySpecEjecta, 2) * 1e10 * All.UnitMass_in_g / All.UnitEnergy_in_cgs;
#endif

  Hyd = LT_NMet - 1;
  All.Time = All.TimeBegin;
  init_SN();
#endif

#if defined(SFR) && !defined(LT_STELLAREVOLUTION)
  set_units_sfr();
#endif


#define cm (All.HubbleParam/All.UnitLength_in_cm)
#define g  (All.HubbleParam/All.UnitMass_in_g)
#define s  (All.HubbleParam/All.UnitTime_in_s)
#define erg (g*cm*cm/(s*s))
#define keV (1.602e-9*erg)
#define deg 1.0
#define m_p (PROTONMASS * g)
#define k_B (BOLTZMANN * erg / deg)

#ifdef NAVIERSTOKES
  /* Braginskii-Spitzer shear viscosity parametrization */
  /* mu = 0.406 * m_p^0.5 * (k_b* T)^(5/2) / e^4 / logLambda  [g/cm/s]*/
  /* eta = frac * mu */

  meanweight = 4.0 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));  /* assuming full ionization */

#ifdef NAVIERSTOKES_CONSTANT  
  All.NavierStokes_ShearViscosity = All.FractionSpitzerViscosity * 
    0.406 * pow(PROTONMASS, 0.5) * 
    pow(BOLTZMANN * All.ShearViscosityTemperature, 5./2.)
    / pow(ELECTRONCHARGE, 4) / LOG_LAMBDA;  /*in cgs units*/
  
  if(ThisTask == 0)
    printf("Constant shear viscosity in cgs units: eta = %g\n", All.NavierStokes_ShearViscosity);
  
  All.NavierStokes_ShearViscosity *= All.UnitTime_in_s * All.UnitLength_in_cm 
    / All.UnitMass_in_g / All.HubbleParam; /* in internal code units */
  
  if(ThisTask == 0)
    printf("Constant shear viscosity in internal code units: eta = %g\n", All.NavierStokes_ShearViscosity);

#else
  All.NavierStokes_ShearViscosity = All.FractionSpitzerViscosity * 
    0.406 * pow(PROTONMASS, 0.5) * 
    pow((meanweight * PROTONMASS * GAMMA_MINUS1), 5./2.)
    / pow(ELECTRONCHARGE, 4) / LOG_LAMBDA;  /*in cgs units*/ 
  /*T = mu*m_p*(gamma-1)/k_b * E * UnitEnergy/UnitMass */

  All.NavierStokes_ShearViscosity *= pow( (All.UnitEnergy_in_cgs/All.UnitMass_in_g), 5./2.); /* now energy can be multiplied later in the internal code units */
  All.NavierStokes_ShearViscosity *= All.UnitTime_in_s * All.UnitLength_in_cm 
    / All.UnitMass_in_g / All.HubbleParam; /* in internal code units */

  if(ThisTask == 0)
    printf("Variable shear viscosity in internal code units: eta = %g\n", All.NavierStokes_ShearViscosity);
  
#endif

#ifdef NAVIERSTOKES_BULK
 if(ThisTask == 0)
    printf("Costant bulk viscosity in internal code units: zeta = %g\n", All.NavierStokes_BulkViscosity);
#endif

#ifdef VISCOSITY_SATURATION
 /* calculate ion mean free path assuming complete ionization: 
    ion mean free path for hydrogen is similar to that of helium, 
    thus we calculate only for hydrogen */
 /* l_i = 3^(3/2)*(k*T)^2 / (4*\pi^(1/2)*ni*(Z*e)^4*lnL) */

 All.IonMeanFreePath = pow(3.0, 1.5) / 
   (4.0 * sqrt(M_PI) * pow(ELECTRONCHARGE, 4) * LOG_LAMBDA);

 All.IonMeanFreePath *= pow(meanweight * PROTONMASS * GAMMA_MINUS1, 2) * 
   pow((All.UnitEnergy_in_cgs/All.UnitMass_in_g), 2); /*kT -> u*/

 All.IonMeanFreePath /= (HYDROGEN_MASSFRAC / PROTONMASS) * 
   (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam); 
 /* n_H = rho * Hfr / mp */ /* now is cgs units */ //changed / to * in front of the unitdensity
 
 All.IonMeanFreePath *= All.HubbleParam / All.UnitLength_in_cm;
 /* in internal code units */
#endif

#endif

#ifdef CONDUCTION
#ifndef CONDUCTION_CONSTANT

  meanweight = m_p * 4.0 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));
  /* assuming full ionization */

  coulomb_log = 37.8;
  /* accordin1g to Sarazin's book */

  All.ConductionCoeff *=
    (1.84e-5 / coulomb_log * pow(meanweight / k_B * GAMMA_MINUS1, 2.5) * erg / (s * deg * cm));
  /* Kappa_Spitzer definition taken from Zakamska & Narayan 2003 
   * ( ApJ 582:162-169, Eq. (5) )
   */

  /* Note: Because we replace \nabla(T) in the conduction equation with
   * \nable(u), our conduction coefficient is not the usual kappa, but
   * rather kappa*(gamma-1)*mu/kB. We therefore need to multiply with 
   * another factor of (meanweight / k_B * GAMMA_MINUS1).
   */
  All.ConductionCoeff *= meanweight / k_B * GAMMA_MINUS1;

  /* The conversion of  ConductionCoeff between internal units and cgs
   * units involves one factor of 'h'. We take care of this here.
   */
  All.ConductionCoeff /= All.HubbleParam;

#ifdef CONDUCTION_SATURATION
  All.ElectronFreePathFactor = 8 * pow(3.0, 1.5) * pow(GAMMA_MINUS1, 2) / pow(3 + 5 * HYDROGEN_MASSFRAC, 2)
    / (1 + HYDROGEN_MASSFRAC) / sqrt(M_PI) / coulomb_log * pow(PROTONMASS, 3) / pow(ELECTRONCHARGE, 4)
    / (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam)
    * pow(All.UnitPressure_in_cgs / All.UnitDensity_in_cgs, 2);

  /* If the above value is multiplied with u^2/rho in code units (with rho being the physical density), then
   * one gets the electrong mean free path in centimeter. Since we want to compare this with another length
   * scale in code units, we now add an additional factor to convert back to code units.
   */
  All.ElectronFreePathFactor *= All.HubbleParam / All.UnitLength_in_cm;
#endif

#endif /* CONDUCTION_CONSTANT */
#endif /* CONDUCTION */

#if defined(CR_DIFFUSION) || defined(CR_DIFFUSION_GREEN)
  if(All.CR_DiffusionDensZero == 0.0)
    {
      /* Set reference density for CR Diffusion to rhocrit at z=0 */

      All.CR_DiffusionDensZero = 3.0 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
    }

  if(All.CR_DiffusionEntropyZero == 0.0)
    {
      All.CR_DiffusionEntropyZero = 1.0e4;
    }

  /* Set reference entropic function to correspond to 
     Reference Temperature @ ReferenceDensity */
  if(ThisTask == 0)
    {
      printf("CR Diffusion: T0 = %g\n", All.CR_DiffusionEntropyZero);
    }

  /* convert Temperature value in Kelvin to thermal energy per unit mass
     in internal units, and then to entropy */
  All.CR_DiffusionEntropyZero *=
    BOLTZMANN / (4.0 * PROTONMASS / (3.0 * HYDROGEN_MASSFRAC + 1.0)) *
    All.UnitMass_in_g / All.UnitEnergy_in_cgs / pow(All.CR_DiffusionDensZero, GAMMA_MINUS1);

  /* Change the density scaling, so that the temp scaling is mapped
     onto an entropy scaling that is numerically less expensive */
  All.CR_DiffusionDensScaling += GAMMA_MINUS1 * All.CR_DiffusionEntropyScaling;

  if(ThisTask == 0)
    {
      printf("CR Diffusion: Rho0 = %g -- A0 = %g\n", All.CR_DiffusionDensZero, All.CR_DiffusionEntropyZero);
    }

  if(All.CR_DiffusionTimeScale <= 0.0)
    {
      All.CR_DiffusionTimeScale = 0.1;
    }

#endif /* CR_DIFFUSION */

#ifdef STATICNFW
  R200 = pow(NFW_M200 * All.G / (100 * All.Hubble * All.Hubble), 1.0 / 3);
  Rs = R200 / NFW_C;
  Dc = 200.0 / 3 * NFW_C * NFW_C * NFW_C / (log(1 + NFW_C) - NFW_C / (1 + NFW_C));
  RhoCrit = 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
  V200 = 10 * All.Hubble * R200;
  if(ThisTask == 0)
    printf("V200= %g\n", V200);

  fac = 1.0;
  Mtot = enclosed_mass(R200);
  if(ThisTask == 0)
    printf("M200= %g\n", Mtot);
  fac = V200 * V200 * V200 / (10 * All.G * All.Hubble) / Mtot;
  Mtot = enclosed_mass(R200);
  if(ThisTask == 0)
    printf("M200= %g\n", Mtot);
#endif
}

#ifdef STATICNFW
/*! auxiliary function for static NFW potential
 */
double enclosed_mass(double R)
{
  /* Eps is in units of Rs !!!! */

  if(R > Rs * NFW_C)
    R = Rs * NFW_C;

  return fac * 4 * M_PI * RhoCrit * Dc *
    (-(Rs * Rs * Rs * (1 - NFW_Eps + log(Rs) - 2 * NFW_Eps * log(Rs) + NFW_Eps * NFW_Eps * log(NFW_Eps * Rs)))
     / ((NFW_Eps - 1) * (NFW_Eps - 1)) +
     (Rs * Rs * Rs *
      (Rs - NFW_Eps * Rs - (2 * NFW_Eps - 1) * (R + Rs) * log(R + Rs) +
       NFW_Eps * NFW_Eps * (R + Rs) * log(R + NFW_Eps * Rs))) / ((NFW_Eps - 1) * (NFW_Eps - 1) * (R + Rs)));
}
#endif



/*!  This function opens various log-files that report on the status and
 *   performance of the simulstion. On restart from restart-files
 *   (start-option 1), the code will append to these files.
 */
void open_outputfiles(void)
{
#ifdef LT_STELLAREVOLUTION
  int k;
#endif
  char mode[2], buf[200];

  if(RestartFlag == 0)
    strcpy(mode, "w");
  else
    strcpy(mode, "a");

#ifdef BLACK_HOLES
  /* Note: This is done by everyone */
  sprintf(buf, "%sblackhole_details_%d.txt", All.OutputDir, ThisTask);
  if(!(FdBlackHolesDetails = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

#ifdef SFR_PROMOTION
  sprintf(buf, "%spromotion_%d.txt", All.OutputDir, ThisTask);
  if(!(FdPromotion = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

  if(ThisTask != 0)		/* only the root processors writes to the log files */
    return;

//Added by JM
  char my_filename[200] = "Accretion.txt";
  sprintf(buf, "%s%s", All.OutputDir,my_filename);
  if(!(FdAccretion = fopen(buf, mode)))
  {
    printf("error in opening file '%s'\n", buf);
    endrun(1);
  }
//End of Added by JM

  sprintf(buf, "%s%s", All.OutputDir, All.CpuFile);
  if(!(FdCPU = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, All.InfoFile);
  if(!(FdInfo = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, All.EnergyFile);
  if(!(FdEnergy = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, All.TimingsFile);
  if(!(FdTimings = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, "balance.txt");
  if(!(FdBalance = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  fprintf(FdBalance, "\n");
  fprintf(FdBalance, "Treewalk       = '%c' / '%c'\n", CPU_Symbol[CPU_TREEWALK],
	  CPU_SymbolImbalance[CPU_TREEWALK]);
  fprintf(FdBalance, "Treebuild      = '%c' / '%c'\n", CPU_Symbol[CPU_TREEBUILD],
	  CPU_SymbolImbalance[CPU_TREEBUILD]);
  fprintf(FdBalance, "Treeupdate     = '%c' / '%c'\n", CPU_Symbol[CPU_TREEUPDATE],
	  CPU_SymbolImbalance[CPU_TREEUPDATE]);
  fprintf(FdBalance, "Treehmaxupdate = '%c' / '%c'\n", CPU_Symbol[CPU_TREEHMAXUPDATE],
	  CPU_SymbolImbalance[CPU_TREEHMAXUPDATE]);
  fprintf(FdBalance, "Treecomm       = '%c' / '%c'\n", CPU_Symbol[CPU_TREECOMM],
	  CPU_SymbolImbalance[CPU_TREECOMM]);
  fprintf(FdBalance, "Domain decomp  = '%c' / '%c'\n", CPU_Symbol[CPU_DOMAIN],
	  CPU_SymbolImbalance[CPU_DOMAIN]);
  fprintf(FdBalance, "Density estim  = '%c' / '%c'\n", CPU_Symbol[CPU_DENSITY],
	  CPU_SymbolImbalance[CPU_DENSITY]);
  fprintf(FdBalance, "Hydro forces   = '%c' / '%c'\n", CPU_Symbol[CPU_HYDRA], CPU_SymbolImbalance[CPU_HYDRA]);
  fprintf(FdBalance, "Drifts         = '%c' / '%c'\n", CPU_Symbol[CPU_DRIFT], CPU_SymbolImbalance[CPU_DRIFT]);
  fprintf(FdBalance, "Moves          = '%c' / '%c'\n", CPU_Symbol[CPU_MOVE], CPU_SymbolImbalance[CPU_MOVE]);
  fprintf(FdBalance, "Kicks          = '%c' / '%c'\n", CPU_Symbol[CPU_TIMELINE],
	  CPU_SymbolImbalance[CPU_TIMELINE]);
  fprintf(FdBalance, "Potential      = '%c' / '%c'\n", CPU_Symbol[CPU_POTENTIAL],
	  CPU_SymbolImbalance[CPU_POTENTIAL]);
  fprintf(FdBalance, "PM             = '%c' / '%c'\n", CPU_Symbol[CPU_MESH], CPU_SymbolImbalance[CPU_MESH]);
  fprintf(FdBalance, "Peano-Hilbert  = '%c' / '%c'\n", CPU_Symbol[CPU_PEANO], CPU_SymbolImbalance[CPU_PEANO]);
  fprintf(FdBalance, "Cooling & SFR  = '%c' / '%c'\n", CPU_Symbol[CPU_COOLINGSFR],
	  CPU_SymbolImbalance[CPU_COOLINGSFR]);
  fprintf(FdBalance, "Snapshot dump  = '%c' / '%c'\n", CPU_Symbol[CPU_SNAPSHOT],
	  CPU_SymbolImbalance[CPU_SNAPSHOT]);
  fprintf(FdBalance, "FoF            = '%c' / '%c'\n", CPU_Symbol[CPU_FOF], CPU_SymbolImbalance[CPU_FOF]);
  fprintf(FdBalance, "Comm. Grav.    = '%c' / '%c'\n", CPU_Symbol[CPU_GRAVCOMM],
	  CPU_SymbolImbalance[CPU_GRAVCOMM]);
  fprintf(FdBalance, "Comm. Dens.    = '%c' / '%c'\n", CPU_Symbol[CPU_DENSCOMM],
	  CPU_SymbolImbalance[CPU_DENSCOMM]);
  fprintf(FdBalance, "Comm. Hydro    = '%c' / '%c'\n", CPU_Symbol[CPU_HYDCOMM],
	  CPU_SymbolImbalance[CPU_HYDCOMM]);
  fprintf(FdBalance, "\n");

#ifdef SFR
  sprintf(buf, "%s%s", All.OutputDir, "sfr.txt");
  if(!(FdSfr = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

#ifdef BG_STELLAR_EVOLUTION
  sprintf(buf, "%s%s", All.OutputDir, "metals_tot.txt");
  if(!(FdMetTot = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
  sprintf(buf, "%s%s", All.OutputDir, "metals_gas.txt");
  if(!(FdMetGas = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
  sprintf(buf, "%s%s", All.OutputDir, "metals_stars.txt");
  if(!(FdMetStars = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
  sprintf(buf, "%s%s", All.OutputDir, "metals_sf.txt");
  if(!(FdMetSF = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

#ifdef BLACK_HOLES
  sprintf(buf, "%s%s", All.OutputDir, "blackholes.txt");
  if(!(FdBlackHoles = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

#ifdef LT_EXTEGY_INFO
  sprintf(buf, "%s%s", All.OutputDir, "extegy.txt");
  if(!(FdExtEgy = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

#ifdef LT_CharT_INFO
  sprintf(buf, "%s%s", All.OutputDir, "chart.txt");
  if(!(FdCharT = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

#ifdef LT_STELLAREVOLUTION

#ifdef LT_SEv_INFO

#ifdef LT_SEvDbg
  sprintf(buf, "%s%s", All.OutputDir, "met_sumcheck.txt");
  if(!(FdMetSumCheck = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

  sprintf(buf, "%s%s", All.OutputDir, "metals.txt");
  if(!(FdMetals = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, "sn.txt");
  if(!(FdSn = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, "sn_lost.txt");
  if(!(FdSnLost = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

#ifdef WINDS
  sprintf(buf, "%s%s", All.OutputDir, "winds.txt");
  if(!(FdWinds = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

#endif
#endif


#ifdef SFR_METALS
/*  sprintf(buf, "%s%s", All.OutputDir, "mphase.txt"); */
  sprintf(buf, "%s%s", All.OutputDir, "energy_test.txt");
  if(!(FdMphase = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

  sprintf(buf, "%s%s", All.OutputDir, "SNenergy.txt");
  if(!(FdSNE = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }

#if defined(SFR_SNI) || defined(SFR_SNII)
  sprintf(buf, "%s%s", All.OutputDir, "SN.txt");
  if(!(FdSN = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif
#ifdef SFR_PROMOTION
  sprintf(buf, "%spromotion_%d.txt", All.OutputDir, ThisTask);
  if(!(FdPromotion = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
#endif

#endif


#ifdef FORCETEST
  if(RestartFlag == 0)
    {
      sprintf(buf, "%s%s", All.OutputDir, "forcetest.txt");
      if(!(FdForceTest = fopen(buf, "w")))
	{
	  printf("error in opening file '%s'\n", buf);
	  endrun(1);
	}
      fclose(FdForceTest);
    }
#endif

#ifdef XXLINFO
  sprintf(buf, "%s%s", All.OutputDir, "xxl.txt");
  if(!(FdXXL = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
  else
    {
      if(RestartFlag == 0)
	{
	  fprintf(FdXXL, "nstep time ");
#ifdef MAGNETIC
	  fprintf(FdXXL, "<|B|> ");
#ifdef TRACEDIVB
	  fprintf(FdXXL, "max(divB) ");
#endif
#endif
#ifdef TIME_DEP_ART_VISC
	  fprintf(FdXXL, "<alpha> ");
#endif
	  fprintf(FdXXL, "\n");
	  fflush(FdXXL);
	}
    }
#endif

#ifdef DARKENERGY
  sprintf(buf, "%s%s", All.OutputDir, "darkenergy.txt");
  if(!(FdDE = fopen(buf, mode)))
    {
      printf("error in opening file '%s'\n", buf);
      endrun(1);
    }
  else
    {
      if(RestartFlag == 0)
	{
	  fprintf(FdDE, "nstep time H(a) ");
#ifndef TIMEDEPDE
	  fprintf(FdDE, "w0 Omega_L ");
#else
	  fprintf(FdDE, "w(a) Omega_L ");
#endif
#ifdef TIMEDEPGRAV
	  fprintf(FdDE, "dH dG ");
#endif
	  fprintf(FdDE, "\n");
	  fflush(FdDE);
	}
    }
#endif

}


void write_outputfiles_header(void)
{
#ifdef BG_STELLAR_EVOLUTION
  int i;

  if(RestartFlag == 0)
    {
      fprintf(FdMetGas, "%12s", "# Time");
      fprintf(FdMetStars, "%12s", "# Time");
      fprintf(FdMetSF, "%12s", "# Time");
      fprintf(FdMetTot, "%12s", "# Time");

      for(i = 0; i < BG_NELEMENTS; i++)
	{
	  fprintf(FdMetGas, "%13s", SPH_Element_Names[i]);
	  fprintf(FdMetStars, "%13s", SPH_Element_Names[i]);
	  fprintf(FdMetSF, "%13s", SPH_Element_Names[i]);
	  fprintf(FdMetTot, "%13s", SPH_Element_Names[i]);
	}

      fprintf(FdMetGas, "%13s\n", "Total Mass");
      fprintf(FdMetStars, "%13s\n", "Total Mass");
      fprintf(FdMetSF, "%13s\n", "Total Mass");
      fprintf(FdMetTot, "%13s\n", "Total Mass");
    }
#endif
}


/*!  This function closes the global log-files.
 */
void close_outputfiles(void)
{
#ifdef BLACK_HOLES
  fclose(FdBlackHolesDetails);	/* needs to be done by everyone */
#endif

#ifdef SFR_PROMOTION
  fclose(FdPromotion);
#endif

  if(ThisTask != 0)		/* only the root processors writes to the log files */
    return;

//Added by JM
  fclose(FdAccretion);
//End of added by JM
  fclose(FdCPU);
  fclose(FdInfo);
  fclose(FdEnergy);
  fclose(FdTimings);
  fclose(FdBalance);

#ifdef LT_EXTEGY_INFO
  fclose(FdExtEgy);
#endif

#ifdef LT_CharT_INFO
  fclose(FdCharT);
#endif

#ifdef LT_EXTEGY_INFO
  fclose(FdExtEgy);
#endif

#ifdef LT_CharT_INFO
  fclose(FdCharT);
#endif

#ifdef LT_STELLAREVOLUTION

#ifdef LT_SEv_INFO
#ifdef LT_SEvDbg
  fclose(FdMetSumCheck);
#endif
  fclose(FdSn);
  fclose(FdSnLost);
  fclose(FdMetals);
#ifdef WINDS
  fclose(FdWinds);
#endif
#endif
#endif

#ifdef SFR
  fclose(FdSfr);
#endif

#ifdef BG_STELLAR_EVOLUTION
  fclose(FdMetGas);
  fclose(FdMetStars);
  fclose(FdMetSF);
  fclose(FdMetTot);
#endif

#ifdef BLACK_HOLES
  fclose(FdBlackHoles);
#endif
#ifdef XXLINFO
  fclose(FdXXL);
#endif
#ifdef SFR_METALS
  fclose(FdMphase);
  fclose(FdSNE);
#if defined(SFR_SNI) || defined(SFR_SNII)
  fclose(FdSN);
#endif
#ifdef SFR_PROMOTION
  fclose(FdPromotion);
#endif
#endif
}





/*! This function parses the parameterfile in a simple way.  Each paramater is
 *  defined by a keyword (`tag'), and can be either of type douple, int, or
 *  character string.  The routine makes sure that each parameter appears
 *  exactly once in the parameterfile, otherwise error messages are
 *  produced that complain about the missing parameters.
 */
void read_parameter_file(char *fname)
{
#define REAL 1
#define STRING 2
#define INT 3
#define MAXTAGS 300

  FILE *fd, *fdout;
  char buf[200], buf1[200], buf2[200], buf3[400];
  int i, j, nt;
  int id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  int pnum, errorFlag = 0;

  All.StarformationOn = 0;	/* defaults */


  if(sizeof(long long) != 8)
    {
      if(ThisTask == 0)
	printf("\nType `long long' is not 64 bit on this platform. Stopping.\n\n");
      endrun(0);
    }

  if(sizeof(int) != 4)
    {
      if(ThisTask == 0)
	printf("\nType `int' is not 32 bit on this platform. Stopping.\n\n");
      endrun(0);
    }

  if(sizeof(float) != 4)
    {
      if(ThisTask == 0)
	printf("\nType `float' is not 32 bit on this platform. Stopping.\n\n");
      endrun(0);
    }

  if(sizeof(double) != 8)
    {
      if(ThisTask == 0)
	printf("\nType `double' is not 64 bit on this platform. Stopping.\n\n");
      endrun(0);
    }


  if(ThisTask == 0)		/* read parameter file on process 0 */
    {
      nt = 0;

      strcpy(tag[nt], "InitCondFile");
      addr[nt] = All.InitCondFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputDir");
      addr[nt] = All.OutputDir;
      id[nt++] = STRING;

      strcpy(tag[nt], "SnapshotFileBase");
      addr[nt] = All.SnapshotFileBase;
      id[nt++] = STRING;

      strcpy(tag[nt], "EnergyFile");
      addr[nt] = All.EnergyFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "CpuFile");
      addr[nt] = All.CpuFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "InfoFile");
      addr[nt] = All.InfoFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "TimingsFile");
      addr[nt] = All.TimingsFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "RestartFile");
      addr[nt] = All.RestartFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "ResubmitCommand");
      addr[nt] = All.ResubmitCommand;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputListFilename");
      addr[nt] = All.OutputListFilename;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputListOn");
      addr[nt] = &All.OutputListOn;
      id[nt++] = INT;

      strcpy(tag[nt], "Omega0");
      addr[nt] = &All.Omega0;
      id[nt++] = REAL;

      strcpy(tag[nt], "OmegaBaryon");
      addr[nt] = &All.OmegaBaryon;
      id[nt++] = REAL;

      strcpy(tag[nt], "OmegaLambda");
      addr[nt] = &All.OmegaLambda;
      id[nt++] = REAL;

      strcpy(tag[nt], "HubbleParam");
      addr[nt] = &All.HubbleParam;
      id[nt++] = REAL;

      strcpy(tag[nt], "BoxSize");
      addr[nt] = &All.BoxSize;
      id[nt++] = REAL;

      strcpy(tag[nt], "PeriodicBoundariesOn");
      addr[nt] = &All.PeriodicBoundariesOn;
      id[nt++] = INT;

      strcpy(tag[nt], "TimeOfFirstSnapshot");
      addr[nt] = &All.TimeOfFirstSnapshot;
      id[nt++] = REAL;

      strcpy(tag[nt], "CpuTimeBetRestartFile");
      addr[nt] = &All.CpuTimeBetRestartFile;
      id[nt++] = REAL;

      strcpy(tag[nt], "TimeBetStatistics");
      addr[nt] = &All.TimeBetStatistics;
      id[nt++] = REAL;

      strcpy(tag[nt], "TimeBegin");
      addr[nt] = &All.TimeBegin;
      id[nt++] = REAL;

      strcpy(tag[nt], "TimeMax");
      addr[nt] = &All.TimeMax;
      id[nt++] = REAL;

      strcpy(tag[nt], "TimeBetSnapshot");
      addr[nt] = &All.TimeBetSnapshot;
      id[nt++] = REAL;

      strcpy(tag[nt], "UnitVelocity_in_cm_per_s");
      addr[nt] = &All.UnitVelocity_in_cm_per_s;
      id[nt++] = REAL;

      strcpy(tag[nt], "UnitLength_in_cm");
      addr[nt] = &All.UnitLength_in_cm;
      id[nt++] = REAL;

      strcpy(tag[nt], "UnitMass_in_g");
      addr[nt] = &All.UnitMass_in_g;
      id[nt++] = REAL;

      strcpy(tag[nt], "TreeDomainUpdateFrequency");
      addr[nt] = &All.TreeDomainUpdateFrequency;
      id[nt++] = REAL;

      strcpy(tag[nt], "ErrTolIntAccuracy");
      addr[nt] = &All.ErrTolIntAccuracy;
      id[nt++] = REAL;

      strcpy(tag[nt], "ErrTolTheta");
      addr[nt] = &All.ErrTolTheta;
      id[nt++] = REAL;

      strcpy(tag[nt], "ErrTolForceAcc");
      addr[nt] = &All.ErrTolForceAcc;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinGasHsmlFractional");
      addr[nt] = &All.MinGasHsmlFractional;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxSizeTimestep");
      addr[nt] = &All.MaxSizeTimestep;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinSizeTimestep");
      addr[nt] = &All.MinSizeTimestep;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxRMSDisplacementFac");
      addr[nt] = &All.MaxRMSDisplacementFac;
      id[nt++] = REAL;

      strcpy(tag[nt], "ArtBulkViscConst");
      addr[nt] = &All.ArtBulkViscConst;
      id[nt++] = REAL;

      strcpy(tag[nt], "CourantFac");
      addr[nt] = &All.CourantFac;
      id[nt++] = REAL;

      strcpy(tag[nt], "DesNumNgb");
      addr[nt] = &All.DesNumNgb;
      id[nt++] = INT;

      strcpy(tag[nt], "MaxNumNgbDeviation");
      addr[nt] = &All.MaxNumNgbDeviation;
      id[nt++] = REAL;

      strcpy(tag[nt], "ComovingIntegrationOn");
      addr[nt] = &All.ComovingIntegrationOn;
      id[nt++] = INT;

      strcpy(tag[nt], "ICFormat");
      addr[nt] = &All.ICFormat;
      id[nt++] = INT;

      strcpy(tag[nt], "SnapFormat");
      addr[nt] = &All.SnapFormat;
      id[nt++] = INT;

      strcpy(tag[nt], "NumFilesPerSnapshot");
      addr[nt] = &All.NumFilesPerSnapshot;
      id[nt++] = INT;

      strcpy(tag[nt], "NumFilesWrittenInParallel");
      addr[nt] = &All.NumFilesWrittenInParallel;
      id[nt++] = INT;

      strcpy(tag[nt], "ResubmitOn");
      addr[nt] = &All.ResubmitOn;
      id[nt++] = INT;

      strcpy(tag[nt], "CoolingOn");
      addr[nt] = &All.CoolingOn;
      id[nt++] = INT;

      strcpy(tag[nt], "StarformationOn");
      addr[nt] = &All.StarformationOn;
      id[nt++] = INT;

      strcpy(tag[nt], "TypeOfTimestepCriterion");
      addr[nt] = &All.TypeOfTimestepCriterion;
      id[nt++] = INT;

      strcpy(tag[nt], "TypeOfOpeningCriterion");
      addr[nt] = &All.TypeOfOpeningCriterion;
      id[nt++] = INT;

      strcpy(tag[nt], "TimeLimitCPU");
      addr[nt] = &All.TimeLimitCPU;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningHalo");
      addr[nt] = &All.SofteningHalo;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningDisk");
      addr[nt] = &All.SofteningDisk;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningBulge");
      addr[nt] = &All.SofteningBulge;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningGas");
      addr[nt] = &All.SofteningGas;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningStars");
      addr[nt] = &All.SofteningStars;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningBndry");
      addr[nt] = &All.SofteningBndry;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningHaloMaxPhys");
      addr[nt] = &All.SofteningHaloMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningDiskMaxPhys");
      addr[nt] = &All.SofteningDiskMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningBulgeMaxPhys");
      addr[nt] = &All.SofteningBulgeMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningGasMaxPhys");
      addr[nt] = &All.SofteningGasMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningStarsMaxPhys");
      addr[nt] = &All.SofteningStarsMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "SofteningBndryMaxPhys");
      addr[nt] = &All.SofteningBndryMaxPhys;
      id[nt++] = REAL;

      strcpy(tag[nt], "BufferSize");
      addr[nt] = &All.BufferSize;
      id[nt++] = INT;

      strcpy(tag[nt], "PartAllocFactor");
      addr[nt] = &All.PartAllocFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "TreeAllocFactor");
      addr[nt] = &All.TreeAllocFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "GravityConstantInternal");
      addr[nt] = &All.GravityConstantInternal;
      id[nt++] = REAL;

      strcpy(tag[nt], "InitGasTemp");
      addr[nt] = &All.InitGasTemp;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinGasTemp");
      addr[nt] = &All.MinGasTemp;
      id[nt++] = REAL;

#ifdef OUTPUTLINEOFSIGHT
      strcpy(tag[nt], "TimeFirstLineOfSight");
      addr[nt] = &All.TimeFirstLineOfSight;
      id[nt++] = REAL;
#endif

#ifdef BG_SFR
      strcpy(tag[nt], "GammaEffective");
      addr[nt] = &All.GammaEffective;
      id[nt++] = REAL;
#endif

#if defined(BG_SFR) && defined(BG_STELLAR_EVOLUTION) && defined(BG_COOLING)
      strcpy(tag[nt], "MetalDependentDensThresh");
      addr[nt] = &All.MetalDependentDensThresh;
      id[nt++] = INT;

      strcpy(tag[nt], "EnergySNIa");
      addr[nt] = &All.EnergySNIa;
      id[nt++] = REAL;
#endif

#ifdef BG_WINDS
      strcpy(tag[nt], "WindSpeed");
      addr[nt] = &All.WindSpeed;
      id[nt++] = REAL;

      strcpy(tag[nt], "WindMassLoading");
      addr[nt] = &All.WindMassLoading;
      id[nt++] = REAL;
#endif

#ifdef BG_COOLING
      strcpy(tag[nt], "CoolTablePath");
      addr[nt] = &All.CoolTablePath;
      id[nt++] = STRING;
#endif

#ifdef BG_STELLAR_EVOLUTION
      strcpy(tag[nt], "AGBOn");
      addr[nt] = &All.AGBOn;
      id[nt++] = INT;

      strcpy(tag[nt], "SNIaOn");
      addr[nt] = &All.SNIaOn;
      id[nt++] = INT;

      strcpy(tag[nt], "SNIIOn");
      addr[nt] = &All.SNIIOn;
      id[nt++] = INT;

      strcpy(tag[nt], "YieldTablePath");
      addr[nt] = &All.YieldTablePath;
      id[nt++] = STRING;
#endif

#if defined(ADAPTIVE_GRAVSOFT_FORGAS) && !defined(ADAPTIVE_GRAVSOFT_FORGAS_HSML)
      strcpy(tag[nt], "ReferenceGasMass");
      addr[nt] = &All.ReferenceGasMass;
      id[nt++] = REAL;
#endif

#ifdef NAVIERSTOKES
      strcpy(tag[nt], "FractionSpitzerViscosity");
      addr[nt] = &All.FractionSpitzerViscosity;
      id[nt++] = REAL;
#endif

#ifdef NAVIERSTOKES_CONSTANT
      strcpy(tag[nt], "ShearViscosityTemperature");
      addr[nt] = &All.ShearViscosityTemperature;
      id[nt++] = REAL;
#endif

#ifdef NAVIERSTOKES_BULK
      strcpy(tag[nt], "NavierStokes_BulkViscosity");
      addr[nt] = &All.NavierStokes_BulkViscosity;
      id[nt++] = REAL;
#endif

#ifdef CHEMISTRY
      strcpy(tag[nt], "Epsilon");
      addr[nt] = &All.Epsilon;
      id[nt++] = REAL;
#endif


#ifdef CONDUCTION
      strcpy(tag[nt], "ConductionEfficiency");
      addr[nt] = &All.ConductionCoeff;
      id[nt++] = REAL;
#endif

#if defined(BUBBLES) || defined(MULTI_BUBBLES)
      strcpy(tag[nt], "BubbleDistance");
      addr[nt] = &All.BubbleDistance;
      id[nt++] = REAL;

      strcpy(tag[nt], "BubbleRadius");
      addr[nt] = &All.BubbleRadius;
      id[nt++] = REAL;

      strcpy(tag[nt], "BubbleTimeInterval");
      addr[nt] = &All.BubbleTimeInterval;
      id[nt++] = REAL;

      strcpy(tag[nt], "BubbleEnergy");
      addr[nt] = &All.BubbleEnergy;
      id[nt++] = REAL;

      strcpy(tag[nt], "FirstBubbleRedshift");
      addr[nt] = &All.FirstBubbleRedshift;
      id[nt++] = REAL;
#endif

#ifdef MULTI_BUBBLES
      strcpy(tag[nt], "MinFoFMassForNewSeed");
      addr[nt] = &All.MinFoFMassForNewSeed;
      id[nt++] = REAL;

      strcpy(tag[nt], "ClusterMass200");
      addr[nt] = &All.ClusterMass200;
      id[nt++] = REAL;

      strcpy(tag[nt], "massDMpart");
      addr[nt] = &All.massDMpart;
      id[nt++] = REAL;

#endif


#ifdef COSMIC_RAYS
      strcpy(tag[nt], "CR_SpectralIndex");
      addr[nt] = &All.CR_Alpha;
      id[nt++] = REAL;

      strcpy(tag[nt], "CR_SupernovaEfficiency");
      addr[nt] = &All.CR_SNEff;
      id[nt++] = REAL;

      strcpy(tag[nt], "CR_SupernovaSpectralIndex");
      addr[nt] = &All.CR_SNAlpha;
      id[nt++] = REAL;

#if defined(CR_DIFFUSION) || defined(CR_DIFFUSION_GREEN)
      strcpy(tag[nt], "CR_DiffusionCoeff");
      addr[nt] = &All.CR_DiffusionCoeff;
      id[nt++] = REAL;

      strcpy(tag[nt], "CR_DiffusionDensityScaling");
      addr[nt] = &All.CR_DiffusionDensScaling;
      id[nt++] = REAL;

      /* CR Diffusion scaling: reference density rho_0.
         If value is 0, then rho_crit @ z=0 is used. */
      strcpy(tag[nt], "CR_DiffusionReferenceDensity");
      addr[nt] = &All.CR_DiffusionDensZero;
      id[nt++] = REAL;

      strcpy(tag[nt], "CR_DiffusionTemperatureScaling");
      addr[nt] = &All.CR_DiffusionEntropyScaling;
      id[nt++] = REAL;

      strcpy(tag[nt], "CR_DiffusionReferenceTemperature");
      addr[nt] = &All.CR_DiffusionEntropyZero;
      id[nt++] = REAL;

      strcpy(tag[nt], "CR_DiffusionTimeScale");
      addr[nt] = &All.CR_DiffusionTimeScale;
      id[nt++] = REAL;
#endif /* CR_DIFFUSION */

#ifdef CR_SHOCK
      strcpy(tag[nt], "CR_ShockEfficiency");
      addr[nt] = &All.CR_ShockEfficiency;
      id[nt++] = REAL;
#if ( CR_SHOCK == 1 )		/* Constant Spectral Index method */
      strcpy(tag[nt], "CR_ShockSpectralIndex");
      addr[nt] = &All.CR_ShockAlpha;
      id[nt++] = REAL;
#else /* Mach-Number - Dependent method */
      strcpy(tag[nt], "CR_ShockCutoffFac");
      addr[nt] = &All.CR_ShockCutoff;
      id[nt++] = REAL;
#endif
#endif /* CR_SHOCK */

#ifdef FIX_QINJ
      strcpy(tag[nt], "Shock_Fix_Qinj");
      addr[nt] = &All.Shock_Fix_Qinj;
      id[nt++] = REAL;
#endif

#endif /* COSMIC_RAYS */


#ifdef MACHNUM
      strcpy(tag[nt], "Shock_LengthScale");
      addr[nt] = &All.Shock_Length;
      id[nt++] = REAL;

      strcpy(tag[nt], "Shock_DeltaDecayTimeMax");
      addr[nt] = &All.Shock_DeltaDecayTimeMax;
      id[nt++] = REAL;
#endif


#ifdef BLACK_HOLES
      strcpy(tag[nt], "TimeBetBlackHoleSearch");
      addr[nt] = &All.TimeBetBlackHoleSearch;
      id[nt++] = REAL;

      strcpy(tag[nt], "BlackHoleAccretionFactor");
      addr[nt] = &All.BlackHoleAccretionFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "BlackHoleFeedbackFactor");
      addr[nt] = &All.BlackHoleFeedbackFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "BlackHoleEddingtonFactor");
      addr[nt] = &All.BlackHoleEddingtonFactor;
      id[nt++] = REAL;


      strcpy(tag[nt], "SeedBlackHoleMass");
      addr[nt] = &All.SeedBlackHoleMass;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinFoFMassForNewSeed");
      addr[nt] = &All.MinFoFMassForNewSeed;
      id[nt++] = REAL;

      strcpy(tag[nt], "BlackHoleNgbFactor");
      addr[nt] = &All.BlackHoleNgbFactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "BlackHoleActiveTime");
      addr[nt] = &All.BlackHoleActiveTime;
      id[nt++] = REAL;

#ifdef FOF
      strcpy(tag[nt], "massDMpart");
      addr[nt] = &All.massDMpart;
      id[nt++] = REAL;
#endif

#endif

#ifdef LT_STELLAREVOLUTION
      /* the minimum number of neighbours used to spread
       * metals from Stars.
       */
      strcpy(tag[nt], "InfNeighNum");
      addr[nt] = &All.NeighInfNum;
      id[nt++] = INT;

      /* the desired number of neighbours used to spread
       * metals from Stars.
       */
      strcpy(tag[nt], "DesNumNgbSN");
      addr[nt] = &All.DesNumNgbSN;
      id[nt++] = INT;

      /* the allowed deviation on number of neighbours used
       * to spread metals from Stars.
       */
      strcpy(tag[nt], "SpreadNumNgbDev");
      addr[nt] = &All.SpreadNumNgbDev;
      id[nt++] = INT;

      /* note: this is used just to calculate All.MaxPartMet in 
       * read_file(); then, it should give a reasonable "average"
       * of the generations of each specified IMFs.
       */
      strcpy(tag[nt], "Generations");
      addr[nt] = &All.Generations;
      id[nt++] = INT;

      /* the baryon fraction that you expect to end in stars
       */
      strcpy(tag[nt], "SFfactor");
      addr[nt] = &All.SFfactor;
      id[nt++] = REAL;

      strcpy(tag[nt], "IMFsFileName");
      addr[nt] = All.IMFfilename;
      id[nt++] = STRING;

      /* the minimum mass for CC supernovae (usually 8 Msun) */
      strcpy(tag[nt], "SnII_InfMass");
      addr[nt] = &All.Mup;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinBinSystemMass");
      addr[nt] = &All.MBm;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxBinSystemMass");
      addr[nt] = &All.MBM;
      id[nt++] = REAL;

      strcpy(tag[nt], "BinarySystemFrac");
      addr[nt] = &All.BinFrac;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinStarMass_inBinSystem");
      addr[nt] = &All.MBms;
      id[nt++] = REAL;

#ifdef LT_AVOID_ENRICH_SFGAS
      /* this value is set in Msun/yr ;
       * Gas Particles having a sfr larger that this
       * value will not receive metals from stars.
       */
      strcpy(tag[nt], "SFTh_for_enrich");
      addr[nt] = &All.Enrich_SFGas_Th;
      id[nt++] = REAL;
#endif

#ifdef LT_MOD_EFFM
      strcpy(tag[nt], "ModSEffCrit");
      addr[nt] = &All.Mod_SEff;
      id[nt++] = INT;
#endif

#ifdef LT_HOT_EJECTA
      strcpy(tag[nt], "EgySpecEjecta");
      addr[nt] = &All.EgySpecEjecta;
      id[nt++] = REAL;
#endif

#ifdef LT_SNIa
      strcpy(tag[nt], "SnIaDataFileName");
      addr[nt] = All.SnIaDataFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "SnIaDataSetNum");
      addr[nt] = &All.Ia_Nset_ofYields;
      id[nt++] = INT;

      strcpy(tag[nt], "SnIaEnergy");
      addr[nt] = &All.SnIaEgy;
      id[nt++] = REAL;

      strcpy(tag[nt], "LongLiving_Step_Prec");
      addr[nt] = &All.LLv_Step_Prec;
      id[nt++] = REAL;

      /*
         strcpy(tag[nt], "SnIa_Remnant");
         addr[nt] = &All.SnIaRemn;
         id[nt++] = REAL;
       */
#endif

#ifdef LT_SNII
      strcpy(tag[nt], "SnIIdataFileName");
      addr[nt] = All.SnIIDataFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "SnIIDataSetNum");
      addr[nt] = &All.II_Nset_ofYields;
      id[nt++] = INT;

      strcpy(tag[nt], "SnIIEnergy");
      addr[nt] = &All.SnIIEgy;
      id[nt++] = REAL;

      /*
         strcpy(tag[nt], "ChemTimeStepII");
         addr[nt] = &All.ChemTimeStepII;
         id[nt++] = REAL;

         strcpy(tag[nt], "LongCTStepII");
         addr[nt] = &All.LongChemTimeStepII;
         id[nt++] = REAL;
       */

      strcpy(tag[nt], "metIRAThMass");
      addr[nt] = &All.metIRA_ThMass;
      id[nt++] = REAL;

      strcpy(tag[nt], "egyIRAThMass");
      addr[nt] = &All.egyIRA_ThMass;
      id[nt++] = REAL;

      strcpy(tag[nt], "SnII_Step_Prec");
      addr[nt] = &All.SnII_Step_Prec;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "Z_toset_SF_DensTh");
      addr[nt] = &All.referenceZ_toset_SF_DensTh;
      id[nt++] = REAL;

      strcpy(tag[nt], "Zdependent_SFTh");
      addr[nt] = &All.SFTh_Zdep;
      id[nt++] = INT;

#ifdef LT_AGB
      strcpy(tag[nt], "AGBdataFileName");
      addr[nt] = All.AGBDataFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "AGBDataSetNum");
      addr[nt] = &All.AGB_Nset_ofYields;
      id[nt++] = INT;
#endif

      strcpy(tag[nt], "MinChemTimeStep");
      addr[nt] = &All.MinChemTimeStep;
      id[nt++] = REAL;


      strcpy(tag[nt], "MinSpreadLength");
      addr[nt] = &All.MinChemSpreadL;
      id[nt++] = REAL;

      /*
         strcpy(tag[nt], "MaxSpreadLength");
         addr[nt] = &All.MaxChemSpreadL;
         id[nt++] = REAL;
       */
#endif

#ifdef TIME_DEP_ART_VISC
      strcpy(tag[nt], "ViscositySourceScaling");
      addr[nt] = &All.ViscSource0;
      id[nt++] = REAL;

      strcpy(tag[nt], "ViscosityDecayLength");
      addr[nt] = &All.DecayLength;
      id[nt++] = REAL;

      strcpy(tag[nt], "ViscosityAlphaMin");
      addr[nt] = &All.AlphaMin;
      id[nt++] = REAL;
#endif


#ifdef MOREPARAMS

#ifdef SFR_METALS
      strcpy(tag[nt], "FactorSFR");
      addr[nt] = &All.FactorSFR;
      id[nt++] = REAL;


      strcpy(tag[nt], "TlifeSNII");
      addr[nt] = &All.TlifeSNII;
      id[nt++] = REAL;

      strcpy(tag[nt], "MinTlifeSNI");
      addr[nt] = &All.MinTlifeSNI;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxTlifeSNI");
      addr[nt] = &All.MaxTlifeSNI;
      id[nt++] = REAL;

      strcpy(tag[nt], "RateSNI");
      addr[nt] = &All.RateSNI;
      id[nt++] = REAL;

      strcpy(tag[nt], "FactorSN_Phase");
      addr[nt] = &All.FactorSN_Phase;
      id[nt++] = REAL;

      strcpy(tag[nt], "Tcrit_Phase");
      addr[nt] = &All.Tcrit_Phase;
      id[nt++] = REAL;

      strcpy(tag[nt], "DensFrac_Phase");
      addr[nt] = &All.DensFrac_Phase;
      id[nt++] = REAL;

      strcpy(tag[nt], "FracEnergySN_Phase");
      addr[nt] = &All.FracEnergySN_Phase;
      id[nt++] = REAL;


#endif

      strcpy(tag[nt], "CritOverDensity");
      addr[nt] = &All.CritOverDensity;
      id[nt++] = REAL;

      strcpy(tag[nt], "CritPhysDensity");
      addr[nt] = &All.CritPhysDensity;
      id[nt++] = REAL;

#ifndef BG_SFR
#ifndef LT_STELLAREVOLUTION
      strcpy(tag[nt], "FactorSN");
      addr[nt] = &All.FactorSN;
      id[nt++] = REAL;
#endif

      strcpy(tag[nt], "FactorEVP");
      addr[nt] = &All.FactorEVP;
      id[nt++] = REAL;

      strcpy(tag[nt], "TempSupernova");
      addr[nt] = &All.TempSupernova;
      id[nt++] = REAL;

      strcpy(tag[nt], "TempClouds");
      addr[nt] = &All.TempClouds;
      id[nt++] = REAL;

      strcpy(tag[nt], "MaxSfrTimescale");
      addr[nt] = &All.MaxSfrTimescale;
      id[nt++] = REAL;

      strcpy(tag[nt], "WindEfficiency");
      addr[nt] = &All.WindEfficiency;
      id[nt++] = REAL;

      strcpy(tag[nt], "WindEnergyFraction");
      addr[nt] = &All.WindEnergyFraction;
      id[nt++] = REAL;

      strcpy(tag[nt], "WindFreeTravelLength");
      addr[nt] = &All.WindFreeTravelLength;
      id[nt++] = REAL;

      strcpy(tag[nt], "WindFreeTravelDensFac");
      addr[nt] = &All.WindFreeTravelDensFac;
      id[nt++] = REAL;
#endif

#ifdef SOFTEREQS
      strcpy(tag[nt], "FactorForSofterEQS");
      addr[nt] = &All.FactorForSofterEQS;
      id[nt++] = REAL;
#endif
#ifdef DARKENERGY
#ifndef TIMEDEPDE
      strcpy(tag[nt], "DarkEnergyParam");
      addr[nt] = &All.DarkEnergyParam;
      id[nt++] = REAL;
#endif
#endif

#ifdef RESCALEVINI
      strcpy(tag[nt], "VelIniScale");
      addr[nt] = &All.VelIniScale;
      id[nt++] = REAL;
#endif

#ifdef DARKENERGY
#ifdef TIMEDEPDE
      strcpy(tag[nt], "DarkEnergyFile");
      addr[nt] = All.DarkEnergyFile;
      id[nt++] = STRING;
#endif
#endif

#ifdef MAGNETIC_DISSIPATION
      strcpy(tag[nt], "ArtificialMagneticDissipationConstant");
      addr[nt] = &All.ArtMagDispConst;
      id[nt++] = REAL;

#ifdef TIME_DEP_MAGN_DISP
      strcpy(tag[nt], "ArtificialMagneticDissipationMin");
      addr[nt] = &All.ArtMagDispMin;
      id[nt++] = REAL;

      strcpy(tag[nt], "ArtificialMagneticDissipationSource");
      addr[nt] = &All.ArtMagDispSource;
      id[nt++] = REAL;

      strcpy(tag[nt], "ArtificialMagneticDissipationDecaytime");
      addr[nt] = &All.ArtMagDispTime;
      id[nt++] = REAL;
#endif
#endif

#ifdef DIVBCLEANING_DEDNER
      strcpy(tag[nt], "DivBcleaningParabolicSigma");
      addr[nt] = &All.DivBcleanParabolicSigma;
      id[nt++] = REAL;

      strcpy(tag[nt], "DivBcleaningHyperbolicSigma");
      addr[nt] = &All.DivBcleanHyperbolicSigma;
      id[nt++] = REAL;
#endif

#ifdef MAGNETIC
#ifdef BINISET
      strcpy(tag[nt], "BiniX");
      addr[nt] = &All.BiniX;
      id[nt++] = REAL;

      strcpy(tag[nt], "BiniY");
      addr[nt] = &All.BiniY;
      id[nt++] = REAL;

      strcpy(tag[nt], "BiniZ");
      addr[nt] = &All.BiniZ;
      id[nt++] = REAL;
#endif

#ifdef BSMOOTH
      strcpy(tag[nt], "BSmoothInt");
      addr[nt] = &All.BSmoothInt;
      id[nt++] = INT;

      strcpy(tag[nt], "BSmoothFrac");
      addr[nt] = &All.BSmoothFrac;
      id[nt++] = REAL;
#endif
#endif

#endif

#ifdef VOLUME_CORRECTION
      strcpy(tag[nt], "VolCorrectFile");
      addr[nt] = All.VolCorrectFile;
      id[nt++] = STRING;
#endif

      if((fd = fopen(fname, "r")))
	{
	  sprintf(buf, "%s%s", fname, "-usedvalues");
	  if(!(fdout = fopen(buf, "w")))
	    {
	      printf("error opening file '%s' \n", buf);
	      errorFlag = 1;
	    }
	  else
	    {
	      while(!feof(fd))
		{
		  *buf = 0;
		  fgets(buf, 200, fd);
		  if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
		    continue;

		  if(buf1[0] == '%')
		    continue;

		  for(i = 0, j = -1; i < nt; i++)
		    if(strcmp(buf1, tag[i]) == 0)
		      {
			j = i;
			tag[i][0] = 0;
			break;
		      }

		  if(j >= 0)
		    {
		      switch (id[j])
			{
			case REAL:
			  *((double *) addr[j]) = atof(buf2);
			  fprintf(fdout, "%-35s%g\n", buf1, *((double *) addr[j]));
			  break;
			case STRING:
			  strcpy((char *) addr[j], buf2);
			  fprintf(fdout, "%-35s%s\n", buf1, buf2);
			  break;
			case INT:
			  *((int *) addr[j]) = atoi(buf2);
			  fprintf(fdout, "%-35s%d\n", buf1, *((int *) addr[j]));
			  break;
			}
		    }
		  else
		    {
#ifdef ALLOWEXTRAPARAMS
		      fprintf(stdout, "WARNING from file %s:   Tag '%s' ignored !\n", fname, buf1);
#else
		      fprintf(stdout, "Error in file %s:   Tag '%s' not allowed or multiple defined.\n",
			      fname, buf1);
		      errorFlag = 1;
#endif
		    }
		}
	      fclose(fd);
	      fclose(fdout);

	      i = strlen(All.OutputDir);
	      if(i > 0)
		if(All.OutputDir[i - 1] != '/')
		  strcat(All.OutputDir, "/");

	      sprintf(buf1, "%s%s", fname, "-usedvalues");
	      sprintf(buf2, "%s%s", All.OutputDir, "parameters-usedvalues");
	      sprintf(buf3, "cp %s %s", buf1, buf2);
#ifndef NOCALLSOFSYSTEM
	      system(buf3);
#endif
	    }
	}
      else
	{
	  printf("Parameter file %s not found.\n", fname);
	  errorFlag = 1;
	}


      for(i = 0; i < nt; i++)
	{
#ifdef LT_STELLAREVOLUTION
	  if(*tag[i] &&
	     (strcmp(tag[i], "metIRAThMass") != 0 && strcmp(tag[i], "egyIRAThMass") != 0) &&
	     strcmp(tag[i], "WindEfficiency") != 0 && strcmp(tag[i], "WindEnergyFraction") != 0)
#else
	  if(*tag[i])
#endif
	    {
	      printf("Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
	      errorFlag = 1;
	    }
	}

      if(All.OutputListOn && errorFlag == 0)
	errorFlag += read_outputlist(All.OutputListFilename);
      else
	All.OutputListLength = 0;
    }

  MPI_Bcast(&errorFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if(errorFlag)
    {
      MPI_Finalize();
      exit(0);
    }

  /* now communicate the relevant parameters to the other processes */
  MPI_Bcast(&All, sizeof(struct global_data_all_processes), MPI_BYTE, 0, MPI_COMM_WORLD);



  for(pnum = 0; All.NumFilesWrittenInParallel > (1 << pnum); pnum++);

  if(All.NumFilesWrittenInParallel != (1 << pnum))
    {
      if(ThisTask == 0)
	printf("NumFilesWrittenInParallel MUST be a power of 2\n");
      endrun(0);
    }

  if(All.NumFilesWrittenInParallel > NTask)
    {
      if(ThisTask == 0)
	printf("NumFilesWrittenInParallel MUST be smaller than number of processors\n");
      endrun(0);
    }

#ifdef PERIODIC
  if(All.PeriodicBoundariesOn == 0)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with periodic boundary conditions switched on.\n");
	  printf("You must set `PeriodicBoundariesOn=1', or recompile the code.\n");
	}
      endrun(0);
    }
#else
  if(All.PeriodicBoundariesOn == 1)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with periodic boundary conditions switched off.\n");
	  printf("You must set `PeriodicBoundariesOn=0', or recompile the code.\n");
	}
      endrun(0);
    }
#endif

#ifdef COOLING
  if(All.CoolingOn == 0)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with cooling switched on.\n");
	  printf("You must set `CoolingOn=1', or recompile the code.\n");
	}
      endrun(0);
    }
#else
  if(All.CoolingOn == 1)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with cooling switched off.\n");
	  printf("You must set `CoolingOn=0', or recompile the code.\n");
	}
      endrun(0);
    }
#endif

  if(All.TypeOfTimestepCriterion >= 3)
    {
      if(ThisTask == 0)
	{
	  printf("The specified timestep criterion\n");
	  printf("is not valid\n");
	}
      endrun(0);
    }

#if defined(LONG_X) ||  defined(LONG_Y) || defined(LONG_Z)
#ifndef NOGRAVITY
  if(ThisTask == 0)
    {
      printf("Code was compiled with LONG_X/Y/Z, but not with NOGRAVITY.\n");
      printf("Stretched periodic boxes are not implemented for gravity yet.\n");
    }
  endrun(0);
#endif
#endif

#ifdef SFR

#ifndef MOREPARAMS
  if(ThisTask == 0)
    {
      printf("Code was compiled with SFR, but not with MOREPARAMS.\n");
      printf("This is not allowed.\n");
    }
  endrun(0);
#endif

  if(All.StarformationOn == 0)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with star formation switched on.\n");
	  printf("You must set `StarformationOn=1', or recompile the code.\n");
	}
      endrun(0);
    }
  if(All.CoolingOn == 0)
    {
      if(ThisTask == 0)
	{
	  printf("You try to use the code with star formation enabled,\n");
	  printf("but you did not switch on cooling.\nThis mode is not supported.\n");
	}
      endrun(0);
    }
#else
  if(All.StarformationOn == 1)
    {
      if(ThisTask == 0)
	{
	  printf("Code was compiled with star formation switched off.\n");
	  printf("You must set `StarformationOn=0', or recompile the code.\n");
	}
      endrun(0);
    }
#endif


#ifdef METALS
#ifndef SFR
  if(ThisTask == 0)
    {
      printf("Code was compiled with METALS, but not with SFR.\n");
      printf("This is not allowed.\n");
    }
  endrun(0);
#endif
#ifdef SFR_METALS
#ifdef SFR_SNI
#ifndef SFR_ENRICH
  if(ThisTask == 0)
    {
      printf("Code was compiled with SNI, but not with ENRICH.\n");
      printf("This is not allowed.\n");
    }
  endrun(0);
#endif
#endif
#ifdef SFR_SNII
#ifndef SFR_ENRICH
  if(ThisTask == 0)
    {
      printf("Code was compiled with SNII, but not with ENRICH.\n");
      printf("This is not allowed.\n");
    }
  endrun(0);
#endif
#endif
#endif
#endif

#ifdef LT_STELLAREVOLUTION	/* ======= */

#ifndef LT_SNII
#error >>>LT_STELLAREVOLUTION without LT_SNII is not allowed
#endif
#if !defined(LT_SNII) && !defined(LT_SNIa) && !defined(LT_Neabulae)
#error >>> STELLAREVOLUTION enabled without LT_SNII, LT_SNIa, LT_Nebulae does not make sense
#endif
#if defined(LT_METAL_EFF_MODEL) && !defined(LT_METAL_COOLING)
#error >>> LT_METAL_EFF_MODEL enabled without LT_METAL_COOLING does not make sense
#endif
#else /* === */

#if defined(LT_SNII) || defined(LT_SNIa) || defined(LT_Neabulae)
#error >>>LT_SNII, LT_SNIa, LT_Nebulae enabled without STELLAREVOLUTION does not make sense
#endif
#ifdef LT_SEv_INFO
#error >>> LT_SEv_INFO enabled without LT_STELLAREVOLUTION does not make sense
#endif
#endif /* ======= */

#ifndef MOREPARAMS
#ifdef DARKENERGY
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with DARKENERGY, but not with MOREPARAMS.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif

#ifdef TIMEDEPDE
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with TIMEDEPDE, but not with MOREPARAMS.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif
#endif

#ifdef TIMEDEPDE
#ifndef DARKENERGY
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with TIMEDEPDE, but not with DARKENERGY.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif
#endif

#ifndef MAGNETIC
#ifdef TRACEDIVB
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with TRACEDIVB, but not with MAGNETIC.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif

#ifdef DBOUTPUT
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with DBOUTPUT, but not with MAGNETIC.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif

#ifdef MAGFORCE
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with MAGFORCE, but not with MAGNETIC.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif

#ifdef BSMOOTH
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with BSMOOTH, but not with MAGNETIC.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif

#ifdef BFROMROTA
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with BFROMROTA, but not with MAGNETIC.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif

#ifdef MU0_UNITY
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with MU0_UNITY, but not with MAGNETIC.\n");
      fprintf(stdout, "This makes no sense.\n");
    }
  endrun(0);
#endif

#ifdef MAGNETIC_DISSIPATION
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with MAGNETIC_DISSIPATION, but not with MAGNETIC.\n");
      fprintf(stdout, "This makes no sense.\n");
    }
  endrun(0);
#endif

#ifdef TIME_DEP_MAGN_DISP
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with TIME_DEP_MAGN_DISP, but not with MAGNETIC.\n");
      fprintf(stdout, "This makes no sense.\n");
    }
  endrun(0);
#endif

#ifdef DIVBCLEANING_DEDNER
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with DIVBCLEANING_DEDNER, but not with MAGNETIC.\n");
      fprintf(stdout, "This makes no sense.\n");
    }
  endrun(0);
#endif

#endif

#if defined(NOWINDTIMESTEPPING) && defined(MAGNETIC)
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with NOWINDTIMESTEPPING and with MAGNETIC.\n");
      fprintf(stdout, "This is not allowed, as it leads to inconsitent MHD for wind particles.\n");
    }
  endrun(0);
#endif

#ifndef MAGFORCE
#ifdef ARTBPRES
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with ARTBPRES, but not with MAGFORCE.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif

#ifdef DIVBFORCE
  if(ThisTask == 0)
    {
      fprintf(stdout, "Code was compiled with DIVBFORCE, but not with MAGFORCE.\n");
      fprintf(stdout, "This is not allowed.\n");
    }
  endrun(0);
#endif
#endif


#undef REAL
#undef STRING
#undef INT
#undef MAXTAGS


#ifdef COSMIC_RAYS
  if(ThisTask == 0)
    {
      printf("CR SN Efficiency: %g\n", All.CR_SNEff);
    }
#endif
}


/*! this function reads a table with a list of desired output times. The table
 *  does not have to be ordered in any way, but may not contain more than
 *  MAXLEN_OUTPUTLIST entries.
 */
int read_outputlist(char *fname)
{
  FILE *fd;

  if(!(fd = fopen(fname, "r")))
    {
      printf("can't read output list in file '%s'\n", fname);
      return 1;
    }

  All.OutputListLength = 0;
  do
    {
      if(fscanf(fd, " %lg ", &All.OutputListTimes[All.OutputListLength]) == 1)
	All.OutputListLength++;
      else
	break;
    }
  while(All.OutputListLength < MAXLEN_OUTPUTLIST);

  fclose(fd);

  printf("\nfound %d times in output-list.\n", All.OutputListLength);

  return 0;
}


/*! If a restart from restart-files is carried out where the TimeMax variable
 * is increased, then the integer timeline needs to be adjusted. The approach
 * taken here is to reduce the resolution of the integer timeline by factors
 * of 2 until the new final time can be reached within TIMEBASE.
 */
void readjust_timebase(double TimeMax_old, double TimeMax_new)
{
  int i;
  long long ti_end;

  if(sizeof(long long) != 8)
    {
      if(ThisTask == 0)
	printf("\nType 'long long' is not 64 bit on this platform\n\n");
      endrun(555);
    }

  if(ThisTask == 0)
    {
      printf("\nAll.TimeMax has been changed in the parameterfile\n");
      printf("Need to adjust integer timeline\n\n\n");
    }

  if(TimeMax_new < TimeMax_old)
    {
      if(ThisTask == 0)
	printf("\nIt is not allowed to reduce All.TimeMax\n\n");
      endrun(556);
    }

  if(All.ComovingIntegrationOn)
    ti_end = (long long) (log(TimeMax_new / All.TimeBegin) / All.Timebase_interval);
  else
    ti_end = (long long) ((TimeMax_new - All.TimeBegin) / All.Timebase_interval);

  while(ti_end > TIMEBASE)
    {
      All.Timebase_interval *= 2.0;

      ti_end /= 2;
      All.Ti_Current /= 2;

#ifdef PMGRID
      All.PM_Ti_begstep /= 2;
      All.PM_Ti_endstep /= 2;
#endif

      for(i = 0; i < NumPart; i++)
	{
	  P[i].Ti_begstep /= 2;
	  P[i].Ti_endstep /= 2;
	}

      All.Ti_nextlineofsight /= 2;
    }

  All.TimeMax = TimeMax_new;
}
