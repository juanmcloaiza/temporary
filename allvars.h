/*! \file allvars.h
 *  \brief declares global variables.
 *
 *  This file declares all global variables. Further variables should be added here, and declared as
 *  'extern'. The actual existence of these variables is provided by the file 'allvars.c'. To produce
 *  'allvars.c' from 'allvars.h', do the following:
 *
 *     - Erase all #define statements
 *     - add #include "allvars.h" 
 *     - delete all keywords 'extern'
 *     - delete all struct definitions enclosed in {...}, e.g.
 *        "extern struct global_data_all_processes {....} All;"
 *        becomes "struct global_data_all_processes All;"
 */

#ifndef ALLVARS_H
#define ALLVARS_H

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include "tags.h"
#ifdef CHEMISTRY
#include "chemistry.h"
#endif
#include "assert.h"
#ifdef LT_STELLAREVOLUTION
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#endif

#define  GADGETVERSION   "2.0"   /*!< code version string */

#define  GENERATIONS     2       /*!< Number of star particles that may be created per gas particle
                                  */
#ifdef LT_STELLAREVOLUTION
#undef GENERATIONS
#endif

#define  TIMEBASE        (1<<28) /*!< The simulated timespan is mapped onto the integer interval [0,TIMESPAN],
                                  *   where TIMESPAN needs to be a power of 2. Note that (1<<28) corresponds
                                  *   to 2^29
                                  */
#define MAXTOPNODES    100000


typedef  long long  peanokey;

#define  BITS_PER_DIMENSION 18	/* for Peano-Hilbert order. Note: Maximum is 10 to fit in 32-bit integer ! */
#define  PEANOCELLS (((peanokey)1)<<(3*BITS_PER_DIMENSION))


#ifdef   ISOTHERM_EQS
#define  GAMMA         (1.0)     /*!< index for isothermal gas */
#else
#define  GAMMA         (5.0/3.0)   /*!< adiabatic index of simulated gas */
#endif

#ifdef   MYSWITCH
#define my_a 1.0 //1kpc
#define my_Ma 1.0 //1e10 Msun
#define my_Mbh 1.0e-2 //Units: 1e10 Msun
#define my_Mcore 2.0e-2
#define my_rcore 2.0e-2
#define my_racc 1.0e-3
#endif

#define  GAMMA_MINUS1  (GAMMA-1)

#define  HYDROGEN_MASSFRAC 0.76   /*!< mass fraction of hydrogen, relevant only for radiative cooling */

#define  METAL_YIELD       0.02   /*!< effective metal yield for star formation */

#define  MAX_REAL_NUMBER  1e37
#define  MIN_REAL_NUMBER  1e-37

#define  RNDTABLE 8192

/* ... often used physical constants (cgs units) */

#define  GRAVITY     6.672e-8
#define  SOLAR_MASS  1.989e33
#define  SOLAR_LUM   3.826e33
#define  RAD_CONST   7.565e-15
#define  AVOGADRO    6.0222e23
#define  BOLTZMANN   1.3806e-16
#define  GAS_CONST   8.31425e7
#define  C           2.9979e10
#define  PLANCK      6.6262e-27
#define  CM_PER_MPC  3.085678e24
#define  PROTONMASS  1.6726e-24
#define  ELECTRONMASS 9.10953e-28
#define  THOMPSON     6.65245e-25
#define  ELECTRONCHARGE  4.8032e-10
#define  HUBBLE          3.2407789e-18	/* in h/sec */
#define  LYMAN_ALPHA      1215.6e-8      /* 1215.6 Angstroem */
#define  LYMAN_ALPHA_HeII  303.8e-8       /* 303.8 Angstroem */
#define  OSCILLATOR_STRENGTH       0.41615
#define  OSCILLATOR_STRENGTH_HeII  0.41615

#ifdef NAVIERSTOKES
#define  LOG_LAMBDA      37.8	/* logarithmic Coulomb factor */
#endif

#ifdef CHEMISTRY
#define  T_CMB0      2.728	/* present-day CMB temperature */
#endif

#define  SEC_PER_MEGAYEAR   3.155e13
#define  SEC_PER_YEAR       3.155e7

#ifndef ASMTH
/*! ASMTH gives the scale of the short-range/long-range force split in units of FFT-mesh cells */
#define ASMTH 1.25		
#endif
#ifndef RCUT
/*! RCUT gives the maximum distance (in units of the scale used for the force split) out to which short-range
 * forces are evaluated in the short-range tree walk.
 */
#define RCUT  4.5		
#endif

#ifndef HPM_SMTH
#define HPM_SMTH ASMTH
#endif

#define COND_TIMESTEP_PARAMETER 0.25
#define VISC_TIMESTEP_PARAMETER 0.25

#define MAX_NGB  20000		/*!< defines maximum length of neighbour list */

#define MAXLEN_OUTPUTLIST 350	/*!< maxmimum number of entries in output list */

#define DRIFT_TABLE_LENGTH  1000  /*!< length of the lookup table used to hold the drift and kick factors */ 


#define MAXITER 150

#define LINKLENGTH 0.2
#define GROUP_MIN_LEN 32

#define MINRESTFAC 0.05



#define CPU_ALL        0
#define CPU_TREEWALK   1
#define CPU_TREEBUILD  2
#define CPU_TREEUPDATE 3
#define CPU_TREEHMAXUPDATE 4
#define CPU_TREECOMM   5
#define CPU_DOMAIN     6
#define CPU_DENSITY    7
#define CPU_HYDRA      8
#define CPU_DRIFT      9
#define CPU_MOVE       10
#define CPU_TIMELINE   11
#define CPU_POTENTIAL  12
#define CPU_MESH       13
#define CPU_PEANO      14
#define CPU_COOLINGSFR 15        
#define CPU_SNAPSHOT   16
#define CPU_FOF        17
#define CPU_GRAVCOMM   18
#define CPU_DENSCOMM   19
#define CPU_HYDCOMM    20
#define CPU_PARTS      21  /* this gives the number of parts above (must be last) */

#define CPU_STRING_LEN 140

extern double CPU_Step[CPU_PARTS];
extern char CPU_Symbol[CPU_PARTS];
extern char CPU_SymbolImbalance[CPU_PARTS];
extern char CPU_String[CPU_STRING_LEN+1];


#ifdef DOUBLEPRECISION
#define FLOAT double
#else
#define FLOAT float
#endif

#define MyFloat FLOAT

#ifdef FLTROUNDOFFREDUCTION

#define FLT(x) ((FLOAT)(x))

#ifdef SOFTDOUBLEDOUBLE  /* this requires a C++ compilation */
#include "dd.h"
#define DOUBLE dd
#else
#define DOUBLE long double
#endif

#else

#define FLT(x) (x)
#define DOUBLE FLOAT

#endif  /* end FLTROUNDOFFREDUCTION */



#ifndef  TWODIMS
#define  NUMDIMS 3                                      /*!< For 3D-normalized kernel */
#define  KERNEL_COEFF_1  2.546479089470                 /*!< Coefficients for SPH spline kernel and its derivative */ 
#define  KERNEL_COEFF_2  15.278874536822
#define  KERNEL_COEFF_3  45.836623610466
#define  KERNEL_COEFF_4  30.557749073644
#define  KERNEL_COEFF_5  5.092958178941
#define  KERNEL_COEFF_6  (-15.278874536822)
#define  NORM_COEFF      4.188790204786                 /*!< Coefficient for kernel normalization. Note:  4.0/3 * PI = 4.188790204786 */ 
#define  NORM_COEFF_2    1.0444543
#define  NORM_COEFF_4    0.015542474
#define  NORM_COEFF_5    0.071619720
#else
#define  NUMDIMS 2                                      /*!< For 2D-normalized kernel */
#define  KERNEL_COEFF_1  (5.0/7*2.546479089470)         /*!< Coefficients for SPH spline kernel and its derivative */ 
#define  KERNEL_COEFF_2  (5.0/7*15.278874536822)
#define  KERNEL_COEFF_3  (5.0/7*45.836623610466)
#define  KERNEL_COEFF_4  (5.0/7*30.557749073644)
#define  KERNEL_COEFF_5  (5.0/7*5.092958178941)
#define  KERNEL_COEFF_6  (5.0/7*(-15.278874536822))
#define  NORM_COEFF      M_PI                           /*!< Coefficient for kernel normalization. */
#endif






#ifdef BG_SFR
/* KENNICUTT_COEFF = 2.343e-16 is derived for CGS units assuming
   KENNICUTT_EXP = 1.4. If one would like to change the value of
   KENNICUTT_EXP, KENNICUTT_COEFF must be rescaled as below */ 
#define KENNICUTT_EXP   1.4
#define KENNICUTT_COEFF (2.343e-16*pow(SOLAR_MASS/pow(CM_PER_MPC/1.0e6, 2), 1.4 - KENNICUTT_EXP))
#define KENNICUTT_NORM  1.0
#define TEMP_THRESH     1.0e4
#endif

#if defined(BG_SFR) || defined(BG_COOLING) || defined(BG_STELLAR_EVOLUTION)
#define EL_NAME_LENGTH 20
#define BG_NELEMENTS    8
#endif


#ifdef LT_STELLAREVOLUTION
/* ********************************
 * 
 *                STELLAR EVOLUTION
 *
 ********************************** */
#ifndef LT_NMet
>>> ERROR: you MUST define LT_NMet at compilation time
#endif

#define LT_NMetP (LT_NMet-1)

#define REMN 0

#define EVOLVE_SN 0

#define NO_METAL -5.0

extern FILE *FdSnInit, *FdWarn;

extern gsl_error_handler_t *old_error_handler;
extern int my_gslstatus;

#ifndef LT_NIMFs
#define LT_NIMFs 1
#endif

/* ------------------- *
 *   gsl integration   *
 * ------------------- */

#define gsl_ws_dim 10000
extern gsl_function F;
extern int gsl_status;
extern gsl_integration_workspace *w;

extern gsl_interp_accel *accel;
extern gsl_spline *spline;

/* ----------------- *
 *   metal species   *
 * ----------------- */

extern char *MetNames[LT_NMet];
extern double MetSolarValues[LT_NMet];
extern int Hyd, Hel, Iron, Oxygen, FillEl;

extern float xclouds;
extern int search_for_metalspread;

#ifdef LT_TRACK_CONTRIBUTES
  /*
   * LT_N_INT_PACKING is calculated once you declare at
   * compile time both LT_NBits and LT_NIMFs, that are
   * defined here below.
   * LT_N_INT_PACKING_1IMF is also calculated, that is
   * defined as LT_N_INT_PACKING / Number_of_IMFs_Used
   *
   *
   * LT_N_INT_PACKING can be calculated as follows:
   *
   * contributes run from 0 to 1. say that you want to track them
   * with a precision o nbits: then, you can map those number
   * using (number / 2^-nbits) between 0 and 2^nbits.
   * if you multiply (number / 2^-nbits) by 10^exp_max, you will
   * increase the number if significant digits in the integer
   * representation, then reducing the error when reconstructing
   * the original fraction.
   *
   * if you want to use nbits for the significant part and 2 bits
   * for the exponent, the number of bytes that you need is:
   *
   * ((nbits + 2) * N_species * N_elements * N_imf) / 8
   *
   * because you can infer the last species from
   *
   * using N_species = 3 (snIa, CC sn, agb stars) this becomes
   *
   *   >>>> LT_N_INT_PACKING = ((nbits + 2) * 3 * N_elements * N_imf) / 8  <<<<
   *
   * simply using float will put the memory requirement larger
   * by a factor 32 / (nbits+2);
   */
#ifndef LT_Nbits
#define LT_Nbits 20
#endif
#ifndef LT_NIMFs
#define LT_NIMFs 2
#endif
#define LT_power10_Nbits 2
typedef struct
{
  unsigned     II_el0_imf0 : LT_Nbits;
  unsigned  IIexp_el0_imf0 : LT_power10_Nbits;
  unsigned     II_el1_imf0 : LT_Nbits;
  unsigned  IIexp_el1_imf0 : LT_power10_Nbits;
  unsigned     II_el2_imf0 : LT_Nbits;
  unsigned  IIexp_el2_imf0 : LT_power10_Nbits;
  unsigned     II_el3_imf0 : LT_Nbits;
  unsigned  IIexp_el3_imf0 : LT_power10_Nbits;
  unsigned     II_el4_imf0 : LT_Nbits;
  unsigned  IIexp_el4_imf0 : LT_power10_Nbits;
  unsigned     II_el5_imf0 : LT_Nbits;
  unsigned  IIexp_el5_imf0 : LT_power10_Nbits;
  unsigned     II_el6_imf0 : LT_Nbits;
  unsigned  IIexp_el6_imf0 : LT_power10_Nbits;
  unsigned     II_el7_imf0 : LT_Nbits;
  unsigned  IIexp_el7_imf0 : LT_power10_Nbits;
  unsigned     II_el8_imf0 : LT_Nbits;
  unsigned  IIexp_el8_imf0 : LT_power10_Nbits;

  unsigned     Ia_el0_imf0 : LT_Nbits;
  unsigned  Iaexp_el0_imf0 : LT_power10_Nbits;
  unsigned     Ia_el1_imf0 : LT_Nbits;
  unsigned  Iaexp_el1_imf0 : LT_power10_Nbits;
  unsigned     Ia_el2_imf0 : LT_Nbits;
  unsigned  Iaexp_el2_imf0 : LT_power10_Nbits;
  unsigned     Ia_el3_imf0 : LT_Nbits;
  unsigned  Iaexp_el3_imf0 : LT_power10_Nbits;
  unsigned     Ia_el4_imf0 : LT_Nbits;
  unsigned  Iaexp_el4_imf0 : LT_power10_Nbits;
  unsigned     Ia_el5_imf0 : LT_Nbits;
  unsigned  Iaexp_el5_imf0 : LT_power10_Nbits;
  unsigned     Ia_el6_imf0 : LT_Nbits;
  unsigned  Iaexp_el6_imf0 : LT_power10_Nbits;
  unsigned     Ia_el7_imf0 : LT_Nbits;
  unsigned  Iaexp_el7_imf0 : LT_power10_Nbits;
  unsigned     Ia_el8_imf0 : LT_Nbits;
  unsigned  Iaexp_el8_imf0 : LT_power10_Nbits;

  unsigned    AGB_el0_imf0 : LT_Nbits;
  unsigned AGBexp_el0_imf0 : LT_power10_Nbits;
  unsigned    AGB_el1_imf0 : LT_Nbits;
  unsigned AGBexp_el1_imf0 : LT_power10_Nbits;
  unsigned    AGB_el2_imf0 : LT_Nbits;
  unsigned AGBexp_el2_imf0 : LT_power10_Nbits;
  unsigned    AGB_el3_imf0 : LT_Nbits;
  unsigned AGBexp_el3_imf0 : LT_power10_Nbits;
  unsigned    AGB_el4_imf0 : LT_Nbits;
  unsigned AGBexp_el4_imf0 : LT_power10_Nbits;
  unsigned    AGB_el5_imf0 : LT_Nbits;
  unsigned AGBexp_el5_imf0 : LT_power10_Nbits;
  unsigned    AGB_el6_imf0 : LT_Nbits;
  unsigned AGBexp_el6_imf0 : LT_power10_Nbits;
  unsigned    AGB_el7_imf0 : LT_Nbits;
  unsigned AGBexp_el7_imf0 : LT_power10_Nbits;
  unsigned    AGB_el8_imf0 : LT_Nbits;
  unsigned AGBexp_el8_imf0 : LT_power10_Nbits;

  /* add here below more blocks if more ifms are
   * being used; just change imfX appropriately.
   */

  unsigned     II_el0_imf1 : LT_Nbits;
  unsigned  IIexp_el0_imf1 : LT_power10_Nbits;
  unsigned     II_el1_imf1 : LT_Nbits;
  unsigned  IIexp_el1_imf1 : LT_power10_Nbits;
  unsigned     II_el2_imf1 : LT_Nbits;
  unsigned  IIexp_el2_imf1 : LT_power10_Nbits;
  unsigned     II_el3_imf1 : LT_Nbits;
  unsigned  IIexp_el3_imf1 : LT_power10_Nbits;
  unsigned     II_el4_imf1 : LT_Nbits;
  unsigned  IIexp_el4_imf1 : LT_power10_Nbits;
  unsigned     II_el5_imf1 : LT_Nbits;
  unsigned  IIexp_el5_imf1 : LT_power10_Nbits;
  unsigned     II_el6_imf1 : LT_Nbits;
  unsigned  IIexp_el6_imf1 : LT_power10_Nbits;
  unsigned     II_el7_imf1 : LT_Nbits;
  unsigned  IIexp_el7_imf1 : LT_power10_Nbits;
  unsigned     II_el8_imf1 : LT_Nbits;
  unsigned  IIexp_el8_imf1 : LT_power10_Nbits;
  
} Contrib;

extern unsigned int   Packing_Factor;
extern unsigned int   *Power10_Factors, Max_Power10;
extern unsigned int   TrackMask;
extern float Max_Packed_Int, UnPacking_Factor;
#endif

/* -------------- *
 *   SF related   *
 * -------------- */

extern double *PhysDensTh, *FEVP;
extern int sfrrate_filenum;

/* -------------- *
 *   Sn related   *
 * -------------- */

extern struct SDtype
{
  int Zbin, Zdim, Mdim, Yset;
  double Zstar;
  double *ZArray, *MArray;
  double *Y;
  int ExtrDir;
  int ExtrZbin, ExtrZdim, ExtrMdim, ExtrYset;
  double *ExtrZArray, *ExtrMArray;
  double *ExtrY;
} SD;

#ifdef LT_SNIa
/*
 * SnIaData defines the yields for SnIa.
 * we define it as a pointer so that in the future we can easily extend the
 * number of parameters on which they depend (e.g. mass of the system, time).
 * currently they are fixed
 */
extern double ***SnIaYields;
/*
 * it is common that yields tables use the same mass array and Z array for
 * all elements. in this case, we simply use IaZbins, IaMbins to store these 
 * common data. otherwise (different mass/Z array for each element), you can
 * just use the already existent Yield structure, just allocating as twice as
 * the amount of memory needed for the yields' value and using half of the
 * memory to store also the array's value. You can also use SnIaY_dim to store
 * the dimensions.
 */
extern double **IaZbins, **IaMbins;
extern int *IaZbins_dim, *IaMbins_dim, *SnIaY_dim[LT_NMet];
extern int *SnIaY_Give_Produced;
#endif

#ifdef LT_SNII
/*
 * the same as for SnIa.
 */
extern double ***SnIIYields, **SnIIEj;
extern double ***SnII_ShortLiv_Yields;
extern double **IIZbins, **IIMbins;
extern int *IIZbins_dim, *IIMbins_dim, *SnIIY_dim[LT_NMet];
extern double ***SnII_steps;
extern int *SnII_Nsteps;
extern int *SnIIY_Give_Produced;
#endif

#ifdef LT_AGB
extern double ***AGBYields, **AGBEj;
extern double **AGBZbins, **AGBMbins;
extern int *AGBZbins_dim, *AGBMbins_dim, *AGB_dim[LT_NMet];
extern int *AGBY_Give_Produced;
#endif

#if defined(LT_AGB) || defined(LT_SnIa)
extern double ***LLv_steps;
extern int *LLv_Nsteps;

#endif

/* ------- *
 *   IMF   *
 * ------- */

typedef enum {power_law, whatever} IMF_SPEC;
extern char *IMF_Spec_Labels[];
extern IMF_SPEC IMF_Spec;
extern int IsThere_TimeDep_IMF;
extern int IsThere_ZDep_IMF;

#define IMF_NSPEC 15 /* defines how many fields must be specified in IMFs' format file for each IMF */
#define INC_BH 0
#define EXC_BH 1

typedef struct
{
  IMF_SPEC type;
  double (* IMFfunc_byMass)(double, void*);
  double (* IMFfunc_byNum)(double, void*);
  double (* IMFfunc_byEgy)(double, void*);
  double (* getp)(int, double*, double);

  double A, MU, Mm;
  double inf_lifetime, sup_lifetime;
  double referenceZ_toset_SF_DensTh;
  double egyShortLiv_MassTh, egyShortLiv_TimeTh;
  double metShortLiv_MassTh, metShortLiv_TimeTh;
  double ShortLiv_MassTh, ShortLiv_TimeTh;
  double *Params;
  double *PhysDensThresh, *FEVP, EgySpecSN, FactorSN, totFactorSN, metFactorSN;
  double MassFrac_inIRA, NumFrac_inIRA, totResFrac;
  struct
  {
    double *masses, *slopes;
  }Slopes;
  struct
  {
    double *masses, *ekin;
  }EKin;
  struct
  {
    double *sup, *inf;
  }notBH_ranges;
  double IRA_erg_per_g, TOT_erg_per_g;
  double WindEnergy, WindEnergyFraction, WindEfficiency;
  int NSlopes, NParams, YSet, NEKin, N_notBH_ranges;
  int SFTh_Zdep, timedep, SF_DensTh_Criterion;
  int Generations;
  int nonZeroIRA;
}IMF_Type;

extern IMF_Type *IMFs, *IMFp, IMFu;
extern int IMFs_dim;

extern double IMF_func_arg;
extern double *IMF_CommBuff;

extern FILE *FdIMFin, *FdIMF;


/* ---------- *
 *   SEvDbg   *
 * ---------- */

#ifdef LT_SEvDbg
extern unsigned int FirstID;
extern int checkFirstID;
#endif


/* ------------ *
 *   SEv_INFO   *
 * ------------ */

#ifdef LT_SEv_INFO

/*   temporary   */

/* extern double delta_Trun, *Tpop_num, *Tpop_mass; */
/* extern double *tot_Tpop_num, *tot_Tpop_mass; */
/* extern double delta_TCrun, *TCpop_num, *TCpop_mass; */
/* extern double *tot_TCpop_num, *tot_TCpop_mass; */
/* extern int Trun_size, TCrun_size; */
/* extern FILE *FdTrun, *FdTCrun; */

/* ............. */


#define SEvInfo_GRAIN 15

extern FILE *FdSn, *FdSnLost, *FdMetals, *FdSnDetails;

#ifdef LT_SEvDbg
extern FILE *FdMetSumCheck;
#endif

#ifdef WINDS
extern FILE *FdWinds;
#endif

#ifdef LT_EXTEGY_INFO
extern FILE *FdExtEgy;
#endif

#endif /* close SEv_INFO */

#ifdef LT_SEv_INFO_DETAILS
extern double DetailsW[LT_NMet], DetailsWo[LT_NMet];
#endif

#ifdef LT_CharT_INFO

#define CharT_Info_GRAIN 20
extern FILE *FdCharT;
#endif

/* ------------------------------------------------ */
/* ********************************
 *             END STELLAR EVOLUTION
 ********************************** */
#endif


#ifdef LT_METAL_COOLING
/* ********************************
 * 
 *                    METAL COOLING
 *
 ********************************** */

#define Zlength 91
#define ZMin -4.0
#define ZMax 0.5
#define TMin 4.0
#define TMax 8.5
#define ZBins 8
extern double ZTemp[Zlength], Zvalue[8], Lsuthdop[8][Zlength];
extern double ThInst_onset[ZBins];


#else

#define Zlength 1
#define ZMin -4.0
#define ZMax -4.0
#define TMin 4.0
#define TMax 8.5
#define ZBins 1
extern double ZTemp[Zlength], Zvalue[8], Lsuthdop[8][Zlength];
extern double ThInst_onset[ZBins];

#endif
/* ********************************
 *                end METAL COOLING
 ********************************** */

#if defined(SFR_METALS) || defined (BLACK_HOLES) || defined (BG_SFR)
#define PPP P
#ifdef SFR_FEEDBACK
#define EgySNcgs 1e52
#endif
extern int Flag_promotion;   
extern int Flag_phase; 
#else
#define PPP SphP
#endif

#define DMAX(a,b) (dmax1=(a),dmax2=(b),(dmax1>dmax2)?dmax1:dmax2)
#define DMIN(a,b) (dmin1=(a),dmin2=(b),(dmin1<dmin2)?dmin1:dmin2)
#define IMAX(a,b) (imax1=(a),imax2=(b),(imax1>imax2)?imax1:imax2)
#define IMIN(a,b) (imin1=(a),imin2=(b),(imin1<imin2)?imin1:imin2)



extern int ThisTask;		/*!< the number of the local processor  */
extern int NTask;               /*!< number of processors */
extern int PTask;	        /*!< note: NTask = 2^PTask */

extern double CPUThisRun;	/*!< Sums CPU time of current process */

extern int NumForceUpdate;      /*!< number of active particles on local processor in current timestep  */
extern int NumSphUpdate;        /*!< number of active SPH particles on local processor in current timestep  */

extern int RestartFlag;         /*!< taken from command line used to start code. 0 is normal start-up from
                                  initial conditions, 1 is resuming a run from a set of restart files, while 2
                                  marks a restart from a snapshot file. */

extern char *Exportflag;


extern int Flag_FullStep;       /*!< Flag used to signal that the current step involves all particles */


extern int TreeReconstructFlag;  
                             

extern int NumPart;		/*!< number of particles on the LOCAL processor */
extern int N_gas;		/*!< number of gas particles on the LOCAL processor  */
//Added by JM
extern int flag_accretion; /*!< LOCAL flag. If set, domain decomposition will take place */
//Added by JM

extern long long Ntype[6];      /*!< total number of particles of each type */
extern int NtypeLocal[6];       /*!< local number of particles of each type */
#ifdef LT_STELLAREVOLUTION
extern int N_star;             /*!< number of star particles on the LOCAL processor  */
#endif

extern gsl_rng *random_generator; /*!< the random number generator used */


#ifdef SFR
extern int Stars_converted;	/*!< current number of star particles in gas particle block */

#ifdef SFR_METALS
  extern double Ucrit;
  extern double TotalEnergy;
  extern double DEnergy_spawned, DEnergy_converted;
  extern double DEnergy_radiation, DEnergy_promotion;
  extern double DEnergy_feedback, TotalReservoir;
#ifdef SFR_FEEDBACK
  extern double ESN;
  extern int nhot, ncold;
#endif  
#endif  

#endif


extern double TimeOfLastTreeConstruction; /*!< holds what it says */

extern int *Ngblist;            /*!< Buffer to hold indices of neighbours retrieved by the neighbour search
                                  routines */

#ifdef PERIODIC
extern FLOAT boxSize_X, boxHalf_X;
extern FLOAT boxSize_Y, boxHalf_Y;
extern FLOAT boxSize_Z, boxHalf_Z;
#endif

extern double DomainCorner[3], DomainCenter[3], DomainLen, DomainFac;
extern int DomainMyStart, DomainMyLast;
extern int *DomainStartList, *DomainEndList;



extern double *DomainWork;
extern int *DomainCount;
extern int *DomainCountSph;
extern int *DomainTask;
extern int *DomainNodeIndex;
extern FLOAT *DomainHmax;
#ifdef LT_STELLAREVOLUTION
extern int *DomainCountStar;
#endif


extern struct DomainNODE
{
  FLOAT s[3];
  FLOAT mass;
  FLOAT len;
  FLOAT minbound[3], maxbound[3];
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
  FLOAT maxsoft;
#endif
  int   bitflags;
}
*DomainMoment;


extern peanokey *Key, *KeySorted;

extern struct topnode_data
{
  peanokey Size;
  peanokey StartKey;
  long long Count;
  FLOAT GravCost;
  int Daughter;
  int Pstart;
  int Blocks;
  int Leaf;
} *TopNodes;

extern int NTopnodes, NTopleaves;




extern double RndTable[RNDTABLE];


/* variables for input/output , usually only used on process 0 */


extern char ParameterFile[100];  /*!< file name of parameterfile used for starting the simulation */

extern FILE *FdInfo,   /*!< file handle for info.txt log-file. */
  *FdEnergy,           /*!< file handle for energy.txt log-file. */
  *FdTimings,          /*!< file handle for timings.txt log-file. */
  *FdBalance,          /*!< file handle for balance.txt log-file. */
  *FdCPU;              /*!< file handle for cpu.txt log-file. */
//Added by JM
extern FILE *FdAccretion;
//End of Added by JM

#ifdef SFR
extern FILE *FdSfr;    /*!< file handle for sfr.txt log-file. */
#endif

#ifdef BG_STELLAR_EVOLUTION
extern FILE *FdMetGas;    /*!< file handle for metals_gas.txt log-file. */
extern FILE *FdMetStars;  /*!< file handle for metals_star.txt log-file. */
extern FILE *FdMetSF;     /*!< file handle for metals_sf.txt log-file. */
extern FILE *FdMetTot;    /*!< file handle for metals_tot.txt log-file. */
#endif

#ifdef BLACK_HOLES
extern FILE *FdBlackHoles;    /*!< file handle for blackholes.txt log-file. */
extern FILE *FdBlackHolesDetails;
#endif



#ifdef SFR_METALS
extern FILE *FdMphase; 
extern FILE *FdSNE; 
#if defined(SFR_SNI) || defined(SFR_SNII)
extern FILE *FdSN;
#endif    
#ifdef SFR_PROMOTION
extern FILE *FdPromotion;
#endif
#endif

#ifdef FORCETEST
extern FILE *FdForceTest;  /*!< file handle for forcetest.txt log-file. */
#endif

#ifdef DARKENERGY
extern FILE *FdDE;  /*!< file handle for darkenergy.txt log-file. */
#endif

#ifdef XXLINFO
extern FILE *FdXXL;  /*!< file handle for xxl.txt log-file. */
#ifdef MAGNETIC
extern double MeanB;
#ifdef TRACEDIVB
extern double MaxDivB;
#endif
#endif
#ifdef TIME_DEP_ART_VISC
extern double MeanAlpha;
#endif 
#endif

/*! table for the cosmological drift factors */
extern double DriftTable[DRIFT_TABLE_LENGTH];  

/*! table for the cosmological kick factor for gravitational forces */
extern double GravKickTable[DRIFT_TABLE_LENGTH];

/*! table for the cosmological kick factor for hydrodynmical forces */
extern double HydroKickTable[DRIFT_TABLE_LENGTH];

extern void *CommBuffer;   /*!< points to communication buffer, which is used in the domain decomposition, the
                             parallel tree-force computation, and the SPH routines. */

/*! This structure contains data which is the SAME for all tasks (mostly code parameters read from the
 * parameter file).  Holding this data in a structure is convenient for writing/reading the restart file, and
 * it allows the introduction of new global variables in a simple way. The only thing to do is to introduce
 * them into this structure.
 */
extern struct global_data_all_processes
{
  long long TotNumPart;		/*!<  total particle numbers (global value) */
  long long TotN_gas;		/*!<  total gas particle number (global value) */

//Added by JM
  int dom_dec_flag;         /*!< if set, domain decomposition takes place */
//Added by JM

  int MaxPart;			/*!< This gives the maxmimum number of particles that can be stored on one
				     processor. */
  int MaxPartSph;		/*!< This gives the maxmimum number of SPH particles that can be stored on one
				     processor. */

  int ICFormat;			/*!< selects different versions of IC file-format */

  int SnapFormat;		/*!< selects different versions of snapshot file-formats */

  int NumFilesPerSnapshot;      /*!< number of files in multi-file snapshot dumps */
  int NumFilesWrittenInParallel; /*!< maximum number of files that may be written simultaneously when
                                   writing/reading restart-files, or when writing snapshot files */ 

  int BufferSize;		/*!< size of communication buffer in MB */
  int BunchSizeForce;		/*!< number of particles fitting into the buffer in the parallel tree-force
				   algorithm  */
  int BunchSizeDensity;         /*!< number of particles fitting into the communication buffer in the density
                                  computation */
#ifdef FOF
  int BunchSizeFoF;
#endif

#ifdef BG_SFR
  int BunchSizeSolidAngle;
#endif

#if defined(SFR_METALS) || defined(BG_SFR)
  int BunchSizeMetal;
#ifdef SFR_PROMOTION
  int BunchSizeHotNgbs;
#endif
#endif  

#ifdef BG_SFR
  double GammaEffective;
#endif

#if defined(BG_SFR) && defined(BG_STELLAR_EVOLUTION) && defined(BG_COOLING)
  double MetalDependentDensThresh;
  double EnergySNIa;
#endif

#ifdef BG_WINDS
  double WindSpeed;
  double WindMassLoading;
#endif

#ifdef BG_COOLING
  char CoolTablePath[100];
#endif

#ifdef BG_STELLAR_EVOLUTION
  int AGBOn;
  int SNIaOn;
  int SNIIOn;
  char YieldTablePath[100];
#endif

#ifdef LT_STELLAREVOLUTION
  int BunchSizeMetal_Ngb,
    BunchSizeMetal_Spread;
#endif

#ifdef BLACK_HOLES
  int BunchSizeBlackhole;
#endif

#ifdef MHM
  int BunchSizeKinetic;
#endif

  int BunchSizeHydro;           /*!< number of particles fitting into the communication buffer in the SPH
                                  hydrodynamical force computation */
  int BunchSizeDomain;          /*!< number of particles fitting into the communication buffer in the domain
                                  decomposition */

  double PartAllocFactor;	/*!< in order to maintain work-load balance, the particle load will usually
				   NOT be balanced.  Each processor allocates memory for PartAllocFactor times
				   the average number of particles to allow for that */

  double TreeAllocFactor;	/*!< Each processor allocates a number of nodes which is TreeAllocFactor times
				   the maximum(!) number of particles.  Note: A typical local tree for N
				   particles needs usually about ~0.65*N nodes. */

  /* some SPH parameters */

  int DesNumNgb;                /*!< Desired number of SPH neighbours */
  double MaxNumNgbDeviation;    /*!< Maximum allowed deviation neighbour number */

  double ArtBulkViscConst;      /*!< Sets the parameter \f$\alpha\f$ of the artificial viscosity */
  double InitGasTemp;		/*!< may be used to set the temperature in the IC's */
  double InitGasU;		/*!< the same, but converted to thermal energy per unit mass */
  double MinGasTemp;		/*!< may be used to set a floor for the gas temperature */
  double MinEgySpec;            /*!< the minimum allowed temperature expressed as energy per unit mass */



  /* some force counters  */

  long long TotNumOfForces;	/*!< counts total number of force computations  */

  long long NumForcesSinceLastDomainDecomp;  /*!< count particle updates since last domain decomposition */

  /* some variable for dynamic work-load adjustment based on CPU measurements */

  double Cadj_Cost;
  double Cadj_Cpu;

  /* system of units  */

  double UnitTime_in_s,   	/*!< factor to convert internal time unit to seconds/h */
    UnitMass_in_g,             	/*!< factor to convert internal mass unit to grams/h */
    UnitVelocity_in_cm_per_s,   /*!< factor to convert intqernal velocity unit to cm/sec */
    UnitLength_in_cm,           /*!< factor to convert internal length unit to cm/h */
    UnitPressure_in_cgs,        /*!< factor to convert internal pressure unit to cgs units (little 'h' still
                                     around!) */
    UnitDensity_in_cgs,         /*!< factor to convert internal length unit to g/cm^3*h^2 */
    UnitCoolingRate_in_cgs,     /*!< factor to convert internal cooling rate to cgs units */
    UnitEnergy_in_cgs,          /*!< factor to convert internal energy to cgs units */
    UnitTime_in_Megayears,      /*!< factor to convert internal time to megayears/h */
    GravityConstantInternal,    /*!< If set to zero in the parameterfile, the internal value of the
                                  gravitational constant is set to the Newtonian value based on the system of
                                  units specified. Otherwise the value provided is taken as internal gravity
                                  constant G. */
    G;                          /*!< Gravity-constant in internal units */

  /* Cosmology */

  double Hubble;  /*!< Hubble-constant in internal units */
  double Omega0,  /*!< matter density in units of the critical density (at z=0)*/
    OmegaLambda,  /*!< vaccum energy density relative to crictical density (at z=0) */
    OmegaBaryon,  /*!< baryon density in units of the critical density (at z=0)*/
    HubbleParam;  /*!< little `h', i.e. Hubble constant in units of 100 km/s/Mpc.  Only needed to get absolute
		   * physical values for cooling physics
                   */

  double BoxSize; /*!< Boxsize in case periodic boundary conditions are used */

  /* Code options */

  int ComovingIntegrationOn;	/*!< flags that comoving integration is enabled */
  int PeriodicBoundariesOn;     /*!< flags that periodic boundaries are enabled */
  int ResubmitOn;               /*!< flags that automatic resubmission of job to queue system is enabled */
  int TypeOfOpeningCriterion;   /*!< determines tree cell-opening criterion: 0 for Barnes-Hut, 1 for relative
                                  criterion */
  int TypeOfTimestepCriterion;  /*!< gives type of timestep criterion (only 0 supported right now - unlike
                                  gadget-1.1) */
  int OutputListOn;             /*!< flags that output times are listed in a specified file */
  int CoolingOn;                /*!< flags that cooling is enabled */
  int StarformationOn;          /*!< flags that star formation is enabled */


  /* parameters determining output frequency */

  int SnapshotFileCount;     /*!< number of snapshot that is written next */
  double TimeBetSnapshot,    /*!< simulation time interval between snapshot files */
    TimeOfFirstSnapshot,     /*!< simulation time of first snapshot files */
    CpuTimeBetRestartFile,   /*!< cpu-time between regularly generated restart files */
    TimeLastRestartFile,     /*!< cpu-time when last restart-file was written */
    TimeBetStatistics,       /*!< simulation time interval between computations of energy statistics */
    TimeLastStatistics;      /*!< simulation time when the energy statistics was computed the last time */
  int NumCurrentTiStep;      /*!< counts the number of system steps taken up to this point */

  /* Current time of the simulation, global step, and end of simulation */

  double Time,  /*!< current time of the simulation */
    TimeBegin,  /*!< time of initial conditions of the simulation */
    TimeStep,   /*!< difference between current times of previous and current timestep */
    TimeMax;	/*!< marks the point of time until the simulation is to be evolved */

  /* variables for organizing discrete timeline */

  double Timebase_interval; /*!< factor to convert from floating point time interval to integer timeline */
  int Ti_Current;           /*!< current time on integer timeline */ 
  int Ti_nextoutput;        /*!< next output time on integer timeline */

#ifdef FLEXSTEPS
  int PresentMinStep;		/*!< If FLEXSTEPS is used, particle timesteps are chosen as multiples of the present minimum timestep. */
  int PresentMaxStep;		/*!< If FLEXSTEPS is used, this is the maximum timestep in timeline units, rounded down to the next power 2 division */
#endif


#ifdef PMGRID
  int PM_Ti_endstep, PM_Ti_begstep;
  double Asmth[2], Rcut[2];
  double Corner[2][3], UpperCorner[2][3], Xmintot[2][3], Xmaxtot[2][3];
  double TotalMeshSize[2];
#endif

#ifdef CHEMISTRY
  double Epsilon;
#endif 

  int Ti_nextlineofsight;
#ifdef OUTPUTLINEOFSIGHT
  double TimeFirstLineOfSight;
#endif

  /* variables that keep track of cumulative CPU consumption */

  double TimeLimitCPU;
  double CPU_TreeConstruction;
  double CPU_TreeWalk;
  double CPU_Gravity;
  double CPU_Potential;
  double CPU_Domain;
  double CPU_Snapshot;
  double CPU_Total;
  double CPU_CommSum;
  double CPU_Imbalance;
  double CPU_HydCompWalk;
  double CPU_HydCommSumm;
  double CPU_HydImbalance;
  double CPU_Hydro;
  double CPU_EnsureNgb;
  double CPU_Predict;
  double CPU_TimeLine;
  double CPU_PM;
  double CPU_Peano;
#ifdef COOLING
  double CPU_SfrCool;
#endif
#ifdef LT_STELLAREVOLUTION
  double CPU_Eff_Iter, CPU_EE_info;
  double CPU_SEv;
  double CPU_SN_info;
  double CPU_SN_find;
  double CPU_SN_NeighFind;
  double CPU_SN_NeighCheck;
  double CPU_SN_Spread;
  double CPU_SN_Calc;
  double CPU_SN_Comm;
  double CPU_SN_Imbalance;
#endif

  /* tree code opening criterion */

  double ErrTolTheta;		/*!< BH tree opening angle */
  double ErrTolForceAcc;	/*!< parameter for relative opening criterion in tree walk */


  /* adjusts accuracy of time-integration */

  double ErrTolIntAccuracy;	/*!< accuracy tolerance parameter \f$ \eta \f$ for timestep criterion. The
                                  timesteps is \f$ \Delta t = \sqrt{\frac{2 \eta eps}{a}} \f$ */

  double MinSizeTimestep,       /*!< minimum allowed timestep. Normally, the simulation terminates if the
                                  timestep determined by the timestep criteria falls below this limit. */ 
         MaxSizeTimestep;       /*!< maximum allowed timestep */

  double MaxRMSDisplacementFac; /*!< this determines a global timestep criterion for cosmological simulations
                                     in comoving coordinates.  To this end, the code computes the rms velocity
                                     of all particles, and limits the timestep such that the rms displacement
                                     is a fraction of the mean particle separation (determined from the
                                     particle mass and the cosmological parameters). This parameter specifies
                                     this fraction. */



  double CourantFac;		/*!< SPH-Courant factor */


  /* frequency of tree reconstruction/domain decomposition */


  double TreeDomainUpdateFrequency; /*!< controls frequency of domain decompositions  */


  /* gravitational and hydrodynamical softening lengths (given in terms of an `equivalent' Plummer softening
   * length)
   *
   * five groups of particles are supported 0=gas,1=halo,2=disk,3=bulge,4=stars
   */
  double MinGasHsmlFractional, /*!< minimum allowed SPH smoothing length in units of SPH gravitational
                                  softening length */
    MinGasHsml;                /*!< minimum allowed SPH smoothing length */


  double SofteningGas,    /*!< for type 0 */ 
    SofteningHalo,        /*!< for type 1 */ 
    SofteningDisk,        /*!< for type 2 */ 
    SofteningBulge,       /*!< for type 3 */ 
    SofteningStars,       /*!< for type 4 */ 
    SofteningBndry;       /*!< for type 5 */ 

  double SofteningGasMaxPhys,   /*!< for type 0 */ 
    SofteningHaloMaxPhys,       /*!< for type 1 */ 
    SofteningDiskMaxPhys,       /*!< for type 2 */ 
    SofteningBulgeMaxPhys,      /*!< for type 3 */ 
    SofteningStarsMaxPhys,      /*!< for type 4 */ 
    SofteningBndryMaxPhys;      /*!< for type 5 */ 

  double SofteningTable[6];  /*!< current (comoving) gravitational softening lengths for each particle type */
  double ForceSoftening[6];  /*!< the same, but multiplied by a factor 2.8 - at that scale the force is Newtonian */


  /*! If particle masses are all equal for one type, the corresponding entry in MassTable is set to this
   *  value, * allowing the size of the snapshot files to be reduced
   */
  double MassTable[6];


  /* some filenames */
  char InitCondFile[100],
    OutputDir[100],
    SnapshotFileBase[100],
    EnergyFile[100],
    CpuFile[100],
    InfoFile[100], 
    TimingsFile[100], 
    RestartFile[100], 
    ResubmitCommand[100], 
    OutputListFilename[100];

  /*! table with desired output times */
  double OutputListTimes[MAXLEN_OUTPUTLIST];

  int OutputListLength; /*!< number of times stored in table of desired output times */



#if defined(ADAPTIVE_GRAVSOFT_FORGAS) && !defined(ADAPTIVE_GRAVSOFT_FORGAS_HSML)
  double ReferenceGasMass;
#endif

#ifdef MOREPARAMS		/* star formation and feedback sector */
  double CritOverDensity;
  double CritPhysDensity;
  double OverDensThresh;
  double PhysDensThresh;
  double OrigGasMass; 
#ifndef BG_SFR
  double EgySpecCold;
#ifndef LT_STELLAREVOLUTION
  double EgySpecSN;
  double FactorSN;
#endif
    double FactorEVP;
  double FeedbackEnergy;
  double TempSupernova;
  double TempClouds;
  double MaxSfrTimescale;
  double WindEfficiency;
  double WindEnergyFraction;
  double WindFreeTravelLength;
  double WindFreeTravelDensFac;
  double FactorForSofterEQS;
#endif

#ifdef SFR_METALS
  double FactorSFR;
  double MinTlifeSNI;
  double MaxTlifeSNI;
  double TlifeSNII;
  double RateSNI;
  double FactorSN_Phase;
  double Tcrit_Phase;
  double DensFrac_Phase;
  double FracEnergySN_Phase;
#ifdef SFR_DECOUPLING
  double DensityTailThreshold; 
#endif
#endif

#endif

#ifdef DARKENERGY
  double DarkEnergyParam;	/*!< fixed w for equation of state */
#ifdef TIMEDEPDE
  char DarkEnergyFile[100];	/*!< tabelized w for equation of state */
#ifdef TIMEDEPGRAV
  double Gini;
#endif
#endif
#endif

#ifdef RESCALEVINI
  double VelIniScale;		/*!< Scale the initial velocities by this amount */
#endif

#ifdef TIME_DEP_ART_VISC
  double ViscSource0;		/*!< Given sourceterm in viscosity evolution */
  double DecayLength;		/*!< Number of h for the viscosity decay */
  double ViscSource;		/*!< Reduced sourceterm in viscosity evolution*/
  double DecayTime;		/*!< Calculated decaytimescale */
  double AlphaMin;		/*!< Minimum of allowed viscosity parameter */
#endif

#ifdef CONDUCTION
  double ConductionCoeff;         /*!< Thermal Conductivity */
#ifdef CONDUCTION_SATURATION
  double ElectronFreePathFactor;  /*!< Factor to get electron mean free path */
#endif
#endif

#ifdef BINISET                 
  double BiniX,BiniY,BiniZ;       /*!< Initial values for B */
#endif

#ifdef BSMOOTH
  int BSmoothInt;
  double BSmoothFrac;
  int MainTimestepCounts;
#endif

#ifdef MAGNETIC_DISSIPATION
  double ArtMagDispConst;      /*!< Sets the parameter \f$\alpha\f$ of the artificial magnetic disipation */
#ifdef TIME_DEP_MAGN_DISP
  double ArtMagDispMin;
  double ArtMagDispSource;
  double ArtMagDispTime;
#endif
#endif

#ifdef DIVBCLEANING_DEDNER
  double DivBcleanParabolicSigma;
  double DivBcleanHyperbolicSigma;
#endif

#ifdef VOLUME_CORRECTION
  char VolCorrectFile[100];	/*!< tabelized volume correction function */
#endif

#ifdef BLACK_HOLES
  double TimeNextBlackHoleCheck;
  double TimeBetBlackHoleSearch;
  double BlackHoleAccretionFactor;  /*!< Fraction of BH bondi accretion rate */
  double BlackHoleFeedbackFactor;   /*!< Fraction of the black luminosity feed into thermal feedback */
    double SeedBlackHoleMass;         /*!< Seed black hole mass */
  double MinFoFMassForNewSeed;      /*!< Halo mass required before new seed is put in */
  double BlackHoleNgbFactor;        /*!< Factor by which the normal SPH neighbour should be increased/decreased */
  double BlackHoleActiveTime;
  double BlackHoleEddingtonFactor; /*! Factor above Eddington */
#ifdef FOF
  double massDMpart;
#endif
#ifdef MODIFIEDBONDI
  double BlackHoleRefDensity;
  double BlackHoleRefSoundspeed;
#endif
#endif

#ifdef COSMIC_RAYS
  double CR_Alpha;               /*!< Cosmic ray spectral index [2..3]*/
  double CR_SNEff;               /*!< SN injection efficiency [0..1] */
  double CR_SNAlpha;             /*!< SN injection spectral index [2..3] */
  int bDebugFlag;                /*!< enables debug outputs after triggered */

#if defined(CR_DIFFUSION) || defined (CR_DIFFUSION_GREEN)
  double CR_DiffusionCoeff;      /*!< (temporary) fixed value for CR diffusivity */

  double CR_DiffusionDensScaling; /*!< grade of density dependence of diffusivity */
  double CR_DiffusionDensZero;    /*!< Reference point density for diffusivity */

  double CR_DiffusionEntropyScaling; /*!< grade of specific energy dependence of diffusivity */

  double CR_DiffusionEntropyZero; /*!< Reference Entropic function for diffusivity */

  double CR_DiffusionTimeScale; /*!< Parameter for Diffusion Timestep Criterion */

  double TimeOfLastDiffusion;
#endif /* CR_DIFFUSION */

#if defined(CR_SHOCK)
#if (CR_SHOCK == 1)
  double CR_ShockAlpha;          /*!< spectral index to be used in shock injection */
#else
  double CR_ShockCutoff;         /*!< Cutoff factor x_inj for CR accel */
#endif
  double CR_ShockEfficiency;     /*!< energy fraction of shock energy fed into CR */
#endif /* CR_SHOCK */

#ifdef FIX_QINJ
  double Shock_Fix_Qinj;         /*!< inject only CRps with threshold cutoff Shock_Fix_Qinj */
#endif

#endif /* COSMIC_RAYS */

#ifdef MACHNUM
  double Shock_Length;            /*!< length scale on which the shock is smoothed out */
  double Shock_DeltaDecayTimeMax; /*!< maximum time interval (Dloga) for which the 
				       Mach number is kept at its maximum */
#endif

#ifdef REIONIZATION
  int not_yet_reionized;          /*!< flag that makes sure that there is only one reionization */
#endif


#ifdef HPM
  double HPM_entr0, HPM_entr1;
  double HPM_ne0, HPM_ne1;
  double HPM_rho0, HPM_rho1;
  double HPM_P0, HPM_P1, HPM_alpha;
#endif

#ifdef LT_STELLAREVOLUTION
  long long TotN_star;            /*!<  total gas particle number (global value) */
  int MaxPartMet;                 /*!< This gives the maxmimum number of star particles that can be stored on one
				     processor. */
  int Generations;
  double SFfactor;                 /*!< the expected maximum factor of proportionality between TotN_gas and TotN_star */

  double Time_Age;
  int NeighInfNum;                /* minimum number of neighbours for metal and egy spreading */
  int DesNumNgbSN,                /* desired number of neighbours */
    SpreadNumNgbDev,              /* maximum deviation */
    LeftNumNgbSN, RightNumNgbSN;  /* range of neighbours */
  double SpreadNeighCoeff;        /* store DesNumNgbSN / DesNumNgb */

  double MinChemSpreadL;
  double Enrich_SFGas_Th;
  /*
  double MaxChemSpreadL;
    - so far we use different timesteps for different epochs of a star life
  double ChemTimeStepIa, ChemTimeStepII, LongChemTimeStepII;
  */

  double MinChemTimeStep;
  double SnII_Step_Prec, LLv_Step_Prec;
  int Sn_Step_Criterion;
  double referenceZ_toset_SF_DensTh;
  int SFTh_Zdep, referenceZbin_SFTh;

  /* mass function exponent x for IMF phi(m) = m^ -(1 + x) */
  double Xphi;
  /* normalization for the IMF */
  double Aphi;
  /* stellar definition */
  double Mup,   /* SnII threshold (e.g. 8 Msun) */
    MBm,        /* min bin. system mass */
    MBM,        /* max bin. system mass */
    BinFrac,    /* bin. system fraction */
    MBms;       /* min star mass in SnIa binary systems */
  double SnIaRemn; /* SnIa Remn (1.4Msun ?) */
  /* lifetimes */
  double inf_lifetime, mean_lifetime, sup_lifetime;
  /* energy provided by sn explosions */
  double SnIaEgy, SnIIEgy;
  /* define IRA range */
  double metIRA_ThMass, egyIRA_ThMass;

  double WindEnergy;
#ifdef LT_HOT_EJECTA
  double EgySpecEjecta;
#endif
  int Mod_SEff;
  char IMFfilename[300], SnIaDataFile[300], SnIIDataFile[300], AGBDataFile[300];
  int Ia_Nset_ofYields, II_Nset_ofYields, AGB_Nset_ofYields;
#endif

#ifdef BUBBLES
  double BubbleDistance;
  double BubbleRadius;
  double BubbleTimeInterval;
  double BubbleEnergy;
  double TimeOfNextBubble;
  double FirstBubbleRedshift;
#ifdef FOF
  int    BiggestGroupLen;
  float  BiggestGroupCM[3];
  double BiggestGroupMass;   
#endif
#endif

#if defined(MULTI_BUBBLES) && defined(FOF)
  double MinFoFMassForNewSeed;      /*!< Halo mass required before new seed is put in */
  double BubbleDistance;
  double BubbleRadius;
  double BubbleTimeInterval;
  double BubbleEnergy;
  double TimeOfNextBubble;
  double ClusterMass200;
  double massDMpart;
  double FirstBubbleRedshift;
#endif

#ifdef NAVIERSTOKES
  double NavierStokes_ShearViscosity;
  double FractionSpitzerViscosity;
  double ShearViscosityTemperature;
#endif
#ifdef NAVIERSTOKES_BULK
  double NavierStokes_BulkViscosity;
#endif
#ifdef VISCOSITY_SATURATION
  double IonMeanFreePath;
#endif
}
All;




/*! This structure holds all the information that is
 * stored for each particle of the simulation.
 */
extern struct particle_data
{
  FLOAT Pos[3];			/*!< particle position at its current time */
  FLOAT Mass;			/*!< particle mass */
  FLOAT Vel[3];			/*!< particle velocity at its current time */
  union
  {  
    FLOAT GravAccel[3];		/*!< particle acceleration due to gravity */
    DOUBLE dGravAccel[3];
  } g;
#ifdef PMGRID
  FLOAT GravPM[3];		/*!< particle acceleration due to long-range PM gravity force*/
#endif
#ifdef FORCETEST
  FLOAT GravAccelDirect[3];	/* particle acceleration */
#endif
  union
  {
    FLOAT Potential;		/*!< gravitational potential */
    DOUBLE dPotential;
  } p;
  FLOAT OldAcc;			/*!< magnitude of old gravitational force. Used in relative opening
                                  criterion */

#if defined(EVALPOTENTIAL) && defined(PMGRID)
  FLOAT PM_Potential;
#endif

#ifdef STELLARAGE
  FLOAT StellarAge;		/*!< formation time of star particle */
#endif

#ifdef BG_SFR
  FLOAT SolidAngleWeightSum;
  FLOAT Metals[BG_NELEMENTS];
  FLOAT InitialMass;
  FLOAT Metallicity;
#endif

#ifdef METALS
#ifdef SFR_METALS
  FLOAT Zm[12];
  FLOAT ZmReservoir[12];
#ifdef SFR_FEEDBACK
  FLOAT EnergySN;
  FLOAT EnergySNCold;    
#endif    
#else
  FLOAT Metallicity;		/*!< metallicity of gas or star particle */
#endif
#endif  /* closes METALS */

#if defined(SFR_METALS) || defined (BLACK_HOLES) || defined (BG_SFR)
  FLOAT Hsml;	
  FLOAT Left,                   /*!< lower bound in iterative smoothing length search */  
    Right;                      /*!< upper bound in iterative smoothing length search */ 
#ifndef  NOFIXEDMASSINKERNEL
 union
  {
    FLOAT NumNgb;
    DOUBLE dNumNgb;
  } a1;
#else
  int NumNgb;			/*!< (effective) number of SPH neighbours */
#endif
#endif

#ifndef LONGIDS
  unsigned int ID;			/*!< particle identifier */
#else
  unsigned long long ID;
#endif
  int Type;		   /*!< flags particle type.  0=gas, 1=halo, 2=disk, 3=bulge, 4=stars, 5=bndry */
  int Ti_endstep;          /*!< marks start of current timestep of particle on integer timeline */ 
  int Ti_begstep;          /*!< marks end of current timestep of particle on integer timeline */
#ifdef FLEXSTEPS
  int FlexStepGrp;		/*!< a random 'offset' on the timeline to create a smooth groouping of particles */
#endif
  float GravCost;		/*!< weight factor used for balancing the work-load */
#ifdef PSEUDOSYMMETRIC
  float AphysOld;               /*!< magnitude of acceleration in last timestep. Used to make a first order
                                  prediction of the change of acceleration expected in the future, thereby
                                  allowing to guess whether a decrease/increase of the timestep should occur
                                  in the timestep that is started. */
#endif
#ifdef LT_STELLAREVOLUTION
  unsigned int MetID;
#endif

#ifdef BLACK_HOLES
  int   SwallowID;                                  
  FLOAT BH_Mass;
  FLOAT BH_Mdot;
  FLOAT BH_MdotEddington;   /* in units of the Eddington accretion rate */
  union
  {
    FLOAT BH_Density;
    DOUBLE dBH_Density;
  } b1;
  union
  {
    FLOAT BH_Entropy;
    DOUBLE dBH_Entropy;
  } b2;
  union
  {
    FLOAT BH_SurroundingGasVel[3];
    DOUBLE dBH_SurroundingGasVel[3];
  } b3;
  union 
  {
    FLOAT BH_accreted_Mass;
    DOUBLE dBH_accreted_Mass;
  } b4;
  union 
  {
    FLOAT BH_accreted_BHMass;
    DOUBLE dBH_accreted_BHMass;
  } b5;
  union
  {
    FLOAT BH_accreted_momentum[3];
    DOUBLE dBH_accreted_momentum[3];
  } b6;
#ifdef REPOSITION_ON_POTMIN
  FLOAT BH_MinPotPos[3];
  FLOAT BH_MinPot;
#endif
#ifdef BH_KINETICFEEDBACK
  FLOAT ActiveTime;
  FLOAT ActiveEnergy;
#endif
#endif
}
 *P,              /*!< holds particle data on local processor */
 *DomainPartBuf;  /*!< buffer for particle data used in domain decomposition */

#ifdef LT_STELLAREVOLUTION
extern struct met_particle_data
{
  int NextChemStepIa, NextChemStepII;   /* the next time of evolution for Ia and II */
  float LastIaTime, LastIITime;         /* the last time of evolution in units of Gyr */
  float iMass;                          /* initial mass of this SSP */
  float Metals[LT_NMetP];              /* the metal array (H is not stored here) */
  FLOAT SLguess;                        /* the spreading lenght */
  double weight;                        /* used in spreading */
  float Left, Right;                    /* bounds in iterative spread lenght search */
  int NumNgb;                           /* number of sph neigbours */
  long int PID;
#ifdef LT_SMOOTH_Z
  float Zsmooth;
#endif
#ifdef LT_TRACK_CONTRIBUTES
  Contrib contrib;
#endif
}
 *MetP,                                 /*!< holds metal particle data on local processor */          
 *DomainMetBuf; 			/*!< buffer for metal data in domain decomposition */
#endif


/* the following struture holds data that is stored for each SPH particle in addition to the collisionless
 * variables.
 */
extern struct sph_particle_data
{
  FLOAT Entropy;                /*!< current value of entropy (actually entropic function) of particle */
  union
  {
    FLOAT Density;		/*!< current baryonic mass density of particle */
    DOUBLE dDensity;
  } a2;

#if !defined(SFR_METALS) && !defined(BLACK_HOLES) && !defined(BG_SFR)
  FLOAT Hsml;			/*!< current smoothing length */
  FLOAT Left,                   /*!< lower bound in iterative smoothing length search */  
    Right;                      /*!< upper bound in iterative smoothing length search */ 
#ifndef  NOFIXEDMASSINKERNEL
  union
  {
    FLOAT NumNgb;
    DOUBLE dNumNgb;
  } a1;
#else
  int NumNgb;                   /*!< (effective) number of SPH neighbours */
#endif
#endif
  FLOAT Pressure;		/*!< current pressure */
  union
  {
    FLOAT DtEntropy;              /*!< rate of change of entropy */
    DOUBLE dDtEntropy;
  } e;
  union
  {
    FLOAT HydroAccel[3];		/*!< acceleration due to hydrodynamical force */
    DOUBLE dHydroAccel[3];
  } a;
  FLOAT VelPred[3];		/*!< predicted SPH particle velocity at the current time */

  union
  {
#ifdef NAVIERSTOKES
    FLOAT DV[3][3];
#endif
    struct
    {
      union
      {
        FLOAT DivVel;			/*!< local velocity divergence */
        DOUBLE dDivVel;
      } a4;
      FLOAT CurlVel;	 	        /*!< local velocity curl */
#ifdef NAVIERSTOKES
      FLOAT StressDiag[3];
      FLOAT StressOffDiag[3];
#ifdef NAVIERSTOKES_BULK
      FLOAT StressBulk;
#endif
#else
      union
      {
        FLOAT Rot[3];		        /*!< local velocity curl */
        DOUBLE dRot[3];
      } a5;
#endif
    }
    s;
  }
  u;

#ifndef NOGRADHSML
  union
  {
    FLOAT DhsmlDensityFactor;     /*!< correction factor needed in the equation of motion of the conservative
                                    entropy formulation of SPH */
    DOUBLE dDhsmlDensityFactor; 
  } a3;
#endif
  FLOAT MaxSignalVel;           /*!< maximum signal velocity */
#ifdef COOLING
  FLOAT Ne;			/*!< electron fraction, expressed as local electron number density normalized
                                  to the hydrogen number density. Gives indirectly ionization state and mean
                                  molecular weight. */
#ifdef SFR
  FLOAT Sfr;
#endif
#endif
#ifdef WINDS
  FLOAT DelayTime;              /*!< remaining maximum decoupling time of wind particle */
#endif

#ifdef MAGNETIC
  FLOAT B[3], BPred[3];
  FLOAT DtB[3];
#if defined(TRACEDIVB) || defined(TIME_DEP_MAGN_DISP)
  FLOAT divB;
#endif
#if defined(BSMOOTH) || defined(BFROMROTA)
  FLOAT BSmooth[3];
  FLOAT DensityNorm;
#endif
#ifdef TIME_DEP_MAGN_DISP
  FLOAT Balpha, DtBalpha;
#endif
#ifdef DIVBCLEANING_DEDNER
  FLOAT Phi,PhiPred,DtPhi;
#ifdef SMOOTH_PHI
  FLOAT SmoothPhi;
#endif
#endif
#endif
#ifdef TIME_DEP_ART_VISC
  FLOAT alpha, Dtalpha;
#endif
#ifdef NS_TIMESTEP
  FLOAT ViscEntropyChange;
#endif
#ifdef CONDUCTION
  FLOAT CondEnergyChange;
  FLOAT SmoothedEntr;
#ifdef CONDUCTION_SATURATION
  FLOAT GradEntr[3];
#endif
#ifdef OUTPUTCOOLRATE
  FLOAT CondRate;
#endif
#endif

#if defined(LT_SMOOTH_Z) || defined(BG_SFR)
  FLOAT Zsmooth;
#endif
#ifdef LT_STELLAREVOLUTION
  double EgyRes;                 /* external (Sn) energy resorvoir */
  float Metals[LT_NMetP];       /* H is not stored here */
#ifdef LT_EJECTA_IN_HOTPHASE
  float x;                       /* stores the last value of the cold clouds mass fraction */
#endif
#ifndef LT_LOCAL_IRA
  double mstar;
#endif
#if !defined(LT_LOCAL_IRA) && defined(NOFIXEDMASSINKERNEL)
  float weight;
#endif
#ifdef LT_TRACK_CONTRIBUTES
  Contrib contrib;
#endif
#endif

#if defined(BH_THERMALFEEDBACK) || defined(BH_KINETICFEEDBACK)
   union
      {
      FLOAT Injected_BH_Energy;
      DOUBLE dInjected_BH_Energy;
      } i;
#endif


#if defined (SFR_METALS) && defined (SFR_DECOUPLING)
  FLOAT DensityOld;     
#endif
#ifdef SFR_PROMOTION
  FLOAT DensityAvg;
  FLOAT EntropyAvg;
  FLOAT HotHsml;
  int   HotNgbNum;
  FLOAT DensPromotion;
  FLOAT TempPromotion;
#endif
#ifdef MHM
  FLOAT FeedbackEnergy;
#endif

#ifdef COSMIC_RAYS
  FLOAT CR_C0;                  /*!< Cosmic ray amplitude adiabatic invariable */
  FLOAT CR_q0;                  /*!< Cosmic ray cutoff adiabatic invariable */
  FLOAT CR_E0;                  /*!< Specific Energy at Rho0 */
  FLOAT CR_n0;                  /*!< baryon fraction in cosmic rays */

  FLOAT CR_DeltaE;              /*!< Specific Energy growth during timestep */
  FLOAT CR_DeltaN;              /*!< baryon fraction growth during timestep */
#ifdef MACHNUM
  FLOAT CR_Gamma0;  
#endif

#ifdef CR_OUTPUT_INJECTION
  FLOAT CR_Specific_SupernovaHeatingRate;
#endif

#ifdef CR_DIFFUSION
  FLOAT CR_SmoothE0;         /*!< SPH-smoothed interpolant of diffusion term */
  FLOAT CR_Smoothn0;         /*!< SPH-smoothed interpolant for diffusion term */
#endif /* CR_DIFFUSION */
#ifdef CR_DIFFUSION_GREEN
  FLOAT CR_Kappa;
  FLOAT CR_Kappa_egy;
  FLOAT CR_WeightSum;
  FLOAT CR_WeightSum_egy;
#endif
#endif /* COSMIC_RAYS */

#ifdef MACHNUM
  FLOAT Shock_MachNumber;           /*!< Mach number */
  FLOAT Shock_DecayTime;            /*!< Shock decay time */
#ifdef COSMIC_RAYS		
  FLOAT Shock_DensityJump;          /*!< Density jump at the shock */
  FLOAT Shock_EnergyJump;           /*!< Energy jump at the shock */
  FLOAT PreShock_PhysicalDensity;   /*!< Specific energy in the preshock regime */
  FLOAT PreShock_PhysicalEnergy;    /*!< Density in the preshock regime */
  FLOAT PreShock_XCR;               /*!< XCR = PCR / Pth in the preshock regime */
#endif
#ifdef MACHSTATISTIC
  FLOAT Shock_DtEnergy;             /*!< Change of thermal specific energy at Shocks */
#endif
#endif /* Mach number estimate */


#ifdef CHEMISTRY
  FLOAT elec;            
  FLOAT HI;
  FLOAT HII;

  FLOAT HeI;
  FLOAT HeII;
  FLOAT HeIII;

  FLOAT H2I;
  FLOAT H2II;

  FLOAT HM;

  FLOAT Gamma;
  FLOAT t_elec, t_cool;
#endif

#ifdef VOLUME_CORRECTION
  FLOAT GradDensity[3];
  FLOAT NormDensity;
#endif

}
 *SphP,                        	/*!< holds SPH particle data on local processor */
 *DomainSphBuf;                 /*!< buffer for SPH particle data in domain decomposition */


extern peanokey *DomainKeyBuf;

/* global state of system 
*/
extern struct state_of_system
{
  double Mass,
    EnergyKin,
    EnergyPot,
    EnergyInt,
    EnergyTot,
    Momentum[4],
    AngMomentum[4],
    CenterOfMass[4],
    MassComp[6],
    EnergyKinComp[6],
    EnergyPotComp[6],
    EnergyIntComp[6], 
    EnergyTotComp[6], 
    MomentumComp[6][4], 
    AngMomentumComp[6][4], 
    CenterOfMassComp[6][4];
}
SysState, SysStateAtStart, SysStateAtEnd;


/* Various structures for communication during the gravity computation.
 */
extern struct gravdata_in
{
  union
  {
    FLOAT Pos[3];
    DOUBLE Acc[3];
  }
  u;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
  FLOAT Soft;
#endif
  union
  {
    FLOAT OldAcc;
    int Ninteractions;
  }
  w;  
#if defined(UNEQUALSOFTENINGS) || defined(EVALPOTENTIAL) || defined(OUTPUTPOTENTIAL) || defined(COMPUTE_POTENTIAL_ENERGY) || defined(OUTPUTPOTENTIAL)
  union 
  {
    DOUBLE Potential;
    int Type;
  }
  v;
#endif
}
*GravDataIn, *GravDataGet, *GravDataResult, *GravDataOut;


extern struct gravdata_index
{
  int Task;
  int Index;
  int SortIndex;
}
 *GravDataIndexTable;



/*! Structure for communication during the density computation. Holds data that is sent to other processors.
 */
extern struct densdata_in
{
  FLOAT Pos[3];
  FLOAT Vel[3];
  FLOAT Hsml;
#ifdef WINDS
  FLOAT DelayTime;
#endif
  int Index;
  int Task;
#ifdef SFR_DECOUPLING
  FLOAT DensityOld;
  FLOAT Entropy;
#endif
#ifdef CR_DIFFUSION_GREEN
  FLOAT CR_Kappa;
  FLOAT CR_Kappa_egy;
  FLOAT CR_WeightSum;
  FLOAT CR_WeightSum_egy;
  FLOAT CR_E0;
  FLOAT CR_n0;
#endif
#ifdef VOLUME_CORRECTION
  FLOAT Density;
#endif
}
 *DensDataIn, *DensDataGet;

#ifdef BG_SFR
extern struct solidangledata_in
{
  FLOAT Pos[3];
  FLOAT Hsml; 
  int Index;
  int Task;
}
 *SolidAngleDataIn, *SolidAngleDataGet;
#endif



#ifdef LT_STELLAREVOLUTION
/* structures for computation of Sn products spreading */

extern struct metaldata_in
{
  FLOAT Pos[3];
  float Metals[LT_NMet];
#if defined(LT_EJECTA_IN_HOTPHASE) || defined(LT_HOT_EJECTA) || defined(LT_SNegy_IN_HOTPHASE)
  float LMMass;
#endif
  FLOAT L;
  double energy;
  double weight;
  int index, Task, Type, IMFi;
#ifdef LT_SEvDbg
  unsigned int ID;
#endif
#ifdef LT_TRACK_CONTRIBUTES
  Contrib contrib;
#endif
} *MetalDataIn, *MetalDataGet;

extern struct metal_ngbfindingdata_in
{
  FLOAT Pos[3];
  FLOAT L;
  int index, Task;
} *MetalNgbIn, *MetalNgbGet;

extern struct metal_ngbfindingdata_out
{
  double weight;
  int num_ngb;
} *MetalNgbResultIn, *MetalNgbResultGet;

#endif


#ifdef BG_SFR
extern struct metaldata_in
{
  FLOAT Pos[3];
  FLOAT Metals[BG_NELEMENTS + 1];
  FLOAT Hsml;
  FLOAT SolidAngleWeightSum;
  FLOAT StarMass;
  FLOAT Energy;
  unsigned int ID;
  int   Index;
  int   Task;
} *MetalDataIn, *MetalDataGet;
#endif



#ifdef SFR_METALS
extern struct metaldata_in
{
  float Pos[3];
  float ZmReservoir[12];
  float Hsml;
#ifdef SFR_FEEDBACK
  float EnergySN;
  FLOAT EnergySNCold;    
#endif  
#ifndef  NOFIXEDMASSINKERNEL
  FLOAT NumNgb;
#else	
  int NumNgb;
#endif
  int Index;
  int Task;
} *MetalDataIn, *MetalDataGet;
#endif


#ifdef MHM
extern struct kindata_in
{
  FLOAT Pos[3];
  FLOAT Hsml;
  FLOAT Density;
  FLOAT Energy;
  int Index;
  int Task;
}
 *KinDataIn, *KinDataGet;
#endif


#ifdef BLACK_HOLES
extern struct blackholedata_in
{
  FLOAT Pos[3];
  FLOAT Density;
  FLOAT Mdot;
  FLOAT Dt;
  FLOAT Hsml;
  FLOAT Mass;
  FLOAT BH_Mass;
  FLOAT Vel[3];
  FLOAT Csnd;
#ifdef BH_KINETICFEEDBACK
  FLOAT ActiveTime;
  FLOAT ActiveEnergy;
#endif
  int ID;
  int Index;
  int Task;
}
 *BlackholeDataIn, *BlackholeDataGet;

extern struct blackholedata_out
{
  DOUBLE Mass;
  DOUBLE BH_Mass;
  DOUBLE AccretedMomentum[3];
#ifdef REPOSITION_ON_POTMIN
  FLOAT BH_MinPotPos[3];
  FLOAT BH_MinPot;
#endif
}
 *BlackholeDataResult, *BlackholeDataPartialResult;
#endif






/*! Structure for communication during the density computation. Holds data that is received from other
 * processors.
 */
extern struct densdata_out
{
  DOUBLE Rho;
#if defined(LT_SMOOTH_Z) || defined(BG_SFR)
  DOUBLE ZRho;
#endif
#ifndef NAVIERSTOKES
  DOUBLE Div, Rot[3];
#else
  FLOAT DV[3][3];
#endif
#ifndef NOGRADHSML
  DOUBLE DhsmlDensity;
#endif
#ifndef NOFIXEDMASSINKERNEL
  DOUBLE Ngb;
#else
  int Ngb;
#endif
#ifdef MAGNETIC 
#if defined(BSMOOTH) || defined(BFROMROTA)
  FLOAT BSmooth[3];
  FLOAT DensityNorm;
#endif
#if defined(SMOOTH_PHI) && defined(DIVBCLEANING_DEDNER)
  FLOAT SmoothPhi;
#endif
#endif

#if defined(CONDUCTION) || defined(BLACK_HOLES)
  DOUBLE SmoothedEntr;
#ifdef CONDUCTION_SATURATION
  FLOAT GradEntr[3];
#endif
#endif

#ifdef CR_DIFFUSION
  FLOAT CR_SmoothE0;
  FLOAT CR_Smoothn0;
#endif

#ifdef BLACK_HOLES
  DOUBLE GasVel[3];
#endif
#ifdef CR_DIFFUSION_GREEN
  FLOAT CR_WeightSum;
  FLOAT CR_WeightSum_egy;
#endif

#ifdef VOLUME_CORRECTION
  FLOAT GradDensity[3];
  FLOAT NormDensity;
#endif
}
 *DensDataResult, *DensDataPartialResult;

#ifdef BG_SFR
extern struct solidangledata_out
{
  FLOAT SolidAngleWeightSum;
}
 *SolidAngleDataResult, *SolidAngleDataPartialResult;
#endif


#ifdef SFR_PROMOTION
extern struct hotngbs_in
{
  FLOAT Pos[3];
  FLOAT HotHsml;
  FLOAT Entropy;
  int Index;
  int Task;
}
 *HotNgbsIn, *HotNgbsGet;

extern struct hotngbs_out
{
  FLOAT DensitySum;
  FLOAT EntropySum;
  int   HotNgbNum;
}
*HotNgbsResult, *HotNgbsPartialResult;
#endif


#ifdef FOF
extern struct fofdata_in
{
  FLOAT Pos[3];
  FLOAT Hsml;
  int MinID;
  int MinIDTask;
  int Index;
  int Task;
}
 *FoFDataIn, *FoFDataGet;

extern struct fofdata_out
{
  FLOAT Distance;
  int   MinID;
  int   MinIDTask;
}
*FoFDataResult, *FoFDataPartialResult;

#endif











/* Various structures for communication during the computation of hydrodynamical forces.
 */
extern struct hydrodata_in
{
  FLOAT Pos[3];
  FLOAT Vel[3];
  FLOAT Hsml;
  FLOAT Mass;
  FLOAT Density;
  FLOAT Pressure;
  FLOAT F1;
#ifndef NOGRADHSML
  FLOAT DhsmlDensityFactor;
#endif
#ifdef MAGNETIC
  FLOAT BPred[3];
#ifdef TIME_DEP_MAGN_DISP
  FLOAT Balpha;
#endif
#ifdef DIVBCLEANING_DEDNER
  FLOAT PhiPred;
#endif
#endif
#ifdef TIME_DEP_ART_VISC
  FLOAT alpha;
#endif
  int   Timestep;
  int   Task;
  int   Index;

#if defined(CONDUCTION) || defined(CR_DIFFUSION) || defined(NAVIERSTOKES)
  FLOAT Entropy;
#endif

#ifdef CONDUCTION
  FLOAT SmoothedEntr;
#ifdef CONDUCTION_SATURATION
  FLOAT GradEntr[3];
#endif
#endif

#ifdef SFR_DECOUPLING
  FLOAT DensityOld;
  FLOAT Entropy;
#endif

#ifdef PARTICLE_DEBUG
#ifndef LONGIDS
  unsigned int ID;			/*!< particle identifier */
#else
  unsigned long long ID;
#endif
#endif 

#ifdef CR_DIFFUSION
  FLOAT CR_E0;       /*!< diffusion term for cosmic ray energy */
  FLOAT CR_n0;       /*!< diffusion term for cosmic ray baryon fraction */
  FLOAT CR_q0;
#endif

#ifdef NAVIERSTOKES
  FLOAT stressoffdiag[3];
  FLOAT stressdiag[3];
  FLOAT shear_viscosity;
#endif

#ifdef NAVIERSTOKES_BULK
  FLOAT divvel;
#endif

}
 *HydroDataIn, *HydroDataGet;

extern struct hydrodata_out
{
  DOUBLE Acc[3];
  DOUBLE DtEntropy;
  FLOAT MaxSignalVel;
#ifdef MAGNETIC
  FLOAT DtB[3];
#if defined(TRACEDIVB) || defined(TIME_DEP_MAGN_DISP)
  FLOAT divB;
#endif
#ifdef DIVBCLEANING_DEDNER
  FLOAT DtPhi;
#endif
#endif
#ifdef CONDUCTION
  FLOAT CondEnergyChange;
#ifdef OUTPUTCOOLRATE
  FLOAT CondRate;
#endif
#endif

#if  defined(CR_SHOCK) || defined(CR_DIFFUSION) 
  FLOAT CR_EnergyChange;
  FLOAT CR_BaryonFractionChange;
#endif

#ifdef LT_SMOOTH_Z
    FLOAT ZRho;
#endif
}
 *HydroDataResult, *HydroDataPartialResult;



/*! Header for the standard file format.
 */
extern struct io_header
{
  int npart[6];           /*!< number of particles of each type in this file */
  double mass[6];              /*!< mass of particles of each type. If 0, then the masses are explicitly
                                 stored in the mass-block of the snapshot file, otherwise they are omitted */
  double time;                 /*!< time of snapshot file */
  double redshift;             /*!< redshift of snapshot file */
  int flag_sfr;           /*!< flags whether the simulation was including star formation */
  int flag_feedback;      /*!< flags whether feedback was included (obsolete) */
  unsigned int npartTotal[6];      /*!< total number of particles of each type in this snapshot. This can be
                                        different from npart if one is dealing with a multi-file snapshot. */
  int flag_cooling;       /*!< flags whether cooling was included  */
  int num_files;          /*!< number of files in multi-file snapshot */
  double BoxSize;              /*!< box-size of simulation in case periodic boundaries were used */
  double Omega0;               /*!< matter density in units of critical density */
  double OmegaLambda;          /*!< cosmological constant parameter */
  double HubbleParam;          /*!< Hubble parameter in units of 100 km/sec/Mpc */
  int flag_stellarage;    /*!< flags whether the file contains formation times of star particles */
  int flag_metals;        /*!< flags whether the file contains metallicity values for gas and star
                                 particles */
  unsigned int npartTotalHighWord[6];   /*!< High word of the total number of particles of each type */
  int  flag_entropy_instead_u;         /*!< flags that IC-file contains entropy instead of u */
#ifdef LT_STELLAREVOLUTION
  int flag_metalcooling;
  int flag_stellarevolution;
  char fill[52];               /*!< fills to 256 Bytes */
#else
  char fill[60];	       /*!< fills to 256 Bytes */
#endif
}
header;  /*!< holds header for snapshot files */



enum iofields
{ IO_POS,
  IO_VEL,
  IO_ID,
  IO_MASS,
  IO_U,
  IO_CR_C0,
  IO_CR_Q0,
  IO_CR_P0,
  IO_CR_E0,
  IO_CR_n0,
  IO_CR_ThermalizationTime,
  IO_CR_DissipationTime,
  IO_RHO,
  IO_NE,
  IO_NH,
  IO_ELECT,
  IO_HI,
  IO_HII,
  IO_HeI,
  IO_HeII,
  IO_HeIII,
  IO_H2I,
  IO_H2II,
  IO_HM,
  IO_iMass,
  IO_Zs,
  IO_CLDX,
  IO_CONTRIB,
  IO_HSML,
  IO_SFR,
  IO_AGE,
  IO_Z,
  IO_POT,
  IO_ACCEL,
  IO_DTENTR,
  IO_STRESSDIAG,
  IO_STRESSOFFDIAG,
  IO_STRESSBULK,
  IO_SHEARCOEFF,
  IO_TSTP,
  IO_BFLD,
  IO_DBDT,
  IO_DIVB,
  IO_ABVC,
  IO_AMDC,
  IO_PHI,
  IO_COOLRATE,
  IO_CONDRATE,
  IO_BSMTH,
  IO_DENN,
  IO_EGYPROM,
  IO_EGYCOLD,
  IO_BHMASS,
  IO_BHMDOT,
  IO_MACH,
  IO_DTENERGY,
  IO_PRESHOCK_DENSITY,
  IO_PRESHOCK_ENERGY,
  IO_PRESHOCK_XCR,
  IO_DENSITY_JUMP,
  IO_ENERGY_JUMP,
  IO_CRINJECT,
  IO_BG_METALS,
  IO_BG_INITIAL_MASS,
  IO_BG_METALLICITY,

  IO_LASTENTRY        /* This should be kept - it signals the end of the list */
};



/*
 * Variables for Tree
 * ------------------
 */

extern struct NODE
{
  FLOAT len;			/*!< sidelength of treenode */
  FLOAT center[3];		/*!< geometrical center of node */
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
  FLOAT maxsoft;                /*!< hold the maximum gravitational softening of particles in the 
                                     node if the ADAPTIVE_GRAVSOFT_FORGAS option is selected */
#endif
  union
  {
    int suns[8];		/*!< temporary pointers to daughter nodes */
    struct
    {
      FLOAT s[3];               /*!< center of mass of node */
      FLOAT mass;               /*!< mass of node */
      int bitflags;        /*!< flags certain node properties */
      int sibling;         /*!< this gives the next node in the walk in case the current node can be used */
      int nextnode;        /*!< this gives the next node in case the current node needs to be opened */
      int father;          /*!< this gives the parent node of each node (or -1 if we have the root node) */
    }
    d;
  }
  u;
}
*Nodes_base,                    /*!< points to the actual memory allocted for the nodes */
*Nodes;                         /*!< this is a pointer used to access the nodes which is shifted such that Nodes[All.MaxPart] 
				  gives the first allocated node */


extern struct extNODE
{
  FLOAT hmax;			/*!< maximum SPH smoothing length in node. Only used for gas particles */
}
 *Extnodes, *Extnodes_base;


extern int MaxNodes;		/*!< maximum allowed number of internal nodes */
extern int Numnodestree;	/*!< number of (internal) nodes in each tree */


extern int *Nextnode;	/*!< gives next node in tree walk  (nodes array) */
extern int *Father;	/*!< gives parent node in tree (Prenodes array) */

#ifdef STATICNFW
extern double Rs, R200;
extern double Dc;
extern double RhoCrit, V200;
extern double fac;
#endif


#ifdef CHEMISTRY
/* ----- chemistry part ------- */

#define H_number_fraction 0.76
#define He_number_fraction 0.06

/* ----- Tables ------- */
extern double  T[N_T],J0_nu[N_nu],J_nu[N_nu],nu[N_nu];
extern double  k1a[N_T],k2a[N_T],k3a[N_T],k4a[N_T],k5a[N_T],k6a[N_T],k7a[N_T],k8a[N_T],k9a[N_T],k10a[N_T],k11a[N_T];
extern double  k12a[N_T],k13a[N_T],k14a[N_T],k15a[N_T],k16a[N_T],k17a[N_T],k18a[N_T],k19a[N_T],k20a[N_T],k21a[N_T];
extern double  ciHIa[N_T],ciHeIa[N_T],ciHeIIa[N_T],ciHeISa[N_T],reHIIa[N_T],brema[N_T];
extern double  ceHIa[N_T],ceHeIa[N_T],ceHeIIa[N_T],reHeII1a[N_T],reHeII2a[N_T],reHeIIIa[N_T];

/* cross-sections */
#ifdef RADIATION
extern double sigma24[N_nu],sigma25[N_nu],sigma26[N_nu],sigma27[N_nu],sigma28[N_nu],sigma29[N_nu],sigma30[N_nu],sigma31[N_nu];
#endif
#endif

#endif
