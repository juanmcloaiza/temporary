#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"

//#include "lt_error_codes.h"

#ifdef LT_STELLAREVOLUTION

#define Z_ROUND_OFF_ERROR (1e-7)
#define qag_INT_KEY 4

#ifndef INLINE_FUNC
#ifdef INLINE
#define INLINE_FUNC inline
#else
#define INLINE_FUNC
#endif
#endif

/* ...........: useful quantities :...................... */

static double BoxHalf, BoxSize;
static double UnitMassFact;
static int SEvInfo_grain;

#ifdef LT_TRACK_CONTRIBUTES
Contrib contrib;
float IIcontrib[LT_NMet], Iacontrib[LT_NMet], AGBcontrib[LT_NMet];
#endif

/* ...........: integration :...................... */

/*
#define gsl_ws_dim 10000
static gsl_integration_workspace *w;

static gsl_interp_accel *accel;
static gsl_spline *spline;
*/
/* :::............................................................... */


#define TIME_INTERP_SIZE 10000
static double aarray[TIME_INTERP_SIZE], tarray[TIME_INTERP_SIZE];
static double cosmic_time;

static void *imf_pointer;

#ifdef PERIODIC
#define NGB_PERIODIC_LONG_X(x) (xtmp=fabs(x),(xtmp>boxHalf_X)?(boxSize_X-xtmp):xtmp)
#define NGB_PERIODIC_LONG_Y(x) (xtmp=fabs(x),(xtmp>boxHalf_Y)?(boxSize_Y-xtmp):xtmp)
#define NGB_PERIODIC_LONG_Z(x) (xtmp=fabs(x),(xtmp>boxHalf_Z)?(boxSize_Z-xtmp):xtmp)
#else
#define NGB_PERIODIC_LONG_X(x) fabs(x)
#define NGB_PERIODIC_LONG_Y(x) fabs(x)
#define NGB_PERIODIC_LONG_Z(x) fabs(x)
#endif

#define agefact 1000000000
#define max(x, y) (((x) > (y)) ? (x) : (y))
#define min(x, y) (((x) < (y)) ? (x) : (y))


/* :::............................................................... */

#define IaY(set,element,metallicity,mass) (SnIaYields[set][element][metallicity*IaMbins_dim[set] + mass])
#define IIY(set,element,metallicity,mass) (SnIIYields[set][element][metallicity*IIMbins_dim[set] + mass])
#define AGBY(set,element,metallicity,mass) (AGBYields[set][element][metallicity*AGBMbins_dim[set] + mass])

/* ...........: neighbours finding :...................... */

/*
   T r i c k s  for neighbours finding

   -DLT_HASH_STARS       save export flag and ngblist for the local processor for each star
   -DLT_USE_NEAREST_HSML find the hsml of the nearest gas neighbour. to be used as guess for
			  the spreading length.
   -DLT_H_GUESS          every time a star particle is grav-active, while walking the tree the number
			  of gas neighbours included in a sphere with the current spreading length
			  is estimated. this estimate is used to update smoothly the spreading
			  length.
   -DLT_NEIGHB_GUESS     searching neighbours return the ratio #gas/#stars within the SL-sphere.
			  this is used coupled to NEAREST_HSML to correctly guess the SL (the
			  estimate will be different if we are in a star-crowded or gas-crowded
			  region).
*/

#define FIND_NEIGHB      8
#define EVOLVE           16
#define SNIAflag         32
#define SNIIflag         64
#define ALL_SnFLAGS      (FIND_NEIGHB + EVOLVE + SNIAflag + SNIIflag)

#define MAX_ITER         4

/*
     As long as we use the upper part (> bit 8 so far ) of P[i].Type to store the ordinal number of
     the star, we don't need any hash function.
     NOTE: by this way we are assuming that no more than INT_MAX / 2^8 stars will be active at the
	   same time for each processor ; for 32-bits u-integer this results in ~16.7e6 stars per proc,
	   that will be surely sufficient up to quite big runs.
	   We insert in count_evolving_stars a check on the number of active stars, so that the
	   flagging is stopped if this limit is crossed. The stars excluded will evolve in the next
	   timestep.
*/

#define Gf               1.26	/* 2^(1/3) */
#define Sf               (1.0/Gf)

#define LOCAL            0
#define BUFFER           1

FLOAT *NOlds;

/* :::............................................................... */

gsl_function F;

#ifdef LT_SEvDbg
int do_spread_dbg, *do_spread_dbg_list;
double weight_sum, tot_weightsum;
double metals_spread[LT_NMet], tot_metalsspread[LT_NMet];
#endif


/* ...........: stats and log files related quantities :...................... */

#define INUM          8
#define SN_NeighFind  0
#define SN_NeighCheck 1
#define SN_Comm       2
#define SN_Calc       3
#define SN_Imbalance  4
#define SN_info       5
#define SN_Spread     6

double infos[INUM], sum_infos[INUM];
double SN_Find, sumSN_Find;

int tot_starsnum;

#ifdef LT_SEv_INFO

#define MIN_INIT_VALUE 1e6

/* Chemical Enrichment history */

#define SPECIES     4
#define SPEC_snIa   0
#define SPEC_snII   2
#define SPEC_agb    4
#define SPEC_IRA    6		/* if LT_LOAL_IRA is not defined */

#define SN_INUM     2
#define SN_num      0
#define SN_egy      1

#define AL_INUM_sum 4
#define MEAN_ngb    0
#define MEAN_sl     1
#define NUM_uspread 2
#define NUM_nspread 3

#define AL_INUM_min 2
#define MIN_ngb     0
#define MIN_sl      1

#define AL_INUM_max 2
#define MAX_ngb     0
#define MAX_sl      1

/* log the ejected mass for each metal by SnIa, SnII and nebulae */
/* in [SPEC_XX + 1] there will be the losses (in the case that losses are possible) */
double Zmass[LT_NIMFs][SPECIES * 2][LT_NMet], tot_Zmass[LT_NIMFs][SPECIES * 2][LT_NMet];

/* log number and energy released for each specie */
double SNdata[LT_NIMFs][SPECIES * 2][SN_INUM], tot_SNdata[LT_NIMFs][SPECIES * 2][SN_INUM];

/* log the number of neighbours and the spreading lenght */
/* the min,max and mean will be logged in ALdata_min _max and _sum */
double ALdata_sum[AL_INUM_sum], tot_ALdata_sum[AL_INUM_sum];
double ALdata_min[AL_INUM_min], tot_ALdata_min[AL_INUM_min];
double ALdata_max[AL_INUM_max], tot_ALdata_max[AL_INUM_max];


#define S_INUM_sum (2 * LT_NMet) + 5
#define MEAN_Zs     0		/* mean abundance for each element in gas */
#define MEAN_Z      LT_NMet	/* mean metallicity in gas */
#define MEAN_egyf   MEAN_Z + 1	/* mean (internal egy) / (sn egy) value */
#define NUM_egyf    MEAN_egyf + 1
#define MEAN_Zsstar NUM_egyf + 1	/* mean abundance for each element in stars */
#define MEAN_Zstar  MEAN_Zsstar + LT_NMet	/* mean metallicity in stars */
#define NUM_star    MEAN_Zstar + 1

#define S_INUM_min  3
#define MIN_Z       0		/* min metallicity in gas */
#define MIN_Zstar   1		/* min metallicity in stars */
#define MIN_egyf    2		/* min (internal egy) / (sn egy) value */

#define S_INUM_max  3
#define MAX_Z       0		/* max metallicity in gas */
#define MAX_Zstar   1		/* min metallicity in stars */
#define MAX_egyf    2		/* max (internal egy) / (sn egy) value */

double Stat_sum[S_INUM_sum], tot_Stat_sum[S_INUM_sum];
double Stat_min[S_INUM_min], tot_Stat_min[S_INUM_min];
double Stat_max[S_INUM_max], tot_Stat_max[S_INUM_max];

#if defined(LT_EJECTA_IN_HOTPHASE) && !defined(LT_SNegy_IN_HOTPHASE)
double SpreadEgy[2], tot_SpreadEgy[2];
double SpreadMinMaxEgy[2][2], tot_SpreadMinMaxEgy[2][2];
double AgbFrac[3], tot_AgbFrac[3];
double SpecEgyChange[3], tot_SpecEgyChange[3];
double CFracChange[4], tot_CFracChange[4];
#endif

int spreading_on;
#endif

#ifdef LT_SEv_INFO_DETAILS
struct details
{
  unsigned int ID;
  char type;
  float Data[LT_NMet + 4];	/* stores all metals + the iMass + start_time + end_time */
} *Details;

#define DetailsZ 4
int DetailsPos;

#endif

double dmax1, dmax2;

#ifdef UM_COOLING
double Fill_mol_num;
#endif

/* :::............................................................... */


/* ...........: declare internal functions :...................... */

int INLINE_FUNC sev_compare_key_ngb(const void *, const void *);
int INLINE_FUNC sev_compare_key_spread(const void *, const void *);

/*- neighbours finding -*/
#ifdef DOUBLEPRECISION
double INLINE_FUNC myfloor(double);
#else
float INLINE_FUNC myfloor(double);
#endif
int count_evolving_stars(int);
void build_SN_Stepping(int, int);
double INLINE_FUNC lifetime(double);
void neighbours_count(int, int);
int INLINE_FUNC iterate(int, int, FLOAT, double *);
double INLINE_FUNC inner_integrand(double, void *);
double INLINE_FUNC get_da_dota(double, void *);

/*- stellar evolution -*/
#if !defined(LT_EJECTA_IN_HOTPHASE) && !defined(LT_HOT_EJECTA) && !defined(LT_SNegy_IN_HOTPHASE)
void stellarevolution(int, int, float *, double *);
void spread(int, int, FLOAT *, float *, double);
#else
float stellarevolution(int, int, float *, double *);
void spread(int, int, FLOAT *, float *, float, double);
#endif

/*
void spread(int, int, float*, double);
*/

double INLINE_FUNC dying_mass(double);
double INLINE_FUNC dm_dt(double, double);
double INLINE_FUNC phi(double, void *);
double INLINE_FUNC integrate_imf(double, double, int, void *);
double INLINE_FUNC integrand0(double, void *);

void get_metallicity_stat(void);

/*- snIa -*/
double INLINE_FUNC sec_dist(double, double);
double get_SnIa_product(int, int, double *, double *, double, double);
double INLINE_FUNC nRSnIa(double, void *);
double INLINE_FUNC mRSnIa(double, void *);

/*- snII and nebulae -*/
double get_AGB_product(int, int, double *, double, double);
double get_SnII_product(int, int, double *, double *, double, double);
double INLINE_FUNC nRSnII(double, void *);
double INLINE_FUNC mRSnII(double time, void *p);
double INLINE_FUNC zmRSnII(double m, void *p);
double INLINE_FUNC ejectaSnII(double, void *);
double INLINE_FUNC ztRSnII(double, void *);

void get_Egy_and_Beta(double, double *, double *, IMF_Type *);
void calculate_ShortLiving_related(IMF_Type *, int);

/* :::............................................................... */


/* ********************************************************************************************************** */
/* ********************************************************************************************************** */


/* .........................................
   * ------------------------------------- *
   |                                       |
   |   Chemical Stepping                   |
   * ------------------------------------- *
*/


/* double INLINE_FUNC get_cost_SE(int i) */
/* { */
/*   int index; */
/*   double costIa = 0, costII = 0; */
/*   double lifetime; */

/*   index = get_IMF_index(i); */

/*   lifetime = get_age(P[i].StellarAge) - All.Time_Age; */
/* #if defined(LT_SNIa) || defined(LT_AGB) */
/*   if((costIa = LLv_Nsteps[index] - get_current_chem_bin(lifetime, 1, index)) < 0) */
/*     costIa = LLv_Nsteps[index]; */
/*   costIa = 1 - costIa / LLv_Nsteps[index]; */
/* #endif */
/* #ifdef LT_SNII */
/*   if((costII = SnII_Nsteps[index] - get_current_chem_bin(lifetime, 2, index)) < 0) */
/*     costII = SnII_Nsteps[index]; */
/*   costII = 1 - costII / SnII_Nsteps[index]; */
/* #endif */

/*   return (costIa + costII) / 2; */
/* } */


double INLINE_FUNC get_timing(int type, double lifetime, int *bin, int IMFi)
     /* returns how much Gyr are needed to reach the end
        of the chemical timestep you live in, accordingly
        to the chemical timesteps calculated in
        build_SN_Stepping; *bin will contain the ordinal
        number of the array's timestep in which the
        star currently lives.
      */
{
  double *timings;
  int Ntimings;


  if(type == 2)
    {
      timings = &SnII_steps[IMFi][0][0];
      Ntimings = SnII_Nsteps[IMFi];
      if(Ntimings == 0)
	return 0;

      if(lifetime > timings[Ntimings - 1])
	lifetime = timings[Ntimings - 1];
      /* a rough initial estimate */
      *bin = lifetime / All.mean_lifetime * (Ntimings - 1);
    }
#if defined(LT_SNIa) || defined(LT_AGB)
  else if(type == 1 && lifetime < IMFs[IMFi].sup_lifetime)
    {
      timings = &LLv_steps[IMFi][0][0];
      Ntimings = LLv_Nsteps[IMFi];
      if(Ntimings == 0)
	return 0;

      if(lifetime > timings[Ntimings - 1])
	lifetime = timings[Ntimings - 1];
      /* a rough initial estimate */
      *bin = lifetime / IMFs[IMFi].sup_lifetime * (Ntimings - 1);
    }
#endif
  else if(lifetime == All.mean_lifetime && type == 2)
    {
      *bin = SnII_Nsteps[IMFi] - 1;
      return 0;
    }
#if defined(LT_SNIa) || defined(LT_AGB)
  else if(lifetime == IMFs[IMFi].sup_lifetime)
    {
      *bin = LLv_Nsteps[IMFi] - 1;
      return 0;
    }
#endif
  /* arrays are small enough to avoid binary search; it is
   * very likely that the overhead for routine calling and
   * variables' set-up takes as much time as getting through
   * the whole array. */

  if(timings[*bin] >= lifetime)
    while(timings[*bin] > lifetime)
      (*bin)--;
  else
    while(timings[*bin + 1] <= lifetime)
      (*bin)++;

  return (timings[*bin + 1] - lifetime);
}

int INLINE_FUNC convert_time_to_timesteps(double start_a, double start_time, double delta_time)
{
  double delta_a;

  /* get delta_expansion_factor when moving from start_time to start_time + delta_time   *
   * start_a is the expansion factor corresponding to start_time                         */
  delta_a = gsl_spline_eval(spline, cosmic_time - start_time + delta_time * 1.01, accel) - start_a;

  /* converts in code steps */
  return (int) (log(delta_a / start_a + 1) / All.Timebase_interval);
}

int INLINE_FUNC get_chemstep(int type, double lifetime, double atime, int base_start, int base_end, int IMFi)
     /* returns the number of base intervals needed to reach *
      * the next evolution time.
      * lifetime is the age of the Star in Gyr, atime is the
      * expansion factor from which to calculate the step (usually
      * All.Time). */
{
  double *timings, delta;
  int Ntimings, pos, steps;

  if((lifetime == All.mean_lifetime && type == 2) || (lifetime == IMFs[IMFi].sup_lifetime && type == 1))
    return 0;

  if(type == 2)
    {
      timings = &SnII_steps[IMFi][0][0];
      Ntimings = SnII_Nsteps[IMFi];
      if(Ntimings == 0)
	return 0;
      if(lifetime > timings[Ntimings - 1])
	lifetime = timings[Ntimings - 1];
      /* a rough initial estimate */
      pos = lifetime / All.mean_lifetime * (Ntimings - 1);
    }
#if defined(LT_SNIa) || defined(LT_AGB)
  else
    {
      timings = &LLv_steps[IMFi][0][0];
      Ntimings = LLv_Nsteps[IMFi];
      if(Ntimings == 0)
	return 0;
      /* a rough initial estimate */
      if(lifetime > timings[Ntimings - 1])
	lifetime = timings[Ntimings - 1];
      pos = lifetime / IMFs[IMFi].sup_lifetime * (Ntimings - 1);
    }
#endif

  if(timings[pos] >= lifetime)
    while(timings[pos] > lifetime)
      pos--;
  else
    while(timings[pos + 1] <= lifetime)
      pos++;

  delta = timings[pos + 1] - lifetime;
  steps = convert_time_to_timesteps(atime, All.Time_Age, delta);

  while((base_start + steps < base_end) && (base_start + steps <= TIMEBASE) && (++pos < Ntimings - 1))
    {
      delta = timings[pos + 1] - lifetime;
      steps = convert_time_to_timesteps(atime, All.Time_Age, delta);
    }

  return steps;
}


double INLINE_FUNC get_da_dota(double y, void *param)
     /* gives da/(da/dt) */
{
  return 1 / (y * sqrt(All.Omega0 * y * y * y + All.OmegaLambda +
		       (1 - All.Omega0 - All.OmegaLambda) * y * y));
}

double INLINE_FUNC get_age(double a)
     /* gives the look-back time corresponding to the expansion factor a */
{
  static double sec_in_Gyr = (86400.0 * 365.0 * 1e9);
  double result, error;

  F.function = &get_da_dota;
  F.params = NULL;

  if((my_gslstatus =
      gsl_integration_qag(&F, 1, 1 / a, 1e-4, 1e-5, gsl_ws_dim, qag_INT_KEY, w, &result, &error)))
    {
      printf(">>>>> [%3d] qag integration error %d in get_age\n", ThisTask, my_gslstatus);
      endrun(LT_ERROR_INTEGRATION_ERROR);
    }

  return 1.0 / (HUBBLE * All.HubbleParam) * result / sec_in_Gyr;
}


/* ********************************************************************************************************** */
/* ********************************************************************************************************** */

#ifdef DOUBLEPRECISION
double INLINE_FUNC myfloor(double v)
{
  return floor(v);
}
#else
float INLINE_FUNC myfloor(double v)
{
  return floorf((float) v);
}
#endif


/* ...........................................................
   * ------------------------------------------------------- *
   |                                                         |
   |   Searching for Evolving Stars                          |
   * ------------------------------------------------------- *
*/

int count_evolving_stars(int mode)
{
  int I;
  int flag;
  long int i, num_of_stars, num_tobe_done;

  num_tobe_done = 0;

  for(i = num_of_stars = 0; i < NumPart; i++)
    {
      if(P[i].Type == 4)
	{
	  I = P[i].MetID;
	  flag = 0;

#if defined(LT_SNIa) || defined(LT_AGB)
	  flag += SNIAflag * (MetP[I].NextChemStepIa <= All.Ti_Current);
#endif
#ifdef LT_SNII
	  flag += SNIIflag * (MetP[I].NextChemStepII <= All.Ti_Current);
#endif

	  if(flag)
	    {
	      num_of_stars++;
	      P[i].Type |= (flag |= FIND_NEIGHB);
	      MetP[I].Left = MetP[I].Right = 0;
	    }
	}
#ifndef LT_LOCAL_IRA
      else if((P[i].Type == 0) && (SphP[i].mstar > 0))
	{
	  num_of_stars++;
	  P[i].Type |= EVOLVE;	/* no need here to find neighbours for gas particles ! */
	}
#endif
    }

  return num_of_stars;
}


/* ********************************************************************************************************** */
/* ********************************************************************************************************** */



/*
   =======================================================
     M A I N   R O U T I N E
     .........................................
     * ------------------------------------- *
     |                                       |
     |   Evolving Stars                      |
     * ------------------------------------- *
   =======================================================
*/

static int *noffset, *nbuffer, *nsend, *nsend_local;
static int *ndone, *ntotdone, *noff, *ntotoff;
int ndoneoff[2], ntotdoneoff[2];

int evolve_SN(int mode)
{
  /*- general variables -*/
  unsigned int I;
  int starsnum, ntotleft, tot_tobe_done;
  int i, j, k, source;		/* general useful indexes */
  unsigned int place;
  int iter, maxiter;
  int ncount;
  int nexport, ntotexport;
  int ndoneoff[2], ntotdoneoff[2];
  int sendTask, recvTask, level, ngrp, maxfill;
  double Factor, FactorExponent;
  MPI_Status status;

  /*- bufferize physics -*/
  float Metals[LT_NMet];

#if defined(LT_EJECTA_IN_HOTPHASE) || defined(LT_HOT_EJECTA) || defined(LT_SNegy_IN_HOTPHASE)
  float LMMass;
#endif
  double energy;

  /*- infos -*/
  double tstart, tend;

  int few, many, bounce, min_sl, unres, tot_few, tot_many, tot_bounce, tot_minsl, tot_unres;
  int Lguess, NNgb;

  int Yset, YZbin;
  double Z;
  float int_part, frac_part;


  /* debug useful variables */
  /* step, task, buffer position, private position */
  /*
     #define dStep 21067
     #define dTask 13
     #define dSL 12.172611
     #define dPart_in_buff 38
   */
  /* - */

  /* . -------------------------------- .
     .                                  .
     . are there evolving stars ?       .
     . ________________________________ .
   */

  tstart = second();

  starsnum = count_evolving_stars(mode);
  tot_tobe_done = starsnum;

  tend = second();
  SN_Find = timediff(tstart, tend);
  MPI_Reduce(&SN_Find, &sumSN_Find, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  All.CPU_SN_find += sumSN_Find / NTask;

  MPI_Allreduce(&starsnum, &tot_starsnum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(tot_starsnum == 0)
    {
#if defined(LT_SEv_INFO) && defined(LT_SEvDbg)
      if(SEvInfo_GRAIN == 0)
	get_metals_sumcheck(9);
#endif
      return 0;
    }

  for(i = 0; i < INUM; i++)
    infos[i] = 0;

#ifdef LT_SEv_INFO_DETAILS
  if(ThisTask == 0)
    Details = (struct details *) malloc(tot_starsnum * 3 * sizeof(struct details));
  else
    Details = (struct details *) malloc(starsnum * 3 * sizeof(struct details));
  DetailsPos = 0;
#endif


  /* . ---------------------------------- .
     .                                    .
     . ok, let's initialize some vars     .
     . __________________________________ .
   */

  BoxSize = All.BoxSize;
  BoxHalf = 0.5 * All.BoxSize;


#ifdef MAX_CHEM_SPREAD_LENGTH
  /* assume the length will grow by a 26% each iteration - the case of non-guessed growth */
  maxiter = (int) (ceil(log(All.MaxChemSpreadL) / log(Gf)) * 1.5);
#else
  maxiter = MAX_ITER;
#endif

  /* . ---------------------------------- .
     .                                    .
     . ok, let's initialize info vars     .
     . __________________________________ .
   */

#ifdef LT_SEv_INFO

  for(k = 0; k < LT_NIMFs; k++)
    for(j = 0; j < SPECIES * 2; j++)
      for(i = 0; i < LT_NMet; i++)
	{
	  Zmass[k][j][i] = tot_Zmass[k][j][i] = 0;
	  SNdata[k][j][i] = tot_SNdata[k][j][i] = 0;
	}

#if defined(LT_EJECTA_IN_HOTPHASE) && !defined(LT_SNegy_IN_HOTPHASE)
  for(i = 0; i < 2; i++)
    {
      SpreadEgy[i] = tot_SpreadEgy[i] = 0;
      SpreadMinMaxEgy[i][0] = SpreadMinMaxEgy[i][1] = 0;
      tot_SpreadMinMaxEgy[i][0] = tot_SpreadMinMaxEgy[i][1] = 0;
    }
  for(i = 0; i < 3; i++)
    {
      AgbFrac[i] = 0;
      SpecEgyChange[i] = 0;
      CFracChange[i] = 0;
    }
  SpreadMinMaxEgy[0][0] = MIN_INIT_VALUE;
  AgbFrac[1] = MIN_INIT_VALUE;
  SpecEgyChange[1] = MIN_INIT_VALUE;
  CFracChange[2] = MIN_INIT_VALUE;
  CFracChange[3] = 0;
#endif

  if(++SEvInfo_grain > SEvInfo_GRAIN)
    {
      for(i = 0; i < AL_INUM_sum; i++)
	ALdata_sum[i] = 0;
      for(i = 0; i < AL_INUM_max; i++)
	ALdata_max[i] = 0;
      for(i = 0; i < AL_INUM_min; i++)
	ALdata_min[i] = MIN_INIT_VALUE;

      for(i = 0; i < S_INUM_sum; i++)
	Stat_sum[i] = 0;
      for(i = 0; i < S_INUM_max; i++)
	Stat_max[i] = 0;
      for(i = 0; i < S_INUM_min; i++)
	Stat_min[i] = MIN_INIT_VALUE;
    }
#endif

  /* ******************************************* ****** *** **  * */



  /* . -------------------------------- .
     .                                  .
     . NEIGHBOURS FINDING               .
     . ________________________________ .
   */


  iter = 1;
  ntotleft = tot_starsnum;

  search_for_metalspread = 1;

  /* MAIN (outer) cycle                 */
  /* .................................. */
  while(ntotleft > 0 && iter <= maxiter)
    {

      for(j = 0; j < NTask; j++)
	nsend_local[j] = 0;
      few = many = bounce = min_sl = unres = 0;
      tot_few = tot_many = tot_bounce = tot_minsl = tot_unres = 0;
      tstart = second();

      /* LOCAL cycle */
      for(nexport = ntotexport = i = ncount = 0; i < NumPart && ncount < All.BunchSizeMetal_Ngb - NTask; i++)
	if(P[i].Type & FIND_NEIGHB)
	  {
	    I = P[i].MetID;
	    ncount++;
	    for(j = 0; j < NTask; j++)
	      Exportflag[j] = 0;

	    /* check for local neighbours
	       check for Exportflag entries
	       set .NumNgb and .weight fields of the star particle i
	     */
	    neighbours_count(i, LOCAL);


	    for(j = 0; j < NTask; j++)
	      if(Exportflag[j] > 0)
		{
		  MetalNgbIn[nexport].Pos[0] = P[i].Pos[0];
		  MetalNgbIn[nexport].Pos[1] = P[i].Pos[1];
		  MetalNgbIn[nexport].Pos[2] = P[i].Pos[2];
		  MetalNgbIn[nexport].Task = j;
		  MetalNgbIn[nexport].index = i;
		  MetalNgbIn[nexport].L = MetP[I].SLguess;

		  nexport++;
		  nsend_local[j]++;
		}
	  }
      /* END of LOCAL cycle */
      tend = second();
      /* sum-up time */
      infos[SN_NeighFind] += timediff(tstart, tend);

      qsort(MetalNgbIn, nexport, sizeof(struct metal_ngbfindingdata_in), sev_compare_key_ngb);

      for(j = 1, noffset[0] = 0; j < NTask; j++)
	noffset[j] = noffset[j - 1] + nsend_local[j - 1];

      tstart = second();
      MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);
      MPI_Allreduce(&nexport, &ntotexport, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      tend = second();
      infos[SN_Imbalance] += timediff(tstart, tend);

      /* now do the particles that need to be exported */
      /* NON-LOCAL cycle */
      if(ntotexport > 0)
	for(level = 1; level < (1 << PTask); level++)
	  {

	    tstart = second();
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
		if(maxfill >= All.BunchSizeMetal_Ngb)
		  break;

		sendTask = ThisTask;
		recvTask = ThisTask ^ ngrp;

		if(recvTask < NTask)
		  {
		    if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		      {
			/* get the particles */
			MPI_Sendrecv(&MetalNgbIn[noffset[recvTask]],
				     nsend_local[recvTask] * sizeof(struct metal_ngbfindingdata_in), MPI_BYTE,
				     recvTask, TAG_SE,
				     &MetalNgbGet[nbuffer[ThisTask]],
				     nsend[recvTask * NTask +
					   ThisTask] * sizeof(struct metal_ngbfindingdata_in), MPI_BYTE,
				     recvTask, TAG_SE, MPI_COMM_WORLD, &status);
		      }
		  }

		for(j = 0; j < NTask; j++)
		  if((j ^ ngrp) < NTask)
		    nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];

	      }
	    tend = second();
	    infos[SN_Comm] += timediff(tstart, tend);

	    tstart = second();

	    /* find the number of neighbours */
	    for(j = 0; j < nbuffer[ThisTask]; j++)
	      neighbours_count(j, BUFFER);

	    tend = second();
	    infos[SN_NeighFind] += timediff(tstart, tend);

	    /* get the number of neighbours */
	    tstart = second();
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
		if(maxfill >= All.BunchSizeMetal_Ngb)
		  break;

		sendTask = ThisTask;
		recvTask = ThisTask ^ ngrp;

		if(recvTask < NTask)
		  {
		    if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		      {
			/* send results */
			MPI_Sendrecv(&MetalNgbResultIn[nbuffer[ThisTask]],
				     nsend[recvTask * NTask +
					   ThisTask] * sizeof(struct metal_ngbfindingdata_out), MPI_BYTE,
				     recvTask, TAG_SE, &MetalNgbResultGet[noffset[recvTask]],
				     nsend_local[recvTask] * sizeof(struct metal_ngbfindingdata_out),
				     MPI_BYTE, recvTask, TAG_SE, MPI_COMM_WORLD, &status);

			/* upload results */
			for(j = 0; j < nsend_local[recvTask]; j++)
			  {
			    source = j + noffset[recvTask];
			    place = MetalNgbIn[source].index;
			    I = P[place].MetID;
			    MetP[I].NumNgb += MetalNgbResultGet[source].num_ngb;
			    MetP[I].weight += MetalNgbResultGet[source].weight;
			  }
		      }
		  }

		for(j = 0; j < NTask; j++)
		  if((j ^ ngrp) < NTask)
		    nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	      }
	    tend = second();
	    infos[SN_Comm] += timediff(tstart, tend);

	    level = ngrp - 1;
	  }
      /* END of NON-LOCAL cycle */


      /* ok, let's check wheter each particle has the right number of neighbs */
      /* NGB-NUMCHECK cycle */

      tstart = second();

      for(*ndone = *noff = ncount = nexport = i = 0; i < NumPart && ncount < All.BunchSizeMetal_Ngb - NTask;
	  i++)
	if(P[i].Type & FIND_NEIGHB)
	  {
	    I = P[i].MetID;

	    ncount++;
	    if((MetP[I].NumNgb >= All.LeftNumNgbSN && MetP[I].NumNgb <= All.RightNumNgbSN) ||
	       MetP[I].SLguess == All.MinChemSpreadL)
	      /* ok, the NumNgb is right */
	      {
		*ndone++;
		P[i].Type &= ~FIND_NEIGHB;
		P[i].Type |= EVOLVE;
	      }
	    else
	      /* NumNgb is either too large or too small */
	      {
		if(MetP[I].Left > 0 && MetP[I].Right > 0)
		  if(((MetP[I].Right - MetP[I].Left) < 0.001 * MetP[I].Left))
		    /* no way to further modify the spreading length. take it as it comes! */
		    {
		      *ndone++;
		      unres++;
		      P[i].Type &= ~FIND_NEIGHB;

		      if(MetP[I].NumNgb > All.NeighInfNum)
			P[i].Type |= EVOLVE;
		      else
			{
			  *noff++;
			  P[i].Type = 4;
			}
		      continue;
		    }

		if(iter == maxiter)
		  /* oops, anyway too much iterarions done */
		  {
		    P[i].Type &= ~FIND_NEIGHB;

		    if(MetP[I].NumNgb < All.LeftNumNgbSN)
		      few++;
		    else
		      many++;
		    if(MetP[I].Right > 0 && MetP[I].Left > 0)
		      bounce++;

		    *ndone++;
		    *noff++;
		    P[i].Type = 4;

		    continue;
		  }

		/* set the right and left boundaries */
		if(MetP[I].NumNgb < All.LeftNumNgbSN)
		  {
		    MetP[I].Left = max(MetP[I].SLguess, MetP[I].Left);
		    few++;
		  }
		else
		  {
		    many++;
		    if(MetP[I].Right != 0)
		      {
			if(MetP[I].SLguess < MetP[I].Right)
			  MetP[I].Right = MetP[I].SLguess;
		      }
		    else
		      MetP[I].Right = MetP[I].SLguess;
		  }

		if(MetP[I].Right > 0 && MetP[I].Left > 0)
		  /* SL is bouncing accross the permitted interval: set it to the average of the
		   *    last two values
		   */
		  {
		    MetP[I].SLguess = pow(0.5 * (pow(MetP[I].Left, 3) + pow(MetP[I].Right, 3)), 0.33333333);
		    bounce++;
		  }
		else
		  {
		    if(MetP[I].Right == 0 && MetP[I].Left == 0)
		      {
			printf(" !!! Very Bad Error at ThisTask = %i for particle.ID %i.%ui [%ui]\n%s",
			       ThisTask, i, P[i].ID,
			       P[MetP[I].PID].MetID - I,
			       ((MetP[I].SLguess ==
				 0) ? "the reason is likely to be that this particle got a 0 SL\n" : '\0'));
			fflush(stdout);
			endrun(777);
		      }

		    /* calculate the multiplying factor to estimate the best hint for SL */
		    FactorExponent = 1.0 / 3;
		    Factor = pow((double) All.DesNumNgbSN / MetP[I].NumNgb, FactorExponent);

		    /* grow */
		    if(MetP[I].Right == 0 && MetP[I].Left > 0)
		      {
			if(MetP[I].NumNgb == 0)
			  Factor = 2;

			MetP[I].SLguess *= Factor;
		      }
		    /* shrink */
		    else
		      MetP[I].SLguess *= Factor;

		  }

		if(MetP[I].SLguess < All.MinChemSpreadL)
		  /* gulp, SL has got a value below the allowed minimum */
		  {
		    P[i].Type &= ~FIND_NEIGHB;

		    min_sl++;
		    *ndone++;

		    if(MetP[I].NumNgb > All.NeighInfNum)
		      P[i].Type |= EVOLVE;
		    else
		      {
			*noff++;
			P[i].Type = 4;
		      }

		    MetP[I].SLguess = All.MinChemSpreadL;
		  }

	      }
	  }
      /* END of NGB-NUMCHECK cycle */
      tend = second();
      infos[SN_NeighCheck] += timediff(tstart, tend);

      tstart = second();
      MPI_Allreduce(&ndoneoff[0], &ntotdoneoff[0], 2, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      tend = second();
      infos[SN_Imbalance] += timediff(tstart, tend);

      ntotleft -= *ntotdone;
      tot_starsnum -= *ntotoff;
      if(iter == maxiter)
	{
	  MPI_Reduce(&few, &tot_few, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&many, &tot_many, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&bounce, &tot_bounce, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&min_sl, &tot_minsl, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&unres, &tot_unres, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	}
      iter++;
    }

  /* END of MAIN NEIGHBOURS cycle    */
  /* ................................ */

  fflush(stdout);
  if(ThisTask == 0)
    {
      printf("%i iteration%c done to ensure neighbours to evolve %i stars (tot stars: %lli)\n ",
	     iter, (iter > 1) ? 's' : '\0', tot_starsnum, Ntype[4]);


      if(iter > maxiter)
	printf("  [%i tf; %i tm; %i b; %i msl; %i u]\n ",
	       tot_few, tot_many, tot_bounce, tot_minsl, tot_unres);

      fflush(stdout);
    }


  /* ******************************************** ****** *** **  * */



  /* . ----------------------------- .
     .                               .
     . SPREADING                     .
     . _____________________________ .
   */



  /* now, all the evolving stars have the right spreading lenght and
   *  have stored the total weight of neighbours.
   *
   * the nex step is to communicate the stellar evolution details and
   *  the total weight in order to spread over the neighbours
   */

  i = 0;
  ntotleft = tot_starsnum;

  /* MAIN (outer) cycle                 */
  /* .................................. */

  while(ntotleft > 0)
    {

      for(j = 0; j < NTask; j++)
	nsend_local[j] = 0;

      /* LOCAL cycle */
      for(*ndone = nexport = 0; i < NumPart && ncount < All.BunchSizeMetal_Spread - NTask; i++)
	if(P[i].Type & EVOLVE)
	  {
	    ncount++;
	    *ndone++;

	    if(P[i].Type & 4)
	      I = P[i].MetID;

	    /* deal with those parts which did not reach the desired num of neighb */

#ifdef LT_SEv_INFO

	    spreading_on = 1;

	    /* collect some information */
	    if(SEvInfo_grain > SEvInfo_GRAIN && (P[i].Type & 4))
	      {
		tstart = second();

		/* ..about the minimum, maximum and mean spreading lenght */
		if(MetP[I].SLguess < ALdata_min[MIN_sl])
		  ALdata_min[MIN_sl] = MetP[I].SLguess;
		if(MetP[I].SLguess > ALdata_max[MAX_sl])
		  ALdata_max[MAX_sl] = MetP[I].SLguess;
		ALdata_sum[MEAN_sl] += MetP[I].SLguess;

		/* ..about the minimum, maximum and mean number of neighbours and the associate spreading legnth */
		if((MetP[I].NumNgb >= All.NeighInfNum) && (MetP[I].NumNgb < (int) ALdata_min[MIN_ngb]))
		  ALdata_min[MIN_ngb] = MetP[I].NumNgb;
		if(MetP[I].NumNgb < All.LeftNumNgbSN)
		  ALdata_sum[NUM_uspread]++;
		if(MetP[I].NumNgb > (int) ALdata_max[MAX_ngb])
		  ALdata_max[MAX_ngb] = MetP[I].NumNgb;
		ALdata_sum[MEAN_ngb] += MetP[I].NumNgb;

		tend = second();
		infos[SN_info] += timediff(tstart, tend);
	      }
#endif

	    for(j = 0; j < LT_NMet; j++)
	      Metals[j] = 0;
	    energy = 0;
#if defined(LT_EJECTA_IN_HOTPHASE) || defined(LT_HOT_EJECTA) || defined(LT_SNegy_IN_HOTPHASE)
	    LMMass = 0;
#endif
	    k = get_IMF_index(i);
	    IMFp = &IMFs[k];

	    /* if this is an active star, calculate its evolution and update fields */
	    tstart = second();
	    if(P[i].Type > 4 + EVOLVE)
	      {
#if !defined(LT_EJECTA_IN_HOTPHASE) && !defined(LT_HOT_EJECTA) && !defined(LT_SNegy_IN_HOTPHASE)
		stellarevolution(i, mode, &Metals[0], &energy);
#else
		LMMass = stellarevolution(i, mode, &Metals[0], &energy);
#endif
	      }
#ifndef LT_LOCAL_IRA
	    else
	      {
		if(IMFp->nonZeroIRA == 0)
		  {
		    P[i].Type &= ~EVOLVE;
		    continue;
		  }
		Yset = IMFp->YSet;
		Z = get_metallicity_solarunits(get_metallicity(i, 2));
		for(YZbin = IIZbins_dim[Yset] - 1; Z < IIZbins[Yset][YZbin] && YZbin > 0; YZbin--)
		  ;
		for(Factor = 0, j = 0; j < LT_NMetP; j++)
		  {
		    Factor += (Metals[j] = SphP[i].mstar * SnII_ShortLiv_Yields[Yset][j][YZbin]);
#ifdef LT_SEv_INFO
		    if(spreading_on && j < LT_NMet)
		      Zmass[k][SPEC_IRA][j] += Metals[j];
		    else
		      Zmass[k][SPEC_IRA + 1][j] += Metals[j];
#endif
#ifdef LT_SEv_INFO_DETAILS
		    DetailsWo[j] += Metals[j];
#endif
		  }
#ifdef UM_COOLING
		Fill_mol_num = SnII_ShortLiv_FillMolNum[YSet][YZbin];
#endif

#ifdef LT_SEv_INFO
		if(spreading_on)
		  {
		    SNdata[k][SPEC_IRA][SN_num] += SphP[i].mstar;
		    SNdata[k][SPEC_IRA][SN_num] += IMFp->IRA_erg_per_g * SphP[i].mstar;
		  }
		else
		  {
		    SNdata[k][SPEC_IRA + 1][SN_num] += SphP[i].mstar;
		    SNdata[k][SPEC_IRA][SN_num] += IMFp->IRA_erg_per_g * SphP[i].mstar;
		  }
#endif
#ifndef LT_SNegy_IN_HOTPHASE
		energy = 0;	/* already used in effective model */
#endif
#ifndef LT_EJECTA_IN_HOTPHASE
#ifndef LT_AVOID_SELFENRICHMENT
		P[i].Mass -= (1 - (2.546479089470 * 4.188790204786) / PPP[i].a1.dNumNgb) * Factor;
#else
		P[i].Mass -= Factor;
#endif
#endif

#if defined(LT_EJECTA_IN_HOTPHASE) || defined(LT_HOT_EJECTA) || defined(LT_SNegy_IN_HOTPHASE)
		LMMass = 0;
#endif

#ifdef LT_TRACK_CONTRIBUTES
		for(j = 0; j < LT_NMet; j++)
		  {
		    IIcontrib[j] = 1;
		    Iacontrib[j] = AGBcontrib[j] = 0;
		  }
#endif
	      }
#endif

#ifdef LT_TRACK_CONTRIBUTES
	    pack_contrib(&contrib, k, IIcontrib, Iacontrib, AGBcontrib);
#endif

	    tend = second();
	    infos[SN_Calc] += timediff(tstart, tend);

	    tstart = second();

	    for(j = 0; j < NTask; j++)
	      Exportflag[j] = 0;
	    if(P[i].Type & 4)
	      ngb_treefind_flagexport(&P[i].Pos[0], MetP[I].SLguess);
	    else
	      ngb_treefind_flagexport(&P[i].Pos[0], PPP[i].Hsml);

#ifdef LT_SEvDbg
	    if(FirstID > 0 && P[i].ID == FirstID)
	      {
		do_spread_dbg = 1;
		printf("[SEvDBG] spread: spreading particle on %i neighbours :: ", MetP[I].NumNgb);
		for(j = 0; j < LT_NMet; j++)
		  printf("%g ", Metals[j]);
		printf("\n");
		fflush(stdout);
	      }
#endif
	    for(j = 0; j < NTask; j++)

	      if(Exportflag[j] > 0)
		{
		  MetalDataIn[nexport].Type = P[i].Type & 4;
		  MetalDataIn[nexport].Pos[0] = P[i].Pos[0];
		  MetalDataIn[nexport].Pos[1] = P[i].Pos[1];
		  MetalDataIn[nexport].Pos[2] = P[i].Pos[2];
		  MetalDataIn[nexport].Task = j;
		  MetalDataIn[nexport].IMFi = k;
		  if(P[i].Type & 4)
		    MetalDataIn[nexport].L = MetP[I].SLguess;
		  else
		    MetalDataIn[nexport].L = PPP[i].Hsml;
		  if(P[i].Type & 4)
		    MetalDataIn[nexport].weight = MetP[I].weight;
		  else
		    {
#ifndef NOFIXEDMASSINKERNEL
		      MetalDataIn[nexport].weight = PPP[i].a1.dNumNgb;
#else
		      MetalDataIn[nexport].weight = SphP[i].weight;
#endif
		    }
#ifdef UM_COOLING
		  MetalDataIn[nexport].FillMolNum = Fill_mol_num;
#endif
		  if(MetalDataIn[nexport].weight == 0)
		    {
		      printf
			("Task %i > something strange distributing metals:\n particle %i (type %i / %i) has 0 weight with %i neighbours\n",
			 ThisTask, i, P[i].Type, P[i].Type & 4,
			 (P[i].Type & 4 > 0) ? MetP[I].NumNgb : PPP[i].a1.dNumNgb);
		      fflush(stdout);
		      endrun(666);
		    }

		  MetalDataIn[nexport].energy = energy;
		  for(k = 0; k < LT_NMet; k++)
		    MetalDataIn[nexport].Metals[k] = Metals[k];
#if defined(LT_EJECTA_IN_HOTPHASE) || defined(LT_HOT_EJECTA) || defined(LT_SNegy_IN_HOTPHASE)
		  MetalDataIn[nexport].LMMass = LMMass;
#endif


#ifdef LT_SEvDbg
		  if(FirstID > 0 && P[i].ID == FirstID)
		    MetalDataIn[nexport].ID = P[i].ID;
		  else
		    MetalDataIn[nexport].ID = 0;
#endif
#ifdef LT_TRACK_CONTRIBUTES
		  MetalDataIn[nexport].contrib = contrib;
#endif

		  nexport++;
		  nsend_local[j]++;
		}

#if !defined(LT_EJECTA_IN_HOTPHASE) && !defined(LT_HOT_EJECTA) && !defined(LT_SNegy_IN_HOTPHASE)
	    spread(i, LOCAL, &P[i].Pos[0], &Metals[0], energy);
#else
	    spread(i, LOCAL, &P[i].Pos[0], &Metals[0], LMMass, energy);
#endif

	    tend = second();
	    infos[SN_Spread] += timediff(tstart, tend);

	    P[i].Type &= ~EVOLVE;
	  }
      /* END of LOCAL cycle */

      qsort(MetalDataIn, nexport, sizeof(struct metaldata_in), sev_compare_key_spread);

      for(j = 1, noffset[0] = 0; j < NTask; j++)
	noffset[j] = noffset[j - 1] + nsend_local[j - 1];

      tstart = second();
      MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);
      tend = second();
      infos[SN_Imbalance] += timediff(tstart, tend);

#ifdef LT_SEvDbg
      MPI_Allgather(&do_spread_dbg, 1, MPI_INT, &do_spread_dbg_list[0], 1, MPI_INT, MPI_COMM_WORLD);
      for(j = 0; j < NTask; j++)
	if((do_spread_dbg = do_spread_dbg_list[j]))
	  break;
#endif

      /* now do the particles that need to be exported */

      /* NON-LOCAL cycle */
      for(level = 1; level < (1 << PTask); level++)
	{

	  tstart = second();
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
	      if(maxfill >= All.BunchSizeMetal_Spread)
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
				   recvTask, TAG_SE,
				   &MetalDataGet[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct metaldata_in),
				   MPI_BYTE, recvTask, TAG_SE, MPI_COMM_WORLD, &status);
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];

	    }
	  tend = second();
	  infos[SN_Comm] += timediff(tstart, tend);

	  tstart = second();

	  /* do the spreading */
	  for(j = 0; j < nbuffer[ThisTask]; j++)
#if !defined(LT_EJECTA_IN_HOTPHASE) && !defined(LT_HOT_EJECTA) && !defined(LT_SNegy_IN_HOTPHASE)
	    spread(j, BUFFER, &MetalDataGet[j].Pos[0], &MetalDataGet[j].Metals[0], MetalDataGet[j].energy);
#else
	    spread(j, BUFFER, &MetalDataGet[j].Pos[0], &MetalDataGet[j].Metals[0], MetalDataGet[j].LMMass,
		   MetalDataGet[j].energy);
#endif

	  tend = second();
	  infos[SN_Spread] += timediff(tstart, tend);

	  level = ngrp - 1;
	}

      tstart = second();
      MPI_Allreduce(&ndoneoff[0], &ntotdoneoff[0], 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      tend = second();
      infos[SN_Imbalance] += timediff(tstart, tend);
      ntotleft -= *ntotdone;
    }
  /* END of NON-LOCAL cycle */

  /* ******************************************** ****** *** **  * */

#ifdef LT_SEvDbg
  if(do_spread_dbg)
    {
      MPI_Reduce(&weight_sum, &tot_weightsum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      if(ThisTask == 0)
	{
	  printf("[SEvDBG] spread sum weight: %f\n", tot_weightsum);
	  fflush(stdout);
	}
      do_spread_dbg = 0;
      weight_sum = tot_weightsum = 0;
      for(j = 0; j < NTask; j++)
	do_spread_dbg_list[j] = 0;
    }
#endif


  /* free memory */
  myfree(nsend);
  myfree(nsend_local);
  myfree(nbuffer);
  myfree(noffset);

#ifdef LT_SEv_INFO

  tstart = second();

  MPI_Reduce(&SNdata[0][0][0], &tot_SNdata[0][0][0], SPECIES * 2 * SN_INUM * LT_NIMFs, MPI_DOUBLE, MPI_SUM, 0,
	     MPI_COMM_WORLD);

  MPI_Reduce(&Zmass[0][0][0], &tot_Zmass[0][0][0], SPECIES * 2 * LT_NMet * LT_NIMFs,
	     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

#if defined(LT_EJECTA_IN_HOTPHASE) && !defined(LT_SNegy_IN_HOTPHASE)
  MPI_Reduce(&SpreadEgy[0], &tot_SpreadEgy[0], 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&SpreadMinMaxEgy[0][0], &tot_SpreadMinMaxEgy[0][0], 1, MPI_2DOUBLE_PRECISION,
	     MPI_MINLOC, 0, MPI_COMM_WORLD);
  MPI_Reduce(&SpreadMinMaxEgy[1][0], &tot_SpreadMinMaxEgy[1][0], 1, MPI_2DOUBLE_PRECISION,
	     MPI_MAXLOC, 0, MPI_COMM_WORLD);
  MPI_Reduce(&AgbFrac[0], &tot_AgbFrac[0], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&AgbFrac[1], &tot_AgbFrac[1], 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&AgbFrac[2], &tot_AgbFrac[2], 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  MPI_Reduce(&SpecEgyChange[0], &tot_SpecEgyChange[0], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&SpecEgyChange[1], &tot_SpecEgyChange[1], 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&SpecEgyChange[2], &tot_SpecEgyChange[2], 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  MPI_Reduce(&CFracChange[0], &tot_CFracChange[0], 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&CFracChange[2], &tot_CFracChange[2], 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&CFracChange[3], &tot_CFracChange[3], 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

#endif

  if(ThisTask == 0)
    {

#if defined(LT_EJECTA_IN_HOTPHASE) && !defined(LT_SNegy_IN_HOTPHASE)
      if(tot_SpreadEgy[0] > 0)
	{
	  fprintf(FdExtEgy, "> %g %g min %g %g max %g %g - %g %g %g"
		  "- %g %g %g - %g %g %g \n", All.Time,
		  tot_SpreadEgy[1] / tot_SpreadEgy[0],
		  tot_SpreadMinMaxEgy[0][0], tot_SpreadMinMaxEgy[0][1],
		  tot_SpreadMinMaxEgy[1][0], tot_SpreadMinMaxEgy[1][1],
		  tot_AgbFrac[0] / tot_SpreadEgy[0],
		  tot_AgbFrac[1], tot_AgbFrac[2],
		  tot_SpecEgyChange[0] / tot_SpreadEgy[0],
		  tot_SpecEgyChange[1], tot_SpecEgyChange[2],
		  tot_CFracChange[1] / tot_CFracChange[0], tot_CFracChange[2], tot_CFracChange[3]);
	  fflush(FdExtEgy);
	}
#endif

      for(k = 0; k < LT_NIMFs; k++)
	{
#ifdef LT_SNIa
	  if(tot_SNdata[k][SPEC_snIa][SN_num])
	    {
	      fprintf(FdSn, "%3s %1d %9.7g %9i %9.7g %9.7g ", "Ia", k,
		      All.Time, tot_starsnum, tot_SNdata[k][SPEC_snIa][SN_num],
		      tot_SNdata[k][SPEC_snIa][SN_egy]);
	      for(j = 0; j < LT_NMet; j++)
		fprintf(FdSn, "%9.7g ", tot_Zmass[k][SPEC_snIa][j]);
	      fprintf(FdSn, "\n");
	      fflush(FdSn);
	    }

	  /* write losses */
	  if(tot_SNdata[k][SPEC_snIa + 1][SN_num])
	    {
	      fprintf(FdSnLost, "%3s %1d %9.7g %9.7g %9.7g ", "Ia", k,
		      All.Time, tot_SNdata[k][SPEC_snIa + 1][SN_num], tot_SNdata[k][SPEC_snIa + 1][SN_egy]);
	      for(j = 0; j < LT_NMet; j++)
		fprintf(FdSnLost, "%9.7g ", tot_Zmass[k][SPEC_snIa + 1][j]);
	      fprintf(FdSnLost, "\n");
	      fflush(FdSnLost);
	    }
#endif
#ifdef LT_AGB
	  if(tot_SNdata[k][SPEC_agb][SN_num])
	    {
	      fprintf(FdSn, "%3s %1d %9.7g %9i %9.7g 0 ", "AGB", k, All.Time,
		      tot_starsnum, tot_SNdata[k][SPEC_agb][SN_num]);
	      for(j = 0; j < LT_NMet; j++)
		fprintf(FdSn, "%9.7g ", tot_Zmass[k][SPEC_agb][j]);
	      fprintf(FdSn, "\n");
	      fflush(FdSn);
	    }

	  /* write losses */
	  if(tot_SNdata[k][SPEC_agb + 1][SN_num])
	    {
	      fprintf(FdSnLost, "%3s %1d %9.7g %9.7g ", "AGB", k, All.Time,
		      tot_SNdata[k][SPEC_agb + 1][SN_num]);
	      for(j = 0; j < LT_NMet; j++)
		fprintf(FdSnLost, "%9.7g ", tot_Zmass[k][SPEC_agb + 1][j]);
	      fprintf(FdSnLost, "\n");
	      fflush(FdSnLost);
	    }
#endif
#ifdef LT_SNII
	  if(tot_SNdata[k][SPEC_snII][SN_num] > 0)
	    {
	      fprintf(FdSn, "%3s %1d %9.7g %9i %9.7g %9.7g ", "II", k,
		      All.Time, tot_starsnum, tot_SNdata[k][SPEC_snII][SN_num],
		      tot_SNdata[k][SPEC_snII][SN_egy]);
	      for(j = 0; j < LT_NMet; j++)
		fprintf(FdSn, "%9.7g ", tot_Zmass[k][SPEC_snII][j]);
	      fprintf(FdSn, "\n");
	      fflush(FdSn);
	    }

	  /* write losses */
	  if(tot_SNdata[k][SPEC_snII + 1][SN_num] > 0)
	    {
	      fprintf(FdSnLost, "%3s %1d %9.7g %9.7g %9.7g ", "II", k,
		      All.Time, tot_SNdata[k][SPEC_snII + 1][SN_num], tot_SNdata[k][SPEC_snII + 1][SN_egy]);
	      for(j = 0; j < LT_NMet; j++)
		fprintf(FdSnLost, "%9.7g ", tot_Zmass[k][SPEC_snII + 1][j]);
	      fprintf(FdSnLost, "\n");
	      fflush(FdSnLost);
	    }

	  if(tot_SNdata[k][SPEC_IRA][SN_num] > 0)
	    {
	      fprintf(FdSn, "%3s %1d %9.7g %9i %9.7g %9.7g ", "IRA", k,
		      All.Time, tot_starsnum, tot_SNdata[k][SPEC_IRA][SN_num],
		      tot_SNdata[k][SPEC_snII][SN_egy]);
	      for(j = 0; j < LT_NMet; j++)
		fprintf(FdSn, "%9.7g ", tot_Zmass[k][SPEC_snII][j]);
	      fprintf(FdSn, "\n");
	      fflush(FdSn);
	    }

	  /* write losses */
	  if(tot_SNdata[k][SPEC_IRA + 1][SN_num] > 0)
	    {
	      fprintf(FdSnLost, "%3s %1d %9.7g %9.7g %9.7g ", "IRA", k,
		      All.Time, tot_SNdata[k][SPEC_snII + 1][SN_num], tot_SNdata[k][SPEC_snII + 1][SN_egy]);
	      for(j = 0; j < LT_NMet; j++)
		fprintf(FdSnLost, "%9.7g ", tot_Zmass[k][SPEC_snII + 1][j]);
	      fprintf(FdSnLost, "\n");
	      fflush(FdSnLost);
	    }
	}
    }
#endif

  if(SEvInfo_grain > SEvInfo_GRAIN)
    {
      SEvInfo_grain = 0;
      write_metallicity_stat();
    }

  tend = second();
  infos[SN_info] += timediff(tstart, tend);
#endif

#ifdef LT_SEv_INFO_DETAILS

  for(j = 0; j < NTask; j++)
    /* consider to use Gatherv.. */
    {
      if(ThisTask == j && ThisTask > 0)
	{
	  MPI_Send(&DetailsPos, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	  MPI_Send(&Details[0], DetailsPos * sizeof(struct details), MPI_BYTE, 0, 1, MPI_COMM_WORLD);
	}
      else if(ThisTask == 0)
	{
	  if(j > 0)
	    {
	      MPI_Recv(&starsnum, 1, MPI_INT, j, 0, MPI_COMM_WORLD, &status);
	      MPI_Recv(&Details[DetailsPos], starsnum * sizeof(struct details), MPI_BYTE, j, 1,
		       MPI_COMM_WORLD, &status);
	      DetailsPos += starsnum;
	    }
	}
      MPI_Barrier(MPI_COMM_WORLD);
    }

  if(ThisTask == 0)
    fwrite(&Details[0], sizeof(struct details), DetailsPos, FdSnDetails);

  free(Details);

#endif

  MPI_Reduce(&infos[0], &sum_infos[0], INUM, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      All.CPU_SN_info += sum_infos[SN_info] / NTask;
      All.CPU_SN_Comm += sum_infos[SN_Comm] / NTask;
      All.CPU_SN_Calc += sum_infos[SN_Calc] / NTask;
      All.CPU_SN_NeighFind += sum_infos[SN_NeighFind] / NTask;
      All.CPU_SN_NeighCheck += sum_infos[SN_NeighCheck] / NTask;
      All.CPU_SN_Imbalance += sum_infos[SN_Imbalance] / NTask;
      All.CPU_SN_Spread += sum_infos[SN_Spread] / NTask;
    }

  search_for_metalspread = 0;

  return tot_starsnum;
}

/*
   =======================================================
     END of   M a i n   R o u t i n e
     .........................................
   =======================================================
*/


/* ********************************************************************************************************** */
/* ********************************************************************************************************** */


/* NOTE: not used now */
int iterate(int n, int ntrue, FLOAT grad, double *deltax)
{
  FLOAT old, Lold, L, L3, rho;
  int i = 0;

  Lold = *deltax;		/* at the begin *deltax contains the guess for L */
  rho = n / (*deltax * *deltax * *deltax);
  L = *deltax * pow(ntrue / n, 0.33333333);
  L3 = L * L * L;
  *deltax = L - *deltax;
  old = *deltax / 2;

  while(fabs(*deltax - old) / (*deltax) > 0.001)
    {
      i++;
      old = *deltax;
      L3 = L * L * L;
      *deltax = (ntrue - rho * L3) / (grad * L3);
      L = Lold + *deltax;
    }
  return i;
}


/* ...........................................................
   * ------------------------------------------------------- *
   |                                                         |
   |   Searching for Neighbours                              |
   * ------------------------------------------------------- *
*/

void neighbours_count(int i, int mode)
{
  int I;
  int numngb, realngb, *pnumngb;
  int j, n, startnode;
  double linv, linv3, tot_weight, myweight;
  double dx, dy, dz, dist, L, L2, u;
  FLOAT Pos[3];
  double *w, wk;

#ifdef PERIODIC
  FLOAT xtmp;
#endif

  if(mode)
    {
      for(j = 0; j < 3; j++)
	Pos[j] = MetalNgbGet[i].Pos[j];
      L = MetalNgbGet[i].L;
      pnumngb = &MetalNgbResultIn[i].num_ngb;
      w = &MetalNgbResultIn[i].weight;
    }
  else
    {
      I = P[i].MetID;
      for(j = 0; j < 3; j++)
	Pos[j] = P[i].Pos[j];
      L = MetP[I].SLguess;
      pnumngb = &MetP[I].NumNgb;
      w = &MetP[I].weight;
    }

  realngb = 0;
  tot_weight = 0;
  startnode = All.MaxPart;
  linv = 1.0 / (L);
  linv3 = linv * linv * linv;
  L2 = L * L;

  do
    {
      numngb = ngb_treefind_variable(&Pos[0], L, &startnode);

      if(numngb > 0)
	{
	  for(n = 0; n < numngb; n++)
	    {
	      j = Ngblist[n];

	      dx = NGB_PERIODIC_LONG_X(Pos[0] - P[j].Pos[0]);
	      dy = NGB_PERIODIC_LONG_Y(Pos[1] - P[j].Pos[1]);
	      dz = NGB_PERIODIC_LONG_Z(Pos[2] - P[j].Pos[2]);

	      if((dist = dx * dx + dy * dy + dz * dz) <= L2)
		{
		  if(dist == 0)
		    /* this is a newborn star not yet kicked along with its gas parent */
		    continue;

		  realngb++;

#if defined(LT_USE_TOP_HAT_WEIGHT)
		  myweight = 1;
#else
#ifndef LT_USE_SOLIDANGLE_WEIGHT
		  dist = sqrt(dist);
		  u = dist * linv;

		  if(u < 0.5)
		    myweight = (wk = linv3 * (2.546479089470 + 15.278874536822 * (u - 1.0) * u * u));
		  else
		    myweight = (wk = linv3 * 5.092958178941 * (1.0 - u) * (1.0 - u) * (1.0 - u));
#else
		  myweight = SphP[j].Hsml * SphP[j].Hsml / (dist * dist);
#endif
#endif
		  tot_weight += myweight * P[j].Mass;
		}
	    }
	}
    }
  while(startnode >= 0);

  *w = tot_weight;
  *pnumngb = realngb;

  return;
}

/*! This routine is a comparison kernel used in a sort routine to group
 *  particles that are exported to the same processor.
 */

int INLINE_FUNC sev_compare_key_ngb(const void *a, const void *b)
{
  if(((struct metal_ngbfindingdata_in *) a)->Task < (((struct metal_ngbfindingdata_in *) b)->Task))
    return -1;

  if(((struct metal_ngbfindingdata_in *) a)->Task > (((struct metal_ngbfindingdata_in *) b)->Task))
    return +1;

  return 0;
}

int INLINE_FUNC sev_compare_key_spread(const void *a, const void *b)
{
  if(((struct metaldata_in *) a)->Task < (((struct metaldata_in *) b)->Task))
    return -1;

  if(((struct metaldata_in *) a)->Task > (((struct metaldata_in *) b)->Task))
    return +1;

  return 0;
}


/* ...........................................................
   * ------------------------------------------------------- *
   |                                                         |
   |   Evolving Stars                                        |
   * ------------------------------------------------------- *
*/

/* This routine return the stellar evolution products for a given star
 * it also upgrade the mass of the star and the last chemical times to the
 * current values.
 */
#if !defined(LT_EJECTA_IN_HOTPHASE) && !defined(LT_HOT_EJECTA) && !defined(LT_SNegy_IN_HOTPHASE)
void stellarevolution(int i, int mode, float *metals, double *energy)
#else
float stellarevolution(int i, int mode, float *metals, double *energy)
#endif
{
  int I, IMFi, Yset, Zset;
  double mymetals[LT_NMet], sum_mymetals[LT_NMet];
  double myenergy;
  double Z, mass, nonproc_metals, nonproc_metal;
  double lifetime, delta_evolve_time, prec_evolve_time;
  double inf_mass, sup_mass;

  int current_chem_bin;

#if defined(LT_EJECTA_IN_HOTPHASE) || defined(LT_HOT_EJECTA) || defined(LT_SNegy_IN_HOTPHASE)
  float LMmass = 0;
#endif

  int j, ti_step;

  I = P[i].MetID;
  /* find the IMF associated to this particle */
  IMFi = get_IMF_index(i);
  IMFp = &IMFs[IMFi];
  /* find the yield set associated to this IMF */
  Yset = IMFs[IMFi].YSet;

  for(j = 0; j < LT_NMet; j++)
    mymetals[j] = sum_mymetals[j] = 0;
  *energy = 0;

  /* calculate the lifetime of the star */
  if(P[i].StellarAge < (float) All.Time)
    lifetime = get_age(P[i].StellarAge) - All.Time_Age;
  else
    /* due to float-double conversion possible round off */
    lifetime = 0;

#ifdef UM_COOLING
  Fill_mol_num = 0;
#endif

#ifdef LT_SNII
  /* SnII */
  if(P[i].Type & SNIIflag)
    {

#ifdef LT_TRACK_CONTRIBUTES
      for(j = nonproc_metals = 0; j < LT_NMet; j++)
	{
	  IIcontrib[j] = 1;
	  Iacontrib[j] = AGBcontrib[j] = 0;
	}
#endif

      get_timing(2, lifetime, &current_chem_bin, IMFi);

      if((lifetime >= IMFs[IMFi].ShortLiv_TimeTh) &&
	 (current_chem_bin > 0) && (lifetime - MetP[I].LastIITime >= All.MinChemTimeStep))
	{
	  prec_evolve_time = MetP[I].LastIITime;
	  if(lifetime > All.mean_lifetime)
	    lifetime = All.mean_lifetime;

	  inf_mass = dying_mass(lifetime);
	  sup_mass = dying_mass(prec_evolve_time);

	  mass = get_SnII_product(i, Yset, &mymetals[0], &myenergy, inf_mass, sup_mass);

	  if(mass >= 0)
	    {
	      P[i].Mass -= mass;

#ifdef LT_ACCOUNT_NONPROC_METALS
	      for(j = nonproc_metals = 0; j < LT_NMetP; j++)
		if(j != Hel)
		  {
		    mymetals[j] += (nonproc_metal = mymetals[Hyd] * (MetP[I].Metals[j] / MetP[I].iMass));
		    nonproc_metals += nonproc_metal;
		  }
	      mymetals[Hyd] -= nonproc_metals;
	      for(j = nonproc_metals = 0; j < LT_NMetP; j++)
		if(j != Hel)
		  {
		    mymetals[j] += (nonproc_metal = mymetals[Hel] * (MetP[I].Metals[j] / MetP[I].iMass));
		    nonproc_metals += nonproc_metal;
		  }
	      mymetals[Hel] -= nonproc_metals;
#endif

#ifdef LT_SEv_INFO_DETAILS
	      Details[DetailsPos].ID = P[i].ID;
	      Details[DetailsPos].type = (char) 2;
	      Details[DetailsPos].Data[0] = All.Time;
	      Details[DetailsPos].Data[1] = MetP[I].iMass;
	      Details[DetailsPos].Data[2] = prec_evolve_time;
	      Details[DetailsPos].Data[3] = lifetime;
	      for(j = 0; j < LT_NMet; j++)
		Details[DetailsPos].Data[DetailsZ + j] = mymetals[j];
	      DetailsPos++;
#endif
	      for(j = 0; j < LT_NMet; j++)
		sum_mymetals[j] += mymetals[j];
	      *energy += myenergy;
	    }

	  MetP[I].LastIITime = lifetime;
	}

      if(mode == EVOLVE_SN)
	{
	  if(All.Ti_Current < TIMEBASE && lifetime < All.mean_lifetime)
	    {
	      ti_step = get_chemstep(2, lifetime, All.Time, All.Ti_Current, All.Ti_Current, IMFi);
	      if(ti_step > 0)
		{
		  if((MetP[I].NextChemStepII = All.Ti_Current + ti_step) > TIMEBASE)
		    MetP[I].NextChemStepII = TIMEBASE - 1;
		}
	      else
		MetP[I].NextChemStepII = TIMEBASE + 1;
	    }
	  else
	    MetP[I].NextChemStepII = TIMEBASE + 1;
	}


#ifdef LT_SEvDbg
      if(FirstID != 0 && P[i].ID == FirstID)
	{
	  printf("[SEvDBG] II %g %g %g %g %g ", MetP[I].iMass, prec_evolve_time, lifetime, mass, myenergy);
	  for(j = 0; j < LT_NMet; j++)
	    printf("%g ", mymetals[j]);
	  printf("\n   %i :: %i %g\n", MetP[I].NextChemStepII, current_chem_bin,
		 SnII_steps[IMFi][1][current_chem_bin]);
	  fflush(stdout);
	}
#endif


      P[i].Type &= ~SNIIflag;
    }
#endif


#if defined(LT_SNIa) || defined(LT_AGB)
  if(P[i].Type & SNIAflag)
    {
#ifdef LT_SEvDbg
      get_timing(1, lifetime, &current_chem_bin, IMFi);
#endif
      prec_evolve_time = MetP[I].LastIaTime;
      delta_evolve_time = lifetime - prec_evolve_time;

      if(delta_evolve_time >= All.MinChemTimeStep)
	{
	  inf_mass = dying_mass(prec_evolve_time + delta_evolve_time);
	  sup_mass = dying_mass(prec_evolve_time);

#ifdef LT_SNIa

	  mass = get_SnIa_product(i, Yset, &mymetals[0], &myenergy, prec_evolve_time, delta_evolve_time);

	  if(mass >= 0)
	    {
	      P[i].Mass -= mass;

	      *energy += myenergy;
	      for(j = 0; j < LT_NMet; j++)
		{
		  sum_mymetals[j] += mymetals[j];
#ifdef LT_TRACK_CONTRIBUTES
		  Iacontrib[j] = (float) mymetals[j];
		  IIcontrib[j] = 0;
#endif
		}

#ifdef LT_SEv_INFO_DETAILS
	      Details[DetailsPos].ID = P[i].ID;
	      Details[DetailsPos].type = (char) 1;
	      Details[DetailsPos].Data[0] = All.Time;
	      Details[DetailsPos].Data[1] = MetP[I].iMass;
	      Details[DetailsPos].Data[2] = prec_evolve_time;
	      Details[DetailsPos].Data[3] = prec_evolve_time + delta_evolve_time;
	      for(j = 0; j < LT_NMet; j++)
		Details[DetailsPos].Data[DetailsZ + j] = mymetals[j];
	      DetailsPos++;
#endif
	    }

#ifdef LT_SEvDbg
	  if(FirstID > 0 && P[i].ID == FirstID)
	    {
	      printf("[SEvDBG] I %g %g %g %g %g ", MetP[I].iMass, prec_evolve_time, lifetime, mass, myenergy);
	      for(j = 0; j < LT_NMet; j++)
		printf("%g ", mymetals[j]);
	      printf("\n");
	      fflush(stdout);
	    }
#endif

#endif
#ifdef LT_AGB
	  mass = get_AGB_product(i, Yset, &mymetals[0], inf_mass, sup_mass);

	  if(mass >= 0)
	    {
	      P[i].Mass -= mass;

#ifdef LT_ACCOUNT_NONPROC_METALS
	      for(j = nonproc_metals = 0; j < LT_NMetP; j++)
		if(j != Hel)
		  {
		    mymetals[j] += (nonproc_metal = mymetals[Hyd] * (MetP[I].Metals[j] / MetP[I].iMass));
		    nonproc_metals += nonproc_metal;
		  }
	      mymetals[Hyd] -= nonproc_metals;
	      for(j = nonproc_metals = 0; j < LT_NMetP; j++)
		if(j != Hel)
		  {
		    mymetals[j] += (nonproc_metal = mymetals[Hel] * (MetP[I].Metals[j] / MetP[I].iMass));
		    nonproc_metals += nonproc_metal;
		  }
	      mymetals[Hel] -= nonproc_metals;
	      /* this second term try to account for the enrichment due to pristine metals carried by ejected H and He:
	       * we assume that this non processed gas is enriched as much as the gas which formed the star */
#endif

	      for(j = 0; j < LT_NMet; j++)
		{
		  sum_mymetals[j] += mymetals[j];
#if defined(LT_EJECTA_IN_HOTPHASE) || defined(LT_HOT_EJECTA) || defined(LT_SNegy_IN_HOTPHASE)
		  /* keep track of ejecta that are not injected explosively */
		  LMmass += mymetals[j];
#endif
#ifdef LT_TRACK_CONTRIBUTES
		  if(sum_mymetals[j] > 0)
		    {
		      AGBcontrib[j] = (float) (mymetals[j] / sum_mymetals[j]);
		      Iacontrib[j] /= (float) sum_mymetals[j];
		    }
		  else
		    {
		      AGBcontrib[j] = 0;
		      Iacontrib[j] = 0;
		    }
#endif
		}

#ifdef LT_SEv_INFO_DETAILS
	      Details[DetailsPos].ID = P[i].ID;
	      Details[DetailsPos].type = (char) 3;
	      Details[DetailsPos].Data[0] = All.Time;
	      Details[DetailsPos].Data[1] = MetP[I].iMass;
	      Details[DetailsPos].Data[2] = prec_evolve_time;
	      Details[DetailsPos].Data[3] = lifetime;
	      for(j = 0; j < LT_NMet; j++)
		Details[DetailsPos].Data[DetailsZ + j] = mymetals[j];
	      DetailsPos++;
#endif

	    }
#ifdef LT_SEvDbg
	  if(FirstID != 0 && P[i].ID == FirstID)
	    {
	      printf("[SEvDBG] AGB %g %g %g %g 0 ", MetP[I].iMass, prec_evolve_time, lifetime, mass);
	      for(j = 0; j < LT_NMet; j++)
		printf("%g ", mymetals[j]);
	      printf("\n");
	      fflush(stdout);
	    }
#endif
	  MetP[I].LastIaTime = lifetime;
#endif
	}

      if(mode == EVOLVE_SN)
	{
	  if(All.Ti_Current < TIMEBASE && lifetime < IMFs[IMFi].sup_lifetime)
	    {
	      ti_step = get_chemstep(1, lifetime, All.Time, All.Ti_Current, All.Ti_Current, IMFi);
	      if(ti_step > 0)
		{
		  if((MetP[I].NextChemStepIa = All.Ti_Current + ti_step) > TIMEBASE)
		    MetP[I].NextChemStepIa = TIMEBASE - 1;
		}
	      else
		MetP[I].NextChemStepIa = TIMEBASE + 1;
	    }
	  else
	    MetP[I].NextChemStepIa = TIMEBASE + 1;
	}


#ifdef LT_SEvDbg
      if(FirstID != 0 && P[i].ID == FirstID)
	{
	  printf("[SEvDBG] IN %i :: %i %g \n", MetP[I].NextChemStepIa, current_chem_bin,
		 LLv_steps[get_IMF_index(i)][1][current_chem_bin]);
	}
      fflush(stdout);
#endif

      P[i].Type &= ~SNIAflag;
    }
#endif


  if(P[i].Mass <= 0)
    {
      /*.. shouldn't occour (oh, really??) */
      printf("  @@@@@@@@@ %i %i %u %i %g %g %g\n", ThisTask, i, P[i].ID, All.NumCurrentTiStep, P[i].Mass,
	     P[i].StellarAge, lifetime);
      fflush(stdout);
    }

  for(j = 0; j < LT_NMet; j++)
    {
      metals[j] = (float) sum_mymetals[j];
#ifdef LT_SEv_INFO_DETAILS
      if(j < LT_NMetP)
	DetailsWo[j] += sum_mymetals[j];
#endif
    }

#if !defined(LT_EJECTA_IN_HOTPHASE) && !defined(LT_HOT_EJECTA) && !defined(LT_SNegy_IN_HOTPHASE)
  return;
#else
  return (float) LMmass;
#endif
}



/* ...........................................................
   * ------------------------------------------------------- *
   |                                                         |
   |   Spreading                                             |
   * ------------------------------------------------------- *
*/


#if !defined(LT_EJECTA_IN_HOTPHASE) && !defined(LT_HOT_EJECTA) && !defined(LT_SNegy_IN_HOTPHASE)
void spread(int i, int mode, FLOAT * Pos, float *metals, double energy)
#else
void spread(int i, int mode, FLOAT * Pos, float *metals, float LMmass, double energy)
#endif
{
  int I;
  int index, j, k, n, startnode, numngb, Type;
  int *ngblist;
  double linv, linv3, myweight, wfac;
  double dx, dy, dz, dist, L2, u;
  double add_mass;
  double egyfrac;

  /*FLOAT Pos[3], L; */
  FLOAT L;
  double weight;

  double a3inv;

#ifdef LT_EJECTA_IN_HOTPHASE
  double dt, current_egy, current_egyhot, hotmass;

#ifdef LT_SNegy_IN_HOTPHASE
  double x_hotejecta, sn_spec_egy;
#endif
#if defined(LT_HOT_EJECTA) || defined(LT_SNegy_IN_HOTPHASE)
  double v;
#endif
#endif

#ifdef LT_SEvDbg
  unsigned int ID = 0;
#endif

#ifdef PERIODIC
  FLOAT xtmp;
#endif

#ifdef LT_TRACK_CONTRIBUTES
  float contrib_metals[LT_NMet];
#endif

#if defined(LT_SEv_INFO) && defined(LT_EJECTA_IN_HOTPHASE) && !defined(LT_SNegy_IN_HOTPHASE)
  double tot_spreadegy = 0, spreadegy_ratio, tot_agbspreadegy = 0;
  double agb_frac;
  double egy_ratio = 0, x_ratio = 0;
#endif

#ifdef UM_COOLING
  double myFillMolNum;
#endif

  a3inv = pow(All.Time, -3.0);

  startnode = All.MaxPart;

  for(k = 0, add_mass = 0; k < LT_NMet; k++)
    add_mass += metals[k];

  if(mode)
    {
      Type = MetalDataGet[i].Type;
      for(j = 0; j < 3; j++)
	Pos[j] = MetalDataGet[i].Pos[j];
      L = MetalDataGet[i].L;
      weight = MetalDataGet[i].weight;
#ifdef LT_SEvDbg
      ID = MetalDataGet[i].ID;
#endif
#ifdef LT_TRACK_CONTRIBUTES
      contrib = MetalDataGet[i].contrib;
#endif
      IMFp = &IMFs[MetalDataGet[i].IMFi];
#ifdef UM_COOLING
      myFillMolNum = MetalDataGet[i].FillMolNum;
#endif
    }
  else
    {
      Type = P[i].Type & 4;
      if(Type)
	{
	  I = P[i].MetID;
	  L = MetP[I].SLguess;
	  weight = MetP[I].weight;
	}
      else
	/* if we enter this else then LT_LOCAL_IRA is not defined;
	   then it is needless to use #ifndef directives.
	 */
	{
	  L = PPP[i].Hsml;
#ifdef NOFIXEDMASSINKERNEL
	  weight = SphP[i].weight;
#else
	  weight = PPP[i].a1.dNumNgb;
#endif
	}
      for(j = 0; j < 3; j++)
	Pos[j] = P[i].Pos[j];

#ifdef LT_SEvDbg
      if(FirstID > 0 && P[i].ID == FirstID)
	ID = P[i].ID;
#endif
#ifdef UM_COOLING
      myFillMolNum = FillMolNum;
#endif
    }

#ifdef LT_EJECTA_IN_HOTPHASE
  agb_frac = LMmass / add_mass;
#endif

#ifdef LT_SNegy_IN_HOTPHASE
  sn_spec_egy = energy / add_mass;
#endif

#ifdef LT_AVOID_SELF_ENRICHMENT
  if(Type == 0)
    {
      weight -= (2.546479089470 * 4.188790204786) / weight;
      if(mode == 0)
	P[i].Mass -= add_mass;
    }
#endif

  linv = 1.0 / (L);
  linv3 = linv * linv * linv;
  L2 = L * L;

  do
    {
      numngb = ngb_treefind_variable(&Pos[0], L, &startnode);
      ngblist = Ngblist;

      if(numngb > 0)
	{
#ifndef LT_AVOID_SELFENRICHMENT
	  if(Type == 0 && mode == 0)
	    /* particle didn't find himself because of modified .Type */
	    ngblist[numngb++] = i;
#endif

	  for(n = 0; n < numngb; n++)
	    {
	      j = ngblist[n];

	      dx = NGB_PERIODIC_LONG_X(Pos[0] - P[j].Pos[0]);
	      dy = NGB_PERIODIC_LONG_Y(Pos[1] - P[j].Pos[1]);
	      dz = NGB_PERIODIC_LONG_Z(Pos[2] - P[j].Pos[2]);

	      if((dist = dx * dx + dy * dy + dz * dz) <= L2)
		{
		  if(Type)
		    /* we are spreading from a star */
		    {
#ifndef LT_USE_TOP_HAT_WEIGHT
#ifndef LT_USE_SOLIDANGLE_WEIGHT
		      dist = sqrt(dist);
		      u = dist * linv;

		      if(u < 0.5)
			myweight = linv3 * (2.546479089470 + 15.278874536822 * (u - 1.0) * u * u);
		      else
			myweight = linv3 * 5.092958178941 * (1.0 - u) * (1.0 - u) * (1.0 - u);
#else
		      myweight = SphP[j].Hsml * SphP[j].Hsml / (dist * dist);
#endif
#else
		      myweight = 1;
#endif
		    }
		  else
		    /* we are spreading from a gas particle; this means
		       that LT_LOCAL_IRA is not defined.
		     */
		    {
		      dist = sqrt(dist);
		      u = dist * linv;

		      if(u < 0.5)
			myweight = linv3 * (2.546479089470 + 15.278874536822 * (u - 1.0) * u * u);
		      else
			myweight = linv3 * 5.092958178941 * (1.0 - u) * (1.0 - u) * (1.0 - u);

#ifndef NOFIXEDMASSINKERNEL
		      myweight = 4.188790204786 * myweight / linv3;	/* note: 4.188790204786 = 4.0/3 * PI */
#endif
		    }

		  if(Type)
		    myweight *= P[j].Mass;
		  else
		    {
#ifdef NOFIXEDMASSINKERNEL
		      myweight *= P[j].Mass;
#endif
		    }


		  if((wfac = myweight / weight) > (1 + 1e-3))
		    {
		      printf("Something strange arose when calculating weight in metal distribution:\n"
			     "weightsum: %.8g <  wfac : %.8g\n  *  spreading length : %.8g\n"
			     "Task : %i, neighbour : %i, particle : %i, index [%i] : %i, numngb: %i\n",
			     weight, wfac, L, ThisTask, j, i, mode, ((mode == LOCAL) ? index : -1), numngb);
		      fflush(stdout);
		      endrun(101010);
		    }


		  if(wfac < 0 && wfac > -5e-3)
		    {
		      printf
			("warning: particle has got a negative weight factor! possibly it is a round-off error.\n"
			 "         [%d][%d][%d] %g %g %g %g\n", ThisTask, mode, j, u, myweight, weight, wfac);
		      wfac = 0;
		    }

#ifdef LT_TRACK_CONTRIBUTES
		  for(k = 0; k < LT_NMet; k++)
		    contrib_metals[k] = metals[k] * wfac;
		  update_contrib(&SphP[j].contrib, &SphP[j].Metals[0], &contrib, &contrib_metals[0]);
#endif

#if defined(LT_SEvDbg)
		  if(ID)
		    weight_sum += wfac;
#endif
		  for(k = 0; k < LT_NMetP; k++)
		    {
		      /* Hydrogen is the last element, we don't store it in the Metals array */
		      if((SphP[j].Metals[k] += metals[k] * wfac) < 0)
			{
			  printf(" \n ooops... it's a shame! %i %i %i %i %i %i %g %g %g \n",
				 ThisTask, mode, j, i, k, Type, metals[k], wfac, weight);
			  endrun(333333);
			}
#ifdef LT_SEv_INFO_DETAILS
		      DetailsW[k] += metals[k] * wfac;
#endif
		    }
#ifdef UM_COOLING
		  SphP[j].FillMolNum += myFillMolNum * wfac;
#endif


		  /* ===============================================================
		   *
		   * a feedback form
		   * ejecta from sn (not from agb!) are put in the hot phase of the
		   * gas; thie means:
		   *  (1) the current intrinsic energy is update supposing that
		   *      the added mass has the same erg/g than the hot phase
		   *  (2) the entropy or the entropy change rate is update accordingly
		   *
		   * also, you can choose to put the ejecta into the hot phase with
		   * some their own specific energy, either using the specific energy
		   * of the supernovae (1^51 erg / ejecta mass for each SN) or a
		   * specific energy that you specify in paramfile. These two options
		   * correspond to switching on either LT_SNegy_IN_HOTPHASE or
		   * LT_HOT_EJECTA.
		   *
		   * =============================================================== */

#ifdef LT_EJECTA_IN_HOTPHASE
		  /* calculate the current specific energy */
		  if(P[j].Ti_endstep == All.Ti_Current)
		    {
		      /* if this is an active particle, calculate from the end-step quantities */
		      dt = (P[j].Ti_endstep - P[j].Ti_begstep) * All.Timebase_interval;
		      current_egy = DMAX(All.MinEgySpec, (SphP[j].Entropy + SphP[j].e.DtEntropy * dt) /
					 GAMMA_MINUS1 * pow(SphP[j].a2.Density * a3inv, GAMMA_MINUS1));
		    }
		  else
		    {
		      dt = 0;
		      current_egy =
			SphP[j].Entropy / GAMMA_MINUS1 * pow(SphP[j].a2.Density * a3inv, GAMMA_MINUS1);
		    }

		  /* calculate the current specific energy of the hot phase */
		  current_egyhot = (current_egy - All.EgySpecCold * SphP[j].x) / (1 - SphP[j].x);
		  hotmass = P[j].Mass * (1 - SphP[j].x);

		  /* update the mass and the cold fraction */
		  if(SphP[j].x > 0)
		    {
		      x_ratio = SphP[j].x;
#ifndef LT_LOCAL_IRA
		      if(dist == 0 && !Type && mode == 0)
			/* the particle itself; should be sufficient dist == 0 */
			{
			  SphP[j].x = (SphP[j].x * P[j].Mass - add_mass) /
			    (P[j].Mass - (1 - wfac) * add_mass);
			  P[j].Mass -= (1 - wfac) * add_mass;
			}
		      else
#endif
			{
			  SphP[j].x *= P[j].Mass / (P[j].Mass + add_mass * wfac);
			  P[j].Mass += add_mass * wfac;
			}

		      x_ratio /= SphP[j].x;
		    }
		  else
		    x_ratio = 0;

#if defined(LT_HOT_EJECTA) || defined(LT_SNegy_IN_HOTPHASE)
		  v = (add_mass - LMmass) * wfac / (hotmass + add_mass * wfac);
#ifdef LT_HOT_EJECTA
		  current_egyhot = current_egyhot * (1 - v) + All.EgySpecEjecta * v;
#endif
#ifdef LT_SNegy_IN_HOTPHASE
		  if(All.NumCurrentTiStep == 1474)
		    printf("  +++ %g %g %g %g ", v, add_mass, LMmass, current_egyhot);
		  current_egyhot = current_egyhot * (1 - v) + sn_spec_egy * v;
		  if(All.NumCurrentTiStep == 1474)
		    printf("  +++ %g\n", current_egyhot);
		  fflush(stdout);

#endif
#endif
		  egy_ratio = current_egy;
		  current_egy = current_egyhot * (1 - SphP[j].x) + All.EgySpecCold * SphP[j].x;
		  egy_ratio /= current_egy;
#if defined(LT_SEv_INFO) &&  defined(LT_EJECTA_IN_HOTPHASE) && !defined(LT_SNegy_IN_HOTPHASE)
		  tot_spreadegy += current_egyhot * wfac;
#endif



		  if(dt > 0 && SphP[j].e.DtEntropy != 0)
		    {
		      SphP[j].e.DtEntropy =
			(current_egy * GAMMA_MINUS1 / pow(SphP[j].a2.Density * a3inv, GAMMA_MINUS1) -
			 SphP[j].Entropy) / dt;
		      if(SphP[j].e.DtEntropy < -0.5 * SphP[j].Entropy / dt)
			SphP[j].e.DtEntropy = -0.5 * SphP[j].Entropy / dt;
		    }
		  else
		    SphP[j].Entropy =
		      current_egy * GAMMA_MINUS1 / pow(SphP[j].a2.Density * a3inv, GAMMA_MINUS1);
#endif

		  /* ===============================================================
		   *
		   * end of feedback
		   * =============================================================== */


#ifndef LT_SNegy_IN_HOTPHASE
		  SphP[j].EgyRes += energy * wfac;
#endif


#ifdef LT_SEv_INFO
#if defined(LT_EJECTA_IN_HOTPHASE) && !defined(LT_SNegy_IN_HOTPHASE)
		  if(energy > 0)
		    {
		      tot_spreadegy *= add_mass / energy;
		      tot_agbspreadegy = tot_spreadegy * agb_frac;

		      SpreadEgy[0] += 1;
		      SpreadEgy[1] += tot_spreadegy;

		      if(SpreadMinMaxEgy[0][0] > tot_spreadegy)
			{
			  SpreadMinMaxEgy[0][0] = tot_spreadegy;
			  SpreadMinMaxEgy[0][1] = tot_agbspreadegy;
			}
		      if(SpreadMinMaxEgy[1][0] < tot_spreadegy)
			{
			  SpreadMinMaxEgy[1][0] = tot_spreadegy;
			  SpreadMinMaxEgy[1][1] = tot_agbspreadegy;
			}

		      AgbFrac[0] += tot_agbspreadegy;
		      if(AgbFrac[1] > tot_agbspreadegy)
			AgbFrac[1] = tot_agbspreadegy;
		      if(AgbFrac[2] < tot_agbspreadegy)
			AgbFrac[2] = tot_agbspreadegy;

		      SpecEgyChange[0] += egy_ratio;
		      if(SpecEgyChange[1] > egy_ratio)
			SpecEgyChange[1] = egy_ratio;
		      if(SpecEgyChange[2] < egy_ratio)
			SpecEgyChange[2] = egy_ratio;

		      if(x_ratio > 0)
			{
			  CFracChange[0] += 1;
			  CFracChange[1] += x_ratio;
			  if(CFracChange[2] > x_ratio)
			    CFracChange[2] = x_ratio;
			  if(CFracChange[3] < x_ratio)
			    CFracChange[3] = x_ratio;
			}
		    }
#endif

		  if(SEvInfo_grain > SEvInfo_GRAIN)
		    {
		      egyfrac = energy * wfac /
			(SphP[j].Entropy / GAMMA_MINUS1 *
			 pow(SphP[j].a2.Density * All.Time * All.Time * All.Time, GAMMA_MINUS1));

		      if(egyfrac > 0 && egyfrac < Stat_min[MIN_egyf])
			Stat_min[MIN_egyf] = egyfrac;
		      else if(egyfrac > Stat_max[MAX_egyf])
			Stat_max[MAX_egyf] = egyfrac;
		      if(Stat_min[MIN_egyf] == MIN_INIT_VALUE)
			Stat_min[MIN_egyf] = 0;
		      Stat_sum[MEAN_egyf] += egyfrac;
		      Stat_sum[NUM_egyf]++;
		    }

#endif

		}
	    }
	}
    }
  while(startnode >= 0);

  return;
}



/* ******************************************************************************
 *
 * III   III   III
 *  II    II    II
 *  II    II    II
 * IIII  IIII  IIII
 * ..............................................................................
 *
 * SN Ia
 *
 * ****************************************************************************** */


#ifdef LT_SNIa

double get_SnIa_product(int i, int Yset, double *metals, double *energy, double prec_evolve_time,
			double delta_time)
{
  int I, Zbin, Mbin;
  double Zstar, ejecta, result, abserr;
  double next_evolve_time;
  int j;

  /* does not care about the size of the integration time; you must care about that elsewhere */

  I = P[i].MetID;

  next_evolve_time = prec_evolve_time + delta_time;

  F.function = &nRSnIa;
  F.params = NULL;


  if((my_gslstatus =
      gsl_integration_qag(&F, prec_evolve_time, next_evolve_time, 1, 1e-5, 1000, qag_INT_KEY, w, &result,
			  &abserr)))
    {
      printf("  >> Task %i, gsl integration error %i in Sn Ia integration [%.6g - %.6g]: %.6g %.6g\n",
	     ThisTask, my_gslstatus, prec_evolve_time, next_evolve_time, result, abserr);
      fflush(stdout);
      endrun(LT_ERROR_INTEGRATION_ERROR);
    }

  if(result < 0 && result > -Z_ROUND_OFF_ERROR)
    {
      printf(" > probable round-off error in integrating Ia num [%d][%d][%u] %g\n",
	     ThisTask, i, P[i].ID, result);
      fflush(stdout);
      result = 0;
    }

  /* calculate number of SnIa explosions: MetP[I].iMass * UnitMassFact
   * is the initial mass population in unit of Msun
   */

  result *= MetP[I].iMass * UnitMassFact;	/* *agefact */

  if(result == 0)
    {
      for(j = 0; j < LT_NMet; j++)
	metals[j] = 0;
      *energy = 0;
      return 0;
    }

  /* ----------------------------------------------------------------- */

#ifdef LT_SEv_INFO
  if(spreading_on)
    SNdata[IMFp - IMFs][SPEC_snIa][SN_num] += result;	/* sum-up number of SnIa explosions */
  else
    SNdata[IMFp - IMFs][SPEC_snIa + 1][SN_num] += result;
#endif

  Zstar = get_metallicity(i, 2);
  for(Zbin = IaZbins_dim[Yset] - 1; Zstar < IaZbins[Yset][Zbin] && Zbin > 0; Zbin--)
    ;

  if(IaMbins_dim[Yset] > 1)
    {
      /*
       * !!!!!!!!! note: not supported yet !!!!!!!!!
       * the following line is just for completeness
       */
      for(Mbin = IaMbins_dim[Yset] - 1; MetP[I].iMass < IaMbins[Yset][Mbin] && Mbin > 0; Mbin--)
	;
    }
  else
    Mbin = 0;


  for(j = ejecta = 0; j < LT_NMet; j++)
    {
      /* gives mass of jth element in unit of Msun */
      metals[j] = result * IaY(Yset, j, Zbin, Mbin);

#ifdef LT_SEv_INFO
      if(spreading_on && j < LT_NMet)
	Zmass[IMFp - IMFs][SPEC_snIa][j] += metals[j];
      else
	Zmass[IMFp - IMFs][SPEC_snIa + 1][j] += metals[j];
#endif
      metals[j] /= UnitMassFact;	/* normalize to code units */
      ejecta += metals[j];
    }

  *energy = result * All.SnIaEgy;	/* erg in code units */

#ifdef LT_SEv_INFO
  if(spreading_on)
    SNdata[IMFp - IMFs][SPEC_snIa][SN_egy] += *energy * All.UnitEnergy_in_cgs;
  else
    SNdata[IMFp - IMFs][SPEC_snIa + 1][SN_egy] += *energy * All.UnitEnergy_in_cgs;
#endif

  /* return ejected mass (1.4 Msun) */
  /*result /= UnitMassFact;
     return result * SnIaEjectedMass;  return ejected mass */

  return ejecta;
}

double INLINE_FUNC inner_integrand(double m, void *params)
{
  double m2 = m * m;
  double (*phi) (double, void *);

  phi = IMFp->IMFfunc_byNum;

  return (phi(m, 0x0) / (m2 * m2));
}


double INLINE_FUNC nRSnIa(double t, void *params)
{

  /* calculates Sn type Ia explosion rate in #/yr    */
  /* at some time t, for a SSP of 1Msun.             */
  /* in snIa_mass return the total mass in SnIa.     */
  /* in Rejecta return the rate of ejecta in Msun/yr */
  /* t is in Gyr */

  gsl_function localF;
  double result, err;
  size_t n_eval;

  double m2, Mb_inf, Mb_sup, fact1, fact2;
  int i, sbin1, sbin2;
  double (*phi) (double, void *);

  phi = IMFp->IMFfunc_byNum;


  m2 = dying_mass(t);

  if((m2 < All.MBms) || (m2 > All.Mup))
    return 0;

  Mb_inf = max(2 * m2, All.MBm);
  Mb_sup = 0.5 * All.MBM + m2;


  if(IMFp->type == power_law)
    {
      if(IMFp->NSlopes == 1)
	sbin1 = sbin2 = 0;
      else
	{
	  sbin1 = get_IMF_SlopeBin(Mb_sup);
	  sbin2 = get_IMF_SlopeBin(Mb_inf);
	}

      if(sbin1 == sbin2)
	{
	  fact2 = 1.0 / (3 + IMFp->Slopes.slopes[sbin1]);
	  result = ((fact2 * pow(Mb_inf, -(3 + IMFp->Slopes.slopes[sbin1]))) -
		    (fact2 * pow(Mb_sup, -(3 + IMFp->Slopes.slopes[sbin1]))));
	}
      else
	{
	  fact1 = Mb_inf;
	  for(i = sbin1, result = 0; i <= sbin2; i++)
	    {
	      if(i > sbin1)
		Mb_sup = IMFp->Slopes.masses[i];

	      if(i < sbin2)
		Mb_sup = IMFp->Slopes.masses[i + 1];
	      else
		Mb_sup = fact1;

	      fact2 = 1.0 / (3 + IMFp->Slopes.slopes[i]);
	      result += ((fact2 * pow(Mb_inf, -(3 + IMFp->Slopes.slopes[i]))) -
			 (fact2 * pow(Mb_sup, -(3 + IMFp->Slopes.slopes[i]))));
	    }
	}
    }
  else
    {
      localF.function = &inner_integrand;
      localF.params = 0x0;

      if((gsl_status = gsl_integration_qng(&localF, Mb_inf, Mb_sup, 1, 1e-5, &result, &err, &n_eval)))
	{
	  printf("  >> Task %i, gsl integration error %i in Sn (k) [%.6g - %.6g]: %.6g %.6g\n",
		 ThisTask, gsl_status, Mb_inf, Mb_sup, result, err);
	  fflush(stdout);
	}
    }

  fact1 = IMFp->A * All.BinFrac * 24 * m2 * m2 * (-dm_dt(m2, t));
  return (fact1 * result);
}

#endif


/* ******************************************************************************
 *
 * III  III   II
 *  II   II   II
 *  II    II II
 * IIII    III
 * ..............................................................................
 *
 * LT_AGB - SN II
 *
 * ****************************************************************************** */


#ifdef LT_AGB

double get_AGB_product(int i, int Yset, double *metals, double inf_mass, double sup_mass)
{
  int I;
  double Zstar, err, result, ejecta, tot_num;
  int j;

  /* does not care about the size of the integration time; you must care about that elsewhere */

  I = P[i].MetID;

  tot_num = IntegrateIMF_byNum(inf_mass, sup_mass, IMFp, EXC_BH);

  if(tot_num < 0 && tot_num > -Z_ROUND_OFF_ERROR)
    {
      printf(" > probable round-off error in integrating AGB num [%d][%d][%u] %g\n",
	     ThisTask, i, P[i].ID, result);
      fflush(stdout);
      result = 0;
    }
  tot_num *= (MetP[I].iMass * UnitMassFact) / All.HubbleParam * (1 - All.BinFrac);

  if(tot_num == 0)
    {
      for(j = 0; j < LT_NMet; j++)
	metals[j] = 0;
      return 0;
    }

#ifdef LT_SEv_INFO
  if(spreading_on)
    SNdata[IMFp - IMFs][SPEC_agb][SN_num] += tot_num;	/* sum-up number */
  else
    SNdata[IMFp - IMFs][SPEC_agb + 1][SN_num] += tot_num;
#endif

  Zstar = get_metallicity(i, 2);
  for(SD.Zbin = AGBZbins_dim[Yset] - 1; Zstar < AGBZbins[Yset][SD.Zbin] && SD.Zbin > 0; SD.Zbin--)
    ;

  SD.Zstar = Zstar;
  SD.Zdim = AGBZbins_dim[Yset];
  SD.ZArray = AGBZbins[Yset];
  SD.Mdim = AGBMbins_dim[Yset];
  SD.MArray = AGBMbins[Yset];

  ejecta = 0;
  F.function = &zmRSnII;
  F.params = 0x0;

  for(j = 0; j < LT_NMet; j++)
    {
      SD.Y = AGBYields[Yset][j];

      if((my_gslstatus =
	  gsl_integration_qag(&F, inf_mass, sup_mass, 1e-6, 1e-4, gsl_ws_dim, qag_INT_KEY, w, &result, &err)))
	{
	  printf("  >> Task %i, gsl integration error %i in AGB z integration [%.6g - %.6g] : %.6g %.6g\n",
		 ThisTask, my_gslstatus, inf_mass, sup_mass, result, err);
	  fflush(stdout);
	  endrun(LT_ERROR_INTEGRATION_ERROR);
	}

      if(result < 0)
	{
	  printf(" $$$$$ AGB [%d][%d] %g %g %g %d\n", ThisTask, i, inf_mass, sup_mass, result, SD.Zbin);
	  fflush(stdout);
	}

      if(result < 0 && result > -Z_ROUND_OFF_ERROR)
	{
	  printf(" > probable round-off error in integrating AGB [%d][%d][%u][%d] %g\n",
		 ThisTask, i, P[i].ID, j, result);
	  fflush(stdout);
	  result = 0;
	}

      ejecta += (metals[j] = result * MetP[I].iMass * (1 - All.BinFrac));	/* *agefact */

#ifdef LT_SEv_INFO
      if(spreading_on)
	Zmass[IMFp - IMFs][SPEC_agb][j] += metals[j] * UnitMassFact;
      else
	Zmass[IMFp - IMFs][SPEC_agb + 1][j] += metals[j] * UnitMassFact;
#endif

#ifdef UM_COOLING
      if(j == FillEl)
	{
	  SD.Y = AGB_Fill_mol_num[Yset];
	  if((my_gslstatus =
	      gsl_integration_qag(&F, inf_mass, sup_mass, 1e-6, 1e-4, gsl_ws_dim, qag_INT_KEY, w, &result,
				  &err)))
	    {
	      printf
		("  >> Task %i, gsl integration error %i in SnII z integration [%.6g - %.6g] : %.6g %.6g\n",
		 ThisTask, my_gslstatus, inf_mass, sup_mass, result, err);
	      fflush(stdout);
	      endrun(LT_ERROR_INTEGRATION_ERROR);
	    }

	  if(result < 0 && result > -Z_ROUND_OFF_ERROR)
	    {
	      printf(" > probable round-off error in integrating SnII [%d][%d][%u] %g\n",
		     ThisTask, i, P[i].ID, result);
	      fflush(stdout);
	      result = 0;
	    }

	  *Fill_mol_num += result * MetP[i].iMass;
	}
#endif

    }

  return ejecta;
}

#endif

#ifdef LT_SNII

double get_SnII_product(int i, int Yset, double *metals, double *energy, double inf_mass, double sup_mass)
{
  int I;
  double Zstar, err, result, ejecta;
  double sup_mass_store, tot_num;
  int j;

  /* does not care about the size of the integration time; you must care about that elsewhere */

  I = P[i].MetID;

  *energy = 0;
  ejecta = 0;
  result = 0;
  for(j = 0; j < LT_NMet; j++)
    metals[j] = 0;

  if(sup_mass < All.Mup)
    return 0;

  tot_num = IntegrateIMF_byNum(inf_mass, sup_mass, IMFp, EXC_BH);
  tot_num *= (MetP[I].iMass * UnitMassFact) / All.HubbleParam;

  sup_mass_store = sup_mass;

  if(inf_mass < IMFp->egyShortLiv_MassTh)
    {
      if(sup_mass_store > IMFp->egyShortLiv_MassTh)
	sup_mass_store = IMFp->egyShortLiv_MassTh;

      *energy = IntegrateIMF_byEgy(inf_mass, sup_mass_store, IMFp, EXC_BH);

      if(result < 0 && result > -Z_ROUND_OFF_ERROR)
	{
	  printf(" > probable round-off error in integrating SnII num [%d][%d][%u] %g\n",
		 ThisTask, i, P[i].ID, result);
	  fflush(stdout);
	  result = 0;
	}

      *energy *= (MetP[I].iMass * UnitMassFact) / All.HubbleParam;	/* *agefact */
    }

#ifdef LT_SEv_INFO
  if(spreading_on)
    {
      SNdata[IMFp - IMFs][SPEC_snII][SN_num] += tot_num;
      SNdata[IMFp - IMFs][SPEC_snII][SN_egy] += *energy * All.UnitEnergy_in_cgs;
    }
  else
    {
      SNdata[IMFp - IMFs][SPEC_snII + 1][SN_num] += tot_num;
      SNdata[IMFp - IMFs][SPEC_snII + 1][SN_egy] += *energy * All.UnitEnergy_in_cgs;
    }
#endif

  /*
   * integrate metal production and H&He ejection
   *
   *    -> metals
   */

  if(inf_mass < IMFp->metShortLiv_MassTh)
    {
      if(sup_mass > IMFp->metShortLiv_MassTh)
	sup_mass = IMFp->metShortLiv_MassTh;

      Zstar = get_metallicity(i, 2);
      for(SD.Zbin = IIZbins_dim[Yset] - 1; Zstar < IIZbins[Yset][SD.Zbin] && SD.Zbin > 0; SD.Zbin--)
	;

      SD.Zstar = Zstar;
      SD.Zdim = IIZbins_dim[Yset];
      SD.ZArray = IIZbins[Yset];
      SD.Mdim = IIMbins_dim[Yset];
      SD.MArray = IIMbins[Yset];

      F.function = &zmRSnII;
      F.params = 0x0;

      for(j = 0; j < LT_NMet; j++)
	{
	  SD.Y = SnIIYields[Yset][j];

	  if((my_gslstatus =
	      gsl_integration_qag(&F, inf_mass, sup_mass, 1e-6, 1e-4, gsl_ws_dim, qag_INT_KEY, w, &result,
				  &err)))
	    {
	      printf
		("  >> Task %i, gsl integration error %i in SnII z integration [%.6g - %.6g] : %.6g %.6g\n",
		 ThisTask, my_gslstatus, inf_mass, sup_mass, result, err);
	      fflush(stdout);
	      endrun(LT_ERROR_INTEGRATION_ERROR);
	    }

	  if(result < 0 && result > -Z_ROUND_OFF_ERROR)
	    {
	      printf(" > probable round-off error in integrating SnII [%d][%d][%u] %g\n",
		     ThisTask, i, P[i].ID, result);
	      fflush(stdout);
	      result = 0;
	    }

	  /* give metal mass in Msun in code units */
	  ejecta += (metals[j] = result * MetP[I].iMass);	/* *agefact */

#ifdef LT_SEv_INFO
	  if(spreading_on)
	    Zmass[IMFp - IMFs][SPEC_snII][j] += metals[j] * UnitMassFact;
	  else
	    Zmass[IMFp - IMFs][SPEC_snII + 1][j] += metals[j] * UnitMassFact;
#endif

#ifdef UM_COOLING
	  if(j == FillEl)
	    {
	      SD.Y = SnII_Fill_mol_num[Yset];
	      if((my_gslstatus =
		  gsl_integration_qag(&F, inf_mass, sup_mass, 1e-6, 1e-4, gsl_ws_dim, qag_INT_KEY, w, &result,
				      &err)))
		{
		  printf
		    ("  >> Task %i, gsl integration error %i in SnII z integration [%.6g - %.6g] : %.6g %.6g\n",
		     ThisTask, my_gslstatus, inf_mass, sup_mass, result, err);
		  fflush(stdout);
		  endrun(LT_ERROR_INTEGRATION_ERROR);
		}

	      if(result < 0 && result > -Z_ROUND_OFF_ERROR)
		{
		  printf(" > probable round-off error in integrating SnII [%d][%d][%u] %g\n",
			 ThisTask, i, P[i].ID, result);
		  fflush(stdout);
		  result = 0;
		}

	      *Fill_mol_num += result * MetP[i].iMass;
	    }
#endif

	}
    }

  return ejecta;
}
#endif

double INLINE_FUNC nRSnII(double time, void *p)
{
  /* calculates how many stars are dying */
  /* at some time t, for a SSP of 1Msun. */
  /* t is in Gyr */

  double m;
  double (*phi) (double, void *);

  phi = IMFp->IMFfunc_byNum;

  if((time >= IMFp->inf_lifetime) &&	/* Mm Msun < m < MUP Msun */
     (time <= IMFp->sup_lifetime))
    {
      m = dying_mass(time);
      return (phi(m, p) * (-dm_dt(m, time)));
    }
  else
    return 0;
}


double INLINE_FUNC mRSnII(double time, void *p)
{
  /* calculates how many stars are dying */
  /* at some time t, for a SSP of 1Msun. */
  /* t is in Gyr                         */
  /* results is in Msun/Gyr              */

  int BHi;
  double m, fact;
  double (*phi) (double, void *);

  phi = IMFp->IMFfunc_byNum;

  if((time >= IMFp->inf_lifetime) &&	/* Mm Msun < m < MUP Msun */
     (time <= IMFp->sup_lifetime))
    {
      m = dying_mass(time);

      for(BHi = 0; (BHi < IMFp->N_notBH_ranges) && (m <= IMFp->notBH_ranges.inf[BHi]); BHi++)
	;
      if((BHi < IMFp->N_notBH_ranges) && (m > IMFp->notBH_ranges.sup[BHi]))
	return 0;

      fact = phi(m, p) * (-dm_dt(m, time));
      return fact * m;
    }
  else
    return 0;
}


double INLINE_FUNC zmRSnII(double m, void *p)
{
  /* calculates Sn type II explosion rate in Msun/Gyr */
  /* at some time t, for a SSP of 1Msun.              */
  /* t is in Gyr */

  int ir, BHi, zone;
  double y, t, u;
  double (*phi) (double, void *);

  for(BHi = 0; (BHi < IMFp->N_notBH_ranges) && (m <= IMFp->notBH_ranges.inf[BHi]); BHi++)
    ;
  if((BHi < IMFp->N_notBH_ranges) && (m > IMFp->notBH_ranges.sup[BHi]))
    return 0;

  phi = IMFp->IMFfunc_byNum;

  if((m >= IMFp->Mm) &&		/* Mm Msun < m < MU Msun */
     (m <= IMFp->MU))
    {
      for(ir = SD.Mdim - 1; m < SD.MArray[ir] && ir > 0; ir--)
	;

      if(SD.Zdim > 1)
	{
	  zone = 0;
	  zone |= (SD.Zstar >= SD.ZArray[SD.Zdim - 1]);
	  zone |= ((m >= SD.MArray[SD.Mdim - 1]) << 1);
	  zone |= ((SD.Zstar < SD.ZArray[0]) << 2);
	  zone |= ((m < SD.MArray[0]) << 3);
	}
      else
	{
	  zone = 16;
	  zone |= (m >= SD.MArray[SD.Mdim - 1]);
	  zone |= ((m < SD.MArray[0]) << 1);
	}

      switch (zone)
	{
	case 0:
	  t = (m - SD.MArray[ir]) / (SD.MArray[ir + 1] - SD.MArray[ir]);
	  u = (SD.Zstar - SD.ZArray[SD.Zbin]) / (SD.ZArray[SD.Zbin + 1] - SD.ZArray[SD.Zbin]);
	  y = (1 - t) * (1 - u) * SD.Y[SD.Zbin * SD.Mdim + ir] +
	    t * (1 - u) * SD.Y[SD.Zbin * SD.Mdim + ir + 1] +
	    t * u * SD.Y[(SD.Zbin + 1) * SD.Mdim + ir + 1] + (1 - t) * u * SD.Y[(SD.Zbin + 1) * SD.Mdim + ir];
	  break;
	case 2:
	case 8:
	  u = (SD.Zstar - SD.ZArray[SD.Zbin]) / (SD.ZArray[SD.Zbin + 1] - SD.ZArray[SD.Zbin]);
	  y = (1 - u) * SD.Y[SD.Zbin * SD.Mdim + ir] + u * SD.Y[(SD.Zbin + 1) * SD.Mdim + ir];
	  if(zone == 8)
	    y *= (m / SD.MArray[0]);
	  break;
	case 1:
	case 4:
	  t = (m - SD.MArray[ir]) / (SD.MArray[ir + 1] - SD.MArray[ir]);
	  y = (1 - t) * SD.Y[SD.Zbin * SD.Mdim + ir] + t * SD.Y[SD.Zbin * SD.Mdim + ir + 1];
	  break;
	case 9:
	case 12:
	  y = SD.Y[SD.Zbin * SD.Mdim + ir] * (m / SD.MArray[0]);
	  break;
	case 3:
	case 6:
	  y = SD.Y[SD.Zbin * SD.Mdim + ir];
	  break;
	case 16:
	  t = (m - SD.MArray[ir]) / (SD.MArray[ir + 1] - SD.MArray[ir]);
	  y = (1 - t) * SD.Y[SD.Zbin * SD.Mdim + ir] + t * SD.Y[SD.Zbin * SD.Mdim + ir + 1];
	  break;
	case 17:
	  if(SD.ExtrDir & 1 && !SD.ExtrDir & 1 << 31)
	    {
	      t = (m - SD.MArray[SD.Mdim - 1]) / (SD.ExtrMArray[0] - SD.MArray[SD.Mdim - 1]);
	      y = (1 - t) * SD.Y[SD.Zbin * SD.Mdim + SD.Mdim - 1] + t * SD.ExtrY[SD.ExtrZbin * SD.ExtrMdim];
	    }
	  else
	    y = SD.Y[SD.Zbin * SD.Mdim + SD.Mdim - 1];
	  break;
	case 18:
	  if(SD.ExtrDir & 1 && SD.ExtrDir & 1 << 31)
	    {
	      t = (m - SD.ExtrMArray[SD.ExtrMdim - 1]) / (SD.MArray[0] - SD.ExtrMArray[SD.ExtrMdim - 1]);
	      y = (1 - t) * SD.ExtrY[SD.ExtrZbin * SD.ExtrMdim + SD.ExtrMdim - 1] +
		t * SD.Y[SD.Zbin * SD.Mdim];
	    }
	  else
	    y = SD.Y[SD.Zbin * SD.Mdim] * (m / SD.MArray[0]);
	  break;
	}

      return y * phi(m, p);
    }
  else
    return 0;
}


double INLINE_FUNC ejectaSnII(double m, void *p)
{
  /* calculates Sn type II restored mass in Msun */
  /* for a star of mass m */

  int ir, BHi;
  double t, mr;
  double (*phi) (double, void *);

  for(BHi = 0; (BHi < IMFp->N_notBH_ranges) && (m <= IMFp->notBH_ranges.inf[BHi]); BHi++)
    ;
  if((BHi < IMFp->N_notBH_ranges) && (m > IMFp->notBH_ranges.sup[BHi]))
    return 0;

  phi = IMFp->IMFfunc_byNum;

  if(m >= IMFp->Mm && m <= IMFp->MU)
    {
      if(getindex((double *) &SD.MArray[0], 0, SD.Mdim - 1, &m, &ir) < 0)
	{
	  mr = SD.Y[SD.Zbin * SD.Mdim + ir];
	  if(m < SD.MArray[0])
	    mr = SD.Y[SD.Zbin * SD.Mdim] * (m / SD.MArray[0]);
	}
      else
	{
	  t = (m - SD.MArray[ir]) / (SD.MArray[ir + 1] - SD.MArray[ir]);
	  mr = (1 - t) * SD.Y[SD.Zbin * SD.Mdim + ir] + t * SD.Y[SD.Zbin * SD.Mdim + ir + 1];
	}
      return mr * phi(m, p);
    }
  else
    return 0;
}

double INLINE_FUNC ztRSnII(double time, void *p)
{
  /* return the production rate in mass for elements at time t */

  int ir, BHi, zone;
  double m, t, u, y;

  /*
   * linear interpolation in 2 dim:
   *   y1 = y[j][k]
   *   y2 = y[j+1][k]
   *   y3 = y[j+1][k+1]
   *   y4 = y[j][k+1]
   *
   *   y(x1,x2) = (1-t)*(1-u)*y1 + t*(1-u)*y2 + t*u*y3 + (1-t)*u * y4
   *
   */


  if((time >= IMFp->inf_lifetime) &&	/* Mm Msun < m < MUP Msun */
     (time <= IMFp->sup_lifetime))
    {
      m = dying_mass(time);

      for(BHi = 0; (BHi < IMFp->N_notBH_ranges) && (m <= IMFp->notBH_ranges.inf[BHi]); BHi++)
	;
      if((BHi < IMFp->N_notBH_ranges) && (m > IMFp->notBH_ranges.sup[BHi]))
	return 0;


      for(ir = SD.Mdim - 1; m < SD.MArray[ir] && ir > 0; ir--)
	;

      if(SD.Zdim > 1)
	{
	  zone = 0;
	  zone |= (SD.Zstar >= SD.ZArray[SD.Zdim - 1]);
	  zone |= ((m >= SD.MArray[SD.Mdim - 1]) << 1);
	  zone |= ((SD.Zstar < SD.ZArray[0]) << 2);
	  zone |= ((m < SD.MArray[0]) << 3);
	}
      else
	{
	  zone = 16;
	  zone |= (m >= SD.MArray[SD.Mdim - 1]);
	  zone |= ((m < SD.MArray[0]) << 1);
	}

      switch (zone)
	{
	case 0:
	  t = (m - SD.MArray[ir]) / (SD.MArray[ir + 1] - SD.MArray[ir]);
	  u = (SD.Zstar - SD.ZArray[SD.Zbin]) / (SD.ZArray[SD.Zbin + 1] - SD.ZArray[SD.Zbin]);
	  y = (1 - t) * (1 - u) * SD.Y[SD.Zbin * SD.Mdim + ir] +
	    t * (1 - u) * SD.Y[SD.Zbin * SD.Mdim + ir + 1] +
	    t * u * SD.Y[(SD.Zbin + 1) * SD.Mdim + ir + 1] + (1 - t) * u * SD.Y[(SD.Zbin + 1) * SD.Mdim + ir];
	  break;
	case 2:
	case 8:
	  u = (SD.Zstar - SD.ZArray[SD.Zbin]) / (SD.ZArray[SD.Zbin + 1] - SD.ZArray[SD.Zbin]);
	  y = (1 - u) * SD.Y[SD.Zbin * SD.Mdim + ir] + u * SD.Y[(SD.Zbin + 1) * SD.Mdim + ir];
	  if(zone == 8)
	    y *= (m / SD.MArray[0]);
	  break;
	case 1:
	case 4:
	  t = (m - SD.MArray[ir]) / (SD.MArray[ir + 1] - SD.MArray[ir]);
	  y = (1 - t) * SD.Y[SD.Zbin * SD.Mdim + ir] + t * SD.Y[SD.Zbin * SD.Mdim + ir + 1];
	  break;
	case 9:
	case 12:
	  y = SD.Y[SD.Zbin * SD.Mdim + ir] * (m / SD.MArray[0]);
	  break;
	case 3:
	case 6:
	  y = SD.Y[SD.Zbin * SD.Mdim + ir];
	  break;
	case 16:
	  t = (m - SD.MArray[ir]) / (SD.MArray[ir + 1] - SD.MArray[ir]);
	  y = (1 - t) * SD.Y[SD.Zbin * SD.Mdim + ir] + t * SD.Y[SD.Zbin * SD.Mdim + ir + 1];
	  break;
	case 17:
	  if(SD.ExtrDir & 1 && !SD.ExtrDir & 1 << 31)
	    {
	      t = (m - SD.MArray[SD.Mdim - 1]) / (SD.ExtrMArray[0] - SD.MArray[SD.Mdim - 1]);
	      y = (1 - t) * SD.Y[SD.Zbin * SD.Mdim + SD.Mdim - 1] + t * SD.ExtrY[SD.ExtrZbin * SD.ExtrMdim];
	    }
	  else
	    y = SD.Y[SD.Zbin * SD.Mdim + SD.Mdim - 1];
	  break;
	case 18:
	  if(SD.ExtrDir & 1 && SD.ExtrDir & 1 << 31)
	    {
	      t = (m - SD.ExtrMArray[SD.ExtrMdim - 1]) / (SD.MArray[0] - SD.ExtrMArray[SD.ExtrMdim - 1]);
	      y = (1 - t) * SD.ExtrY[SD.ExtrZbin * SD.ExtrMdim + SD.ExtrMdim - 1] +
		t * SD.Y[SD.Zbin * SD.Mdim];
	    }
	  else
	    y = SD.Y[SD.Zbin * SD.Mdim] * (m / SD.MArray[0]);
	  break;
	}
      if(p != 0x0)
	{
	  printf(" [ %g %d %d %d %g %g %g %g %g\n", time, zone, SD.Zbin, ir, t, y, m, SD.MArray[ir],
		 SD.Y[SD.Zbin * SD.Mdim + ir]);
	  fflush(stdout);
	}

      return nRSnII(time, p) * y;
    }
  else
    return 0;
}





/* ******************************************************************************
 *
 * III   II
 *  II   II
 *   II II
 *    III
 * ..............................................................................
 *
 * General Routines
 *
 * ****************************************************************************** */

void initialize_star_lifetimes(void)
     /*- initialize stellar lifetimes for given range of star masses -*/
{
  int i;

  All.mean_lifetime = lifetime(All.Mup);
  All.sup_lifetime = lifetime(All.MBms);
  if(ThisTask == 0)
    {
      printf("\nstellar mean lifetime:   mean (%5.4g Msun) = %g Gyrs\n\n", All.Mup, All.mean_lifetime);
      fflush(stdout);
    }

  for(i = 0; i < IMFs_dim; i++)
    {
      IMFs[i].inf_lifetime = lifetime(IMFs[i].MU);
      if(IMFs[i].Mm > All.MBms)
	IMFs[i].sup_lifetime = lifetime(IMFs[i].Mm);
      else
	IMFs[i].sup_lifetime = All.sup_lifetime;
      if(ThisTask == 0)
	{
	  printf("\n"
		 "   IMF %3d      inf (%5.4g Msun) = %g Gyrs\n"
		 "                sup (%5.4g Msun) = %g Gyrs\n",
		 i, IMFs[i].MU, IMFs[i].inf_lifetime, IMFs[i].Mm, IMFs[i].sup_lifetime);
	  fflush(stdout);
	}
    }

  if(ThisTask == 0)
    {
      printf("\nstellar mean lifetime:   mean (%5.4g Msun) = %g Gyrs\n\n", All.Mup, All.mean_lifetime);
      fflush(stdout);
    }
}


double INLINE_FUNC dm_dt(double m, double t)
{
  /* t is in Gyr */
  /* the last factor 1/agefact normalize in 1/yr: otherwise
     the results would be 1/Gyr */
  /*
     if(t > 0.0302233)
     return -0.37037 * m / (t - 0.012);
     else
     return -0.54054 * m / (t - 0.003);
   */

#ifdef LT_PM_LIFETIMES
  /* padovani & matteucci 1993 */
  if(t > 0.039765318659064693)
    return -m / t * (1.338 - 0.1116 * (9 + log10(t)));
  else
    return -0.54054 * m / (t - 0.003);
#endif

#ifdef LT_MM_LIFETIMES
  /* maeder & meynet 1989 */
  if(m <= 1.3)
    return -m / t / 0.6545;
  if(m > 1.3 && m <= 3)
    return -m / t / 3.7;
  if(m > 3 && m <= 7)
    return -m / t / 2.51;
  if(m > 7 && m <= 15)
    return -m / t / 1.78;
  if(m > 15 && m <= 53.054)
    return -m / t / 0.86;
  if(m > 53.054)
    return -0.54054054054 * m / (t - 0.003);
#endif

}

double INLINE_FUNC lifetime(double mass)
{
  /* calculates lifetime for a given mass  */
  /*                                       */
  /* padovani & matteucci (1993) approach  */
  /* move to gibson (1997) one             */
  /*                                       */
  /* mass is intended in solar units, life */
  /* time in Gyr                           */

  if(mass > 100)
    return 0;
  else if(mass < 0.6)
    return 160;			/* should be INF */


  /* padovani & matteucci 1993 */
#ifdef LT_PM_LIFETIMES
  if(mass <= 6.6)
    return pow(10, ((1.338 - sqrt(1.790 - 0.2232 * (7.764 - log10(mass)))) / 0.1116) - 9);
  else
    return 1.2 * pow(mass, -1.85) + 0.003;
#endif

  /*
     if(mass <= 8)
     return 5 * pow(mass, -2.7) + 0.012;
     else
     return 1.2 * pow(mass, -1.85) + 0.003;
   */

#ifdef LT_MM_LIFETIMES
  /* maeder & meynet 1989 */
  if(mass <= 1.3)
    return pow(10, -0.6545 * log10(mass) + 1);
  if(mass > 1.3 && mass <= 3)
    return pow(10, -3.7 * log10(mass) + 1.35);
  if(mass > 3 && mass <= 7)
    return pow(10, -2.51 * log10(mass) + 0.77);
  if(mass > 7 && mass <= 15)
    return pow(10, -1.78 * log10(mass) + 0.17);
  if(mass > 15 && mass <= 53.054)
    return pow(10, -0.86 * log10(mass) - 0.94);
  if(mass > 53.054)
    return 1.2 * pow(mass, -1.85) + 0.003;
#endif
}

double INLINE_FUNC dying_mass(double time)
{
  /* calculates mass dying at some time */
  /*                                    */
  /* time is time_in_Gyr                */

  if((time < All.inf_lifetime) || (time > All.sup_lifetime))
    return 0;

  /*
     if(time > 0.0302233)
     return pow((time - 0.012) / 5, -0.37037);
     else
     return pow((time - 0.003) / 1.2, -0.54054);
   */

#ifdef LT_PM_LIFETIMES
  /* padovani & matteucci 1993 */
  if(time > 0.039765318659064693)
    return pow(10, 7.764 - (1.79 - pow(1.338 - 0.1116 * (9 + log10(time)), 2)) / 0.2232);
  else
    return pow((time - 0.003) / 1.2, -1 / 1.85);
#endif

#ifdef LT_MM_LIFETIMES
  /* maeder & meynet 1989 */
  if(time >= 8.4221714076)
    return pow(10, (1 - log10(time)) / 0.6545);
  if(time < 8.4221714076 && time >= 0.38428316376)
    return pow(10, (1.35 - log10(time)) / 3.7);
  if(time < 0.38428316376 && time >= 0.044545508363)
    return pow(10, (0.77 - log10(time)) / 2.51);
  if(time < 0.044545508363 && time >= 0.01192772338)
    return pow(10, (0.17 - log10(time)) / 1.78);
  if(time < 0.01192772338 && time >= 0.0037734864318)
    return pow(10, -(0.94 + log10(time)) / 0.86);
  if(time < 0.0037734864318)
    return pow((time - 0.003) / 1.2, -0.54054);
#endif
}


double INLINE_FUNC sec_dist(double gamma, double mu)
{
  /* calculates secondary distribution function */
  /* for Sn type Ia                             */
  /* as far, we take gamma=2, so this function  */
  /* isn't used                                 */

  return pow(2, 1 + gamma) * (1 + gamma) * pow(mu, gamma);
}


/* double INLINE_FUNC integrand0(double m, void *params) */
/*      /\*- used by integrate_imf to integrate imf by mass -*\/ */
/* { */
/*   return m * phi(m, params); */
/* } */

/* double integrate_imf(double inf, double sup, int mode, void *params) */
/*      /\* integrate imf. */
/*       * mode = 0 means integration by mass */
/*       * mode = 1 means integration by number */
/*       *\/ */
/* { */
/*   double result, err; */

/*   if(mode == 0) */
/*     F.function = &integrand0; */
/*   else */
/*     F.function = &phi; */
/*   F.params = params; */

/*   if( (my_gslstatus = gsl_integration_qag(&F, inf, sup, 1e-6, 1e-4, gsl_ws_dim, qag_INT_KEY, w, &result, &err)) ) */
/*     { */
/*       if(ThisTask == 0) */
/*      printf("  >> Task %i, gsl integration error %i in integrate_imf [%.6g - %.6g] : %.6g %.6g\n", */
/*             ThisTask, my_gslstatus, inf, sup, result, err); */
/*       fflush(stdout); */
/*       endrun(999672); */
/*     } */

/*   return result; */
/* } */


/* double INLINE_FUNC normalize_imf(double inf, double sup, void *param) */
/* { */
/*   return (All.Aphi = (All.Xphi - 1) / (1 - pow(inf / sup, All.Xphi - 1)) * pow(inf, All.Xphi - 1)); */
/* } */

/* double INLINE_FUNC phi(double m, void *params) */
/* { */
/*   return All.Aphi * pow(m, -(1 + All.Xphi)); */
/* } */


#ifdef LT_SEv_INFO

void write_metallicity_stat(void)
{
  int i, count, count_b;

#ifdef LT_SEvDbg
  get_metals_sumcheck(1);
#endif

  MPI_Reduce(&ALdata_min[0], &tot_ALdata_min[0], AL_INUM_min, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&ALdata_max[0], &tot_ALdata_max[0], AL_INUM_max, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&ALdata_sum[0], &tot_ALdata_sum[0], AL_INUM_sum, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(tot_starsnum > 1)
    {
      if(tot_starsnum > 3)
	{
	  tot_ALdata_sum[MEAN_sl] -= tot_ALdata_max[MAX_sl] + tot_ALdata_min[MIN_sl];
	  tot_ALdata_sum[MEAN_ngb] -= tot_ALdata_max[MAX_ngb] + tot_ALdata_min[MIN_ngb];
	  count = tot_starsnum - 2;
	}
      else
	count = tot_starsnum;
      tot_ALdata_sum[MEAN_sl] /= count;
      tot_ALdata_sum[MEAN_ngb] /= count;
    }

  get_metallicity_stat();

  MPI_Reduce(&Stat_sum[0], &tot_Stat_sum[0], S_INUM_sum, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&Stat_min[0], &tot_Stat_min[0], S_INUM_min, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&Stat_max[0], &tot_Stat_max[0], S_INUM_max, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      /* metallicities */
      if(tot_Stat_sum[NUM_star] > 3)
	{
	  count = tot_Stat_sum[NUM_star];
	  count_b = tot_Stat_sum[NUM_star] - 2;
	  tot_Stat_sum[MEAN_Zstar] -= tot_Stat_max[MAX_Zstar] + tot_Stat_min[MIN_Zstar];
	}
      else if(tot_Stat_sum[NUM_star] == 0)
	count = count_b = 1;

      for(i = 0; i < LT_NMet; i++)
	{
	  tot_Stat_sum[MEAN_Zs + i] /= All.TotN_gas;
	  tot_Stat_sum[MEAN_Zsstar + i] /= count;
	}
      tot_Stat_sum[MEAN_Z] -= tot_Stat_max[MAX_Z] + tot_Stat_min[MIN_Z];
      tot_Stat_sum[MEAN_Z] /= All.TotN_gas - 2;	/* assume that we have more than 3 gas particles! */
      tot_Stat_sum[MEAN_Zstar] /= count_b;

      /* energy ratios */
      if(tot_Stat_sum[NUM_egyf] > 3)
	{
	  count = tot_Stat_sum[NUM_egyf];
	  count_b = tot_Stat_sum[NUM_egyf] - 2;
	  tot_Stat_sum[MEAN_egyf] -= tot_Stat_max[MAX_egyf] + tot_Stat_min[MIN_egyf];
	}
      else if(tot_Stat_sum[NUM_egyf] == 0)
	count = count_b = 1;
      tot_Stat_sum[MEAN_egyf] /= count_b;

      fprintf(FdMetals, "%9.7g ", All.Time);
      fprintf(FdMetals, "%9.7g %9.7g %9.7g %9.7g %9.7g %9.7g ",
	      tot_Stat_min[MIN_Z], tot_Stat_max[MAX_Z], tot_Stat_sum[MEAN_Z],
	      tot_Stat_min[MIN_Zstar], tot_Stat_max[MAX_Zstar], tot_Stat_sum[MEAN_Zstar]);
      for(i = 0; i < LT_NMet; i++)
	fprintf(FdMetals, "%9.7g %9.7g ", tot_Stat_sum[MEAN_Zs + i], tot_Stat_sum[MEAN_Zsstar + i]);
      fprintf(FdMetals, "%9.7g %9.7g %9.7g ", tot_Stat_min[MIN_egyf], tot_Stat_max[MAX_egyf],
	      tot_Stat_sum[MEAN_egyf]);
      fprintf(FdMetals, "%9.7g %9.7g %9.7g %9.7g %9.7g %9.7g\n", tot_ALdata_min[MIN_sl],
	      tot_ALdata_max[MAX_sl], tot_ALdata_sum[MEAN_sl], tot_ALdata_min[MIN_ngb],
	      tot_ALdata_max[MAX_ngb], tot_ALdata_sum[MEAN_ngb]);
      fflush(FdMetals);
    }
  MPI_Barrier(MPI_COMM_WORLD);
}

void get_metallicity_stat(void)
{
  int i, j;
  double Z, Ztmp, zmass, R;


  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 0)
	{
	  /* find total metal mass (exclude Helium) */
	  zmass = 0;
	  for(j = 0, Z = 0; j < LT_NMetP; j++)
	    {
	      if(j != Hel)
		zmass += SphP[i].Metals[j];
	    }
	  /* metallicities are calculated as metal_mass / hydrogen_mass */
	  R = 1.0 / (P[i].Mass - Z - SphP[i].Metals[Hel]);
	  for(j = 0, Z = 0; j < LT_NMetP; j++)
	    Stat_sum[MEAN_Zs + j] += SphP[i].Metals[j] * R;

	  Stat_sum[MEAN_Z] += (Z = zmass * R);
	  if(Z < Stat_min[MIN_Z])
	    Stat_min[MIN_Z] = Z;
	  else if(Z > Stat_max[MAX_Z])
	    Stat_max[MAX_Z] = Z;
	}
      else
	{
	  /* find total metal mass (exclude Helium) */
	  zmass = 0;
	  for(j = 0, Z = 0; j < LT_NMetP; j++)
	    {
	      if(j != Hel)
		zmass += MetP[P[i].MetID].Metals[j];
	    }
	  /* metallicities are calculated as metal_mass / hydrogen_mass */
	  R = 1.0 / (P[i].Mass - Z - MetP[P[i].MetID].Metals[Hel]);
	  for(j = 0, Z = 0; j < LT_NMetP; j++)
	    Stat_sum[MEAN_Zs + j] += MetP[P[i].MetID].Metals[j] * R;

	  Stat_sum[MEAN_Z] += (Z = zmass * R);
	  if(Z < Stat_min[MIN_Zstar])
	    Stat_min[MIN_Zstar] = Z;
	  else if(Z > Stat_max[MAX_Zstar])
	    Stat_max[MAX_Zstar] = Z;
	}
    }

  return;
}
#endif

#ifdef LT_SEvDbg
void get_metals_sumcheck(int mode)
{
  FILE *outstream;

#define star (LT_NMetP-1)
#define sum_gas (2 * star)
#define sum_star sum_gas + 1
#define sum sum_star + 1

  int i, j;
  double metals[2 * (LT_NMetP - 1) + 3], tot_metals[2 * (LT_NMetP - 1) + 3];

  for(i = 0; i < 2 * (LT_NMetP - 1) + 3; i++)
    metals[i] = tot_metals[i] = 0;

  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 0)
	{
	  for(j = 1; j < LT_NMetP; j++)
	    {
	      metals[j - 1] += SphP[i].Metals[j];
	      metals[sum_gas] += SphP[i].Metals[j];
	      metals[sum] += SphP[i].Metals[j];
	    }
	}
      else if(P[i].Type == 4)
	{
	  for(j = 1; j < LT_NMetP; j++)
	    {
	      metals[star + j - 1] += MetP[P[i].MetID].Metals[j];
	      metals[sum_star] += MetP[P[i].MetID].Metals[j];
	      metals[sum] += MetP[P[i].MetID].Metals[j];
	    }
	}
    }

  MPI_Reduce(&metals[0], &tot_metals[0], 2 * (LT_NMetP - 1) + 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      if(mode & 1)
	{
	  outstream = FdMetSumCheck;
	  mode &= ~1;
	}
      else
	outstream = stdout;

      fprintf(outstream, "%c %7.6g ", (mode) ? ((mode & 2) ? '>' : ((mode & 4) ? '<' : '@')) : ' ', All.Time);
      for(i = 0; i < 2 * (LT_NMetP - 1) + 3; i++)
	fprintf(outstream, "%9.7g ", tot_metals[i] / All.HubbleParam);
      fprintf(outstream, "\n");
    }

}
#endif


int calculate_effective_yields(double inf, double sup, int IMFi)
{
  double abserr, result;
  int i;
  int nonZero = 0;

  F.function = &zmRSnII;

  SD.ZArray = IIZbins[SD.Yset];
  SD.MArray = IIMbins[SD.Yset];
  SD.Mdim = IIMbins_dim[SD.Yset];
  SD.Zdim = IIZbins_dim[SD.Yset];

  for(SD.Zbin = 0; SD.Zbin < IIZbins_dim[SD.Yset]; SD.Zbin++)
    {
      if(IIZbins_dim[SD.Yset] > 1)
	SD.Zstar = IIZbins[SD.Yset][SD.Zbin];
      else
	SD.Zstar = 0;
      for(i = 0; i < LT_NMet; i++)
	{
	  SD.Y = SnIIYields[SD.Yset][i];
	  if((my_gslstatus =
	      gsl_integration_qag(&F, inf, sup, 1e-6, 1e-4, gsl_ws_dim, qag_INT_KEY, w, &result, &abserr)))
	    {
	      if(ThisTask == 0)
		printf("  >> Task %i, gsl integration error %i in calculating effective yields"
		       " [%9.7g - %9.7g] : %9.7g %9.7g\n", ThisTask, my_gslstatus, inf, sup, result, abserr);
	      fflush(stdout);
	      endrun(LT_ERROR_INTEGRATION_ERROR);
	    }

	  if((SnII_ShortLiv_Yields[IMFi][i][SD.Zbin] = result) > 0)
	    nonZero++;
	}

#ifdef UM_COOLING
      if(nonZero)
	{
	  SD.Y = SnII_Fill_mol_num[SD.Yset];
	  if((my_gslstatus =
	      gsl_integration_qag(&F, inf, sup, 1e-6, 1e-4, gsl_ws_dim, qag_INT_KEY, w, &result, &abserr)))
	    {
	      if(ThisTask == 0)
		printf("  >> Task %i, gsl integration error %i in calculating effective yields"
		       " [%9.7g - %9.7g] : %9.7g %9.7g\n", ThisTask, my_gslstatus, inf, sup, result, abserr);
	      fflush(stdout);
	      endrun(LT_ERROR_INTEGRATION_ERROR);
	    }
	  SnII_ShortLiv_FillMolNum[IMFi][SD.Zbin] = result;
	}
#endif
    }

  if(ThisTask == 0)
    printf("\n");

  return nonZero;
}

double calculate_FactorSN(double m_inf, double m_sup, void *params)
{
  double abserr, result;

  F.function = &ejectaSnII;
  F.params = params;

  if((my_gslstatus =
      gsl_integration_qag(&F, m_inf, m_sup, 1e-4, 1e-3, gsl_ws_dim, qag_INT_KEY, w, &result, &abserr)))
    {
      if(ThisTask == 0)
	printf
	  ("  >> Task %i, gsl integration error %i in calculating FactorSN [%9.7g - %9.7g] : %9.7g %9.7g\n",
	   ThisTask, my_gslstatus, m_inf, m_sup, result, abserr);
      fflush(stdout);
      endrun(LT_ERROR_INTEGRATION_ERROR);
    }
  return result;
}



void init_SN(void)
{
  int i, j, guess_on_Nsteps;
  double astep, meanweight;
  char buf[200], mode[2];

  noffset = (int *) mymalloc(sizeof(int) * NTask);
  nbuffer = (int *) mymalloc(sizeof(int) * NTask);
  nsend_local = (int *) mymalloc(sizeof(int) * NTask);
  nsend = (int *) mymalloc(sizeof(int) * NTask * NTask);

  ndone = &ndoneoff[0];
  noff = &ndoneoff[1];
  ntotdone = &ntotdoneoff[0];
  ntotoff = &ntotdoneoff[1];

  if(ThisTask == 0)
    {
      if(RestartFlag == 0)
	strcpy(mode, "w");
      else
	strcpy(mode, "a");
      sprintf(buf, "%s%s", All.OutputDir, "sn_init.txt");
      if(!(FdSnInit = fopen(buf, mode)))
	{
	  printf("error in opening file '%s'\n", buf);
	  endrun(1);
	}
      if(mode == "a")
	fprintf(FdSnInit, "========================================\n" "restarting from a= %g\n\n", All.Time);
      sprintf(buf, "%s%s", All.OutputDir, "sn_details.dat");
      if(!(FdSnDetails = fopen(buf, mode)))
	{
	  printf("error in opening file '%s'\n", buf);
	  endrun(1);
	}
      sprintf(buf, "%s%s", All.OutputDir, "warnings.txt");
      if(!(FdWarn = fopen(buf, mode)))
	{
	  printf("error in opening file '%s'\n", buf);
	  endrun(1);
	}
    }

  /* temporary */
/*   init_T_TC_runs(); */
  /* */

  UnitMassFact = All.UnitMass_in_g / SOLAR_MASS;

#ifdef LT_TRACK_CONTRIBUTES
  init_packing();
#endif

  /* ...........: fix the range for the desired number in neighbours
     finding for metal spreading :...................... */

  if(All.NeighInfNum == 0)
    All.NeighInfNum = 1;
  All.LeftNumNgbSN = All.DesNumNgbSN - All.SpreadNumNgbDev;
  All.RightNumNgbSN = All.DesNumNgbSN + All.SpreadNumNgbDev;
  All.SpreadNeighCoeff = (double) All.DesNumNgbSN / All.DesNumNgb;

  /* ...........: allocate workspace statically :...................... */
  w = gsl_integration_workspace_alloc(gsl_ws_dim);
  old_error_handler = gsl_set_error_handler(&fsolver_error_handler);


  /*  I M F StartUp  */
  /* --------------- */

  read_imfs();

  for(i = 0; i < IMFs_dim; i++)
    /* normalize IMF by mass */
    {
      if(IMFs[i].NParams > 0 || IMFs[i].timedep)
	IMFs[i].getp(i, IMFs[i].Params, All.Time);	/* get parameters for a (time-dependent) IMF */
      IMFs[i].A = 1;
      IMFs[i].A = 1 / IntegrateIMF_byMass(IMFs[i].Mm, IMFs[i].MU, &IMFs[i], INC_BH);
    }

  if(ThisTask == 0)
    {
      for(i = 0; i < IMFs_dim; i++)
	write_IMF_info(i, stdout);
      for(i = 0; i < IMFs_dim; i++)
	write_IMF_info(i, FdSnInit);
    }

  /* ...........: initialize time scales for stellar evolutions :...................... */
  initialize_star_lifetimes();

  if(ThisTask == 0)
    {
      fprintf(FdSnInit, "[stellar evolution initialization - lifetimes]");
      fprintf(FdSnInit, "\nstellar lifetimes: mean (%5.4g Msun) = %g Gyrs\n", All.Mup, All.mean_lifetime);
      for(i = 0; i < IMFs_dim; i++)
	fprintf(FdSnInit, "   IMF %3d         inf (%5.4g Msun) = %g Gyrs\n"
		"                   sup (%5.4g Msun) = %g Gyrs\n",
		i, IMFs[i].MU, IMFs[i].inf_lifetime, min(IMFs[i].Mm, All.MBms), IMFs[i].sup_lifetime);
    }


  /*  metals StartUp  */
  /* ---------------- */

  read_metals();


  /*  Sn StartUp  */
  /* ------------ */

  /* ...........: initialize stepping :...................... */

#ifdef LT_SNII
  SnII_steps = (double ***) malloc(IMFs_dim * sizeof(double **));
  SnII_Nsteps = (int *) malloc(IMFs_dim * sizeof(int));
  guess_on_Nsteps = (int) (1.3 / All.SnII_Step_Prec);

  for(i = 0; i < IMFs_dim; i++)
    if(!IMFs[i].timedep)
      {
	SnII_steps[i] = (double **) calloc(2, sizeof(double *));
	SnII_Nsteps[i] = guess_on_Nsteps;
	for(j = 0; j < 2; j++)
	  SnII_steps[i][j] = (double *) calloc(guess_on_Nsteps, sizeof(double));
      }
#endif

#if defined(LT_SNIa) || defined(LT_AGB)
  LLv_steps = (double ***) malloc(IMFs_dim * sizeof(double **));
  LLv_Nsteps = (int *) malloc(IMFs_dim * sizeof(int));
  guess_on_Nsteps = (int) (1.3 / All.LLv_Step_Prec);

  for(i = 0; i < IMFs_dim; i++)
    if(!IMFs[i].timedep)
      {
	LLv_steps[i] = (double **) calloc(2, sizeof(double *));
	LLv_Nsteps[i] = guess_on_Nsteps;
	for(j = 0; j < 2; j++)
	  LLv_steps[i][j] = (double *) calloc(guess_on_Nsteps, sizeof(double));
	IMFp = &IMFs[i];
	build_SN_Stepping(1, i);
      }
#endif

  /* ...........: reading yields :...................... */

#ifdef LT_SNII
  if(ThisTask == 0)
    printf(":: starting SnII..\n");
  read_SnII_yields();
#endif


#ifdef LT_SNIa
  if(ThisTask == 0)
    printf(":: starting SnIa..\n");
  read_SnIa_yields();
#endif


#ifdef LT_AGB
  if(ThisTask == 0)
    printf(":: starting AGB..\n");
  read_AGB_yields();
#endif

  /*  SF StartUp  */
  /* ------------ */

  All.Time = All.TimeBegin;
  InitCool();

  if(All.SFTh_Zdep > 0)
    for(i = 0; i < ZBins; i++)
      ThInst_onset[i] = pow(10, ThInst_onset[i]);

  /* Set some SF-related Units (old set_units_sfr) */
  meanweight = 4 / (1 + 3 * HYDROGEN_MASSFRAC);	/* note: assuming NEUTRAL GAS */

  All.EgySpecCold = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.TempClouds;
  All.EgySpecCold *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

  meanweight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));	/* note: assuming FULL ionization */

  All.OverDensThresh =
    All.CritOverDensity * All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

#ifdef INTERNAL_CRIT_DENSITY
  All.PhysDensThresh = All.CritPhysDensity;
#else
  All.PhysDensThresh = All.CritPhysDensity * PROTONMASS / HYDROGEN_MASSFRAC / All.UnitDensity_in_cgs /
    (All.HubbleParam * All.HubbleParam);
#endif
  /* ----------- */

  SnII_ShortLiv_Yields = (double ***) malloc(IMFs_dim * sizeof(double **));
#ifdef UM_COOLING
  SnII_ShortLiv_FillMolNum = (double **) malloc(IMFs_dim * sizeof(double *));
#endif

  for(i = 0; i < IMFs_dim; i++)
    {
      if(All.PhysDensThresh > 0)
	IMFs[i].PhysDensThresh[0] = All.PhysDensThresh;
      else
	/* Phys Dens Thresh is supposed to be in hydrogen number density cm^-3 */
	IMFs[i].PhysDensThresh[0] *= PROTONMASS / HYDROGEN_MASSFRAC / All.UnitDensity_in_cgs /
	  (All.HubbleParam * All.HubbleParam);
      if(All.egyIRA_ThMass > 0)
	IMFs[i].egyShortLiv_MassTh = All.egyIRA_ThMass;
      if(All.metIRA_ThMass > 0)
	IMFs[i].metShortLiv_MassTh = All.metIRA_ThMass;
      if(All.WindEfficiency > 0)
	IMFs[i].WindEfficiency = All.WindEfficiency;
      if(All.WindEnergyFraction > 0)
	IMFs[i].WindEnergyFraction = All.WindEnergyFraction;

      IMFp = &IMFs[i];
      setup_SF_related(i);

      /* ...........: build the timestepping for Sn :................... */

      IMFs[i].ShortLiv_MassTh = max(IMFs[i].metShortLiv_MassTh, IMFs[i].egyShortLiv_MassTh);
      IMFs[i].ShortLiv_TimeTh = min(IMFs[i].metShortLiv_TimeTh, IMFs[i].egyShortLiv_TimeTh);


      build_SN_Stepping(2, i);

      /* ...........: write down :................... */

      if(ThisTask == 0)
	{
	  fprintf(FdSnInit, "\n[time stepping]\ntime (Gyr)      timestep (Gyr)\n");
	  for(j = 0; j < SnII_Nsteps[i]; j++)
	    fprintf(FdSnInit, "%8.5lg        %8.6lg\n", SnII_steps[i][0][j], SnII_steps[i][1][j]);
#if defined(LT_SNIa) || defined(LT_AGB)
	  for(j = 0; j < LLv_Nsteps[i]; j++)
	    fprintf(FdSnInit, "%8.5lg        %8.6lg\n", LLv_steps[i][0][j], LLv_steps[i][1][j]);
#endif
	  fflush(FdSnInit);
	}
    }

  /* ...........: build arrays to interpolate a vs t :................... */

  accel = gsl_interp_accel_alloc();
  spline = gsl_spline_alloc(gsl_interp_cspline, TIME_INTERP_SIZE);

  astep = pow(10, log10(1.0 / 0.001) / TIME_INTERP_SIZE);
  aarray[0] = 0.001;
  tarray[0] = 0;
  cosmic_time = get_age(0.001);

  for(i = 1; i < TIME_INTERP_SIZE; i++)
    {
      aarray[i] = aarray[i - 1] * astep;
      tarray[i] = cosmic_time - get_age(aarray[i]);
    }
  gsl_spline_init(spline, tarray, aarray, TIME_INTERP_SIZE);

  /* */

  if(ThisTask == 0)
    fflush(stdout);

  /* ....... */

#ifdef LT_SEvDbg
  do_spread_dbg_list = (int *) malloc(sizeof(int) * NTask);
#endif

  return;
}



void build_SN_Stepping(int SNType, int IMFi)
{
  gsl_function G;

  double f0, Step_Prec;
  double time, end_time, delta_time, timed_frac, Left, Right, m;
  double err;
  double next_notBH_sup, next_notBH_inf;	/* these vars will be set to the sup and inf lifetimes
						 * corresponding to the boundaries of the not-BH mass
						 * ranges in the current IMF.
						 */
  int i, limit, BHi;

  int *Nsteps;
  double *times, *delta_times;

#if defined(LT_SNIa) || defined(LT_AGB)
  if(SNType == 1)
    {
      /* set-up things for SnIa and AGB stars */

      /* calculates the total number of expected stars below 8 Solar masses
         in the adopted model this is related to the total number of SnIa
       */
      f0 = IntegrateIMF_byNum(All.MBms, All.Mup, &IMFs[IMFi], INC_BH);

      if(f0 > 0)
	{
	  /* set the final time */
	  end_time = IMFs[IMFi].sup_lifetime;
	  /* set the precision of the stepping */
	  Step_Prec = All.LLv_Step_Prec;
	  /* guess of the needed steps */
	  Nsteps = &LLv_Nsteps[IMFi];
	  /* set pointers */
	  times = LLv_steps[IMFi][0];
	  delta_times = LLv_steps[IMFi][1];

	  /* first time */
	  times[0] = 0;
	  /* first delta_t */
	  delta_times[0] = All.mean_lifetime;
	  /* start time */
	  time = All.mean_lifetime;
	}
      else
	{
	  LLv_Nsteps[IMFi] = 0;
	  LLv_steps[IMFi][0][0] = 0;
	  LLv_steps[IMFi][1][0] = 0;
	  return;
	}
    }
#endif
  if(SNType == 2)
    {
      /* set-up things for SnII */
      /* note: see comments in the SnIa block above */

      f0 = IntegrateIMF_byNum(max(IMFs[IMFi].Mm, All.Mup), IMFs[IMFi].MU, &IMFs[IMFi], INC_BH);

      if((f0 > 0)
	 && (max(max(IMFs[IMFi].Mm, All.Mup), IMFs[IMFi].notBH_ranges.inf[IMFs[IMFi].N_notBH_ranges - 1]) <
	     IMFs[IMFi].ShortLiv_MassTh))
	{
	  /* note: see comments in the SnIa block above */
	  end_time = min(All.mean_lifetime, IMFs[IMFi].sup_lifetime);
	  Step_Prec = All.SnII_Step_Prec;
	  Nsteps = &SnII_Nsteps[IMFi];
	  times = SnII_steps[IMFi][0];
	  delta_times = SnII_steps[IMFi][1];

	  times[0] = 0;
	  /* note that ShortLiv_TimeTh is set to inf_lifetime in case no threshold is used */
	  delta_times[0] = IMFs[IMFi].ShortLiv_TimeTh;
	  time = IMFs[IMFi].ShortLiv_TimeTh;
	}
      else
	{
	  SnII_Nsteps[IMFi] = 0;
	  SnII_steps[IMFi][0][0] = 0;
	  SnII_steps[IMFi][1][0] = 0;
	  return;
	}
    }

  for(m = dying_mass(time), BHi = 0;
      (BHi < IMFs[IMFi].N_notBH_ranges) && (m <= IMFs[IMFi].notBH_ranges.inf[BHi]); BHi++)
    ;
  next_notBH_sup = lifetime(IMFs[IMFi].notBH_ranges.inf[BHi]);
  next_notBH_inf = lifetime(IMFs[IMFi].notBH_ranges.sup[BHi]);

  i = 1;
  /* initial guess for delta_time */
  modf(log10(time) + 9, &delta_time);
  delta_time = pow(10, delta_time) / 1e9;

  /* cycle to calculate time steps */
  while(time < end_time)
    {
      timed_frac = 1;
      Left = Right = limit = 0;
      while((timed_frac / f0 < Step_Prec * 0.9 || timed_frac / f0 > Step_Prec * 1.1) && !limit)
	{
	  if(time + delta_time > end_time)
	    delta_time = end_time - time;
	  if(time + delta_time > next_notBH_sup)
	    delta_time = next_notBH_sup - time;
	  timed_frac =
	    IntegrateIMF_byNum(dying_mass(time + delta_time), dying_mass(time), &IMFs[IMFi], INC_BH);

	  if(timed_frac / f0 < Step_Prec * 0.9)
	    {
	      if(!(limit = (time + delta_time == end_time)))
		{
		  Left = max(Left, delta_time);
		  if(Right == 0)
		    delta_time *= 2;
		  else
		    delta_time = (Left + Right) / 2;
		}
	    }
	  if(timed_frac / f0 > Step_Prec * 1.1)
	    {
	      if(Right == 0)
		Right = delta_time;
	      else
		Right = min(Right, delta_time);
	      if(Left == 0)
		delta_time /= 2;
	      else
		delta_time = (Left + Right) / 2;
	    }
	}
      if(timed_frac / f0 < Step_Prec / 2 && limit)
	delta_times[i - 1] += delta_time;
      else
	{
	  times[i] = time;
	  delta_times[i] = delta_time;
	  i++;
	}
      if((time += delta_time) == next_notBH_sup)
	{
	  if(BHi < IMFs[IMFi].N_notBH_ranges)
	    {
	      BHi++;
	      delta_time = lifetime(IMFs[IMFi].notBH_ranges.sup[BHi]) - next_notBH_sup;
	      next_notBH_sup = lifetime(IMFs[IMFi].notBH_ranges.inf[BHi]);
	      next_notBH_inf = lifetime(IMFs[IMFi].notBH_ranges.sup[BHi]);
	    }
	  else
	    {
	      limit = 1;
	      time = end_time;
	    }
	}

    }

  *Nsteps = i + 1;
  if(SNType == 1)
    times[i] = IMFs[IMFi].sup_lifetime;
  else
    times[i] = All.mean_lifetime;
  delta_times[i] = 0;
}




void get_Egy_and_Beta(double ThMass, double *Egy, double *Beta, IMF_Type * IMF)
{
  double NumFrac_inIRA;

  *Beta = calculate_FactorSN(ThMass, IMF->MU, imf_pointer);
  NumFrac_inIRA = IntegrateIMF_byMass(ThMass, IMF->MU, IMF, EXC_BH);
  *Egy = IntegrateIMF_byEgy(IMFp->egyShortLiv_MassTh, IMFp->MU, IMFp, EXC_BH) / SOLAR_MASS;
  /**Egy = (All.SnIIEgy * All.UnitEnergy_in_cgs) * NumFrac_inIRA / SOLAR_MASS;*/
  *Egy *= All.UnitMass_in_g;
  *Egy *= (1 - *Beta) / *Beta;

  return;
}

void calculate_ShortLiving_related(IMF_Type * IMFp, int loud)
{
  double num_snII;
  int set;

  if(ThisTask == 0 && loud)
    printf(":: energy IRA is active in the range [%9.7g - %9.7g]Msun <> [%9.7g - %9.7g]Gyr\n",
	   IMFp->MU, IMFp->egyShortLiv_MassTh, IMFp->inf_lifetime, IMFp->egyShortLiv_TimeTh);

  /* mass fraction involved in SnII */
  IMFp->MassFrac_inIRA = IntegrateIMF_byMass(IMFp->egyShortLiv_MassTh, IMFp->MU, IMFp, INC_BH);
  /* expected number of snII per initial stellar population mass in units of solar masses IN IRA RANGE */
  IMFp->NumFrac_inIRA = IntegrateIMF_byNum(IMFp->egyShortLiv_MassTh, IMFp->MU, IMFp, EXC_BH);
  /* expected number of stars per initial stellar population mass in units of solar masses */
  num_snII = IntegrateIMF_byNum(max(All.Mup, IMFp->Mm), IMFp->MU, IMFp, EXC_BH);
  /* erg/g per solar masses of initial stellar population mass, due to IRA snII */
  /*  ~7.4x10^48 erg/(s * Msun) for all Sn (8-100Msun) for a Salp. IMF */
  IMFp->IRA_erg_per_g =
    IntegrateIMF_byEgy(IMFp->egyShortLiv_MassTh, IMFp->MU, IMFp, EXC_BH) * All.UnitEnergy_in_cgs / SOLAR_MASS;
  /*  IMFp->IRA_erg_per_g = (All.SnIIEgy * All.UnitEnergy_in_cgs) * IMFp->NumFrac_inIRA / SOLAR_MASS; */
  /*  erg/g per solar masses of initial stellar population mass, due to ALL IRA snII */
  IMFp->TOT_erg_per_g =
    IntegrateIMF_byEgy(max(All.Mup, IMFp->Mm), IMFp->MU, IMFp, EXC_BH) * All.UnitEnergy_in_cgs / SOLAR_MASS;
  /*  IMFp->TOT_erg_per_g = (All.SnIIEgy * All.UnitEnergy_in_cgs) * num_snII / SOLAR_MASS; */

  set = IMFp->YSet;

/* LT_SNII must be defined for the following lines */
  SD.Zbin = 0;
  SD.MArray = IIMbins[set];
  SD.Mdim = IIMbins_dim[set];
  SD.Y = SnIIEj[set];

  IMFp->FactorSN = calculate_FactorSN(IMFp->egyShortLiv_MassTh, IMFp->MU, IMFp);	/* restored by Short-Living stars */
  IMFp->totFactorSN = calculate_FactorSN(All.Mup, IMFp->MU, IMFp);

#ifdef LT_AGB
  SD.Zbin = 0;
  SD.MArray = AGBMbins[set];
  SD.Mdim = AGBMbins_dim[set];
  SD.Y = AGBEj[set];
  IMFp->totResFrac = IMFp->totFactorSN + calculate_FactorSN(IMFp->Mm, All.Mup, IMFp);
#else
  IMFp->totResFrac = IMFp->totFactorSN;
#endif

  if(ThisTask == 0 && loud)
    {
      printf("    %.3lg%% of mass in egy Short-Living tail\n"
	     "    energy per g due to Short-Living snII is %9.7g erg/g\n",
	     IMFp->MassFrac_inIRA * 100, IMFp->IRA_erg_per_g);

      printf("   restored fraction from short-living stars.. (beta parameter) ");
      printf(" : %.8g\n", IMFp->FactorSN);

      fprintf(FdSnInit,
	      "\n[ Energy from Short-Living Stars]\n"
	      "range is %6.3g - %6.3g Msun [%6.3g - %6.3g Gyr]\n"
	      "Mass fraction in tail is         : %.3lg%%\n"
	      "Number fraction in tail is       : %.3lg%%\n"
	      "energy per g due to S-L  snII is : %.8g erg/g (%.3lg%% of the total budget)\n"
	      "restored mass by Short-Liv stars (`beta' used in SF computation) is    :",
	      IMFp->MU, IMFp->egyShortLiv_MassTh, IMFp->inf_lifetime, IMFp->egyShortLiv_TimeTh,
	      IMFp->MassFrac_inIRA * 100,
	      IMFp->NumFrac_inIRA / num_snII * 100,
	      IMFp->IRA_erg_per_g, IMFp->IRA_erg_per_g / IMFp->TOT_erg_per_g * 100);
      fprintf(FdSnInit, " %.8g", IMFp->FactorSN);
      fprintf(FdSnInit, "\n");
      fflush(FdSnInit);
    }

  /* convert to code units */
  IMFp->IRA_erg_per_g *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;
  IMFp->TOT_erg_per_g *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

  IMFp->EgySpecSN = IMFp->IRA_erg_per_g * (1 - IMFp->FactorSN) / IMFp->FactorSN;
}


void setup_SF_related(int IMFi)
 /*
  *             |         RHOth indep. of Z               |       RHOth dep. on Z
  *  ==========================================================================================================
  *             |                                         |                                                  |
  *             | A calc. RHOth using Mth @ a given       |                                                  |
  *             |   Z using metal cooling                 |                                                  |
  *             |                                         |                                                  |
  *  RHOth from |   ................................      | F calculate RHOth using the given Mth for each Z |
  *  parameters |                                         |                                                  |
  *             | B calc. RHOth using Mth not using       |                                                  |
  *             |   the metal cooling                     |                                                  |
  *  ==========================================================================================================
  *             |                                         |                                                  |
  *             | C calc. Mth @ a given metallicity       | G calculate Mth @                                |
  *             |   keep Mth for all Z                    |   a given Z                                      |
  *  RHOth set  |   ...............................       |         \                                        |
  *  to a given | D calc. Mth not using the metal cool    |           --> keep Mth for all Z --> calc. RHOth |
  *  value      | --------------------------------------  |         /                                        |
  *             |                                         | H set Mth to a  given value                      |
  *             | E set Mth to a given value for all Z    |                                                  |
  * ===========================================================================================================
  *
  */
{
  int j, k, filecount, set;
  double num_snII;
  double myPhysDensTh, myFEVP, m;
  double left, right;
  double feedbackenergyinergs;

  IMFp = &IMFs[IMFi];

  if(IMFs[IMFi].egyShortLiv_MassTh > 0)
    /*
     *  Mth has been specified, check that it stays
     *  within boundaries.
     *  then translate it in a time threshold and calculate
     *  beta and egyspecSN.
     */
    {
      /* checks that the defined mass th does not fall in a BH mass range */
      for(j = 0;
	  (j < IMFs[IMFi].N_notBH_ranges) &&
	  (IMFs[IMFi].egyShortLiv_MassTh <= IMFs[IMFi].notBH_ranges.inf[j]); j++)
	;
      if(j < IMFs[IMFi].N_notBH_ranges)
	{
	  if(IMFs[IMFi].egyShortLiv_MassTh > IMFs[IMFi].notBH_ranges.sup[j])
	    {
	      if(j == 0)
		/* if the mass th falls ini the upper BH range, and that BH range reaches
		 *  MU, then it's impossible to correctly adjust the mass th
		 */
		{
		  if(ThisTask == 0)
		    {
		      printf("The Short-Living mass limit that you defined is not compatible with\n"
			     "the BH ranges that you defined. Pls check all that in param file and\n"
			     "in IMF file\n");
		      endrun(LT_ERR_MISMATCH_ShL_TH_BHranges);
		    }
		}
	      else
		{
		  IMFs[IMFi].egyShortLiv_MassTh = IMFs[IMFi].notBH_ranges.inf[j - 1];
		  if(ThisTask == 0)
		    fprintf(FdWarn, "The Short-Living mass limit has been set to %g in order to match"
			    " with the non-BH range %d\n", IMFs[IMFi].egyShortLiv_MassTh, j);
		}
	    }
	}
      else
	{
	  /* the mass th lives beyond the inf of the last non-BH mass renge */
	  IMFs[IMFi].egyShortLiv_MassTh = IMFs[IMFi].notBH_ranges.inf[j - 1];
	  if(ThisTask == 0)
	    fprintf(FdWarn, "The Short-Living mass limit has been set to %g in order to match"
		    " with the non-BH range %d\n", IMFs[IMFi].egyShortLiv_MassTh, j - 1);
	}

      /* checks that the mass th is not below the SnII inf mass */
      if(IMFs[IMFi].egyShortLiv_MassTh < All.Mup)
	{
	  if(ThisTask == 0)
	    printf("incorrect IRA limit set by your parameters:\n"
		   "  your inf IRA mass is %9.7g vs %9.7g Msun limit\n"
		   "standard %9.7g Msun limit has been set\n", IMFs[IMFi].egyShortLiv_TimeTh, All.Mup,
		   All.Mup);
	  IMFs[IMFi].egyShortLiv_MassTh = All.Mup;
	}

      /* checks that the mass th is not above the upper mass in the imf */
      if(IMFs[IMFi].egyShortLiv_MassTh > IMFs[IMFi].MU)
	{
	  if(ThisTask == 0)
	    printf("incorrect IRA limit set by your parameters:\n"
		   "  your inf IRA mass is %9.7g vs your %9.7g Msun upper limit\n"
		   "%9.7g Msun limit has been set\n", IMFs[IMFi].egyShortLiv_TimeTh, IMFs[IMFi].MU, All.Mup);
	  IMFs[IMFi].egyShortLiv_MassTh = All.Mup;
	}

      /* set the time th accordingly */
      IMFs[IMFi].egyShortLiv_TimeTh = lifetime(IMFs[IMFi].egyShortLiv_MassTh);

      calculate_ShortLiving_related(&IMFs[IMFi], 1);
    }

  if(IMFs[IMFi].SFTh_Zdep == 1)
    /*
     * if metallicity dependence is on, check whether the
     * thresholds' file does exist
     */
    {
      if(ThisTask == 0)
	{
	  if(RestartFlag == 2)
	    {
	      filecount = atoi(All.InitCondFile + strlen(All.InitCondFile) - 3);
	      j = read_eff_model(filecount, IMFi);
	    }
	  else
	    j = read_eff_model(-1, IMFi);
	}
      MPI_Bcast((void *) &IMFs[IMFi].PhysDensThresh[0], ZBins * sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);
      MPI_Bcast((void *) &IMFs[IMFi].FEVP[0], ZBins * sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);
    }

  if(IMFs[IMFi].PhysDensThresh[0] == 0)

    /*
     * we should calculate the value of the density threshold
     * from physical parameters.
     */
    {
      if(IMFs[IMFi].egyShortLiv_MassTh == 0)
	endrun(190000);		/* you must set Mth before! */

      if(IMFs[IMFi].SFTh_Zdep == 1)
	{
	  /*
	   * metal dependence is on but the thresholds' file does not
	   * exist.
	   */
	  init_clouds_cm(IMFs[IMFi].PhysDensThresh, IMFs[IMFi].FEVP,
			 IMFs[IMFi].EgySpecSN, IMFs[IMFi].FactorSN, ZBins, &Zvalue[0]);
	  write_eff_model(-1, IMFi);
	}
      else
	{
	  init_clouds(2, IMFs[IMFi].EgySpecSN, IMFs[IMFi].FactorSN,
		      IMFs[IMFi].referenceZ_toset_SF_DensTh, IMFs[IMFi].PhysDensThresh, FEVP);

	  if(IMFs[IMFi].SFTh_Zdep == 1)
	    {
	      /* note: the dens th at the referenceZ will be calculated once more if referenceZ is equal to one of
	       *       the Z bins in the metallicity array used to calculate the thresholds
	       */
	      init_clouds_cm(IMFs[IMFi].PhysDensThresh, IMFs[IMFi].FEVP,
			     IMFs[IMFi].EgySpecSN, IMFs[IMFi].FactorSN, ZBins, &Zvalue[0]);
	      write_eff_model(-1, IMFi);
	    }
	  else
	    {
	      if(ThisTask == 0)
		{
		  sfrrate_filenum = IMFi;
		  init_clouds(0, IMFs[IMFi].EgySpecSN, IMFs[IMFi].FactorSN,
			      get_metallicity_solarunits(IMFs[IMFi].referenceZ_toset_SF_DensTh),
			      IMFs[IMFi].PhysDensThresh, IMFs[IMFi].FEVP);
		}
	    }
	}
    }
  else
    /*
     * density threshold has been specified.
     * calculate the corresponding mass threshold, if needed.
     */
    {
      if(IMFs[IMFi].egyShortLiv_MassTh == 0)
	/*
	 * we have to calculate Mth
	 */
	{
	  /* initial guess for Mth */
	  IMFs[IMFi].egyShortLiv_MassTh = All.Mup * 1.1;

	  set = IMFs[IMFi].YSet;
	  SD.Zbin = 0;
	  SD.MArray = IIMbins[set];
	  SD.Mdim = IIMbins_dim[set];
	  SD.Y = SnIIEj[set];

	  get_Egy_and_Beta(IMFs[IMFi].egyShortLiv_TimeTh, &IMFs[IMFi].EgySpecSN, &IMFs[IMFi].FactorSN,
			   &IMFs[IMFi]);
	  init_clouds(2, IMFs[IMFi].EgySpecSN, IMFs[IMFi].FactorSN, IMFs[IMFi].referenceZ_toset_SF_DensTh,
		      &myPhysDensTh, &myFEVP);

	  while((fabs(myPhysDensTh - IMFs[IMFi].PhysDensThresh[0]) / IMFs[IMFi].PhysDensThresh[0] > 1e-2) &&
		IMFs[IMFi].egyShortLiv_MassTh > All.Mup &&
		IMFs[IMFi].egyShortLiv_MassTh < IMFs[IMFi].notBH_ranges.sup[0])
	    /* note: because of using get_Egy_an_Beta, this calculations automatically
	     * exclude BH ranges
	     */
	    {
	      if((myPhysDensTh - IMFs[IMFi].PhysDensThresh[0]) < 0)
		{
		  left = IMFs[IMFi].egyShortLiv_TimeTh;
		  if(right == 0)
		    IMFs[IMFi].egyShortLiv_MassTh *= 0.8;
		  else
		    IMFs[IMFi].egyShortLiv_MassTh = (IMFs[IMFi].egyShortLiv_MassTh + right) / 2;
		}
	      else
		{
		  right = IMFs[IMFi].egyShortLiv_MassTh;
		  if(left == 0)
		    IMFs[IMFi].egyShortLiv_MassTh *= 1.2;
		  else
		    IMFs[IMFi].egyShortLiv_MassTh = (IMFs[IMFi].egyShortLiv_MassTh + left) / 2;
		}

	      if(IMFs[IMFi].egyShortLiv_MassTh < All.Mup)
		IMFs[IMFi].egyShortLiv_MassTh = All.Mup;
	      else if(IMFs[IMFi].egyShortLiv_MassTh > IMFs[IMFi].notBH_ranges.sup[0])
		/* this prevent to fall in the last BH range: if it extend up to MU,
		 * we would obtain zero energy. If the energy is still not sufficient,
		 * the cycle will interrupt anyway due to the 3rd head condition */
		IMFs[IMFi].egyShortLiv_MassTh = IMFs[IMFi].notBH_ranges.sup[0];

	      get_Egy_and_Beta(IMFs[IMFi].egyShortLiv_MassTh, &IMFs[IMFi].EgySpecSN, &IMFs[IMFi].FactorSN,
			       &IMFs[IMFi]);

	      init_clouds(2, IMFs[IMFi].EgySpecSN, IMFs[IMFi].FactorSN,
			  IMFs[IMFi].referenceZ_toset_SF_DensTh, &myPhysDensTh, &myFEVP);
	    }

	  /* checks that found mass does not fall in a BH mass range */
	  for(j = 0;
	      (j < IMFs[IMFi].N_notBH_ranges) &&
	      (IMFs[IMFi].egyShortLiv_MassTh <= IMFs[IMFi].notBH_ranges.inf[j]); j++)
	    ;
	  if(j < IMFs[IMFi].N_notBH_ranges)
	    if(IMFs[IMFi].egyShortLiv_MassTh > IMFs[IMFi].notBH_ranges.sup[j])
	      /* this case can arise only when at least 1 BH interval lives
	       * between 2 nom-BH intervals: then, j shold be > 0.
	       */
	      IMFs[IMFi].egyShortLiv_MassTh = IMFs[IMFi].notBH_ranges.sup[j - 1];


	  if((myPhysDensTh - IMFs[IMFi].PhysDensThresh[0]) / IMFs[IMFi].PhysDensThresh[0] > 1e-2 &&
	     ThisTask == 0)
	    {
	      if(IMFs[IMFi].egyShortLiv_MassTh == All.Mup)
		{
		  fprintf(FdWarn, "  . warning: the needed Mass Threshold for a star to be\n"
			  "             short--living would be lower that %5.3f Msun.\n"
			  "             we set it to %5.3f Msun; pls check sfrrate.txt\n", All.Mup, All.Mup);
		  fflush(stdout);
		}
	      else if(IMFs[IMFi].egyShortLiv_MassTh == IMFs[IMFi].MU)
		{
		  fprintf(FdWarn, "  . warning: the needed Mass Threshold for a star to be\n"
			  "             short--living would be larger that %5.3f Msun.\n"
			  "             we set it to %5.3f Msun; pls check sfrrate.txt\n",
			  IMFs[IMFi].MU, IMFs[IMFi].MU);
		  fflush(stdout);
		}
	      else
		{
		  fprintf(FdWarn, "  . warning: Mass Threshold for a star to be short-living is\n"
			  "             %5.3f Msun; this would set a physical density threshold for\n"
			  "             the star formation different than that you set in paramfile.\n"
			  "             pls, check sfrrate.txt\n", IMFs[IMFi].egyShortLiv_MassTh);
		  fflush(stdout);
		}
	    }

	  calculate_ShortLiving_related(&IMFs[IMFi], 1);
	}

      if(IMFs[IMFi].SFTh_Zdep == 0)
	{
	  IMFs[IMFi].FEVP = &All.FactorEVP;
	  if(ThisTask == 0)
	    {
	      printf("determining Kennicutt law..\n");
	      sfrrate_filenum = IMFi;
	      init_clouds(0, IMFs[IMFi].EgySpecSN, IMFs[IMFi].FactorSN,
			  get_metallicity_solarunits(IMFs[IMFi].referenceZ_toset_SF_DensTh),
			  IMFs[IMFi].PhysDensThresh, IMFs[IMFi].FEVP);

	      printf("Energy Fraction in Winds= %g\n", IMFs[IMFi].WindEnergyFraction);
	      printf("PhysDensThresh= %g (internal units)\n", IMFs[IMFi].PhysDensThresh[0]);

	      fprintf(FdSnInit, "\n[Winds]\nEnergy Fraction in Winds= %g \n", IMFs[IMFi].WindEnergyFraction);
	      fprintf(FdSnInit, "\n[SF Rho Threshold]\nPhysDensThresh= %g (internal units)\n",
		      IMFs[IMFi].PhysDensThresh[0]);
	    }

	  MPI_Barrier(MPI_COMM_WORLD);
	}
    }

#ifndef LT_WIND_VELOCITY
#ifdef LT_ALL_SNII_EGY_INWINDS
  feedbackenergyinergs = IMFs[IMFi].TOT_erg_per_g / All.UnitMass_in_g * (All.UnitEnergy_in_cgs * SOLAR_MASS);
  IMFs[IMFi].WindEnergy =
    sqrt(2 * IMFs[IMFi].WindEnergyFraction * IMFs[IMFi].TOT_erg_per_g / IMFs[IMFi].WindEfficiency);
#else
  feedbackenergyinergs = IMFs[IMFi].IRA_erg_per_g / All.UnitMass_in_g * (All.UnitEnergy_in_cgs * SOLAR_MASS);
  IMFs[IMFi].WindEnergy =
    sqrt(2 * IMFs[IMFi].WindEnergyFraction * IMFs[IMFi].IRA_erg_per_g / IMFs[IMFi].WindEfficiency);
#endif
#else
  feedbackenergyinergs = All.FeedbackEnergy / All.UnitMass_in_g * (All.UnitEnergy_in_cgs * SOLAR_MASS);
  IMFs[IMFi].WindEnergy = LT_WIND_VELOCITY;
#ifdef LT_ALL_SNII_EGY_INWINDS
  IMFs[IMFi].WindEnergyFraction = IMFs[IMFi].WindEnergy * IMFs[IMFi].WindEnergy * IMFs[IMFi].WindEfficiency /
    (2 * IMFs[IMFi].TOT_erg_per_g);
#else
  IMFs[IMFi].WindEnergyFraction = IMFs[IMFi].WindEnergy * IMFs[IMFi].WindEnergy * IMFs[IMFi].WindEfficiency /
    (2 * IMFs[IMFi].IRA_erg_per_g);
#endif
#endif

  if(ThisTask == 0)
    {
      printf("Feedback energy per formed solar mass in stars= %g  ergs\n"
	     "Wind Velocity= %g Km/s\n", feedbackenergyinergs, IMFs[IMFi].WindEnergy);

      fprintf(FdSnInit, "Feedback energy per formed solar mass in stars= %g  ergs\n"
	      "Wind Velocity= %g Km/s\n", feedbackenergyinergs, IMFs[IMFi].WindEnergy);
    }

  m = max(IMFs[IMFi].Mm, All.Mup);

  if(IMFs[IMFi].metShortLiv_MassTh <= m)
    IMFs[IMFi].metShortLiv_MassTh = m;
  else if(IMFs[IMFi].metShortLiv_MassTh > IMFs[IMFi].MU)
    IMFs[IMFi].metShortLiv_MassTh = IMFs[IMFi].MU;

  /* checks that the defined mass th does not fall in a BH mass range */
  for(j = 0;
      (j < IMFs[IMFi].N_notBH_ranges) &&
      (IMFs[IMFi].metShortLiv_MassTh <= IMFs[IMFi].notBH_ranges.inf[j]); j++)
    ;

  if(j < IMFs[IMFi].N_notBH_ranges)
    if(IMFs[IMFi].metShortLiv_MassTh > IMFs[IMFi].notBH_ranges.sup[j])
      {
	if(j == 0)
	  /* if the mass th falls ini the upper BH range, and that BH range reaches
	   *  MU, then it's impossible to correctly adjust the mass th
	   */
	  {
	    if(ThisTask == 0)
	      {
		printf("The Short-Living mass limit for metal release that you defined is not\n"
		       "compatible with the BH ranges that you defined. Pls check all that in\n"
		       "param file and in IMF file\n");
		endrun(LT_ERR_MISMATCH_ShL_TH_BHranges);
	      }
	  }
	else
	  {
	    IMFs[IMFi].metShortLiv_MassTh = IMFs[IMFi].notBH_ranges.inf[j - 1];
	    fprintf(FdWarn,
		    "The Short-Living mass limit for metal release has been set to %g in order to match"
		    " with the non-BH range %d\n", IMFs[IMFi].egyShortLiv_MassTh, j);
	  }
      }

  IMFs[IMFi].metShortLiv_TimeTh = lifetime(IMFs[IMFi].metShortLiv_MassTh);

  if(ThisTask == 0)
    {
      if(IMFs[IMFi].metShortLiv_MassTh < IMFs[IMFi].MU)
	printf(":: metal IRA is active in the range [%9.7g - %9.7g]Msun <> [%9.7g - %9.7g]Gyr\n",
	       IMFs[IMFi].MU, IMFs[IMFi].metShortLiv_MassTh, IMFs[IMFi].inf_lifetime,
	       IMFs[IMFi].metShortLiv_TimeTh);
      else
	printf(":: metal IRA is not active\n");
    }

  if(IMFs[IMFi].metShortLiv_MassTh < IMFs[IMFi].MU)
    {
      /* mass fraction involved in SnII */
      IMFs[IMFi].MassFrac_inIRA =
	IntegrateIMF_byMass(IMFs[IMFi].metShortLiv_MassTh, IMFs[IMFi].MU, &IMFs[IMFi], INC_BH);
      /* expected number of SnII per initial stellar population mass in units of solar masses IN IRA RANGE */
      IMFs[IMFi].NumFrac_inIRA =
	IntegrateIMF_byNum(IMFs[IMFi].metShortLiv_MassTh, IMFs[IMFi].MU, &IMFs[IMFi], EXC_BH);
      /* expected number of SnII per initial stellar population mass in units of solar masses */
      num_snII = IntegrateIMF_byNum(max(All.Mup, IMFs[IMFi].Mm), IMFs[IMFi].MU, &IMFs[IMFi], EXC_BH);

      if(ThisTask == 0)
	printf("    %.3lg%% of mass in metal IRA tail\n", IMFs[IMFi].MassFrac_inIRA * 100);

      /* ...........: calculates the restored fraction :................... */
      if(ThisTask == 0)
	printf("   metal IRA, calculating restored fraction..\n");



/* ...........: calculates the effective yields :................... */

      if(ThisTask == 0)
	printf(":: calculating effective yields for Short-Living Stars..\n");


      SD.Yset = IMFs[IMFi].YSet;

      SnII_ShortLiv_Yields[IMFi] = (double **) malloc(LT_NMet * sizeof(double *));
      for(k = 0; k < LT_NMet; k++)
	SnII_ShortLiv_Yields[IMFi][k] = (double *) calloc(IIZbins_dim[SD.Yset], sizeof(double));

#ifdef UM_COOLING
      SnII_ShortLiv_FillMolNum[IMFi] = (double *) calloc(IIZbins_dim[SD.Yset], sizeof(double));
#endif

      IMFs[IMFi].nonZeroIRA = calculate_effective_yields(IMFs[IMFi].metShortLiv_MassTh, IMFs[IMFi].MU, IMFi);

      for(j = 0; j < IIZbins_dim[SD.Yset]; j++)
	{
	  if(j == 0)
	    {
	      SD.Zbin = j;
	      SD.MArray = IIMbins[SD.Yset];
	      SD.Mdim = IIMbins_dim[SD.Yset];
	      SD.Y = SnIIYields[SD.Yset][FillEl];
	      IMFs[IMFi].metFactorSN = calculate_FactorSN(IMFs[IMFi].metShortLiv_MassTh, IMFs[IMFi].MU, &IMFs[IMFi]);	/* restored by IRA SnII */
	    }

	  if(ThisTask == 0)
	    {
	      fprintf(FdSnInit, "\n[Effective Yields for Short-Living stars :: set %03d - Z %8.6g]\n",
		      SD.Yset, (IIZbins_dim[SD.Yset] > 1) ? IIZbins[SD.Yset][j] : 0);
	      fprintf(FdSnInit, " %3s  %6s   %10s %10s\n", "IMF", "Z", "name", "Yield");
	      printf(" %3s  %6s   %10s %10s\n", "IMF", "Z", "name", "Yield");
	      for(k = 0; k < LT_NMet; k++)
		{
		  printf(" [%3d][%8.6g]   %10s %10.7lg\n",
			 IMFi, (IIZbins_dim[SD.Yset] > 1) ? IIZbins[SD.Yset][j] : 0,
			 MetNames[k], SnII_ShortLiv_Yields[IMFi][k][j]);
		  fprintf(FdSnInit, " [%3d][%8.6g] %10s %10.7lg\n",
			  SD.Yset, (IIZbins_dim[SD.Yset] > 1) ? IIZbins[SD.Yset][j] : 0,
			  MetNames[k], SnII_ShortLiv_Yields[IMFi][k][j]);
		}
	    }
	}
    }
  else
    {
      IMFs[IMFi].MassFrac_inIRA = 0;
      IMFs[IMFi].NumFrac_inIRA = 0;
      IMFs[IMFi].nonZeroIRA = 0;
      num_snII = 0;
      SnII_ShortLiv_Yields[IMFi] = 0x0;
      IMFs[IMFi].metFactorSN = 0;
      if(ThisTask == 0)
	{
	  printf(":: No metals are supposed to be promptly ejected..\n");
	  fprintf(FdSnInit, "\n[Effective Yields for IRA part]\n"
		  "No metals are supposed to be promptly ejected..\n");
	}
    }

  if(ThisTask == 0)
    printf("\n\n");


}


void fsolver_error_handler(const char *reason, const char *file, int line, int err)
{
  if(err == GSL_EINVAL)
    {
      my_gslstatus = err;
      return;
    }
  return;
}

#endif

/*
 *
 * SCRATCH.NOTES AREA
 *


 *
 * END of SCRATCH.NOTES AREA
 *
 */
