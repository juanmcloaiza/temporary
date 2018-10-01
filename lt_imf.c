#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"

#ifdef LT_STELLAREVOLUTION

//#include "lt_error_codes.h"

#ifndef INLINE_FUNC
#ifdef INLINE
#define INLINE_FUNC inline
#else
#define INLINE_FUNC
#endif
#endif

/*
 * This file contains all the routines related to IMF.
 * Here you should add your own stuffs.
 */


/* ============================================================
 *
 *     ---=={[  IMF selection section  ]}==---
 * ============================================================ */

/* here place your own code for selecting IMF depending on gas
   physical conditions

   prototype:
   int INLINE_FUNC Get_IMF_index(int i)

   where i is the ordinal number of particle in the P and SphP
   structures

   returned integer is the ordinal number of the IMF
*/


/* no selection, only one IMF active */
/*
int INLINE_FUNC get_IMF_index(int i)
{
  return 0;
}
*/

/* POPIII-POPII selection */
int INLINE_FUNC get_IMF_index(int i)
{
#define INV_SOLAR 50
  double Z;

#ifndef LT_SMOOTH_Z
  int j;
  float *Zs;

  if(P[i].Type & 4)
    Zs = &MetP[P[i].MetID].Metals[0];
  else
    Zs = &SphP[i].Metals[0];

  for(Z = 0, j = 0; j < LT_NMetP; j++)
    if(j != Hel)
      Z += Zs[j];
  Z /= (P[i].Mass - Z - Zs[Hel]);
#else
  if(P[i].Type & 4)
    Z = MetP[P[i].MetID].Zsmooth;
  else
    Z = SphP[i].Zsmooth;
#endif


/*   int j; */
/*   float *Zs; */

/*   if(P[i].Type & 4) */
/*     Zs = &MetP[P[i].MetID].Metals[0]; */
/*   else */
/*     Zs = &SphP[i].Metals[0]; */

/*   for(Z = 0, j = 0; j < LT_NMetP; j++) */
/*     if(j != Hel) */
/*       Z += Zs[j]; */
/*   Z /= (P[i].Mass - Z - Zs[Hel]); */

/* #ifdef LT_SMOOTH_Z */
/*   if(Z > 0) */
/*     { */
/*       if(P[i].Type & 4) */
/*         Z = (Z + MetP[P[i].MetID].Zsmooth) / 2; */
/*       else */
/*         Z = (Z + SphP[i].Zsmooth) / 2; */
/*     } */
/* #endif   */

  if(Z * INV_SOLAR < 1e-4)
    return 1;
  else
    return 0;
}

/* ============================================================
 *
 *     ---=={[  IMF functions section  ]}==---
 * ============================================================ */


double IntegrateIMF_byMass(double Inf, double Sup, IMF_Type * IMFv, int mode)
{
  double inf, sup, partial_result;
  double *vinf, *vsup;
  double Sum, inf_, *params, slope;
  int i, Nint, Int;

  if(Sup < IMFv->Mm || Inf > IMFv->MU)
    return 0;
  if(Inf < IMFv->Mm)
    Inf = IMFv->Mm;
  if(Sup > IMFv->MU)
    Sup = IMFv->MU;

  sup = inf = Sup;
  Nint = IMFv->N_notBH_ranges;
  vinf = &IMFv->notBH_ranges.inf[0];
  vsup = &IMFv->notBH_ranges.sup[0];

  if(mode == INC_BH)
    /* includes also those stars that will end in BH */
    {
      inf = Inf;
      Int = Nint - 1;
    }
  else
    /* does not include those stars.. */
    {
      for(Int = 0; Int < Nint && sup <= vinf[Int]; Int++)
	/* get to the first non-BH mass range */
	;
    }

  partial_result = 0;
  while(inf > Inf || Int < Nint)
    {
      if(mode != INC_BH)
	{
	  if(sup > vsup[Int])
	    sup = vsup[Int];
	  if(Inf < vinf[Int])
	    inf = vinf[Int];
	  else
	    inf = Inf;
	}

      if(sup > inf)
	{
	  switch (IMFv->type)
	    {
	    case power_law:
	      if(IMFv->NParams == 0)
		{
		  /* pure power-law */
		  for(i = 0, Sum = 0; i < IMFv->NSlopes && inf >= IMFv->Slopes.masses[i + 1]; i++)
		    if(sup >= IMFv->Slopes.masses[i + 1])
		      {
			slope = IMFv->Slopes.slopes[i];

			if(inf < IMFv->Slopes.masses[i + 1])
			  inf_ = IMFv->Slopes.masses[i + 1];
			else
			  inf_ = inf;

			Sum += 1 / (1 - slope) * (pow(sup, 1 - slope) - pow(inf_, 1 - slope));

			sup = IMFv->Slopes.masses[i];
		      }
		  partial_result += IMFv->A * Sum;
		}
	      else
		{
		  /* ============================  */
		  /*   place here your own code    */
		  /* i.e. a call to IMFv->func()   */
		  /* ============================  */

		  /* built_in function: Larson IMF, time-dependent with one parameter */

		  params = IMFv->Params;
		  partial_result += IMFv->A * 1 / (1 - slope) *
		    (params[0] * (pow(1 + sup / params[0], 1 - slope) - pow(1 + inf / params[0], 1 - slope)));

		  /*
		   * as an example on how to proceed with your function if not (easily) analitically
		   * integrable, look at integration by number of the Larson IMF
		   */

		}
	      break;
	    case whatever:
	      /* ============================  */
	      /*   place here your own code    */
	      /* i.e. a call to IMFv->func()   */
	      /* ============================  */
	      break;
	    }
	}
      sup = vinf[Int];
      Int++;			/* get to the next non-BH mass range */
    }

  return partial_result;
}

double IntegrateIMF_byNum(double Inf, double Sup, IMF_Type * IMFv, int mode)
{
  double inf, sup, partial_result;
  double *vinf, *vsup;
  double Sum, inf_, slope;
  double result, err;
  int i, Nint, Int;

  if(Sup < IMFv->Mm || Inf > IMFv->MU)
    return 0;
  if(Inf < IMFv->Mm)
    Inf = IMFv->Mm;
  if(Sup > IMFv->MU)
    Sup = IMFv->MU;

  sup = inf = Sup;
  Nint = IMFv->N_notBH_ranges;
  vinf = IMFv->notBH_ranges.inf;
  vsup = IMFv->notBH_ranges.sup;

  if(mode == INC_BH)
    /* includes also those stars that will end in BH */
    {
      inf = Inf;
      Int = Nint - 1;
    }
  else
    /* does not include those stars.. */
    {
      for(Int = 0; Int < Nint && sup <= vinf[Int]; Int++)
	/* get to the first non-BH mass range */
	;
    }

  partial_result = 0;
  while(inf > Inf || Int < Nint)
    {
      if(mode != INC_BH)
	{
	  if(sup > vsup[Int])
	    sup = vsup[Int];
	  if(Inf < vinf[Int])
	    inf = vinf[Int];
	  else
	    inf = Inf;
	}

      if(sup > inf)
	{
	  switch (IMFv->type)
	    {
	    case power_law:
	      slope = IMFv->Slopes.slopes[0];
	      if(IMFv->NParams == 0)
		/* pure power-law */
		{
		  for(i = 0, Sum = 0; i < IMFv->NSlopes && inf >= IMFv->Slopes.masses[i + 1]; i++)
		    if(sup >= IMFv->Slopes.masses[i + 1])
		      {
			slope = IMFv->Slopes.slopes[i];

			if(inf < IMFv->Slopes.masses[i + 1])
			  inf_ = IMFv->Slopes.masses[i + 1];
			else
			  inf_ = inf;

			Sum += 1 / slope * (pow(inf_, -slope) - pow(sup, -slope));

			sup = IMFv->Slopes.masses[i];
		      }
		  partial_result += IMFv->A * Sum;
		}
	      else
		{
		  /* ============================  */
		  /*   place here your own code    */
		  /* i.e. a call to IMFv->func()   */
		  /* ============================  */

		  /* built_in function: Larson IMF, time-dependent with one parameter */

		  IMFp = IMFv;
		  F.function = IMFv->IMFfunc_byNum;
		  F.params = IMFv->Params;
		  if((gsl_status =
		      gsl_integration_qag(&F, inf, sup, 1e-6, 1e-4, gsl_ws_dim, 1, w, &result, &err)))
		    {
		      if(ThisTask == 0)
			printf
			  ("  >> Task %i, gsl integration error %i in Sn (iimf) [%.6g - %.6g] : %.6g %.6g\n",
			   ThisTask, gsl_status, inf, sup, result, err);
		      fflush(stdout);
		      endrun(1010);
		    }
		  partial_result += IMFv->A * 1 / (1 - slope) *
		    (IMFv->Params[0] *
		     (pow(1 + sup / IMFv->Params[0], 1 - slope) - pow(1 + inf / IMFv->Params[0], 1 - slope)));
		}
	      break;
	    case whatever:
	      /* ============================  */
	      /*   place here your own code    */
	      /* i.e. a call to IMFv->func()   */
	      /* ============================  */
	      break;
	    }
	}
      sup = vinf[Int];
      Int++;			/* get to the next non-BH mass range */
    }

  return partial_result;
}


double IntegrateIMF_byEgy(double Inf, double Sup, IMF_Type * IMFv, int mode)
     /*
        Integrates the energy coming form SnII explosions in the
        mass range [Inf, Sup], not accounting for those stars that
        end in BH.
      */
{
  double inf, sup, partial_result;
  double *vinf, *vsup;
  double Sum, inf_, slope;
  double result, err;
  int i, Nint, Int;

  if(Sup < IMFv->Mm || Inf > IMFv->MU)
    return 0;
  if(Inf < IMFv->Mm)
    Inf = IMFv->Mm;
  if(Sup > IMFv->MU)
    Sup = IMFv->MU;

  if(IMFv->NEKin == 1)
    /*
       All SnII explodes with the same energy
     */
    {
      Sum = IntegrateIMF_byNum(Inf, Sup, IMFv, EXC_BH);
      partial_result = Sum * IMFv->EKin.ekin[0];
    }
  else
    /*
       The explosion energy depends on the mass
     */
    {
      sup = inf = Sup;
      Nint = IMFv->N_notBH_ranges;
      vinf = IMFv->notBH_ranges.inf;
      vsup = IMFv->notBH_ranges.sup;

      for(Int = 0; Int < Nint && sup <= vinf[Int]; Int++)
	/* get to the first non-BH mass range */
	;

      partial_result = 0;
      while(inf > Inf || Int < Nint)
	{
	  if(mode != INC_BH)
	    {
	      if(sup > vsup[Int])
		sup = vsup[Int];
	      if(Inf < vinf[Int])
		inf = vinf[Int];
	      else
		inf = Inf;
	    }

	  if(sup > inf)
	    {
	      IMFp = IMFv;
	      F.function = IMFv->IMFfunc_byEgy;
	      F.params = IMFv->Params;
	      if((gsl_status =
		  gsl_integration_qag(&F, inf, sup, 1e-6, 1e-4, gsl_ws_dim, 1, w, &result, &err)))
		{
		  if(ThisTask == 0)
		    printf("  >> Task %i, gsl integration error %i in Sn (iimfe) [%.6g - %.6g] : %.6g %.6g\n",
			   ThisTask, gsl_status, inf, sup, result, err);
		  fflush(stdout);
		  endrun(1010);
		}
	      partial_result += result;
	    }
	  sup = vinf[Int];
	  Int++;		/* get to the next non-BH mass range */
	}
    }


  return partial_result;
}



/* ---------------------
 * by Mass IMF functions
 * --------------------- */

double INLINE_FUNC IMFevaluate_byMass_powerlaw(double arg, void *v)
     /*
        note: this function does not care whether or not the mass
        arg will end in a BH.
      */
{
  int i;

  /* pure power-law */
  for(i = 0; i < IMFp->NSlopes; i++)
    if(arg >= IMFp->Slopes.masses[i + 1])
      break;
  return IMFp->A * pow(arg, -IMFp->Slopes.slopes[i]);
}

double INLINE_FUNC IMFevaluate_byMass_Larson(double arg, void *v)
     /*
        note: this function does not care whether or not the mass
        arg will end in a BH.
      */
{
  /* built_in function: Larson IMF, time-dependent with one parameter */

  return IMFp->A * pow(1 + arg / *(double *) v, -(IMFp->Slopes.slopes[0]));
}



/* -----------------------
 * by Number IMF functions
 * ----------------------- */


double INLINE_FUNC IMFevaluate_byNum_powerlaw(double arg, void *v)
     /*
        note: this function does not care whether or not the mass
        arg will end in a BH.
      */
{
  int i;

  /* pure power-law */
  for(i = 0; i < IMFp->NSlopes; i++)
    if(arg >= IMFp->Slopes.masses[i + 1])
      break;
  return IMFp->A * pow(arg, -(1 + IMFp->Slopes.slopes[i]));
}

double INLINE_FUNC IMFevaluate_byNum_Larson(double arg, void *v)
     /*
        note: this function does not care whether or not the mass
        arg will end in a BH.
      */
{
  /* built_in function: Larson IMF, time-dependent with one parameter */

  return IMFp->A * pow(1 + arg / *(double *) v, -(IMFp->Slopes.slopes[0])) / arg;
}

/* -----------------------
 * get energy
 * ----------------------- */


double INLINE_FUNC IMFevaluate_byEgy_powerlaw(double arg, void *v)
     /*
        note: this function does not care whether the mass arg
        will end in a BH or in a SnII.
      */
{
  int i, j;
  double egy;

  /* pure power-law */
  for(i = 0; i < IMFp->NSlopes; i++)
    if(arg >= IMFp->Slopes.masses[i + 1])
      break;
  for(j = 0; j < IMFp->NEKin; j++)
    if(arg >= IMFp->EKin.masses[j + 1])
      break;
  egy = (IMFp->EKin.masses[j] - arg) / (IMFp->EKin.masses[j] - IMFp->EKin.masses[j + 1]);
  egy = IMFp->EKin.ekin[j] + (IMFp->EKin.ekin[j + 1] - IMFp->EKin.ekin[j]) * egy;
  return IMFp->A * pow(arg, -(1 + IMFp->Slopes.slopes[i])) * egy;
}

double INLINE_FUNC IMFevaluate_byEgy_Larson(double arg, void *v)
     /*
        note: this function does not care whether the mass arg
        will end in a BH or in a SnII.
      */
{
  double egy;
  int j;

  /* built_in function: Larson IMF, time-dependent with one parameter */

  for(j = 0; j < IMFp->NEKin; j++)
    if(arg >= IMFp->EKin.masses[j + 1])
      break;
  egy = (IMFp->EKin.masses[j] - arg) / (IMFp->EKin.masses[j] - IMFp->EKin.masses[j + 1]);
  egy = IMFp->EKin.ekin[j] + (IMFp->EKin.ekin[j + 1] - IMFp->EKin.ekin[j]) * egy;
  return IMFp->A * pow(1 + arg / *(double *) v, -(IMFp->Slopes.slopes[0])) / arg;
}


/* prototypes */

double INLINE_FUNC myIMF_byMass(double arg, void *v)
{
  return 0;
}

double INLINE_FUNC myIMF_byNum(double arg, void *v)
{
  return 0;
}


/* ============================================================
 *
 *     ---=={[  IMF parameters section  ]}==---
 * ============================================================ */

/*                                                           */
/*  you should place here your own code giving parameters of */
/*  imfs                                                     */

/*  this example is for Larson IMF with params taken from    */
/*  Boehringer et al, 2003                                   */

double get_imf_params(int I, double *v, double time)
{
#define z_imf_inf 2
#define time_imf_inf (1/(z_imf_inf + 1))
#define z_imf_sup 10
#define time_imf_sup (1/(z_imf_sup + 1))
#define ms_inf 0.4467
#define ms_sup 10
#define logsup  log10(ms_sup)
#define loginf log10(ms_inf)

  double z = 1.0 / time - 1;

  if(z <= z_imf_inf)
    return 0.35;
  else if(z >= z_imf_sup)
    return 10;
  else
    return pow(10, ((z - z_imf_inf) / (z_imf_sup - z_imf_inf) * (logsup - loginf) + loginf));
}


/* ============================================================
 *
 *     ---=={[  IMF utils  ]}==---
 * ============================================================ */

int INLINE_FUNC get_IMF_SlopeBin(double m)
{
  int i;

  for(i = 0; i < IMFp->NSlopes && m >= IMFp->Slopes.masses[i]; i++)
    ;
  return i;
}


/* ============================================================
 *
 *     ---=={[  IMF I/O section  ]}==---
 * ============================================================ */

void write_IMF_info(int num, FILE * file)
{
  int i;
  double sup, inf;

  switch (IMFs[num].type)
    {
    case power_law:
      i = 0;
      break;
    case whatever:
      i = 1;
      break;
    }

  fprintf(file,
	  "::  IMF %3d\n"
	  "    type: %10s\n"
	  "    depends on time: %3s\n"
	  "    normalization: %8.6f\n",
	  num, IMF_Spec_Labels[i], (IMFs[num].timedep) ? "yes" : "no", IMFs[num].A);
  if(IMFs[num].NParams > 0)
    {
      fprintf(file, "    Parameters: \n");
      for(i = 0; i < IMFs[num].NParams; i++)
	fprintf(file, "      %8.6f ", IMFs[num].Params[i]);
      fprintf(file, "\n");
    }
  else
    fprintf(file, "    No parameters\n");

  fprintf(file, "    Slopes: \n");
  for(i = 0, sup = IMFs[num].MU; i < IMFs[num].NSlopes; i++)
    {
      if(i + 1 < IMFs[num].NSlopes)
	inf = IMFs[num].Slopes.masses[i + 1];
      else
	inf = IMFs[num].Mm;
      fprintf(file, "      [%6.4f:%6.4f] -> %6.4f\n", sup, inf, IMFs[num].Slopes.slopes[i]);
      sup = IMFs[num].Slopes.masses[i];
    }

  return;
}


int read_ekin_file(void *, double **, double **);

void read_imfs(void)
     /*
      * this routine read in imfs' data
      */
{
  FILE *file;
  char *charp, buff[500], param[100], value[100];

  int IMF_index, IMF_Nspec;
  int i, j, p, n, line_num;

  if(ThisTask == 0)
    /* Task 0 read all data and then broadcast them */
    {
      sprintf(buff, "%s%s", All.OutputDir, All.IMFfilename);
      if((file = fopen(buff, "r")) == 0x0)
	{
	  printf("it's impossible to open the imf file <%s>\nwe terminate here!\n", buff);
	  endrun(LT_ERR_IMF_FILE);
	}

      /* now read in the imfs data.
       * the file is organized as follows:
       * first, 4 lines contain respectively the imf type, the time dependence of
       * the imf (0|1), the maximum and minimum mass.
       * then a line contains the number of slopes, wheter the imf is multi-slope, and
       * as much lines as the slopes follow, each with a pair (limiting mass, slope).
       * finally a line contains the number of parameters: each of the following lines
       * contains a parameter.
       * lines beginning with a "#" are ignored.
       */


      fscanf(file, "%d\n", &IMFs_dim);
      if(IMFs_dim != LT_NIMFs)
	{
	  printf
	    ("sever error: the number of IMFs in IMF file <%s> does not match with the compile-time NIMFs variable\n",
	     buff);
	  endrun(LT_ERR_NIMF);
	}
      IMFs = (IMF_Type *) calloc(IMFs_dim, sizeof(IMF_Type));

      if(IMFs == 0x0)
	{
	  printf("[Task 0] memory allocation failed when reading imfs' file\n");
	  endrun(LT_ERR_IMF_ALLOCATE);
	}

      IsThere_TimeDep_IMF = 0;
      IsThere_ZDep_IMF = 0;
      i = 0;
      line_num = 0;
      IMF_index = -1;
      IMF_Nspec = 0;
      do
	{
	  charp = fgets(buff, 500, file);
	  line_num++;
	  if(charp != 0x0 && buff[0] != '#')
	    {
	      ((p = sscanf(buff, "%s %s", &param[0], &value[0])) >= 1) ? n++ : n;

	      if(strcmp(param, "TYPE") == 0)
		{
		  if(IMF_index >= 0)
		    {
		      if(IMF_Nspec < IMF_NSPEC)
			{
			  printf("error in ifms file format at line %d:"
				 " IMF num %d has not been completely specified\n", line_num, IMF_index);
			  endrun(LT_ERR_IMF_FILE_FORMAT);
			}

		      if(IMFs[IMF_index].NSlopes == 1)
			IMFs[IMF_index].Slopes.masses[0] = IMFs[IMF_index].MU;
		    }
		  if(++IMF_index > IMFs_dim)
		    {
		      printf("error in ifms file format at line %d\n", line_num);
		      endrun(LT_ERR_IMF_FILE_FORMAT);
		    }

		  IMF_Nspec = 1;
		  IMFs[IMF_index].PhysDensThresh = (double *) malloc(sizeof(double));
		  IMFs[IMF_index].FEVP = (double *) malloc(sizeof(double));
		  if(strcmp(value, "PowerLaw") == 0)
		    {
		      IMFs[IMF_index].type = power_law;
		      IMFs[IMF_index].IMFfunc_byMass = &IMFevaluate_byMass_powerlaw;
		      IMFs[IMF_index].IMFfunc_byNum = &IMFevaluate_byNum_powerlaw;
		    }
		  else if(strcmp(buff, "Whatever") == 0)
		    IMFs[IMF_index].type = whatever;
		  /* else if: place here your own code */
		}
	      else if(strcmp(param, "TDEP") == 0)
		{
		  IMF_Nspec++;
		  IsThere_TimeDep_IMF += (IMFs[IMF_index].timedep = atoi(value));
		}
	      else if(strcmp(param, "ZDEP") == 0)
		{
		  IMF_Nspec++;
		  if(atoi(value) == 0)
		    IMFs[IMF_index].SFTh_Zdep = 0;
		  else
		    {
		      IsThere_ZDep_IMF += (IMFs[IMF_index].SFTh_Zdep = 1);
		      IMFs[IMF_index].PhysDensThresh = (double *) realloc(IMFs[IMF_index].PhysDensThresh,
									  ZBins * sizeof(double));
		      IMFs[IMF_index].FEVP = (double *) realloc(IMFs[IMF_index].FEVP, ZBins * sizeof(double));
		    }
		}
	      else if(strcmp(param, "MU") == 0)
		{
		  IMF_Nspec++;
		  IMFs[IMF_index].MU = atof(value);
		}
	      else if(strcmp(param, "Mm") == 0)
		{
		  IMF_Nspec++;
		  IMFs[IMF_index].Mm = atof(value);
		}
	      else if(strcmp(param, "Z4Th") == 0)
		{
		  IMF_Nspec++;
		  IMFs[IMF_index].referenceZ_toset_SF_DensTh = atof(value);
		}
	      else if(strcmp(param, "SFTh") == 0)
		{
		  IMF_Nspec++;
		  IMFs[IMF_index].PhysDensThresh[0] = atof(value);
		}
	      else if(strcmp(param, "EgyMTh") == 0)
		{
		  IMF_Nspec++;
		  IMFs[IMF_index].egyShortLiv_MassTh = atof(value);
		}
	      else if(strcmp(param, "MetMTh") == 0)
		{
		  IMF_Nspec++;
		  IMFs[IMF_index].metShortLiv_MassTh = atof(value);
		}
	      else if(strcmp(param, "NPar") == 0)
		{
		  IMF_Nspec++;

		  if((IMFs[IMF_index].NParams = atoi(value)) > 0)
		    {
		      if((IMFs[IMF_index].Params =
			  (double *) malloc(IMFs[IMF_index].NParams * sizeof(double))) == 0x0)
			{
			  printf("[Task 0][a] memory allocation failed when reading imfs' file\n");
			  endrun(LT_ERR_IMF_ALLOCATE);
			}
		      for(j = 0; j < IMFs[IMF_index].NParams; j++)
			{
			  i = fscanf(file, "%s\n", buff);
			  line_num++;
			  if(i == EOF || i < 1)
			    {
			      printf("error in ifms' file format at line %d\n", line_num);
			      endrun(LT_ERR_IMF_FILE_FORMAT);
			    }
			  if(buff[0] != '#')
			    i = sscanf(buff, "%lf\n", &IMFs[IMF_index].Params[j]);
			}
		    }
		  else
		    IMFs[IMF_index].Params = 0x0;
		}
	      else if(strcmp(param, "NSlopes") == 0)
		{
		  IMF_Nspec++;

		  IMFs[IMF_index].NSlopes = atoi(value);

		  IMFs[IMF_index].Slopes.masses =
		    (double *) malloc((IMFs[IMF_index].NSlopes + 1) * sizeof(double));
		  IMFs[IMF_index].Slopes.slopes =
		    (double *) malloc((IMFs[IMF_index].NSlopes + 1) * sizeof(double));
		  if(IMFs[IMF_index].Slopes.masses == 0x0 || IMFs[IMF_index].Slopes.slopes == 0x0)
		    {
		      printf("[Task 0][c] memory allocation failed when reading imfs' file\n");
		      endrun(LT_ERR_IMF_ALLOCATE);
		    }

		  IMFs[IMF_index].Slopes.masses[IMFs[IMF_index].NSlopes] = 0;
		  IMFs[IMF_index].Slopes.slopes[IMFs[IMF_index].NSlopes] = 0;
		  for(j = 0; j < IMFs[IMF_index].NSlopes; j++)
		    {
		      charp = fgets(buff, 500, file);
		      line_num++;
		      if(charp == 0x0)
			{
			  printf("error in IMFs' file format at line %d: not enough slopes specified",
				 line_num);
			  endrun(LT_ERR_IMF_FILE_FORMAT);
			}
		      if(buff[0] != '#')
			/* the slope file format is as follos:
			 *   inf mass, slope
			 * where inf mass is the inf limit of the mass interval to which
			 * the slope refers. then, larger masses come first. */
			i = sscanf(buff, "%lf %lf\n",
				   &IMFs[IMF_index].Slopes.masses[j], &IMFs[IMF_index].Slopes.slopes[j]);
		    }
		  IMFs[IMF_index].Slopes.masses[j] = 0;
		  IMFs[IMF_index].Slopes.slopes[j] = 0;
		}
	      else if(strcmp(param, "N_notBH") == 0)
		{
		  IMF_Nspec++;

		  IMFs[IMF_index].N_notBH_ranges = atoi(value);

		  IMFs[IMF_index].notBH_ranges.sup =
		    (double *) malloc((IMFs[IMF_index].N_notBH_ranges + 1) * sizeof(double));
		  IMFs[IMF_index].notBH_ranges.inf =
		    (double *) malloc((IMFs[IMF_index].N_notBH_ranges + 1) * sizeof(double));
		  if(IMFs[IMF_index].notBH_ranges.sup == 0x0 || IMFs[IMF_index].notBH_ranges.inf == 0x0)
		    {
		      printf("[Task 0][d] memory allocation failed when reading imfs' file\n");
		      endrun(LT_ERR_IMF_ALLOCATE);
		    }

		  for(j = 0; j < IMFs[IMF_index].N_notBH_ranges; j++)
		    {
		      charp = fgets(buff, 500, file);
		      line_num++;
		      if(charp == 0x0)
			{
			  printf("error in IMFs' file format at line %d: not enough slopes specified",
				 line_num);
			  endrun(LT_ERR_IMF_FILE_FORMAT);
			}
		      if(buff[0] != '#')
			/* the slope file format is as follos:
			 *   inf mass, slope
			 * where inf mass is the inf limit of the mass interval to which
			 * the slope refers. then, larger masses come first. */
			i = sscanf(buff, "%lf %lf\n",
				   &IMFs[IMF_index].notBH_ranges.sup[j],
				   &IMFs[IMF_index].notBH_ranges.inf[j]);
		    }
		  IMFs[IMF_index].notBH_ranges.sup[j] = 0;
		  IMFs[IMF_index].notBH_ranges.inf[j] = 0;
		}
	      else if(strcmp(param, "YSet") == 0)
		{
		  IMF_Nspec++;
		  IMFs[IMF_index].YSet = atoi(value);
		}
	      else if(strcmp(param, "NGen") == 0)
		{
		  IMF_Nspec++;
		  IMFs[IMF_index].Generations = atoi(value);
		}
	      else if(strcmp(param, "EKin") == 0)
		{
		  IMF_Nspec++;

		  if(isdigit(*value))
		    {
		      IMFs[IMF_index].EKin.masses = (double *) malloc(sizeof(double));
		      IMFs[IMF_index].EKin.ekin = (double *) malloc(sizeof(double));
		      if(IMFs[IMF_index].MU > 0)
			*(IMFs[IMF_index].EKin.masses) = IMFs[IMF_index].MU;
		      else
			{
			  printf("please, specify single-value EKin after MU in the IMF file\n");
			  endrun(LT_ERR_IMF_FILE_FORMAT);
			}
		      *(IMFs[IMF_index].EKin.ekin) = atof(value) * 1e51 / All.UnitEnergy_in_cgs;
		      IMFs[IMF_index].NEKin = 1;
		      IMFs[IMF_index].IMFfunc_byEgy = 0x0;
		    }
		  else
		    {
		      IMFs[IMF_index].NEKin =
			read_ekin_file(value, &IMFs[IMF_index].EKin.masses, &IMFs[IMF_index].EKin.ekin);
		      if(IMFs[IMF_index].type == power_law)
			IMFs[IMF_index].IMFfunc_byEgy = IMFevaluate_byEgy_powerlaw;
		      /* place your own else if here */
		    }
		}
	      else if(strcmp(param, "WindEff") == 0)
		{
		  IMF_Nspec++;
		  IMFs[IMF_index].WindEfficiency = atof(value);
		}
	      else if(strcmp(param, "WindEF") == 0)
		{
		  IMF_Nspec++;
		  IMFs[IMF_index].WindEnergyFraction = atof(value);
		}
	      else
		{
		  printf("error in IMFs' format file at line %d: unknown parameter\n", line_num);
		  endrun(LT_ERR_IMF_FILE_FORMAT);
		}
	    }
	}			/* close the while */
      while(charp != 0x0);

    }				/* close if(ThisTask == 0) */


  MPI_Bcast(&IMFs_dim, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if(ThisTask != 0)
    if((IMFs = (IMF_Type *) malloc(IMFs_dim * sizeof(IMF_Type))) == 0x0)
      {
	printf("[Task %d] memory allocation failed when reading imfs' file\n", ThisTask);
	endrun(LT_ERR_IMF_ALLOCATE);
      }

  for(i = 0; i < IMFs_dim; i++)
    MPI_Bcast(&IMFs[i], sizeof(IMF_Type), MPI_BYTE, 0, MPI_COMM_WORLD);

  if(ThisTask != 0)
    for(i = 0; i < IMFs_dim; i++)
      {
	if((IMFs[i].PhysDensThresh = (double *) malloc(ZBins * sizeof(double))) == 0x0)
	  {
	    printf("[Task %d][a] memory allocation failed when reading imfs' file\n", ThisTask);
	    endrun(LT_ERR_IMF_ALLOCATE);
	  }
	if((IMFs[i].Params = (double *) malloc(IMFs[i].NParams * sizeof(double))) == 0x0)
	  {
	    printf("[Task %d][b] memory allocation failed when reading imfs' file\n", ThisTask);
	    endrun(LT_ERR_IMF_ALLOCATE);
	  }
	IMFs[i].Slopes.masses = (double *) malloc((IMFs[i].NSlopes + 1) * sizeof(double));
	IMFs[i].Slopes.slopes = (double *) malloc((IMFs[i].NSlopes + 1) * sizeof(double));
	if(IMFs[i].Slopes.masses == 0x0 || IMFs[i].Slopes.slopes == 0x0)
	  {
	    printf("[Task %d][c] memory allocation failed when reading imfs' file\n", ThisTask);
	    endrun(LT_ERR_IMF_ALLOCATE);
	  }
	IMFs[i].notBH_ranges.sup = (double *) malloc((IMFs[i].N_notBH_ranges + 1) * sizeof(double));
	IMFs[i].notBH_ranges.inf = (double *) malloc((IMFs[i].N_notBH_ranges + 1) * sizeof(double));
	if(IMFs[i].notBH_ranges.sup == 0x0 || IMFs[i].notBH_ranges.inf == 0x0)
	  {
	    printf("[Task 0][d] memory allocation failed when reading imfs' file\n");
	    endrun(LT_ERR_IMF_ALLOCATE);
	  }
	IMFs[i].EKin.masses = (double *) malloc((IMFs[i].NEKin) * sizeof(double));
	IMFs[i].EKin.ekin = (double *) malloc((IMFs[i].NEKin) * sizeof(double));
	if(IMFs[i].EKin.masses == 0x0 || IMFs[i].EKin.ekin == 0x0)
	  {
	    printf("[Task %d][e] memory allocation failed when reading imfs' file\n", ThisTask);
	    endrun(LT_ERR_IMF_ALLOCATE);
	  }
      }

  MPI_Barrier(MPI_COMM_WORLD);

  for(i = 0; i < IMFs_dim; i++)
    {
      MPI_Bcast(&IMFs[i].PhysDensThresh[0], sizeof(double) * ZBins, MPI_BYTE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&IMFs[i].Params[0], sizeof(double) * IMFs[i].NParams, MPI_BYTE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&IMFs[i].Slopes.slopes[0], sizeof(double) * (IMFs[i].NSlopes + 1), MPI_BYTE, 0,
		MPI_COMM_WORLD);
      MPI_Bcast(&IMFs[i].Slopes.masses[0], sizeof(double) * (IMFs[i].NSlopes + 1), MPI_BYTE, 0,
		MPI_COMM_WORLD);
      MPI_Bcast(&IMFs[i].notBH_ranges.sup[0], sizeof(double) * (IMFs[i].N_notBH_ranges), MPI_BYTE, 0,
		MPI_COMM_WORLD);
      MPI_Bcast(&IMFs[i].notBH_ranges.inf[0], sizeof(double) * (IMFs[i].N_notBH_ranges), MPI_BYTE, 0,
		MPI_COMM_WORLD);
      MPI_Bcast(&IMFs[i].EKin.masses[0], sizeof(double) * (IMFs[i].NEKin), MPI_BYTE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&IMFs[i].EKin.ekin[0], sizeof(double) * (IMFs[i].NEKin), MPI_BYTE, 0, MPI_COMM_WORLD);
    }

  if((IMF_CommBuff = (double *) malloc(sizeof(double) * 3)) == 0x0)
    {
      printf("[Task %d][e] memory allocation failed when reading imfs' file\n", ThisTask);
      endrun(LT_ERR_IMF_ALLOCATE);
    }

  if(ThisTask == 0)
    {
      printf("%d IMF%c read\n", IMFs_dim, (IMFs_dim > 1) ? 's' : '\0');
      fflush(stdout);
    }

  return;
}


int read_ekin_file(void *string, double **mass_array, double **ekin_array)
{
  double *masses, *ekin;
  FILE *myfile;
  char *filename, *charp, buff[500];
  int i, Num;

  filename = (char *) string;
  if((myfile = fopen(filename, "r")) == 0x0)
    {
      printf("it's impossible to open the kinetic energy file <%s>\nwe terminate here!\n", filename);
      endrun(LT_ERR_IMF_FILE);
    }

  do
    charp = fgets(buff, 500, myfile);
  while(buff[0] == '#' && charp != 0x0);

  sscanf(buff, "%i %*s\n", &Num);
  masses = (double *) malloc(Num * sizeof(double));
  ekin = (double *) malloc(Num * sizeof(double));

  for(i = 0; i < Num; i++)
    {
      do
	charp = fgets(buff, 500, myfile);
      while((charp != 0x0) && (buff[0] == '#'));
      sscanf(buff, "%lg %lg\n", &masses[i], &ekin[i]);
      ekin[i] *= 1e51 / All.UnitEnergy_in_cgs;
    }

  *mass_array = masses;
  *ekin_array = ekin;
  return Num;
}

#endif
