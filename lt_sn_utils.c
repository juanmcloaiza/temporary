#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"

//#include "lt_error_codes.h"

#ifdef LT_STELLAREVOLUTION

#ifndef INLINE_FUNC
#ifdef INLINE
#define INLINE_FUNC inline
#else
#define INLINE_FUNC
#endif
#endif

#ifdef LT_TRACK_CONTRIBUTES
void INLINE_FUNC unpack_fraction(Contrib *, float *, float *, float *);
void INLINE_FUNC pack_fraction(Contrib *, float *, float *, float *);
#endif

double INLINE_FUNC get_metallicity_solarunits(FLOAT Z)
{
  if((Z *= 381.7) < 1e-5)
    return NO_METAL;
  else
    return log10(Z);
}

void INLINE_FUNC get_element_metallicity(int i, double *zs)
{
  int j;
  double MetalMass, Hmass;

  for(j = 0, MetalMass = 0; j < LT_NMetP; j++)
    MetalMass += SphP[i].Metals[j];

  /* ok, you're right, here we should use the real chemical composition.. be patient */
  Hmass = (P[i].Mass - MetalMass);

  for(j = 0; j < LT_NMetP; j++)
    zs[j] = SphP[i].Metals[j] / Hmass;

  return;
}

double INLINE_FUNC get_metallicity(int i, int mode)
     /*
      * mode = 0 :: return metallicity in dex
      * mode = 1 :: return the Hydrogen Mass
      * mode = 2 :: return the ratio between Iron Mass and Hydrogen Mass
      * mode = 3 :: return the Metal Mass
      */
{
  /* Sutherland&Dopita 1993 used Anders & Grevesse 1989 as for solar abundances.
     Then,
     Fe / H = 2.62e-3  by mass
     O / H = 1.36e-2 by mass
   */

  float *Metals;
  double MetalMass = 0, Hmass, r;
  int j, go;

  if(P[i].Type == 0)
    Metals = &SphP[i].Metals[0];
  else
    Metals = &MetP[P[i].MetID].Metals[0];

  for(j = 0; j < LT_NMetP; j++)
    MetalMass += Metals[j];

  Hmass = (P[i].Mass - MetalMass);

  if(mode == 1)
    return Hmass;
  if(mode == 3)
    return MetalMass;

  go = 1;

#ifdef LT_MEAN_Z_ESTIMATE
  /*
     we use a mean over the current iron content and an iron amount inferred from
     oxigen content as for a solar chemical composition
     O * (Iron/O)_solar = inferred iron content using solar composition
   */

  if(Iron >= 0 && Oxygen >= 0)
    {
      r = (Metals[Iron] + (Metals[Oxygen] * 0.193)) / 2 / Hmass * 381.7;
      go = 0;
    }
#endif

  if(go)
    r = Metals[Iron] / Hmass;

  if(mode)
    return r;
  else
    {
      if(r * 381.7 >= 1e-4)
	return log10(r * 381.7);
      else
	return NO_METAL;
    }
}


/* temporary */

/* #ifdef LT_SEv_INFO */
/* void init_T_TC_runs() */
/* { */
/*   int i; */
/*   char mode[2]; */

/*   Trun_size = 100; */
/*   TCrun_size = 100; */

/*   Tpop_num = (double*)calloc(Trun_size, sizeof(double)); */
/*   Tpop_mass = (double*)calloc(Trun_size, sizeof(double)); */
/*   tot_Tpop_num = (double*)calloc(Trun_size, sizeof(double)); */
/*   tot_Tpop_mass = (double*)calloc(Trun_size, sizeof(double)); */


/*   TCpop_num = (double*)calloc(TCrun_size, sizeof(double)); */
/*   TCpop_mass = (double*)calloc(TCrun_size, sizeof(double)); */
/*   tot_TCpop_num = (double*)calloc(TCrun_size, sizeof(double)); */
/*   tot_TCpop_mass = (double*)calloc(TCrun_size, sizeof(double)); */


/*   delta_Trun = (8.0 - 3.0) / Trun_size; */

/*   delta_TCrun = (1.0 + 5.0) / Trun_size; */

/*   if(ThisTask == 0) */
/*     { */
/*       if(RestartFlag == 0) */
/* 	strcpy(mode, "w"); */
/*       else */
/* 	strcpy(mode, "a"); */
      
/*       FdTrun = fopen("Trun.txt", mode); */
/*       FdTCrun = fopen("TCrun.txt", mode); */
/*     } */

/*   return; */
/* } */
/* #endif */
/* */


#ifdef LT_TRACK_CONTRIBUTES

void init_packing()
{
  int i, N;

  Packing_Factor = (1 << LT_Nbits) - 1;
  UnPacking_Factor = 1.0 / Packing_Factor;

  Max_Power10 = (1 << LT_power10_Nbits) - 1;

  Power10_Factors = (int *) calloc(Max_Power10, sizeof(int));
  Power10_Factors[0] = Packing_Factor;
  for(i = 1; i < Max_Power10; i++)
    Power10_Factors[i] = Power10_Factors[i - 1] * 0.1;

  return;
}

int INLINE_FUNC get_packing_power10(float fraction)
{
  int i;

  for(i = Max_Power10 - 1; i > 0; fraction *= 10, i--)
    if(fraction <= (float) Power10_Factors[i])
      break;

  return i;
}

void INLINE_FUNC pack_contrib(Contrib * contrib, int IMFi, float *II, float *Ia, float *AGB)
{
  int Power10;

  /* example using LT_NMet = 9 */

  if(IMFi == 0)
    {
      /* Massive stars */
      contrib->II_el0_imf0 = II[0] * Packing_Factor;
      contrib->IIexp_el0_imf0 = get_packing_power10(contrib->II_el0_imf0);
      contrib->II_el1_imf0 = II[1] * Packing_Factor;
      contrib->IIexp_el1_imf0 = get_packing_power10(contrib->II_el1_imf0);
      contrib->II_el2_imf0 = II[2] * Packing_Factor;
      contrib->IIexp_el2_imf0 = get_packing_power10(contrib->II_el2_imf0);
      contrib->II_el3_imf0 = II[3] * Packing_Factor;
      contrib->IIexp_el3_imf0 = get_packing_power10(contrib->II_el3_imf0);
      contrib->II_el4_imf0 = II[4] * Packing_Factor;
      contrib->IIexp_el4_imf0 = get_packing_power10(contrib->II_el4_imf0);
      contrib->II_el5_imf0 = II[5] * Packing_Factor;
      contrib->IIexp_el5_imf0 = get_packing_power10(contrib->II_el5_imf0);
      contrib->II_el6_imf0 = II[6] * Packing_Factor;
      contrib->IIexp_el6_imf0 = get_packing_power10(contrib->II_el6_imf0);
      contrib->II_el7_imf0 = II[7] * Packing_Factor;
      contrib->IIexp_el7_imf0 = get_packing_power10(contrib->II_el7_imf0);
      contrib->II_el8_imf0 = II[8] * Packing_Factor;
      contrib->IIexp_el8_imf0 = get_packing_power10(contrib->II_el8_imf0);

      /* Ia supernovae */
      contrib->Ia_el0_imf0 = Ia[0] * Packing_Factor;
      contrib->Iaexp_el0_imf0 = get_packing_power10(contrib->Ia_el0_imf0);
      contrib->Ia_el1_imf0 = Ia[1] * Packing_Factor;
      contrib->Iaexp_el1_imf0 = get_packing_power10(contrib->Ia_el1_imf0);
      contrib->Ia_el2_imf0 = Ia[2] * Packing_Factor;
      contrib->Iaexp_el2_imf0 = get_packing_power10(contrib->Ia_el2_imf0);
      contrib->Ia_el3_imf0 = Ia[3] * Packing_Factor;
      contrib->Iaexp_el3_imf0 = get_packing_power10(contrib->Ia_el3_imf0);
      contrib->Ia_el4_imf0 = Ia[4] * Packing_Factor;
      contrib->Iaexp_el4_imf0 = get_packing_power10(contrib->Ia_el4_imf0);
      contrib->Ia_el5_imf0 = Ia[5] * Packing_Factor;
      contrib->Iaexp_el5_imf0 = get_packing_power10(contrib->Ia_el5_imf0);
      contrib->Ia_el6_imf0 = Ia[6] * Packing_Factor;
      contrib->Iaexp_el6_imf0 = get_packing_power10(contrib->Ia_el6_imf0);
      contrib->Ia_el7_imf0 = Ia[7] * Packing_Factor;
      contrib->Iaexp_el7_imf0 = get_packing_power10(contrib->Ia_el7_imf0);
      contrib->Ia_el8_imf0 = Ia[8] * Packing_Factor;
      contrib->Iaexp_el8_imf0 = get_packing_power10(contrib->Ia_el8_imf0);

      /* AGB stars */
      contrib->AGB_el0_imf0 = AGB[0] * Packing_Factor;
      contrib->AGBexp_el0_imf0 = get_packing_power10(contrib->AGB_el0_imf0);
      contrib->AGB_el1_imf0 = AGB[1] * Packing_Factor;
      contrib->AGBexp_el1_imf0 = get_packing_power10(contrib->AGB_el1_imf0);
      contrib->AGB_el2_imf0 = AGB[2] * Packing_Factor;
      contrib->AGBexp_el2_imf0 = get_packing_power10(contrib->AGB_el2_imf0);
      contrib->AGB_el3_imf0 = AGB[3] * Packing_Factor;
      contrib->AGBexp_el3_imf0 = get_packing_power10(contrib->AGB_el3_imf0);
      contrib->AGB_el4_imf0 = AGB[4] * Packing_Factor;
      contrib->AGBexp_el4_imf0 = get_packing_power10(contrib->AGB_el4_imf0);
      contrib->AGB_el5_imf0 = AGB[5] * Packing_Factor;
      contrib->AGBexp_el5_imf0 = get_packing_power10(contrib->AGB_el5_imf0);
      contrib->AGB_el6_imf0 = AGB[6] * Packing_Factor;
      contrib->AGBexp_el6_imf0 = get_packing_power10(contrib->AGB_el6_imf0);
      contrib->AGB_el7_imf0 = AGB[7] * Packing_Factor;
      contrib->AGBexp_el7_imf0 = get_packing_power10(contrib->AGB_el7_imf0);
      contrib->AGB_el8_imf0 = AGB[8] * Packing_Factor;
      contrib->AGBexp_el8_imf0 = get_packing_power10(contrib->AGB_el8_imf0);
    }
  /* Add here below more blocks for additional IMFs
   *
   * else if(IMFi == ..)
   */
  else
    {
      contrib->II_el0_imf1 = II[0] * Packing_Factor;
      contrib->IIexp_el0_imf1 = get_packing_power10(contrib->II_el0_imf1);
      contrib->II_el1_imf1 = II[1] * Packing_Factor;
      contrib->IIexp_el1_imf1 = get_packing_power10(contrib->II_el1_imf1);
      contrib->II_el2_imf1 = II[2] * Packing_Factor;
      contrib->IIexp_el2_imf1 = get_packing_power10(contrib->II_el2_imf1);
      contrib->II_el3_imf1 = II[3] * Packing_Factor;
      contrib->IIexp_el3_imf1 = get_packing_power10(contrib->II_el3_imf1);
      contrib->II_el4_imf1 = II[4] * Packing_Factor;
      contrib->IIexp_el4_imf1 = get_packing_power10(contrib->II_el4_imf1);
      contrib->II_el5_imf1 = II[5] * Packing_Factor;
      contrib->IIexp_el5_imf1 = get_packing_power10(contrib->II_el5_imf1);
      contrib->II_el6_imf1 = II[6] * Packing_Factor;
      contrib->IIexp_el6_imf1 = get_packing_power10(contrib->II_el6_imf1);
      contrib->II_el7_imf1 = II[7] * Packing_Factor;
      contrib->IIexp_el7_imf1 = get_packing_power10(contrib->II_el7_imf1);
      contrib->II_el8_imf1 = II[8] * Packing_Factor;
      contrib->IIexp_el8_imf1 = get_packing_power10(contrib->II_el8_imf1);
    }
}


void INLINE_FUNC unpack_contrib(Contrib * contrib, float *II, float *Ia, float *AGB)
{
  int i, Power10;

  /* example using LT_NMet = 9 */

  /* Massive stars */
  for(Power10 = 1, i = 0; i < (int) contrib->IIexp_el0_imf0; Power10 *= 10, i++);
  II[0] = contrib->II_el0_imf0 * UnPacking_Factor / Power10;
  for(Power10 = 1, i = 0; i < (int) contrib->IIexp_el1_imf0; Power10 *= 10, i++);
  II[1] = contrib->II_el1_imf0 * UnPacking_Factor / Power10;
  for(Power10 = 1, i = 0; i < (int) contrib->IIexp_el2_imf0; Power10 *= 10, i++);
  II[2] = contrib->II_el2_imf0 * UnPacking_Factor / Power10;
  for(Power10 = 1, i = 0; i < (int) contrib->IIexp_el3_imf0; Power10 *= 10, i++);
  II[3] = contrib->II_el3_imf0 * UnPacking_Factor / Power10;
  for(Power10 = 1, i = 0; i < (int) contrib->IIexp_el4_imf0; Power10 *= 10, i++);
  II[4] = contrib->II_el4_imf0 * UnPacking_Factor / Power10;
  for(Power10 = 1, i = 0; i < (int) contrib->IIexp_el5_imf0; Power10 *= 10, i++);
  II[5] = contrib->II_el5_imf0 * UnPacking_Factor / Power10;
  for(Power10 = 1, i = 0; i < (int) contrib->IIexp_el6_imf0; Power10 *= 10, i++);
  II[6] = contrib->II_el6_imf0 * UnPacking_Factor / Power10;
  for(Power10 = 1, i = 0; i < (int) contrib->IIexp_el7_imf0; Power10 *= 10, i++);
  II[7] = contrib->II_el7_imf0 * UnPacking_Factor / Power10;
  for(Power10 = 1, i = 0; i < (int) contrib->IIexp_el8_imf0; Power10 *= 10, i++);
  II[8] = contrib->II_el8_imf0 * UnPacking_Factor / Power10;

  /* Ia supernovae */
  for(Power10 = 1, i = 0; i < (int) contrib->Iaexp_el0_imf0; Power10 *= 10, i++);
  Ia[0] = contrib->Ia_el0_imf0 * UnPacking_Factor / Power10;
  for(Power10 = 1, i = 0; i < (int) contrib->Iaexp_el1_imf0; Power10 *= 10, i++);
  Ia[1] = contrib->Ia_el1_imf0 * UnPacking_Factor / Power10;
  for(Power10 = 1, i = 0; i < (int) contrib->Iaexp_el2_imf0; Power10 *= 10, i++);
  Ia[2] = contrib->Ia_el2_imf0 * UnPacking_Factor / Power10;
  for(Power10 = 1, i = 0; i < (int) contrib->Iaexp_el3_imf0; Power10 *= 10, i++);
  Ia[3] = contrib->Ia_el3_imf0 * UnPacking_Factor / Power10;
  for(Power10 = 1, i = 0; i < (int) contrib->Iaexp_el4_imf0; Power10 *= 10, i++);
  Ia[4] = contrib->Ia_el4_imf0 * UnPacking_Factor / Power10;
  for(Power10 = 1, i = 0; i < (int) contrib->Iaexp_el5_imf0; Power10 *= 10, i++);
  Ia[5] = contrib->Ia_el5_imf0 * UnPacking_Factor / Power10;
  for(Power10 = 1, i = 0; i < (int) contrib->Iaexp_el6_imf0; Power10 *= 10, i++);
  Ia[6] = contrib->Ia_el6_imf0 * UnPacking_Factor / Power10;
  for(Power10 = 1, i = 0; i < (int) contrib->Iaexp_el7_imf0; Power10 *= 10, i++);
  Ia[7] = contrib->Ia_el7_imf0 * UnPacking_Factor / Power10;
  for(Power10 = 1, i = 0; i < (int) contrib->Iaexp_el8_imf0; Power10 *= 10, i++);
  Ia[8] = contrib->Ia_el8_imf0 * UnPacking_Factor / Power10;

  /* AGB stars */
  for(Power10 = 1, i = 0; i < (int) contrib->AGBexp_el0_imf0; Power10 *= 10, i++);
  AGB[0] = contrib->AGB_el0_imf0 * UnPacking_Factor / Power10;
  for(Power10 = 1, i = 0; i < (int) contrib->AGBexp_el1_imf0; Power10 *= 10, i++);
  AGB[1] = contrib->AGB_el1_imf0 * UnPacking_Factor / Power10;
  for(Power10 = 1, i = 0; i < (int) contrib->AGBexp_el2_imf0; Power10 *= 10, i++);
  AGB[2] = contrib->AGB_el2_imf0 * UnPacking_Factor / Power10;
  for(Power10 = 1, i = 0; i < (int) contrib->AGBexp_el3_imf0; Power10 *= 10, i++);
  AGB[3] = contrib->AGB_el3_imf0 * UnPacking_Factor / Power10;
  for(Power10 = 1, i = 0; i < (int) contrib->AGBexp_el4_imf0; Power10 *= 10, i++);
  AGB[4] = contrib->AGB_el4_imf0 * UnPacking_Factor / Power10;
  for(Power10 = 1, i = 0; i < (int) contrib->AGBexp_el5_imf0; Power10 *= 10, i++);
  AGB[5] = contrib->AGB_el5_imf0 * UnPacking_Factor / Power10;
  for(Power10 = 1, i = 0; i < (int) contrib->AGBexp_el6_imf0; Power10 *= 10, i++);
  AGB[6] = contrib->AGB_el6_imf0 * UnPacking_Factor / Power10;
  for(Power10 = 1, i = 0; i < (int) contrib->AGBexp_el7_imf0; Power10 *= 10, i++);
  AGB[7] = contrib->AGB_el7_imf0 * UnPacking_Factor / Power10;
  for(Power10 = 1, i = 0; i < (int) contrib->AGBexp_el8_imf0; Power10 *= 10, i++);
  AGB[8] = contrib->AGB_el8_imf0 * UnPacking_Factor / Power10;

  /* Add here below more blocks for additional IMFs; elements will be addressed
   * as lt_NMetP + N_IMF*#_ELEMENT, e.g. for the 2nd IMF the 1st element will 
   * have address 9 = 8 + 1*0 etc.
   */
  /* Massive stars */

  for(Power10 = 1, i = 0; i < (int) contrib->IIexp_el0_imf1; Power10 *= 10, i++);
  II[LT_NMet + 0] = contrib->II_el0_imf1 * UnPacking_Factor / Power10;
  for(Power10 = 1, i = 0; i < (int) contrib->IIexp_el1_imf1; Power10 *= 10, i++);
  II[LT_NMet + 1] = contrib->II_el1_imf1 * UnPacking_Factor / Power10;
  for(Power10 = 1, i = 0; i < (int) contrib->IIexp_el2_imf1; Power10 *= 10, i++);
  II[LT_NMet + 2] = contrib->II_el2_imf1 * UnPacking_Factor / Power10;
  for(Power10 = 1, i = 0; i < (int) contrib->IIexp_el3_imf1; Power10 *= 10, i++);
  II[LT_NMet + 3] = contrib->II_el3_imf1 * UnPacking_Factor / Power10;
  for(Power10 = 1, i = 0; i < (int) contrib->IIexp_el4_imf1; Power10 *= 10, i++);
  II[LT_NMet + 4] = contrib->II_el4_imf1 * UnPacking_Factor / Power10;
  for(Power10 = 1, i = 0; i < (int) contrib->IIexp_el5_imf1; Power10 *= 10, i++);
  II[LT_NMet + 5] = contrib->II_el5_imf1 * UnPacking_Factor / Power10;
  for(Power10 = 1, i = 0; i < (int) contrib->IIexp_el6_imf1; Power10 *= 10, i++);
  II[LT_NMet + 6] = contrib->II_el6_imf1 * UnPacking_Factor / Power10;
  for(Power10 = 1, i = 0; i < (int) contrib->IIexp_el7_imf1; Power10 *= 10, i++);
  II[LT_NMet + 7] = contrib->II_el7_imf1 * UnPacking_Factor / Power10;
  for(Power10 = 1, i = 0; i < (int) contrib->IIexp_el8_imf1; Power10 *= 10, i++);
  II[LT_NMet + 8] = contrib->II_el8_imf1 * UnPacking_Factor / Power10;


}

void update_contrib(Contrib * current_contrib, float *current_metals, Contrib * new_contrib,
		    float *new_metals)
     /*

      */
{
  float fractionII[LT_NMet * LT_NIMFs], nfractionII[LT_NMet * LT_NIMFs];
  float fractionIa[LT_NMet * LT_NIMFs], nfractionIa[LT_NMet * LT_NIMFs];
  float fractionAGB[LT_NMet * LT_NIMFs], nfractionAGB[LT_NMet * LT_NIMFs];

  float inv_total;
  int power, i, j;

  /* gets current fractional contributes */
  unpack_contrib(current_contrib, fractionII, fractionIa, fractionAGB);
  unpack_contrib(new_contrib, nfractionII, nfractionIa, nfractionAGB);

  /* updates the fractional contributes */
  for(j = 0; j < LT_NIMFs; j++)
    {
      for(i = 0; i < LT_NMet; i++)
	{
	  if(current_metals[i] + new_metals[i] > 0)
	    {
	      inv_total = 1.0 / (current_metals[i] + new_metals[i]);

	      fractionII[j * LT_NMet + i] =
		(fractionII[j * LT_NMet + i] * current_metals[i] +
		 nfractionII[j * LT_NMet + i] * new_metals[i]) * inv_total;

	      fractionIa[j * LT_NMet + i] =
		(fractionIa[j * LT_NMet + i] * current_metals[i] +
		 nfractionIa[j * LT_NMet + i] * new_metals[i]) * inv_total;

	      fractionAGB[j * LT_NMet + i] =
		(fractionAGB[j * LT_NMet + i] * current_metals[i] +
		 nfractionAGB[j * LT_NMet + i] * new_metals[i]) * inv_total;
	    }
	  else
	    fractionII[j * LT_NMet + i] = fractionIa[j * LT_NMet + i] = fractionAGB[j * LT_NMet + i] = 0;
	}
      /* store in a packed forme the fractional contributes */
      pack_contrib(current_contrib, j, fractionII, fractionIa, fractionAGB);
    }


  return;
}
#endif


#endif
