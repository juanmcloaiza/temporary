#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"


/*! \file c_enrichment.c 
 *  
 *  This file contains the routines for the enrichment by SNIa and SNII.
 * 
 */


#ifdef SFR_METALS

#include "c_metals.h"

#ifdef SFR_ENRICH


double sum_reservoir, total_reservoir;
double sum_distributed, total_distributed;
double dist_metals, tot_met_distributed;



/* This function marks stars that have to explode as SNII or SNIa. */

void flag_SN_starparticles(void)
{
  int n, numSNI, ntotI;
  int numSNII, ntotII;
  double time_hubble_a, hubble_a;
  float age;
  float deltaSNI;
  double tlife_SNII = 0, metal;
  int ik;


  if(All.ComovingIntegrationOn)
    {
      hubble_a = All.Hubble * sqrt(All.Omega0 / (All.Time * All.Time * All.Time)
				   + (1 - All.Omega0 - All.OmegaLambda) / (All.Time * All.Time) +
				   All.OmegaLambda);
      time_hubble_a = All.Time * hubble_a;
    }
  else
    time_hubble_a = 1;

  for(n = 0, numSNI = 0, numSNII = 0; n < NumPart; n++)
    if(P[n].Type == 4)
      if(P[n].Ti_endstep == All.Ti_Current)
	{

	  if(All.ComovingIntegrationOn)
	    age = integrated_time(n, time_hubble_a);
	  else
	    age = All.Time - P[n].StellarAge;

#ifdef SFR_SNII
	  tlife_SNII = All.TlifeSNII;

	  if(PPP[n].NumNgb != 0)
	     /*SNII*/
	    {
	      metal = 0;
	      for(ik = 1; ik < 12; ik++)	/* all chemical elements but H & He */
		if(ik != 6)
		  metal += P[n].Zm[ik];

	      metal /= P[n].Mass;

	      /* Raitieri estimations of mean tlife for SNII according to metallicity
	         if(metal <= 1.e-5)
	         tlife_SNII = 0.00682912;
	         else
	         if(metal <= 1.e-4)
	         tlife_SNII =  0.0124270;
	         else
	         if(metal <= 1.e-3)
	         tlife_SNII = 0.0163819;
	         else
	         if(metal <= 1.e-2)
	         tlife_SNII = 0.0152885;

	         if(metal > 1.e-2)
	         tlife_SNII = 0.0142681;
	       */

	      if(age >= tlife_SNII)
		 /*SNII*/
		{
		  numSNII++;
		  P[n].Type |= 3;	/* we mark SNII */

		  if(age >= All.MinTlifeSNI)
		    endrun(32890);
		}
	    }
#endif
#ifdef SFR_SNI
	  deltaSNI = All.MinTlifeSNI + (All.MaxTlifeSNI - All.MinTlifeSNI) * drand48();

	  if(age >= deltaSNI && PPP[n].Hsml != 0)
	     /*SNI*/
	    {
	      numSNI++;
	      P[n].Type |= 2;	/* we mark SNI */
	    }
#endif

	}

#ifdef SFR_SNI
  MPI_Allreduce(&numSNI, &ntotI, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(ThisTask == 0 && ntotI > 0)
    printf("flagging ntot=%d stars as exploding type SNIa\n", ntotI);
#endif
#ifdef SFR_SNII
  MPI_Allreduce(&numSNII, &ntotII, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(ThisTask == 0 && ntotII > 0)
    printf("flagging ntot=%d stars as exploding type SNII\n", ntotII);
#endif

}


/* Routine of integration for calculating time interval from expansion factor interval.
This is relevant only if ComovingIntegration is ON. */
float integrated_time(int indice, double time_hubble_a)
{
  float t1, t2, t3;
  float f1, f2, f3;
  float deltat;

  t1 = P[indice].StellarAge;
  t3 = All.Time;
  t2 = t1 + (t3 - t1) / 2;

  f1 = 1 / (t1 * All.Hubble * sqrt(All.Omega0 / (t1 * t1 * t1)
				   + (1 - All.Omega0 - All.OmegaLambda) / (t1 * t1) + All.OmegaLambda));
  f2 = 1 / (t2 * All.Hubble * sqrt(All.Omega0 / (t2 * t2 * t2)
				   + (1 - All.Omega0 - All.OmegaLambda) / (t2 * t2) + All.OmegaLambda));
  f3 = 1 / time_hubble_a;

  deltat = (t3 - t1) / 2. * (f1 / 3. + 4. / 3. * f2 + f3 / 3.);

  return deltat;
}


/* This function  updates the weights for SN before exploding. This is necessary
due to the fact that gas particles neighbours of a given star could have
been transformed into stars and they need to be taken off the neighbour list
for the exploding star. */
void update_weights(void)
{
  int i;

  for(i = 0; i < NumPart; i++)
    if(!(P[i].Type == 6 || P[i].Type == 7))	/*SNIa or SNII */
      P[i].Ti_endstep = -P[i].Ti_endstep - 1;	/* Mark as inactive if we are NOT 
						   an exploding star */

  density();

  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep < 0)
      P[i].Ti_endstep = -P[i].Ti_endstep - 1;	/* unhide */

}


/* This function contains the main enrichment loop for exploding SNIa or SNII.*/
void enrichment(void)
{
  int *noffset, *nbuffer, *nsend, *nsend_local;
  int i, j, n;
  int ndone, ndonetot, ntot, ntotleft, npleft;
  int maxfill;
  int level, ngrp, sendTask, recvTask;
  int nexport;
  MPI_Status status;
  int numofenrichment, ik;
  double delta_metalsI, sum_SNI = 0, rate_SNI, total_SNI;
  double delta_metalsII, sum_SNII = 0, rate_SNII, total_SNII;
  double frac_phase = 0.;
  double hubble_a, time_hubble_a, a3inv;
  double check = 0.;
  double total_metals = 0;



#ifdef SFR_FEEDBACK
  double sum_energy = 0, total_energy;

  if(Flag_phase == 1)
    {
      sum_reservoir = 0;
      sum_distributed = 0;
      dist_metals = 0;
    }
#endif

  if(All.ComovingIntegrationOn)	/* Factors for comoving integration of hydro */
    {
      a3inv = 1 / (All.Time * All.Time * All.Time);
      hubble_a = All.Hubble * sqrt(All.Omega0 / (All.Time * All.Time * All.Time)
				   + (1 - All.Omega0 - All.OmegaLambda) / (All.Time * All.Time) +
				   All.OmegaLambda);
      time_hubble_a = All.Time * hubble_a;
    }
  else
    a3inv = time_hubble_a = 1;


#ifndef SFR_FEEDBACK
  frac_phase = 1;
#else
  if(Flag_phase == 1)
    frac_phase = 1. - All.FracEnergySN_Phase;
  if(Flag_phase == 2)
    frac_phase = All.FracEnergySN_Phase;
#endif

  noffset = mymalloc(sizeof(int) * NTask);	/* offsets of bunches in common list */
  nbuffer = mymalloc(sizeof(int) * NTask);
  nsend_local = mymalloc(sizeof(int) * NTask);
  nsend = mymalloc(sizeof(int) * NTask * NTask);


  for(n = 0, numofenrichment = 0; n < NumPart; n++)
    {
      if(P[n].Type == 6)	/* exploding SNI  */
	if(P[n].Ti_endstep == All.Ti_Current)
	  numofenrichment++;

      if(P[n].Type == 7)	/* exploding SNII */
	if(P[n].Ti_endstep == All.Ti_Current)
	  numofenrichment++;
    }

  MPI_Allreduce(&numofenrichment, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);


  i = 0;			/* beginn with this index */
  npleft = numofenrichment;	/* particles left for this task */
  ntotleft = ntot;		/* particles left for all tasks together */


  while(ntotleft > 0)
    {
      for(j = 0; j < NTask; j++)
	nsend_local[j] = 0;

      /* do local particles and prepare export list */

      for(nexport = 0, ndone = 0; i < NumPart && nexport < All.BunchSizeMetal - NTask; i++)
	if(P[i].Type == 6 || P[i].Type == 7)
	  {
	    if(P[i].Ti_endstep == All.Ti_Current)
	      {
		ndone++;
		/* NOTE:  First we fill in SNII reservoir and then add contribution of SNI. But a given particle
		   should always explode first as SNII and later as SNIa, so a star will be of type 6 OR 7. */
#ifdef SFR_FEEDBACK
		if(Flag_phase == 1)	/* Only once we should produce the metals */
#endif
		  {
#ifdef SFR_SNII
		    if(P[i].Type == 7)
		      {
			delta_metalsII = SNII_yields(i);
			sum_SNII += delta_metalsII;
		      }
#endif
#ifdef SFR_SNI
		    if(P[i].Type == 6)
		      {
			delta_metalsI = SNI_yields(i);
			sum_SNI += delta_metalsI;
		      }
#endif
		  }
#ifdef SFR_FEEDBACK
		sum_energy += P[i].EnergySN;
		sum_energy += P[i].EnergySNCold;
#endif

		for(j = 0; j < NTask; j++)
		  Exportflag[j] = 0;

		enrichment_evaluate(i, 0);

		for(j = 0; j < NTask; j++)
		  {
		    if(Exportflag[j])
		      {
			MetalDataIn[nexport].Pos[0] = P[i].Pos[0];
			MetalDataIn[nexport].Pos[1] = P[i].Pos[1];
			MetalDataIn[nexport].Pos[2] = P[i].Pos[2];

			MetalDataIn[nexport].Hsml = PPP[i].Hsml;
			MetalDataIn[nexport].NumNgb = PPP[i].NumNgb;
			for(ik = 0; ik < 12; ik++)
			  MetalDataIn[nexport].ZmReservoir[ik] = P[i].ZmReservoir[ik];
#ifdef SFR_FEEDBACK
			MetalDataIn[nexport].EnergySN = P[i].EnergySN;
			MetalDataIn[nexport].EnergySNCold = P[i].EnergySNCold;
#endif
			MetalDataIn[nexport].Index = i;
			MetalDataIn[nexport].Task = j;
			nexport++;
			nsend_local[j]++;
		      }
		  }
	      }
	  }

      qsort(MetalDataIn, nexport, sizeof(struct metaldata_in), metals_compare_key);

      for(j = 1, noffset[0] = 0; j < NTask; j++)
	noffset[j] = noffset[j - 1] + nsend_local[j - 1];


      MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);


      /* now do the particles that need to be exported */

      for(level = 1; level < (1 << PTask); level++)
	{
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
	      if(maxfill >= All.BunchSizeMetal)
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
				   recvTask, TAG_ENRICH_A,
				   &MetalDataGet[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct metaldata_in),
				   MPI_BYTE, recvTask, TAG_ENRICH_A, MPI_COMM_WORLD, &status);
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }

	  for(j = 0; j < nbuffer[ThisTask]; j++)
	    enrichment_evaluate(j, 1);

	  level = ngrp - 1;
	}

      MPI_Allreduce(&ndone, &ndonetot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      npleft -= ndone;
      ntotleft -= ndonetot;
    }

  MPI_Reduce(&sum_SNI, &total_SNI, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sum_SNII, &total_SNII, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      if(All.TimeStep > 0)
	{
	  rate_SNII = total_SNII / (All.TimeStep / time_hubble_a);
	  rate_SNI = total_SNI / (All.TimeStep / time_hubble_a);
	}
      else
	{
	  rate_SNII = 0;
	  rate_SNI = 0;
	}

      if(rate_SNI > 0 || rate_SNII > 0)
	fprintf(FdSN, "%g %g %g %g %g \n", All.Time, total_SNI, rate_SNI, total_SNII, rate_SNII);
      fflush(FdSN);
    }


#ifdef SFR_FEEDBACK
  MPI_Reduce(&sum_energy, &total_energy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    if(total_energy > 0)
      {
	fprintf(FdSNE, "%g %d %g %g %g\n", All.Time, Flag_phase, total_energy, total_SNI, total_SNII);
	fflush(FdSNE);

	total_metals = total_SNI + total_SNII;
	printf("EnergyTest  total_produced=%g total_metals=%g\n", total_energy, total_metals);
	fflush(stdout);
      }



  if(Flag_phase == 2)
    for(i = 0; i < NumPart; i++)
      if(P[i].Type == 0)
	sum_reservoir += P[i].EnergySN;


  if(Flag_phase == 2)
    {
      MPI_Reduce(&sum_reservoir, &total_reservoir, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&sum_distributed, &total_distributed, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&dist_metals, &tot_met_distributed, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    }


  if(ThisTask == 0 && Flag_phase == 2)
    printf("EnergyTest   total_reservoir=%g total_distributed=%g total_distr_metals=%g\n",
	   total_reservoir, total_distributed, tot_met_distributed);

#endif


  /* do final operations */

  if(Flag_phase == 0 || Flag_phase == 2)
    for(i = 0; i < NumPart; i++)
      if(P[i].ZmReservoir[4] > 0)
	if(P[i].Ti_endstep == All.Ti_Current)
	  {

	    check = 0;
	    for(ik = 0; ik < 12; ik++)
	      check += P[i].ZmReservoir[ik];
	    if(check == 0)
	      {
		printf("ERROR Part=%d enters cleaning in enrichment with Reservoir=%g\n", P[i].ID, check);
		endrun(334);
	      }

	    for(ik = 0; ik < 12; ik++)
	      P[i].ZmReservoir[ik] = 0;

	    if(P[i].Type == 6)	/* SNIa */
	      {
		PPP[i].Hsml = 0;	/* to avoid exploding as SNIa again */
		PPP[i].NumNgb = 0;	/* to avoid exploding as SNII again */
	      }

	    if(P[i].Type == 7)	/* SNII */
	      PPP[i].NumNgb = 0;	/* to avoid exploding as SNII again */


	    P[i].Type &= 4;

#ifdef SFR_FEEDBACK
	    P[i].EnergySN = 0;	/* only stars should enter here */
	    P[i].EnergySNCold = 0;	/* only stars should enter here */
#endif
	  }

  myfree(nsend);
  myfree(nsend_local);
  myfree(nbuffer);
  myfree(noffset);

}



/*! This function represents the core of the enrichment calculation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
void enrichment_evaluate(int target, int mode)
{
  int j, n, ik;
  int startnode, numngb_inbox;
  double h, h2, hinv, hinv3;
  double wk;
  double dx, dy, dz, r, r2, u, mass_j;
  double energySN, inv_wk_i, sum_metals = 0;
  double energySNCold;
  FLOAT pos[3], reservoir[12];
  FLOAT numngb;
  double frac_phase = 0;

#ifdef SFR_FEEDBACK
  double a3inv;
  double xhyd, yhel, ne, mu;
  double energy, temp;
  double new_temp;
#endif

#ifdef PERIODIC
  double boxSize, boxHalf;

  boxSize = All.BoxSize;
  boxHalf = 0.5 * All.BoxSize;
#endif

#ifdef SFR_FEEDBACK
  if(All.ComovingIntegrationOn)
    a3inv = 1 / (All.Time * All.Time * All.Time);
  else
    a3inv = 1;

  if(Flag_phase == 1)
    frac_phase = 1. - All.FracEnergySN_Phase;
  if(Flag_phase == 2)
    frac_phase = All.FracEnergySN_Phase;
#else
  frac_phase = 1;
#endif

  if(mode == 0)
    {
      pos[0] = P[target].Pos[0];
      pos[1] = P[target].Pos[1];
      pos[2] = P[target].Pos[2];
      h = PPP[target].Hsml;
      numngb = P[target].NumNgb;
      for(ik = 0; ik < 12; ik++)
	reservoir[ik] = P[target].ZmReservoir[ik];
#ifdef SFR_FEEDBACK
      energySN = P[target].EnergySN;
      energySNCold = P[target].EnergySNCold;
#endif
    }
  else
    {
      pos[0] = MetalDataGet[target].Pos[0];
      pos[1] = MetalDataGet[target].Pos[1];
      pos[2] = MetalDataGet[target].Pos[2];
      h = MetalDataGet[target].Hsml;
      numngb = MetalDataGet[target].NumNgb;
      for(ik = 0; ik < 12; ik++)
	reservoir[ik] = MetalDataGet[target].ZmReservoir[ik];
#ifdef SFR_FEEDBACK
      energySN = MetalDataGet[target].EnergySN;
      energySNCold = MetalDataGet[target].EnergySNCold;
#endif
    }

  h2 = h * h;
  hinv = 1.0 / h;
  hinv3 = hinv * hinv * hinv;

  inv_wk_i = 4 * M_PI / 3.0 / numngb;	/* renormalization of  the kernel */

  startnode = All.MaxPart;
  do
    {
#ifdef SFR_DECOUPLING
      numngb_inbox = ngb_treefind_variable(&pos[0], h, &startnode, 0, 0, 0);
#else
      numngb_inbox = ngb_treefind_variable(&pos[0], h, &startnode);
#endif

      for(n = 0; n < numngb_inbox; n++)
	{
	  j = Ngblist[n];

	  dx = pos[0] - P[j].Pos[0];
	  dy = pos[1] - P[j].Pos[1];
	  dz = pos[2] - P[j].Pos[2];

#ifdef PERIODIC			/*  now find the closest image in the given box size  */
#ifdef LONGBOX
	  if(dx > boxHalf * LONGBOX)
	    dx -= boxSize * LONGBOX;
	  if(dx < -boxHalf * LONGBOX)
	    dx += boxSize * LONGBOX;
#else
	  if(dx > boxHalf)
	    dx -= boxSize;
	  if(dx < -boxHalf)
	    dx += boxSize;
#endif
	  if(dy > boxHalf)
	    dy -= boxSize;
	  if(dy < -boxHalf)
	    dy += boxSize;
	  if(dz > boxHalf)
	    dz -= boxSize;
	  if(dz < -boxHalf)
	    dz += boxSize;
#endif
	  r2 = dx * dx + dy * dy + dz * dz;

	  if(r2 < h2)
	    {
	      r = sqrt(r2);

	      u = r * hinv;

	      if(u < 0.5)
		wk = 2.546479089470 + 15.278874536822 * (u - 1) * u * u;
	      else
		wk = 5.092958178941 * (1.0 - u) * (1.0 - u) * (1.0 - u);

	      mass_j = P[j].Mass;


#ifdef SFR_FEEDBACK

	      xhyd = P[j].Zm[6] / P[j].Mass;
	      yhel = (1 - xhyd) / (4. * xhyd);

	      ne = SphP[j].Ne;
	      mu = (1 + 4 * yhel) / (1 + yhel + ne);
	      energy = SphP[j].Entropy * P[j].Mass / GAMMA_MINUS1 * pow(SphP[j].Density * a3inv, GAMMA_MINUS1);	/* Total Energys */
	      temp = GAMMA_MINUS1 / BOLTZMANN * energy / P[j].Mass * PROTONMASS * mu;
	      temp *= All.UnitEnergy_in_cgs / All.UnitMass_in_g;	/* Temperature in Kelvin */

	      if(Flag_phase == 1)
		if(temp < All.Tcrit_Phase
		   && SphP[j].Density * a3inv > All.PhysDensThresh * All.DensFrac_Phase)
		  {
		    printf
		      (" ERROR  Phase=%d Part=%d  temp=%g K tempcrit=%g rho=%g internal rhocrit=%g ne=%g\n",
		       Flag_phase, P[j].ID, temp, All.Tcrit_Phase, SphP[j].Density * a3inv,
		       All.PhysDensThresh * All.DensFrac_Phase, ne);
		    fflush(stdout);
		    endrun(88912);	/* can't occur */
		  }

	      if(Flag_phase == 2)
		if(!(temp < All.Tcrit_Phase
		     && SphP[j].Density * a3inv > All.PhysDensThresh * All.DensFrac_Phase))
		  {
		    printf(" ERROR  Phase=%d  temp=%g K tempcrit=%g rho=%g internal rhocrit=%g \n",
			   Flag_phase, temp, All.Tcrit_Phase,
			   SphP[j].Density * a3inv, All.PhysDensThresh * All.DensFrac_Phase);
		    fflush(stdout);
		    endrun(888);	/* can't occur */
		  }
#endif


	      sum_metals = 0;

	      for(ik = 0; ik < 12; ik++)
		{
		  P[j].Zm[ik] += frac_phase * reservoir[ik] * wk * inv_wk_i;
		  sum_metals += frac_phase * reservoir[ik] * wk * inv_wk_i;
		  dist_metals += frac_phase * reservoir[ik] * wk * inv_wk_i;
		}


	      P[j].Mass += sum_metals;

#ifdef SFR_FEEDBACK
#ifdef SFR_PROMOTION
	      if(Flag_phase == 1)
		{
		  energy += frac_phase * energySN * wk * inv_wk_i;	/* total energy */

		  sum_distributed += frac_phase * energySN * wk * inv_wk_i;

		  SphP[j].Entropy =
		    energy / P[j].Mass * GAMMA_MINUS1 / pow(SphP[j].Density * a3inv, GAMMA_MINUS1);
		}

	      if(Flag_phase == 2)
		{
		  P[j].EnergySN += frac_phase * energySN * wk * inv_wk_i;	/* accumulate energy */

		  P[j].EnergySN += energySNCold * wk * inv_wk_i;	/* adds up reservoir from converted stars */


		  sum_distributed += frac_phase * energySN * wk * inv_wk_i;
		  sum_distributed += energySNCold * wk * inv_wk_i;


		  xhyd = P[j].Zm[6] / P[j].Mass;
		  yhel = (1 - xhyd) / (4. * xhyd);
		  ne = SphP[j].Ne;
		  mu = (1 + 4 * yhel) / (1 + yhel + ne);
		  new_temp = GAMMA_MINUS1 / BOLTZMANN * P[j].EnergySN / P[j].Mass * PROTONMASS * mu;
		  new_temp *= All.UnitEnergy_in_cgs / All.UnitMass_in_g;	/* Temperature in Kelvin */
/* 
		  if(target == 4382  || target == 4922)
		    {
		      printf("Phase=%d i=%d j=%d New-Res--T=%g New-Res-E=%g \n", Flag_phase, target, j, 
			     new_temp, P[j].EnergySN);
		      fflush(stdout);
		    } */
		}
#else
	      energy += frac_phase * energySN * wk * inv_wk_i;

	      SphP[j].Entropy =
		energy / P[j].Mass * GAMMA_MINUS1 / pow(SphP[j].Density * a3inv, GAMMA_MINUS1);
#endif
#endif
	    }
	}
    }
  while(startnode >= 0);


}

#ifdef SFR_PROMOTION
void promote_particles(void)
{
  int i;
  double a3inv, ne, mu, temp;
  double sum_reservoir_add = 0, total_reservoir_add;
  int n_promoted = 0, total_promoted;
  double sum_reservoir_promotion = 0, total_reservoir_promotion;
  double sum_reservoir_rest = 0, total_reservoir_rest;
  double new_egy, promotion_energy, hot_egy, cold_egy;
  double xhyd, yhel, energy;
  double entropyold, entropynew, critical_entropy;

  if(All.ComovingIntegrationOn)
    a3inv = 1 / (All.Time * All.Time * All.Time);
  else
    a3inv = 1;

  for(i = 0; i < NumPart; i++)
    if(P[i].Type == 0)
      sum_reservoir_add += P[i].EnergySN;

  Flag_promotion = 0;


  for(i = 0; i < N_gas; i++)
    {
      if(P[i].Type == 0 && Flag_promotion == 0)
	if(P[i].Ti_endstep == All.Ti_Current)
	  if(P[i].EnergySN > 0)
	    {
	      xhyd = P[i].Zm[6] / P[i].Mass;
	      yhel = (1 - xhyd) / (4. * xhyd);

	      ne = SphP[i].Ne;
	      mu = (1 + 4 * yhel) / (1 + yhel + ne);
	      energy = SphP[i].Entropy * P[i].Mass / GAMMA_MINUS1 * pow(SphP[i].Density * a3inv, GAMMA_MINUS1);	/* Total Energys */
	      temp = GAMMA_MINUS1 / BOLTZMANN * energy / P[i].Mass * PROTONMASS * mu;
	      temp *= All.UnitEnergy_in_cgs / All.UnitMass_in_g;	/* Temperature in Kelvin */


	      /* consider only cold phase particles for promotion */

	      if(temp < All.Tcrit_Phase && SphP[i].Density * a3inv > All.PhysDensThresh * All.DensFrac_Phase)
		{
		  cold_egy = SphP[i].Entropy / GAMMA_MINUS1 * pow(SphP[i].Density * a3inv, GAMMA_MINUS1);
		  hot_egy = SphP[i].EntropyAvg / GAMMA_MINUS1 * pow(SphP[i].DensityAvg * a3inv, GAMMA_MINUS1);

		  promotion_energy = GAMMA * P[i].Mass * (hot_egy - cold_egy);

		  critical_entropy =
		    All.Tcrit_Phase / pow(All.PhysDensThresh * All.DensFrac_Phase, GAMMA_MINUS1);
		  critical_entropy *= All.UnitMass_in_g / All.UnitEnergy_in_cgs * BOLTZMANN / PROTONMASS / mu;


		  new_egy = cold_egy + P[i].EnergySN / P[i].Mass;
		  entropynew = new_egy * GAMMA_MINUS1 / pow(SphP[i].DensityAvg * a3inv, GAMMA_MINUS1);

		  if(P[i].EnergySN > promotion_energy && SphP[i].HotNgbNum > 0 && SphP[i].EntropyAvg > critical_entropy && entropynew > SphP[i].EntropyAvg)	/* ok, we can promote the particle */
		    {

		      Flag_promotion = 1;

		      sum_reservoir_promotion += P[i].EnergySN;

		      n_promoted++;

		      entropyold = SphP[i].Entropy;

		      DEnergy_promotion = 0;	/* should be < 0 */

		      SphP[i].Entropy =
			new_egy * GAMMA_MINUS1 / pow(SphP[i].DensityAvg * a3inv, GAMMA_MINUS1);

		      PPP[i].Hsml *= pow(SphP[i].Density / SphP[i].DensityAvg, 1.0 / 3);
		      P[i].EnergySN = -1;


		      SphP[i].DensPromotion = SphP[i].Density;
		      SphP[i].TempPromotion = temp;

		      printf("Promoted Part=%d DensityOld=%g DensityAvg=%g EntropyOld=%g EntropyNew=%g\n",
			     P[i].ID, SphP[i].DensityOld, SphP[i].DensityAvg, entropyold, SphP[i].Entropy);
		      fflush(stdout);
		    }
		}
	      else
		{
		  /* thermalize energy in the reservoir */

		  SphP[i].Entropy +=
		    P[i].EnergySN / P[i].Mass * GAMMA_MINUS1 / pow(SphP[i].Density * a3inv, GAMMA_MINUS1);
		  P[i].EnergySN = 0;
		}
	    }
    }


  MPI_Reduce(&n_promoted, &total_promoted, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0 && total_promoted > 0)
    {
      printf("---> Promoted %d particles \n", total_promoted);
      fflush(stdout);
    }

  MPI_Reduce(&sum_reservoir_add, &total_reservoir_add, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    printf("EnergyTest total_reservoir_add=%g\n", total_reservoir_add);


  for(i = 0; i < NumPart; i++)
    if(P[i].Type == 0)
      sum_reservoir_rest += P[i].EnergySN;

  MPI_Reduce(&sum_reservoir_rest, &total_reservoir_rest, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sum_reservoir_promotion, &total_reservoir_promotion, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    printf("EnergyTest total_reservoir_rest=%g total_reservoir_promotion=%g\n", total_reservoir_rest,
	   total_reservoir_promotion);

}

#endif




/* This function computes the SNIa yields from Iwamoto et al. 1999 in units of solar mass */
double SNI_yields(int indice)
{
  float nSNI = 0;
  double delta_metals = 0;
  int ik = 0;
  double check = 0.;

  for(ik = 0; ik < 12; ik++)
    check += P[indice].ZmReservoir[ik];
  if(check > 0)
#ifdef SFR_FEEDBACK
    if(Flag_phase == 1)
#endif
      {
	printf("ERROR Part=%d enters SNI_yields  with Reservoir=%g\n", P[indice].ID, check);
	endrun(3534);
      }

  nSNI = P[indice].Mass * All.RateSNI;

  P[indice].ZmReservoir[0] = 0;	/*  He   */
  P[indice].ZmReservoir[1] = 4.83e-2 * nSNI;	/*  14C  */
  P[indice].ZmReservoir[2] = 8.50e-3 * nSNI;	/*  24Mg */
  P[indice].ZmReservoir[3] = 1.43e-1 * nSNI;	/*  16O  */
  P[indice].ZmReservoir[4] = 6.25e-1 * nSNI;	/*  56Fe */
  P[indice].ZmReservoir[5] = 1.54e-1 * nSNI;	/*  28Si */
  P[indice].ZmReservoir[6] = 0;	/*  H    */
  P[indice].ZmReservoir[7] = 1.16e-6 * nSNI;	/*  14N  */
  P[indice].ZmReservoir[8] = 2.02e-3 * nSNI;	/*  20Ne */
  P[indice].ZmReservoir[9] = 8.46e-2 * nSNI;	/*  32S  */
  P[indice].ZmReservoir[10] = 1.19e-2 * nSNI;	/*  40Ca */
  P[indice].ZmReservoir[11] = 0;	/*  62Zn */

  for(ik = 1; ik < 11; ik++)
    if(ik != 6)
      delta_metals += P[indice].ZmReservoir[ik];


  for(ik = 0; ik < 12; ik++)
    P[indice].Zm[ik] *= (1 - delta_metals / P[indice].Mass);

  P[indice].Mass -= delta_metals;


#ifdef SFR_FEEDBACK
  P[indice].EnergySN += nSNI * ESN / SOLAR_MASS * All.UnitMass_in_g / All.HubbleParam;	/* conversion to internal units */

  /*Energy Test */
  DEnergy_feedback += P[indice].EnergySN;
#endif

  return delta_metals;
}
#endif /* end of SFR_ENRICH */


/* This function calculates the mass in the different phases  */
#ifdef SFR_FEEDBACK
void phase_mass(void)
{
  int i;
  double mcold = 0;
  double mhot = 0;
  double mtot_hot, mtot_cold;
  double a3inv;
  double xhyd, yhel, ne, mu;
  double energy, temp;


  if(All.ComovingIntegrationOn)
    a3inv = 1 / (All.Time * All.Time * All.Time);
  else
    a3inv = 1;

  for(i = 0; i < N_gas; i++)
    {
#ifdef SFR
      if(P[i].Type == 0)
#endif
	{
	  xhyd = P[i].Zm[6] / P[i].Mass;
	  yhel = (1 - xhyd) / (4. * xhyd);

	  ne = SphP[i].Ne;
	  mu = (1 + 4 * yhel) / (1 + yhel + ne);
	  energy = SphP[i].Entropy / GAMMA_MINUS1 * pow(SphP[i].Density * a3inv, GAMMA_MINUS1);	/* Energy per mass unit */
	  temp = GAMMA_MINUS1 / BOLTZMANN * PROTONMASS * mu * energy;
	  temp *= All.UnitEnergy_in_cgs / All.UnitMass_in_g;	/* Temperature in Kelvin */

	  if(temp < All.Tcrit_Phase && SphP[i].Density * a3inv > All.PhysDensThresh * All.DensFrac_Phase)
	    mcold += P[i].Mass;
	  else
	    mhot += P[i].Mass;
	}
    }

  MPI_Reduce(&mhot, &mtot_hot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&mcold, &mtot_cold, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    fprintf(FdMphase, "%g %g %g %g \n", All.Time, mtot_cold + mtot_hot, mtot_cold, mtot_hot);
  fflush(FdMphase);

}
#endif /* SFR_FEEDBACK */



void energy_test(void)
{

  double Tot_TotalEnergy = 0;
  double Tot_DEnergy_spawned = 0;
  double Tot_DEnergy_converted = 0;
  double Tot_DEnergy_radiation = 0;
  double Tot_DEnergy_promotion = 0;
  double Tot_DEnergy_feedback = 0;
  double Tot_energy = 0;
  double TotalReservoir = 0, Tot_TotalReservoir = 0;
  double a3inv;
  int i;


  if(All.ComovingIntegrationOn)	/* Factors for comoving integration of hydro */
    {
      a3inv = 1 / (All.Time * All.Time * All.Time);
    }
  else
    a3inv = 1;


#ifdef SFR_FEEDBACK
  /*Energy Test - Not only for active particles! */
  for(i = 0; i < NumPart; i++)
    if(P[i].Type == 0)
      TotalReservoir += P[i].EnergySN;
#endif

  MPI_Reduce(&TotalEnergy, &Tot_TotalEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&DEnergy_spawned, &Tot_DEnergy_spawned, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&DEnergy_converted, &Tot_DEnergy_converted, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&DEnergy_radiation, &Tot_DEnergy_radiation, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&DEnergy_promotion, &Tot_DEnergy_promotion, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&DEnergy_feedback, &Tot_DEnergy_feedback, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&TotalReservoir, &Tot_TotalReservoir, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  Tot_energy = Tot_TotalEnergy + Tot_DEnergy_spawned + Tot_DEnergy_converted + Tot_DEnergy_radiation
    + Tot_DEnergy_promotion;

  if(ThisTask == 0)
    fprintf(FdMphase, "%g %g %g %g %g %g %g %g %g \n", All.Time, Tot_TotalEnergy, Tot_DEnergy_spawned,
	    Tot_DEnergy_converted, Tot_DEnergy_radiation, Tot_DEnergy_promotion, Tot_DEnergy_feedback,
	    Tot_TotalReservoir, Tot_energy);
  fflush(FdMphase);
}




/*! This routine is a comparison kernel used in a sort routine to group
 *  particles that are exported to the same processor.
 */
#ifdef SFR_ENRICH
int metals_compare_key(const void *a, const void *b)
{
  if(((struct metaldata_in *) a)->Task < (((struct metaldata_in *) b)->Task))
    return -1;

  if(((struct metaldata_in *) a)->Task > (((struct metaldata_in *) b)->Task))
    return +1;

  return 0;
}
#endif



#ifdef SFR_DECOUPLING
void copy_densities(void)
{
  int i;

  for(i = 0; i < N_gas; i++)
    if(P[i].Type == 0)
      SphP[i].DensityOld = SphP[i].Density;
}


void find_low_density_tail(void)
{
  float *rho, *rho_common;
  int i, count;

#define NUM_DENSITY_TAIL (2*All.DesNumNgb)

  rho = mymalloc(N_gas * sizeof(float));
  rho_common = mymalloc(NUM_DENSITY_TAIL * NTask * sizeof(float));

  for(i = 0, count = 0; i < N_gas; i++)
    if(P[i].Type == 0)
      rho[count++] = SphP[i].DensityOld;

  if(count < NUM_DENSITY_TAIL)
    {
      printf("number of gas particles on this cpu is less than %d!\n", NUM_DENSITY_TAIL);
      endrun(13123123);
    }

  qsort(rho, count, sizeof(float), compare_density_values);

  MPI_Allgather(rho, NUM_DENSITY_TAIL, MPI_FLOAT, rho_common, NUM_DENSITY_TAIL, MPI_FLOAT, MPI_COMM_WORLD);

  qsort(rho_common, NUM_DENSITY_TAIL * NTask, sizeof(float), compare_density_values);

  All.DensityTailThreshold = rho_common[NUM_DENSITY_TAIL - 1];

  myfree(rho_common);
  myfree(rho);
}


int compare_density_values(const void *a, const void *b)
{
  if(*((float *) a) < *((float *) b))
    return -1;

  if(*((float *) a) > *((float *) b))
    return +1;

  return 0;
}



#endif



#endif
