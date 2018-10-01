#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"

#ifdef SFR_METALS
#include "c_metals.h"
#endif

#ifdef COSMIC_RAYS
#include "cosmic_rays.h"
#endif

#ifdef BG_STELLAR_EVOLUTION
#include "bg_proto.h"
#endif

/*! \file density.c 
 *  \brief SPH density computation and smoothing length determination
 *
 *  This file contains the "first SPH loop", where the SPH densities and some
 *  auxiliary quantities are computed.  There is also functionality that
 *  corrects the smoothing length if needed.
 */


/*! This function computes the local density for each active SPH particle, the
 * number of neighbours in the current smoothing radius, and the divergence
 * and rotation of the velocity field.  The pressure is updated as well.  If a
 * particle with its smoothing region is fully inside the local domain, it is
 * not exported to the other processors. The function also detects particles
 * that have a number of neighbours outside the allowed tolerance range. For
 * these particles, the smoothing length is adjusted accordingly, and the
 * density() computation is called again.  Note that the smoothing length is
 * not allowed to fall below the lower bound set by MinGasHsml (this may mean
 * that one has to deal with substantially more than normal number of
 * neighbours.)
 */
void density(void)
{
  long long ntot, ntotleft;
  int ndone;
  int *noffset, *nbuffer, *nsend, *nsend_local, *numlist, *ndonelist;
  int i, j, n, NNN;
  int npleft;
  int maxfill, source, iter = 0;
  int level, ngrp, sendTask, recvTask;
  int place, nexport;
  double dt_entr, dmax1, dmax2, fac;
  double tstart, tend, tstart_ngb = 0, tend_ngb = 0;
  double sumt, sumcomm;
  double timecomp = 0, timeimbalance = 0, timecommsumm = 0, sumimbalance;
  double timengb, sumtimengb;
  double desnumngb;

#ifdef NAVIERSTOKES
  int k;
  double dvel[3][3];
  double rotx, roty, rotz;
#endif
#ifdef SFR_PROMOTION
  double xhyd, yhel, ne, mu, energy, temp;
#endif

  MPI_Status status;

#if defined(SOFTEREQS) || defined(SFR_PROMOTION)
  double a3inv;
#endif

#if defined(SOFTEREQS) || defined(SFR_PROMOTION)
  if(All.ComovingIntegrationOn)
    a3inv = 1 / (All.Time * All.Time * All.Time);
  else
    a3inv = 1;
#endif


  noffset = (int *) mymalloc(sizeof(int) * NTask);	/* offsets of bunches in common list */
  nbuffer = (int *) mymalloc(sizeof(int) * NTask);
  nsend_local = (int *) mymalloc(sizeof(int) * NTask);
  nsend = (int *) mymalloc(sizeof(int) * NTask * NTask);
  ndonelist = (int *) mymalloc(sizeof(int) * NTask);

#if defined(SFR_METALS) || defined(BLACK_HOLES) || defined(BG_SFR)
  NNN = NumPart;
#else
  NNN = N_gas;
#endif

  for(n = 0, NumSphUpdate = 0; n < NNN; n++)
    {
#ifdef BLACK_HOLES
      if(P[n].Type == 0 || P[n].Type == 5)
#else
#ifdef SFR
#ifdef SFR_METALS
      if(P[n].Type == 0 || P[n].Type == 6 || P[n].Type == 7)
#else
#ifdef BG_SFR
      if(P[n].Type == 0 || P[n].Type == 4)
#else
      if(P[n].Type == 0)
#endif
#endif
#endif
#endif
	{
	  PPP[n].Left = PPP[n].Right = 0;

	  if(P[n].Ti_endstep == All.Ti_Current)
	    NumSphUpdate++;

#ifdef BLACK_HOLES
	  P[n].SwallowID = 0;
#endif
#if defined(BLACK_HOLES) && defined(FLTROUNDOFFREDUCTION)
	  if(P[n].Type == 0)
	    SphP[n].i.dInjected_BH_Energy = SphP[n].i.Injected_BH_Energy;
#endif
	}
    }

  numlist = (int *) mymalloc(NTask * sizeof(int) * NTask);
  MPI_Allgather(&NumSphUpdate, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 0, ntot = 0; i < NTask; i++)
    ntot += numlist[i];
  myfree(numlist);


  desnumngb = All.DesNumNgb;

  /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
  do
    {
      i = 0;			/* begin with this index */
      ntotleft = ntot;		/* particles left for all tasks together */

      while(ntotleft > 0)
	{
	  for(j = 0; j < NTask; j++)
	    nsend_local[j] = 0;

	  /* do local particles and prepare export list */
	  tstart = second();
	  for(nexport = 0, ndone = 0; i < NNN && nexport < All.BunchSizeDensity - NTask; i++)
#ifdef BLACK_HOLES
	    if(P[i].Type == 0 || P[i].Type == 5)
#else
#ifdef SFR
#ifdef SFR_METALS
	    if(P[i].Type == 0 || P[i].Type == 6 || P[i].Type == 7)
#else
#ifdef BG_SFR
	    if(P[i].Type == 0 || P[i].Type == 4)
#else
	    if(P[i].Type == 0)
#endif
#endif
#endif
#endif
	      if(P[i].Ti_endstep == All.Ti_Current)
		{
		  ndone++;

		  for(j = 0; j < NTask; j++)
		    Exportflag[j] = 0;

		  density_evaluate(i, 0);

		  for(j = 0; j < NTask; j++)
		    {
		      if(Exportflag[j])
			{
			  DensDataIn[nexport].Pos[0] = P[i].Pos[0];
			  DensDataIn[nexport].Pos[1] = P[i].Pos[1];
			  DensDataIn[nexport].Pos[2] = P[i].Pos[2];

#if defined(SFR_METALS) || defined(BLACK_HOLES) || defined (BG_SFR)
			  if(P[i].Type == 0)
#endif
			    {
#ifdef SFR_DECOUPLING
			      DensDataIn[nexport].DensityOld = SphP[i].DensityOld;
			      DensDataIn[nexport].Entropy = SphP[i].Entropy;
#endif
			      DensDataIn[nexport].Vel[0] = SphP[i].VelPred[0];
			      DensDataIn[nexport].Vel[1] = SphP[i].VelPred[1];
			      DensDataIn[nexport].Vel[2] = SphP[i].VelPred[2];
#ifdef WINDS
			      DensDataIn[nexport].DelayTime = SphP[i].DelayTime;
#endif
			    }
#ifdef SFR_DECOUPLING
			  else
			    {
			      DensDataIn[nexport].DensityOld = 0;
			    }
#endif

#if defined(BLACK_HOLES)
			  if(P[i].Type == 5)
			    {
			      DensDataIn[nexport].Vel[0] = 0;
			      DensDataIn[nexport].Vel[1] = 0;
			      DensDataIn[nexport].Vel[2] = 0;
			    }
#endif

			  DensDataIn[nexport].Hsml = PPP[i].Hsml;
			  DensDataIn[nexport].Index = i;
			  DensDataIn[nexport].Task = j;
			  nexport++;
			  nsend_local[j]++;
			}
		    }
		}
	  tend = second();
	  timecomp += timediff(tstart, tend);

	  qsort(DensDataIn, nexport, sizeof(struct densdata_in), dens_compare_key);

	  for(j = 1, noffset[0] = 0; j < NTask; j++)
	    noffset[j] = noffset[j - 1] + nsend_local[j - 1];

	  tstart = second();

	  MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);

	  tend = second();
	  timeimbalance += timediff(tstart, tend);


	  /* now do the particles that need to be exported */

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
		  if(maxfill >= All.BunchSizeDensity)
		    break;

		  sendTask = ThisTask;
		  recvTask = ThisTask ^ ngrp;

		  if(recvTask < NTask)
		    {
		      if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
			{
			  /* get the particles */
			  MPI_Sendrecv(&DensDataIn[noffset[recvTask]],
				       nsend_local[recvTask] * sizeof(struct densdata_in), MPI_BYTE,
				       recvTask, TAG_DENS_A,
				       &DensDataGet[nbuffer[ThisTask]],
				       nsend[recvTask * NTask + ThisTask] * sizeof(struct densdata_in),
				       MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, &status);
			}
		    }

		  for(j = 0; j < NTask; j++)
		    if((j ^ ngrp) < NTask)
		      nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
		}
	      tend = second();
	      timecommsumm += timediff(tstart, tend);


	      tstart = second();
	      for(j = 0; j < nbuffer[ThisTask]; j++)
		{
		  density_evaluate(j, 1);
		}
	      tend = second();
	      timecomp += timediff(tstart, tend);

	      /* do a block to measure imbalance */
	      tstart = second();
	      MPI_Barrier(MPI_COMM_WORLD);
	      tend = second();
	      timeimbalance += timediff(tstart, tend);

	      /* get the result */
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
		  if(maxfill >= All.BunchSizeDensity)
		    break;

		  sendTask = ThisTask;
		  recvTask = ThisTask ^ ngrp;

		  if(recvTask < NTask)
		    {
		      if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
			{
			  /* send the results */
			  MPI_Sendrecv(&DensDataResult[nbuffer[ThisTask]],
				       nsend[recvTask * NTask + ThisTask] * sizeof(struct densdata_out),
				       MPI_BYTE, recvTask, TAG_DENS_B,
				       &DensDataPartialResult[noffset[recvTask]],
				       nsend_local[recvTask] * sizeof(struct densdata_out),
				       MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, &status);

			  /* add the result to the particles */
			  for(j = 0; j < nsend_local[recvTask]; j++)
			    {
			      source = j + noffset[recvTask];
			      place = DensDataIn[source].Index;

			      PPP[place].a1.dNumNgb += DensDataPartialResult[source].Ngb;

#if defined(LT_STELLAREVOLUTION) && !defined(LT_LOCAL_IRA) && defined(NOFIXEDMASSINKERNEL)
			      SphP[place].weight += DensDataPartialResult[source].weight;
#endif

#if defined(SFR_METALS) || defined(BLACK_HOLES) || defined (BG_SFR)
			      if(P[place].Type == 0)
#endif
				{
				  SphP[place].a2.dDensity += DensDataPartialResult[source].Rho;
#ifdef BG_SFR
				  SphP[place].Zsmooth += DensDataPartialResult[source].ZRho;
#endif
#ifndef NOGRADHSML
				  SphP[place].a3.dDhsmlDensityFactor +=
				    DensDataPartialResult[source].DhsmlDensity;
#endif
#ifndef NAVIERSTOKES
				  SphP[place].u.s.a4.dDivVel += DensDataPartialResult[source].Div;
				  SphP[place].u.s.a5.dRot[0] += DensDataPartialResult[source].Rot[0];
				  SphP[place].u.s.a5.dRot[1] += DensDataPartialResult[source].Rot[1];
				  SphP[place].u.s.a5.dRot[2] += DensDataPartialResult[source].Rot[2];
#else
				  for(k = 0; k < 3; k++)
				    {
				      SphP[place].u.DV[k][0] += DensDataPartialResult[source].DV[k][0];
				      SphP[place].u.DV[k][1] += DensDataPartialResult[source].DV[k][1];
				      SphP[place].u.DV[k][2] += DensDataPartialResult[source].DV[k][2];
				    }
#endif
#ifdef CONDUCTION
				  SphP[place].SmoothedEntr += DensDataPartialResult[source].SmoothedEntr;
#ifdef CONDUCTION_SATURATION
				  SphP[place].GradEntr[0] += DensDataPartialResult[source].GradEntr[0];
				  SphP[place].GradEntr[1] += DensDataPartialResult[source].GradEntr[1];
				  SphP[place].GradEntr[2] += DensDataPartialResult[source].GradEntr[2];
#endif
#endif
				}
#ifdef BLACK_HOLES
			      if(P[place].Type == 5)
				{
				  P[place].b1.dBH_Density += DensDataPartialResult[source].Rho;
				  P[place].b2.dBH_Entropy += DensDataPartialResult[source].SmoothedEntr;
				  P[place].b3.dBH_SurroundingGasVel[0] +=
				    DensDataPartialResult[source].GasVel[0];
				  P[place].b3.dBH_SurroundingGasVel[1] +=
				    DensDataPartialResult[source].GasVel[1];
				  P[place].b3.dBH_SurroundingGasVel[2] +=
				    DensDataPartialResult[source].GasVel[2];
				}
#endif
			    }
			}
		    }

		  for(j = 0; j < NTask; j++)
		    if((j ^ ngrp) < NTask)
		      nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
		}
	      tend = second();
	      timecommsumm += timediff(tstart, tend);

	      level = ngrp - 1;
	    }

	  MPI_Allgather(&ndone, 1, MPI_INT, ndonelist, 1, MPI_INT, MPI_COMM_WORLD);
	  for(j = 0; j < NTask; j++)
	    ntotleft -= ndonelist[j];
	}

#ifdef FLTROUNDOFFREDUCTION
      for(i = 0; i < NNN; i++)
	if(P[i].Ti_endstep == All.Ti_Current)
	  {
	    PPP[i].a1.NumNgb = FLT(PPP[i].a1.dNumNgb);

	    if(P[i].Type == 0)
	      {
		SphP[i].a2.Density = FLT(SphP[i].a2.dDensity);
		SphP[i].a3.DhsmlDensityFactor = FLT(SphP[i].a3.dDhsmlDensityFactor);
		SphP[i].u.s.a4.DivVel = FLT(SphP[i].u.s.a4.dDivVel);

		for(j = 0; j < 3; j++)
		  SphP[i].u.s.a5.Rot[j] = FLT(SphP[i].u.s.a5.dRot[j]);
	      }
#ifdef BLACK_HOLES
	    if(P[i].Type == 5)
	      {
		P[i].b1.BH_Density = FLT(P[i].b1.dBH_Density);
		P[i].b2.BH_Entropy = FLT(P[i].b2.dBH_Entropy);
		for(j = 0; j < 3; j++)
		  P[i].b3.BH_SurroundingGasVel[j] = FLT(P[i].b3.dBH_SurroundingGasVel[j]);
	      }
#endif
	  }
#endif


      /* do final operations on results */
      tstart = second();
      for(i = 0, npleft = 0; i < NNN; i++)
	{
#ifdef BLACK_HOLES
	  if(P[i].Type == 0 || P[i].Type == 5)
#else
#ifdef SFR
#ifdef SFR_METALS
	  if(P[i].Type == 0 || P[i].Type == 6 || P[i].Type == 7)
#else
#ifdef BG_SFR
	  if(P[i].Type == 0 || P[i].Type == 4)
#else
	  if(P[i].Type == 0)
#endif
#endif
#endif
#endif
	    if(P[i].Ti_endstep == All.Ti_Current)
	      {
#if defined(SFR_METALS) || defined(BLACK_HOLES) || defined (BG_SFR)
		if(P[i].Type == 0)
#endif
		  {

		    if(SphP[i].a2.Density > 0)
		      {
#ifndef NOGRADHSML
			SphP[i].a3.DhsmlDensityFactor *= PPP[i].Hsml / (NUMDIMS * SphP[i].a2.Density);
			if(SphP[i].a3.DhsmlDensityFactor > -0.9)	/* note: this would be -1 if only a single particle at zero lag is found */
			  SphP[i].a3.DhsmlDensityFactor = 1 / (1 + SphP[i].a3.DhsmlDensityFactor);
			else
			  SphP[i].a3.DhsmlDensityFactor = 1;

			/*
			   SphP[i].a3.DhsmlDensityFactor =
			   1 / (1 + PPP[i].Hsml * SphP[i].a3.DhsmlDensityFactor / (NUMDIMS * SphP[i].a2.Density));
			 */
#endif

#ifdef BG_SFR
			SphP[i].Zsmooth /= SphP[i].a2.Density;
#endif

#ifndef NAVIERSTOKES
			SphP[i].u.s.CurlVel = sqrt(SphP[i].u.s.a5.Rot[0] * SphP[i].u.s.a5.Rot[0] +
						   SphP[i].u.s.a5.Rot[1] * SphP[i].u.s.a5.Rot[1] +
						   SphP[i].u.s.a5.Rot[2] * SphP[i].u.s.a5.Rot[2]) /
			  SphP[i].a2.Density;

			SphP[i].u.s.a4.DivVel /= SphP[i].a2.Density;
#else
			for(k = 0; k < 3; k++)
			  {
			    dvel[k][0] = SphP[i].u.DV[k][0] / SphP[i].a2.Density;
			    dvel[k][1] = SphP[i].u.DV[k][1] / SphP[i].a2.Density;
			    dvel[k][2] = SphP[i].u.DV[k][2] / SphP[i].a2.Density;
			  }
			SphP[i].u.s.a4.DivVel = dvel[0][0] + dvel[1][1] + dvel[2][2];

			SphP[i].u.s.StressDiag[0] = 2 * dvel[0][0] - 2.0 / 3 * SphP[i].u.s.a4.DivVel;
			SphP[i].u.s.StressDiag[1] = 2 * dvel[1][1] - 2.0 / 3 * SphP[i].u.s.a4.DivVel;
			SphP[i].u.s.StressDiag[2] = 2 * dvel[2][2] - 2.0 / 3 * SphP[i].u.s.a4.DivVel;

			SphP[i].u.s.StressOffDiag[0] = dvel[0][1] + dvel[1][0];	/* xy */
			SphP[i].u.s.StressOffDiag[1] = dvel[0][2] + dvel[2][0];	/* xz */
			SphP[i].u.s.StressOffDiag[2] = dvel[1][2] + dvel[2][1];	/* yz */

#ifdef NAVIERSTOKES_BULK
			SphP[i].u.s.StressBulk = All.NavierStokes_BulkViscosity * SphP[i].u.s.a4.DivVel;
#endif

			rotx = dvel[1][2] - dvel[2][1];
			roty = dvel[2][0] - dvel[0][2];
			rotz = dvel[0][1] - dvel[1][0];
			SphP[i].u.s.CurlVel = sqrt(rotx * rotx + roty * roty + rotz * rotz);
#endif


#ifdef CONDUCTION
			SphP[i].SmoothedEntr /= SphP[i].a2.Density;
#ifdef CONDUCTION_SATURATION
			SphP[i].GradEntr[0] /= SphP[i].a2.Density;
			SphP[i].GradEntr[1] /= SphP[i].a2.Density;
			SphP[i].GradEntr[2] /= SphP[i].a2.Density;
#endif
#endif
		      }

		    dt_entr =
		      (All.Ti_Current - (P[i].Ti_begstep + P[i].Ti_endstep) / 2) * All.Timebase_interval;

#ifndef MHM
#ifndef SOFTEREQS
		    SphP[i].Pressure =
		      (SphP[i].Entropy + SphP[i].e.DtEntropy * dt_entr) * pow(SphP[i].a2.Density, GAMMA);
#else
		    /* use an intermediate EQS, between isothermal and the full multiphase model */
		    if(SphP[i].a2.Density * a3inv >= All.PhysDensThresh)
		      SphP[i].Pressure = All.FactorForSofterEQS *
			(SphP[i].Entropy + SphP[i].e.DtEntropy * dt_entr) * pow(SphP[i].a2.Density, GAMMA) +
			(1 - All.FactorForSofterEQS) * GAMMA_MINUS1 * SphP[i].a2.Density * All.InitGasU;
		    else
		      SphP[i].Pressure =
			(SphP[i].Entropy + SphP[i].e.DtEntropy * dt_entr) * pow(SphP[i].a2.Density, GAMMA);
#endif
#else
		    /* Here we use an isothermal equation of state */
		    SphP[i].Pressure = GAMMA_MINUS1 * SphP[i].a2.Density * All.InitGasU;
		    SphP[i].Entropy = SphP[i].Pressure / pow(SphP[i].a2.Density, GAMMA);
#endif

#ifdef COSMIC_RAYS
		    CR_Particle_Update(SphP + i);
#ifndef CR_NOPRESSURE
		    SphP[i].Pressure += CR_Comoving_Pressure(SphP + i);
#endif
#endif
		  }

#ifdef BLACK_HOLES
		if(P[i].Type == 5)
		  {
		    if(P[i].b1.BH_Density > 0)
		      {
			P[i].b2.BH_Entropy /= P[i].b1.BH_Density;
			P[i].b3.BH_SurroundingGasVel[0] /= P[i].b1.BH_Density;
			P[i].b3.BH_SurroundingGasVel[1] /= P[i].b1.BH_Density;
			P[i].b3.BH_SurroundingGasVel[2] /= P[i].b1.BH_Density;
		      }
		  }
#endif

		/* now check whether we had enough neighbours */

#ifdef BLACK_HOLES
		if(P[i].Type == 5)
		  desnumngb = All.DesNumNgb * All.BlackHoleNgbFactor;
		else
		  desnumngb = All.DesNumNgb;
#endif

		if(PPP[i].a1.NumNgb < (desnumngb - All.MaxNumNgbDeviation) ||
		   (PPP[i].a1.NumNgb > (desnumngb + All.MaxNumNgbDeviation)
		    && PPP[i].Hsml > (1.01 * All.MinGasHsml)))
		  {
		    /* need to redo this particle */
		    npleft++;

		    if(PPP[i].Left > 0 && PPP[i].Right > 0)
		      if((PPP[i].Right - PPP[i].Left) < 1.0e-3 * PPP[i].Left)
			{
			  /* this one should be ok */
			  npleft--;
			  P[i].Ti_endstep = -P[i].Ti_endstep - 1;	/* Mark as inactive */
			  continue;
			}

		    if(PPP[i].a1.NumNgb < (desnumngb - All.MaxNumNgbDeviation))
		      PPP[i].Left = DMAX(PPP[i].Hsml, PPP[i].Left);
		    else
		      {
			if(PPP[i].Right != 0)
			  {
			    if(PPP[i].Hsml < PPP[i].Right)
			      PPP[i].Right = PPP[i].Hsml;
			  }
			else
			  PPP[i].Right = PPP[i].Hsml;
		      }

		    if(iter >= MAXITER - 10)
		      {
			printf
			  ("i=%d task=%d ID=%d Hsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g\n   pos=(%g|%g|%g)\n",
			   i, ThisTask, (int) P[i].ID, PPP[i].Hsml, PPP[i].Left, PPP[i].Right,
			   (float) PPP[i].a1.NumNgb, PPP[i].Right - PPP[i].Left, P[i].Pos[0], P[i].Pos[1],
			   P[i].Pos[2]);
			fflush(stdout);
		      }

		    if(PPP[i].Right > 0 && PPP[i].Left > 0)
		      PPP[i].Hsml = pow(0.5 * (pow(PPP[i].Left, 3) + pow(PPP[i].Right, 3)), 1.0 / 3);
		    else
		      {
			if(PPP[i].Right == 0 && PPP[i].Left == 0)
			  endrun(8188);	/* can't occur */

			if(PPP[i].Right == 0 && PPP[i].Left > 0)
			  {
#ifndef NOGRADHSML
			    if(P[i].Type == 0 && fabs(PPP[i].a1.NumNgb - desnumngb) < 0.5 * desnumngb)
			      {
				fac = 1 - (PPP[i].a1.NumNgb -
					   desnumngb) / (NUMDIMS * PPP[i].a1.NumNgb) *
				  SphP[i].a3.DhsmlDensityFactor;

				if(fac < 1.26)
				  PPP[i].Hsml *= fac;
				else
				  PPP[i].Hsml *= 1.26;
			      }
			    else
			      PPP[i].Hsml *= 1.26;
#else
			    PPP[i].Hsml *= 1.26;
#endif
			  }

			if(PPP[i].Right > 0 && PPP[i].Left == 0)
			  {
#ifndef NOGRADHSML
			    if(P[i].Type == 0 && fabs(PPP[i].a1.NumNgb - desnumngb) < 0.5 * desnumngb)
			      {
				fac = 1 - (PPP[i].a1.NumNgb -
					   desnumngb) / (NUMDIMS * PPP[i].a1.NumNgb) *
				  SphP[i].a3.DhsmlDensityFactor;

				if(fac > 1 / 1.26)
				  PPP[i].Hsml *= fac;
				else
				  PPP[i].Hsml /= 1.26;
			      }
			    else
			      PPP[i].Hsml /= 1.26;
#else
			    PPP[i].Hsml /= 1.26;
#endif
			  }
		      }

		    if(PPP[i].Hsml < All.MinGasHsml)
		      PPP[i].Hsml = All.MinGasHsml;
		  }
		else
		  P[i].Ti_endstep = -P[i].Ti_endstep - 1;	/* Mark as inactive */
	      }
	}
      tend = second();
      timecomp += timediff(tstart, tend);


      numlist = (int *) mymalloc(NTask * sizeof(int) * NTask);
      MPI_Allgather(&npleft, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
      for(i = 0, ntot = 0; i < NTask; i++)
	ntot += numlist[i];
      myfree(numlist);

      if(ntot > 0)
	{
	  if(iter == 0)
	    tstart_ngb = second();

	  iter++;

	  if(iter > 0 && ThisTask == 0)
	    {
	      printf("ngb iteration %d: need to repeat for %d%09d particles.\n", iter,
		     (int) (ntot / 1000000000), (int) (ntot % 1000000000));
	      fflush(stdout);
	    }

	  if(iter > MAXITER)
	    {
	      printf("failed to converge in neighbour iteration in density()\n");
	      fflush(stdout);
	      endrun(1155);
	    }
	}
      else
	tend_ngb = second();
    }
  while(ntot > 0);


  for(i = 0; i < NNN; i++)
    {
#if defined(SFR_METALS) && defined(SFR_PROMOTION)
      if(P[i].Type == 0 && P[i].EnergySN < 0)	/* only gas particles should enter here */
	{
	  P[i].EnergySN = 0;
	  SphP[i].Entropy *= pow(SphP[i].a2.DensityAvg / SphP[i].a2.Density, GAMMA_MINUS1);

	  printf("Entropy Conservation: Part=%d factor=%g EnergySN=%g\n", P[i].ID,
		 pow(SphP[i].a2.DensityAvg / SphP[i].a2.Density, GAMMA_MINUS1), P[i].EnergySN);
	  fflush(stdout);
	}

      if(P[i].Type == 0 && (SphP[i].TempPromotion > 0 || SphP[i].DensPromotion > 0))
	{
	  xhyd = P[i].Zm[6] / P[i].Mass;
	  yhel = (1 - xhyd) / (4. * xhyd);
	  ne = SphP[i].Ne;
	  mu = (1 + 4 * yhel) / (1 + yhel + ne);
	  energy = SphP[i].Entropy * P[i].Mass / GAMMA_MINUS1 * pow(SphP[i].a2.Density * a3inv, GAMMA_MINUS1);	/* Total Energys */
	  temp = GAMMA_MINUS1 / BOLTZMANN * energy / P[i].Mass * PROTONMASS * mu;
	  temp *= All.UnitEnergy_in_cgs / All.UnitMass_in_g;	/* Temperature in Kelvin */

	  fprintf(FdPromotion, "%g %d %g %g %g %g %g %g\n", All.Time, P[i].ID, SphP[i].DensPromotion,
		  SphP[i].TempPromotion, SphP[i].a2.Density, temp, SphP[i].a2.DensityAvg, SphP[i].EntropyAvg);
	  fflush(FdPromotion);


	  SphP[i].TempPromotion = 0;
	  SphP[i].DensPromotion = 0;
	}
#endif

#ifdef BLACK_HOLES
      if(P[i].Type == 0 || P[i].Type == 5)
#else
#ifdef SFR
#ifdef SFR_METALS
      if(P[i].Type == 0 || P[i].Type == 6 || P[i].Type == 7)
#else
#ifdef BG_SFR
      if(P[i].Type == 0 || P[i].Type == 4)
#else
      if(P[i].Type == 0)
#endif
#endif
#endif
#endif
	{
	  /* mark as active again */
	  if(P[i].Ti_endstep < 0)
	    P[i].Ti_endstep = -P[i].Ti_endstep - 1;
	}
    }
  myfree(ndonelist);
  myfree(nsend);
  myfree(nsend_local);
  myfree(nbuffer);
  myfree(noffset);


  /* collect some timing information */

  CPU_Step[CPU_DENSITY] += timecomp;
  CPU_Step[CPU_DENSCOMM] += timecommsumm;

  All.Cadj_Cpu += timecomp;


  if(iter > 0)
    timengb = timediff(tstart_ngb, tend_ngb);
  else
    timengb = 0;

  MPI_Reduce(&timengb, &sumtimengb, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timecomp, &sumt, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timecommsumm, &sumcomm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timeimbalance, &sumimbalance, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      All.CPU_HydCompWalk += sumt / NTask;
      All.CPU_HydCommSumm += sumcomm / NTask;
      All.CPU_HydImbalance += sumimbalance / NTask;
      All.CPU_EnsureNgb += sumtimengb / NTask;
    }
}



/*! This function represents the core of the SPH density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
void density_evaluate(int target, int mode)
{
  int j, n;
  int startnode, numngb, numngb_inbox;
  double h, h2, fac, hinv, hinv3, hinv4;
  DOUBLE rho;
  double wk, dwk;
  double dx, dy, dz, r, r2, u, mass_j;
  double dvx, dvy, dvz;

#ifdef BG_SFR
  DOUBLE Zrho = 0;
#endif

#ifndef NAVIERSTOKES
  DOUBLE divv, rotv[3];
#else
  int k;
  double dvel[3][3];
#endif

  FLOAT *pos, *vel;

#if defined(SFR_METALS) || defined(BLACK_HOLES)
  static FLOAT veldummy[3] = { 0, 0, 0 };
#endif

#if defined(CONDUCTION) || defined(BLACK_HOLES)
  DOUBLE smoothentr;
  smoothentr = 0;

#ifdef CONDUCTION_SATURATION
  double gradentr[3];
#endif
#endif
#ifdef WINDS
  double delaytime;
#endif
#ifndef NOFIXEDMASSINKERNEL
  DOUBLE weighted_numngb;
#endif
#ifndef NOGRADHSML
  DOUBLE dhsmlrho;
#endif
#ifdef SFR_DECOUPLING
  double densityold, entropy;
#endif
#ifdef BLACK_HOLES
  DOUBLE gasvel[3];
#endif

#if defined(LT_STELLAREVOLUTION) && !defined(LT_LOCAL_IRA) && defined(NOFIXEDMASSINKERNEL)
  DOUBLE weight_sum = 0;
#endif

#ifndef NAVIERSTOKES
  divv = rotv[0] = rotv[1] = rotv[2] = 0;
#else
  for(k = 0; k < 3; k++)
    dvel[k][0] = dvel[k][1] = dvel[k][2] = 0;
#endif

#if defined(CONDUCTION) || defined(BLACK_HOLES)
#ifdef CONDUCTION_SATURATION
  gradentr[0] = gradentr[1] = gradentr[2] = 0;
#endif
#endif

#ifdef BLACK_HOLES
  gasvel[0] = gasvel[1] = gasvel[2] = 0;
#endif
#ifndef NOFIXEDMASSINKERNEL
  weighted_numngb = 0;
#endif
#ifndef NOGRADHSML
  dhsmlrho = 0;
#endif

  rho = 0;

  if(mode == 0)
    {
      pos = P[target].Pos;
#if defined(SFR_METALS) || defined(BLACK_HOLES)
      vel = veldummy;
      if(P[target].Type == 0)
#endif
	{
	  vel = SphP[target].VelPred;
#ifdef WINDS
	  delaytime = SphP[target].DelayTime;
#endif
	}
      h = PPP[target].Hsml;
#ifdef SFR_DECOUPLING
      if(P[target].Type == 0)
	{
	  densityold = SphP[target].a2.DensityOld;
	  entropy = SphP[target].Entropy;
	}
      else
	densityold = 0;
#endif
    }
  else
    {
      pos = DensDataGet[target].Pos;
      vel = DensDataGet[target].Vel;
      h = DensDataGet[target].Hsml;
#ifdef WINDS
      delaytime = DensDataGet[target].DelayTime;
#endif
#ifdef SFR_DECOUPLING
      densityold = DensDataGet[target].a2.DensityOld;
      entropy = DensDataGet[target].Entropy;
#endif
    }


  h2 = h * h;
  hinv = 1.0 / h;
#ifndef  TWODIMS
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv * hinv / boxSize_Z;
#endif
  hinv4 = hinv3 * hinv;

  startnode = All.MaxPart;
  numngb = 0;
  do
    {
#ifdef SFR_DECOUPLING
      numngb_inbox = ngb_treefind_variable(&pos[0], h, &startnode, densityold, entropy, &vel[0]);
#else
      numngb_inbox = ngb_treefind_variable(&pos[0], h, &startnode);
#endif

      for(n = 0; n < numngb_inbox; n++)
	{
	  j = Ngblist[n];
#ifdef WINDS
	  if(SphP[j].DelayTime > 0)	/* partner is a wind particle */
	    if(!(delaytime > 0))	/* if I'm not wind, then ignore the wind particle */
	      continue;
#endif
#if defined(BLACK_HOLES) || defined(MYSWITCH)
	  if(P[j].Mass == 0)
	    continue;
#endif
	  dx = pos[0] - P[j].Pos[0];
	  dy = pos[1] - P[j].Pos[1];
	  dz = pos[2] - P[j].Pos[2];

#ifdef PERIODIC			/*  now find the closest image in the given box size  */
	  if(dx > boxHalf_X)
	    dx -= boxSize_X;
	  if(dx < -boxHalf_X)
	    dx += boxSize_X;
	  if(dy > boxHalf_Y)
	    dy -= boxSize_Y;
	  if(dy < -boxHalf_Y)
	    dy += boxSize_Y;
	  if(dz > boxHalf_Z)
	    dz -= boxSize_Z;
	  if(dz < -boxHalf_Z)
	    dz += boxSize_Z;
#endif
	  r2 = dx * dx + dy * dy + dz * dz;

	  if(r2 < h2)
	    {
	      numngb++;

	      r = sqrt(r2);

	      u = r * hinv;

	      if(u < 0.5)
		{
		  wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
		  dwk = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
		}
	      else
		{
		  wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
		  dwk = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
		}

	      mass_j = P[j].Mass;

	      rho += FLT(mass_j * wk);

#ifdef BG_SFR
#ifdef BG_STELLAR_EVOLUTION
	      Zrho += bg_get_metallicity(j) * mass_j * wk;
#endif
#endif


#if defined(CONDUCTION) || defined(BLACK_HOLES)
	      smoothentr += FLT(mass_j * wk * SphP[j].Entropy);
#ifdef CONDUCTION_SATURATION
	      if(r > 0)
		{
		  gradentr[0] += mass_j * dwk * dx / r * SphP[j].Entropy;
		  gradentr[1] += mass_j * dwk * dy / r * SphP[j].Entropy;
		  gradentr[2] += mass_j * dwk * dz / r * SphP[j].Entropy;
		}
#endif
#endif

#ifdef BLACK_HOLES
	      gasvel[0] += FLT(mass_j * wk * SphP[j].VelPred[0]);
	      gasvel[1] += FLT(mass_j * wk * SphP[j].VelPred[1]);
	      gasvel[2] += FLT(mass_j * wk * SphP[j].VelPred[2]);
#endif


#ifndef NOFIXEDMASSINKERNEL
	      weighted_numngb += FLT(NORM_COEFF * wk / hinv3);	/* 4.0/3 * PI = 4.188790204786 */
#if defined(LT_STELLAREVOLUTION) && !defined(LT_LOCAL_IRA) && defined(NOFIXEDMASSINKERNEL)
	      weight_sum += wk;
#endif
#endif

#ifndef NOGRADHSML
	      dhsmlrho += FLT(-mass_j * (NUMDIMS * hinv * wk + u * dwk));
#endif
	      if(r > 0)
		{
		  fac = mass_j * dwk / r;

		  dvx = vel[0] - SphP[j].VelPred[0];
		  dvy = vel[1] - SphP[j].VelPred[1];
		  dvz = vel[2] - SphP[j].VelPred[2];
#ifndef NAVIERSTOKES
		  divv += FLT(-fac * (dx * dvx + dy * dvy + dz * dvz));

		  rotv[0] += FLT(fac * (dz * dvy - dy * dvz));
		  rotv[1] += FLT(fac * (dx * dvz - dz * dvx));
		  rotv[2] += FLT(fac * (dy * dvx - dx * dvy));
#else
		  dvel[0][0] -= fac * dx * dvx;
		  dvel[0][1] -= fac * dx * dvy;
		  dvel[0][2] -= fac * dx * dvz;
		  dvel[1][0] -= fac * dy * dvx;
		  dvel[1][1] -= fac * dy * dvy;
		  dvel[1][2] -= fac * dy * dvz;
		  dvel[2][0] -= fac * dz * dvx;
		  dvel[2][1] -= fac * dz * dvy;
		  dvel[2][2] -= fac * dz * dvz;
#endif
		}
	    }
	}
    }
  while(startnode >= 0);


  if(mode == 0)
    {
#ifndef NOFIXEDMASSINKERNEL
      PPP[target].a1.dNumNgb = weighted_numngb;
#else
      PPP[target].a1.dNumNgb = numngb;
#if defined(LT_STELLAREVOLUTION) && !defined(LT_LOCAL_IRA) && defined(NOFIXEDMASSINKERNEL)
      if(P[target].Type == 0)
        SphP[target].weight = weight_sum;
#endif
#endif

#if defined(SFR_METALS) || defined(BLACK_HOLES) || defined (BG_SFR)
      if(P[target].Type == 0)
#endif
	{
	  SphP[target].a2.dDensity = rho;
#ifdef BG_SFR
	  SphP[target].Zsmooth = Zrho;
#endif

#ifdef CONDUCTION
	  SphP[target].SmoothedEntr = smoothentr;
#ifdef CONDUCTION_SATURATION
	  SphP[target].GradEntr[0] = gradentr[0];
	  SphP[target].GradEntr[1] = gradentr[1];
	  SphP[target].GradEntr[2] = gradentr[2];
#endif
#endif


#ifndef NOGRADHSML
	  SphP[target].a3.dDhsmlDensityFactor = dhsmlrho;
#endif

#ifndef NAVIERSTOKES
	  SphP[target].u.s.a4.dDivVel = divv;
	  SphP[target].u.s.a5.dRot[0] = rotv[0];
	  SphP[target].u.s.a5.dRot[1] = rotv[1];
	  SphP[target].u.s.a5.dRot[2] = rotv[2];
#else
	  for(k = 0; k < 3; k++)
	    {
	      SphP[target].u.DV[k][0] = dvel[k][0];
	      SphP[target].u.DV[k][1] = dvel[k][1];
	      SphP[target].u.DV[k][2] = dvel[k][2];
	    }
#endif
	}

#ifdef BLACK_HOLES
      P[target].b1.dBH_Density = rho;
      P[target].b2.dBH_Entropy = smoothentr;
      P[target].b3.dBH_SurroundingGasVel[0] = gasvel[0];
      P[target].b3.dBH_SurroundingGasVel[1] = gasvel[1];
      P[target].b3.dBH_SurroundingGasVel[2] = gasvel[2];
#endif
    }
  else
    {
      DensDataResult[target].Rho = rho;
#ifdef BG_SFR
      DensDataResult[target].ZRho = Zrho;
#endif
#ifndef NAVIERSTOKES
      DensDataResult[target].Div = divv;
      DensDataResult[target].Rot[0] = rotv[0];
      DensDataResult[target].Rot[1] = rotv[1];
      DensDataResult[target].Rot[2] = rotv[2];
#else
      for(k = 0; k < 3; k++)
	{
	  DensDataResult[target].DV[k][0] = dvel[k][0];
	  DensDataResult[target].DV[k][1] = dvel[k][1];
	  DensDataResult[target].DV[k][2] = dvel[k][2];
	}
#endif

#if defined(CONDUCTION) || defined(BLACK_HOLES)
      DensDataResult[target].SmoothedEntr = smoothentr;
#ifdef CONDUCTION_SATURATION
      DensDataResult[target].GradEntr[0] = gradentr[0];
      DensDataResult[target].GradEntr[1] = gradentr[1];
      DensDataResult[target].GradEntr[2] = gradentr[2];
#endif
#endif

#ifndef NOFIXEDMASSINKERNEL
      DensDataResult[target].Ngb = weighted_numngb;
#else
      DensDataResult[target].Ngb = numngb;
#if defined(LT_STELLAREVOLUTION) && !defined(LT_LOCAL_IRA) && defined(NOFIXEDMASSINKERNEL)
      DensDataResult[target].weight = weight_sum;
#endif
#endif

#ifndef NOGRADHSML
      DensDataResult[target].DhsmlDensity = dhsmlrho;
#endif

#ifdef BLACK_HOLES
      DensDataResult[target].GasVel[0] = gasvel[0];
      DensDataResult[target].GasVel[1] = gasvel[1];
      DensDataResult[target].GasVel[2] = gasvel[2];
#endif
    }
}





/*! This routine is a comparison kernel used in a sort routine to group
 *  particles that are exported to the same processor.
 */
int dens_compare_key(const void *a, const void *b)
{
  if(((struct densdata_in *) a)->Task < (((struct densdata_in *) b)->Task))
    return -1;

  if(((struct densdata_in *) a)->Task > (((struct densdata_in *) b)->Task))
    return +1;

  return 0;
}


#ifdef NAVIERSTOKES
double get_shear_viscosity(int i)
{
  return All.NavierStokes_ShearViscosity; 
}
#endif

