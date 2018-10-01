#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>


#include "allvars.h"
#include "proto.h"

#if defined(SFR_METALS) && defined(SFR_DECOUPLING) && defined(SFR_PROMOTION)

#include "c_metals.h"


void find_hot_neighbours(void)
{
  int *noffset, *nbuffer, *nsend, *nsend_local;
  int i, j, n;
  int ndone, ndonetot, ntot, ntotleft, npleft;
  int maxfill, source, iter = 0;
  int level, ngrp, sendTask, recvTask;
  int place, nexport;
  double tstart, tend, tstart_ngb = 0, tend_ngb = 0;
  double sumt, sumcomm, dmax1, dmax2;
  double timecomp = 0, timeimbalance = 0, timecommsumm = 0, sumimbalance;
  double timengb, sumtimengb;
  MPI_Status status;
  double xhyd, yhel, ne, mu, energy, temp;
  double a3inv;


  if(All.ComovingIntegrationOn)
    a3inv = 1 / (All.Time * All.Time * All.Time);
  else
    a3inv = 1;


  noffset = mymalloc(sizeof(int) * NTask);	/* offsets of bunches in common list */
  nbuffer = mymalloc(sizeof(int) * NTask);
  nsend_local = mymalloc(sizeof(int) * NTask);
  nsend = mymalloc(sizeof(int) * NTask * NTask);


  for(n = 0, NumSphUpdate = 0; n < N_gas; n++)
    {
      if(P[n].Ti_endstep == All.Ti_Current)
	if(P[n].Type == 0)
	  {
	    /* select reservoir and cold phase particles */
	    if(P[n].EnergySN > 0 && SphP[n].Density * a3inv > All.PhysDensThresh * All.DensFrac_Phase)
	      {
		xhyd = P[n].Zm[6] / P[n].Mass;
		yhel = (1 - xhyd) / (4. * xhyd);

		ne = SphP[n].Ne;
		mu = (1 + 4 * yhel) / (1 + yhel + ne);
		energy = SphP[n].Entropy * P[n].Mass / GAMMA_MINUS1 * pow(SphP[n].Density * a3inv, GAMMA_MINUS1);	/* Total Energys */
		temp = GAMMA_MINUS1 / BOLTZMANN * energy / P[n].Mass * PROTONMASS * mu;
		temp *= All.UnitEnergy_in_cgs / All.UnitMass_in_g;	/* Temperature in Kelvin */

		if(temp < All.Tcrit_Phase)
		  {
		    PPP[n].Left = PPP[n].Right = 0;
		    SphP[n].HotHsml = PPP[n].Hsml;
		    P[n].Type = 10;	/* temporarily mark particles of interest with this number */
		    NumSphUpdate++;
		  }
	      }
	  }
    }

  MPI_Allreduce(&NumSphUpdate, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);



  /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
  do
    {
      i = 0;			/* beginn with this index */
      ntotleft = ntot;		/* particles left for all tasks together */

      while(ntotleft > 0)
	{
	  for(j = 0; j < NTask; j++)
	    nsend_local[j] = 0;

	  /* do local particles and prepare export list */
	  tstart = second();
	  for(nexport = 0, ndone = 0; i < N_gas && nexport < All.BunchSizeHotNgbs - NTask; i++)
	    if(P[i].Ti_endstep == All.Ti_Current && P[i].Type == 10)
	      {
		ndone++;

		for(j = 0; j < NTask; j++)
		  Exportflag[j] = 0;

		hotngbs_evaluate(i, 0);

		for(j = 0; j < NTask; j++)
		  {
		    if(Exportflag[j])
		      {
			HotNgbsIn[nexport].Pos[0] = P[i].Pos[0];
			HotNgbsIn[nexport].Pos[1] = P[i].Pos[1];
			HotNgbsIn[nexport].Pos[2] = P[i].Pos[2];
			HotNgbsIn[nexport].HotHsml = SphP[i].HotHsml;
			HotNgbsIn[nexport].Entropy = SphP[i].Entropy;

			HotNgbsIn[nexport].Index = i;
			HotNgbsIn[nexport].Task = j;
			nexport++;
			nsend_local[j]++;
		      }
		  }
		tend = second();
		timecomp += timediff(tstart, tend);
	      }

	  qsort(HotNgbsIn, nexport, sizeof(struct hotngbs_in), hotngbs_compare_key);

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
		  if(maxfill >= All.BunchSizeHotNgbs)
		    break;

		  sendTask = ThisTask;
		  recvTask = ThisTask ^ ngrp;

		  if(recvTask < NTask)
		    {
		      if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
			{
			  /* get the particles */
			  MPI_Sendrecv(&HotNgbsIn[noffset[recvTask]],
				       nsend_local[recvTask] * sizeof(struct hotngbs_in), MPI_BYTE,
				       recvTask, TAG_HOTNGB_A,
				       &HotNgbsGet[nbuffer[ThisTask]],
				       nsend[recvTask * NTask + ThisTask] * sizeof(struct hotngbs_in),
				       MPI_BYTE, recvTask, TAG_HOTNGB_A, MPI_COMM_WORLD, &status);
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
		hotngbs_evaluate(j, 1);
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
		  if(maxfill >= All.BunchSizeHotNgbs)
		    break;

		  sendTask = ThisTask;
		  recvTask = ThisTask ^ ngrp;

		  if(recvTask < NTask)
		    {
		      if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
			{
			  /* send the results */
			  MPI_Sendrecv(&HotNgbsResult[nbuffer[ThisTask]],
				       nsend[recvTask * NTask + ThisTask] * sizeof(struct hotngbs_out),
				       MPI_BYTE, recvTask, TAG_HOTNGB_B,
				       &HotNgbsPartialResult[noffset[recvTask]],
				       nsend_local[recvTask] * sizeof(struct hotngbs_out),
				       MPI_BYTE, recvTask, TAG_HOTNGB_B, MPI_COMM_WORLD, &status);

			  /* add the result to the particles */
			  for(j = 0; j < nsend_local[recvTask]; j++)
			    {
			      source = j + noffset[recvTask];
			      place = HotNgbsIn[source].Index;

			      SphP[place].DensityAvg += HotNgbsPartialResult[source].DensitySum;
			      SphP[place].EntropyAvg += HotNgbsPartialResult[source].EntropySum;
			      SphP[place].HotNgbNum += HotNgbsPartialResult[source].HotNgbNum;
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

	  MPI_Allreduce(&ndone, &ndonetot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	  ntotleft -= ndonetot;

	}

      /* do final operations on results */
      tstart = second();
      for(i = 0, npleft = 0; i < N_gas; i++)
	{
	  if(P[i].Ti_endstep == All.Ti_Current && P[i].Type == 10)
	    {
	      if(SphP[i].HotNgbNum > 0)
		{
		  SphP[i].DensityAvg /= SphP[i].HotNgbNum;
		  SphP[i].EntropyAvg /= SphP[i].HotNgbNum;
		}
	      else
		{
		  SphP[i].DensityAvg = 0;
		  SphP[i].EntropyAvg = 0;
		}

	      /* now check whether we had enough neighbours */

	      if(SphP[i].HotNgbNum < (All.DesNumNgb - All.MaxNumNgbDeviation) ||
		 (SphP[i].HotNgbNum > (All.DesNumNgb + All.MaxNumNgbDeviation)))
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

		  if(SphP[i].HotNgbNum < (All.DesNumNgb - All.MaxNumNgbDeviation))
		    PPP[i].Left = DMAX(SphP[i].HotHsml, PPP[i].Left);
		  else
		    {
		      if(PPP[i].Right != 0)
			{
			  if(SphP[i].HotHsml < PPP[i].Right)
			    PPP[i].Right = SphP[i].HotHsml;
			}
		      else
			PPP[i].Right = SphP[i].HotHsml;
		    }

		  if(PPP[i].Left > 10 * PPP[i].Hsml)	/* prevent us from searching too far */
		    {
		      npleft--;
		      P[i].Ti_endstep = -P[i].Ti_endstep - 1;	/* Mark as inactive */

		      if(SphP[i].HotNgbNum == 0)
			{
			  SphP[i].DensityAvg = SphP[i].Density / 100;
			  SphP[i].EntropyAvg = SphP[i].Entropy * 1000;
			}

		      continue;
		    }

		  if(iter >= MAXITER - 10)
		    {
		      printf
			("i=%d task=%d ID=%d Hsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g\n   pos=(%g|%g|%g)\n",
			 i, ThisTask, P[i].ID, SphP[i].HotHsml, PPP[i].Left, PPP[i].Right,
			 (float) SphP[i].HotNgbNum, PPP[i].Right - PPP[i].Left, P[i].Pos[0], P[i].Pos[1],
			 P[i].Pos[2]);
		      fflush(stdout);
		    }

		  if(PPP[i].Right > 0 && PPP[i].Left > 0)
		    SphP[i].HotHsml = pow(0.5 * (pow(PPP[i].Left, 3) + pow(PPP[i].Right, 3)), 1.0 / 3);
		  else
		    {
		      if(PPP[i].Right == 0 && PPP[i].Left == 0)
			endrun(8188);	/* can't occur */

		      if(PPP[i].Right == 0 && PPP[i].Left > 0)
			SphP[i].HotHsml *= 1.26;

		      if(PPP[i].Right > 0 && PPP[i].Left == 0)
			SphP[i].HotHsml /= 1.26;
		    }
		}
	      else
		P[i].Ti_endstep = -P[i].Ti_endstep - 1;	/* Mark as inactive */
	    }
	}
      tend = second();
      timecomp += timediff(tstart, tend);


      MPI_Allreduce(&npleft, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      if(ntot > 0)
	{
	  if(iter == 0)
	    tstart_ngb = second();

	  iter++;

	  if(iter > 0 && ThisTask == 0)
	    {
	      printf("hotngb iteration %d: need to repeat for %d particles.\n", iter, ntot);
	      fflush(stdout);
	    }

	  if(iter > MAXITER)
	    {
	      printf("failed to converge in hot-neighbour iteration\n");
	      fflush(stdout);
	      endrun(1155);
	    }
	}
      else
	tend_ngb = second();
    }
  while(ntot > 0);


  for(i = 0; i < N_gas; i++)
    if(P[i].Type == 10)
      {
	P[i].Type = 0;

	/* mark as active again */
	if(P[i].Ti_endstep < 0)
	  P[i].Ti_endstep = -P[i].Ti_endstep - 1;
      }

  myfree(nsend);
  myfree(nsend_local);
  myfree(nbuffer);
  myfree(noffset);


  /* collect some timing information */

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



void hotngbs_evaluate(int target, int mode)
{
  int j, n;
  int startnode, numngb, numngb_inbox;
  double h, h2;
  double dx, dy, dz, r2;
  FLOAT pos[3];
  FLOAT entropy;
  double density_sum;
  double entropy_sum;
  int hotngb_count;


#ifdef PERIODIC
  double boxSize, boxHalf;

  boxSize = All.BoxSize;
  boxHalf = 0.5 * All.BoxSize;
#endif


  density_sum = 0;
  entropy_sum = 0;
  hotngb_count = 0;

  if(mode == 0)
    {
      pos[0] = P[target].Pos[0];
      pos[1] = P[target].Pos[1];
      pos[2] = P[target].Pos[2];
      entropy = SphP[target].Entropy;
      h = SphP[target].HotHsml;
    }
  else
    {
      pos[0] = HotNgbsGet[target].Pos[0];
      pos[1] = HotNgbsGet[target].Pos[1];
      pos[2] = HotNgbsGet[target].Pos[2];
      entropy = HotNgbsGet[target].Entropy;
      h = HotNgbsGet[target].HotHsml;
    }


  h2 = h * h;

  startnode = All.MaxPart;
  numngb = 0;
  do
    {
      numngb_inbox = ngb_treefind_hotngbs(&pos[0], h, &startnode, entropy);

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
	      density_sum += SphP[j].Density;
	      entropy_sum += SphP[j].Entropy;
	      hotngb_count++;
	    }
	}
    }
  while(startnode >= 0);


  if(mode == 0)
    {
      SphP[target].DensityAvg = density_sum;
      SphP[target].EntropyAvg = entropy_sum;
      SphP[target].HotNgbNum = hotngb_count;
    }
  else
    {
      HotNgbsResult[target].DensitySum = density_sum;
      HotNgbsResult[target].EntropySum = entropy_sum;
      HotNgbsResult[target].HotNgbNum = hotngb_count;
    }
}



/*! This routine is a comparison kernel used in a sort routine to group
 *  particles that are exported to the same processor.
 */
int hotngbs_compare_key(const void *a, const void *b)
{
  if(((struct hotngbs_in *) a)->Task < (((struct hotngbs_in *) b)->Task))
    return -1;

  if(((struct hotngbs_in *) a)->Task > (((struct hotngbs_in *) b)->Task))
    return +1;

  return 0;
}


#endif
