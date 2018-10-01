#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "allvars.h"
#include "proto.h"
#include "domain.h"

#include <gsl/gsl_heapsort.h>


/*! \file domain.c
 *  \brief code for domain decomposition
 *
 *  This file contains the code for the domain decomposition of the simulation
 *  volume.  The domains are constructed from a grid of dimension DOMAINGRID
 *  that covers the volume. These grid-cells correspond to a certain level of
 *  a global BH tree constructed for the full simulation volume. Domain
 *  boundaries hence run along tree-node divisions of a fiducial global BH
 *  tree. As a result of this method, the tree force can be made to be
 *  strictly independent of the way the domains are cut. The domain
 *  decomposition can be carried out for an arbitrary number of
 *  CPUs. Individual domains are not cubical, but are rather an arbitrary
 *  collection of cubical nodes drawn from the domain grid. However, it is
 *  advantageous to have spatial proximity of the domain cells assigned to a
 *  given CPU such that their union has a small surface to volume ratio. To
 *  achieve this, domains are cut along a space-filling Peano-Hilbert curve
 *  that traverses the domain grid.
 */

#define TOPNODEFACTOR  40.0

#define REDUC_FAC      0.98



/*! toGo[task*NTask + partner] gives the number of particles in task 'task'
 *  that have to go to 'partner'
 */
static int *toGo, *toGoSph;
static int *local_toGo, *local_toGoSph;
static int *list_NumPart;
static int *list_N_gas;
static int *list_load;
static int *list_loadsph;
static double *list_work;
static double *list_speedfac;
static double *list_cadj_cpu;
static double *list_cadj_cost;

#ifdef LT_STELLAREVOLUTION
static int *local_toGoStar, *toGoStar;
static int *list_N_star;
static int *list_loadstar;
static int send_countstar, recv_countstar;
#endif

static long long maxload, maxloadsph;

static double totgravcost, gravcost;

#ifdef LT_STELLAREVOLUTION
static long long maxloadstar;
#endif

static struct topnode_exchange
{
  peanokey Startkey;
  int Count;
  float GravCost;
}
 *toplist, *toplist_local;



/*! This is the main routine for the domain decomposition. It will try
 * to balance the work-load as defined by the sum of the P[i]-GravCost
 * in each domain.  The decomposition will respect the maximum
 * memory-imbalance given by the value of PartAllocFactor (it this is
 * possible).  
 */
void DomainDecomposition(void)
{
  double t0, t1;

#ifdef PMGRID
  if(All.PM_Ti_endstep == All.Ti_Current)
    {
      All.NumForcesSinceLastDomainDecomp = (long long) (1 + All.TotNumPart * All.TreeDomainUpdateFrequency);
      /* to make sure that we do a domain decomposition before the PM-force is evaluated.
         this is needed to make sure that the particles are wrapped into the box */
    }
#endif

  /* Check whether it is really time for a new domain decomposition */
  if(All.NumForcesSinceLastDomainDecomp > All.TotNumPart * All.TreeDomainUpdateFrequency)
    {
      t0 = second();

#if defined(SFR) || defined(BLACK_HOLES) || defined(MYSWITCH)
      rearrange_particle_sequence();
#endif

#ifdef PERIODIC
      do_box_wrapping();	/* map the particles back onto the box */
#endif
      All.NumForcesSinceLastDomainDecomp = 0;
      TreeReconstructFlag = 1;	/* ensures that new tree will be constructed */

      if(ThisTask == 0)
	{
	  printf("domain decomposition... \n");
	  fflush(stdout);
	}

      Key = (peanokey *) mymalloc(sizeof(peanokey) * All.MaxPart);
      KeySorted = (peanokey *) mymalloc(sizeof(peanokey) * All.MaxPart);

      toGo = (int *) mymalloc(sizeof(int) * NTask * NTask);
      toGoSph = (int *) mymalloc(sizeof(int) * NTask * NTask);
      local_toGo = (int *) mymalloc(sizeof(int) * NTask);
      local_toGoSph = (int *) mymalloc(sizeof(int) * NTask);
      list_NumPart = (int *) mymalloc(sizeof(int) * NTask);
      list_N_gas = (int *) mymalloc(sizeof(int) * NTask);
      list_load = (int *) mymalloc(sizeof(int) * NTask);
      list_loadsph = (int *) mymalloc(sizeof(int) * NTask);
      list_work = (double *) mymalloc(sizeof(double) * NTask);
      list_speedfac = (double *) mymalloc(sizeof(double) * NTask);
      list_cadj_cpu = (double *) mymalloc(sizeof(double) * NTask);
      list_cadj_cost = (double *) mymalloc(sizeof(double) * NTask);

#ifdef LT_STELLAREVOLUTION
      local_toGoStar = mymalloc(sizeof(int) * NTask * NTask);
      toGoStar = mymalloc(sizeof(int) * NTask * NTask);
      list_N_star = mymalloc(sizeof(int) * NTask);
      list_loadstar = mymalloc(sizeof(int) * NTask);
#endif

      MPI_Allgather(&NumPart, 1, MPI_INT, list_NumPart, 1, MPI_INT, MPI_COMM_WORLD);
      MPI_Allgather(&N_gas, 1, MPI_INT, list_N_gas, 1, MPI_INT, MPI_COMM_WORLD);
#ifdef LT_STELLAREVOLUTION
      MPI_Allgather(&N_star, 1, MPI_INT, list_N_star, 1, MPI_INT, MPI_COMM_WORLD);
#endif

      maxload = (int) (All.MaxPart * REDUC_FAC);
      maxloadsph = (int) (All.MaxPartSph * REDUC_FAC);
#ifdef LT_STELLAREVOLUTION
      maxloadstar = (int) (All.MaxPartMet * REDUC_FAC);
#endif

#ifdef LT_SEvDbg
      get_metals_sumcheck(3);
#endif
      domain_decompose();
#ifdef LT_SEvDbg
      get_metals_sumcheck(5);
#endif

#ifdef LT_STELLAREVOLUTION
      myfree(list_loadstar);
      myfree(list_N_star);
#endif
      myfree(list_cadj_cost);
      myfree(list_cadj_cpu);
      myfree(list_speedfac);
      myfree(list_work);
      myfree(list_loadsph);
      myfree(list_load);
      myfree(list_N_gas);
      myfree(list_NumPart);
      myfree(local_toGoSph);
      myfree(local_toGo);
      myfree(toGoSph);
      myfree(toGo);


      if(ThisTask == 0)
	{
	  printf("domain decomposition done. \n");
	  fflush(stdout);
	}

      t1 = second();
      All.CPU_Domain += timediff(t0, t1);
      CPU_Step[CPU_DOMAIN] += timediff(t0, t1);

#ifdef PEANOHILBERT
      t0 = second();
      peano_hilbert_order();
      t1 = second();
      All.CPU_Peano += timediff(t0, t1);
      CPU_Step[CPU_PEANO] += timediff(t0, t1);
      All.Cadj_Cpu += timediff(t0, t1);

      MPI_Barrier(MPI_COMM_WORLD); /* prevent imbalance in peano to be double-counted */
#endif

      myfree(KeySorted);
      myfree(Key);
    }

}

double domain_particle_costfactor(int i)
{
  if(P[i].Ti_endstep > P[i].Ti_begstep)
    return (1.0 + P[i].GravCost) / (P[i].Ti_endstep - P[i].Ti_begstep);
  else
    return (1.0 + P[i].GravCost) / TIMEBASE;
}

/*!  This function does the domain decomposition for all
 *   particle types
 */
void domain_decompose(void)
{
  int i, j, status;
  int ngrp, task, partner, sendcount, recvcount;
  long long sumtogo, sumload;
  int maxload, *temp;
  double sumwork, sumcpu, sumcost, maxwork, cadj_SpeedFac;

#ifdef CPUSPEEDADJUSTMENT
  double min_load, sum_speedfac;
#endif
#ifdef LT_STELLAREVOLUTION
  int tempcountstar;
#endif

  for(i = 0; i < 6; i++)
    NtypeLocal[i] = 0;

  for(i = 0, gravcost = 0; i < NumPart; i++)
    {
      NtypeLocal[P[i].Type]++;
      gravcost += domain_particle_costfactor(i);
      All.Cadj_Cost += domain_particle_costfactor(i);
    }

  MPI_Allreduce(&gravcost, &totgravcost, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  /* because Ntype[] is of type `long long', we cannot do a simple
   * MPI_Allreduce() to sum the total particle numbers 
   */
  temp = (int *) mymalloc(NTask * 6 * sizeof(int));
  MPI_Allgather(NtypeLocal, 6, MPI_INT, temp, 6, MPI_INT, MPI_COMM_WORLD);
  for(i = 0; i < 6; i++)
    {
      Ntype[i] = 0;
      for(j = 0; j < NTask; j++)
	Ntype[i] += temp[j * 6 + i];
    }
  myfree(temp);

#ifndef UNEQUALSOFTENINGS
  for(i = 0; i < 6; i++)
    if(Ntype[i] > 0)
      break;

  for(ngrp = i + 1; ngrp < 6; ngrp++)
    {
      if(Ntype[ngrp] > 0)
	if(All.SofteningTable[ngrp] != All.SofteningTable[i])
	  {
	    if(ThisTask == 0)
	      {
		fprintf(stdout, "Code was not compiled with UNEQUALSOFTENINGS, but some of the\n");
		fprintf(stdout, "softening lengths are unequal nevertheless.\n");
		fprintf(stdout, "This is not allowed.\n");
	      }
	    endrun(0);
	  }
    }
#endif


  /* determine global dimensions of domain grid */
  domain_findExtent();

  domain_determineTopTree();

  All.Cadj_Cpu *= 0.9;
  All.Cadj_Cost *= 0.9;

  MPI_Allgather(&All.Cadj_Cpu, 1, MPI_DOUBLE, list_cadj_cpu, 1, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgather(&All.Cadj_Cost, 1, MPI_DOUBLE, list_cadj_cost, 1, MPI_DOUBLE, MPI_COMM_WORLD);

#ifdef CPUSPEEDADJUSTMENT
  MPI_Allreduce(&All.Cadj_Cost, &min_load, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  if(min_load > 0)
    {
      cadj_SpeedFac = All.Cadj_Cpu / All.Cadj_Cost;

      MPI_Allreduce(&cadj_SpeedFac, &sum_speedfac, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      cadj_SpeedFac /= (sum_speedfac / NTask);
    }
  else
    cadj_SpeedFac = 1;

  MPI_Allgather(&cadj_SpeedFac, 1, MPI_DOUBLE, list_speedfac, 1, MPI_DOUBLE, MPI_COMM_WORLD);
#else
  cadj_SpeedFac = 1;
  MPI_Allgather(&cadj_SpeedFac, 1, MPI_DOUBLE, list_speedfac, 1, MPI_DOUBLE, MPI_COMM_WORLD);
#endif


  /* determine cost distribution in domain grid */
  domain_sumCost();


  /* find the split of the domain grid recursively */
  status = domain_findSplit_balanced(0, NTask, 0, NTopleaves - 1);
  if(status != 0)
    {
      if(ThisTask == 0)
	printf
	  ("\nNote: the domain decomposition is suboptimum because the ceiling for memory-imbalance is reached\n");

      status = domain_findSplit(0, NTask, 0, NTopleaves - 1);
      if(status != 0)
	{
	  if(ThisTask == 0)
	    printf("\nNo domain decomposition that stays within memory bounds is possible.\n");
	  endrun(0);
	}
      else
	domain_shiftSplit();
    }
  else
    domain_shiftSplit();

  DomainMyStart = DomainStartList[ThisTask];
  DomainMyLast = DomainEndList[ThisTask];

  if(ThisTask == 0)
    {
      sumload = maxload = 0;
      sumwork = sumcpu = sumcost = maxwork = 0;
      for(i = 0; i < NTask; i++)
	{
	  sumload += list_load[i];
	  sumwork += list_speedfac[i] * list_work[i];
	  sumcpu += list_cadj_cpu[i];
          sumcost += list_cadj_cost[i];

	  if(list_load[i] > maxload)
	    maxload = list_load[i];

	  if(list_speedfac[i] * list_work[i] > maxwork)
	    maxwork = list_speedfac[i] * list_work[i];
	}

      printf("work-load balance=%g   memory-balance=%g\n",
	     maxwork / (sumwork / NTask), maxload / (((double) sumload) / NTask));

      printf("\nSpeedfac:\n");
      for(i = 0; i < NTask; i++)
	printf("Speedfac [%3d]  speedfac=%8.4f  work=%8.4f   load=%8.4f   cpu=%8.4f   cost=%8.4f \n", i, list_speedfac[i],
	       list_speedfac[i] * list_work[i] / (sumwork / NTask),
	       list_load[i] / (((double) sumload) / NTask),
	       list_cadj_cpu[i] / (sumcpu / NTask),
	       list_cadj_cost[i] / (sumcost / NTask));
      printf("\n");
    }


  /* determine for each cpu how many particles have to be shifted to other cpus */
  domain_countToGo();

  for(i = 0, sumtogo = 0; i < NTask * NTask; i++)
    sumtogo += toGo[i];

  while(sumtogo > 0)
    {
      if(ThisTask == 0)
	{
	  printf("exchange of %d%09d particles\n", (int) (sumtogo / 1000000000),
		 (int) (sumtogo % 1000000000));
	  fflush(stdout);
	}

      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  for(task = 0; task < NTask; task++)
	    {
	      partner = task ^ ngrp;

	      if(partner < NTask && task < partner)
		{
		  /* treat SPH separately */
		  if(All.TotN_gas > 0)
		    {
		      domain_findExchangeNumbers(task, partner, 1, &sendcount, &recvcount);

		      list_NumPart[task] += recvcount - sendcount;
		      list_NumPart[partner] -= recvcount - sendcount;
		      list_N_gas[task] += recvcount - sendcount;
		      list_N_gas[partner] -= recvcount - sendcount;

		      toGo[task * NTask + partner] -= sendcount;
		      toGo[partner * NTask + task] -= recvcount;
		      toGoSph[task * NTask + partner] -= sendcount;
		      toGoSph[partner * NTask + task] -= recvcount;

		      if(task == ThisTask)	/* actually carry out the exchange */
			domain_exchangeParticles(partner, 1, sendcount, recvcount);
		      if(partner == ThisTask)
			domain_exchangeParticles(task, 1, recvcount, sendcount);
		    }

		  domain_findExchangeNumbers(task, partner, 0, &sendcount, &recvcount);

		  list_NumPart[task] += recvcount - sendcount;
		  list_NumPart[partner] -= recvcount - sendcount;

		  toGo[task * NTask + partner] -= sendcount;
		  toGo[partner * NTask + task] -= recvcount;
#ifdef LT_STELLAREVOLUTION
		  toGoStar[task * NTask + partner] -= send_countstar;
		  toGoStar[partner * NTask + task] -= recv_countstar;
		  list_N_star[task] += recv_countstar - send_countstar;
		  list_N_star[partner] -= recv_countstar - send_countstar;
#endif

		  if(task == ThisTask)	/* actually carry out the exchange */
		    domain_exchangeParticles(partner, 0, sendcount, recvcount);
		  if(partner == ThisTask)
#ifdef LT_STELLAREVOLUTION
		    {
		      tempcountstar = recv_countstar;
		      recv_countstar = send_countstar;
		      send_countstar = tempcountstar;
#endif
		      domain_exchangeParticles(task, 0, recvcount, sendcount);
#ifdef LT_STELLAREVOLUTION
		    }
		  recv_countstar = send_countstar = 0;
#endif
		}
	    }
	}

      for(i = 0, sumtogo = 0; i < NTask * NTask; i++)
	sumtogo += toGo[i];
    }

#ifdef LT_STELLAREVOLUTION
  for(i = 0; i < NumPart; i++)
    if(P[i].Type == 4)
      {
	if((i != MetP[P[i].MetID].PID) || (P[i].MetID >= N_star))
	  {
	    printf("************\n  Warning [a] on Task %i at particle %i\n************\n"
		   "P.MetID = %u  MetP.ID = %u N_star = %i\n",
		   ThisTask, i, P[i].MetID, MetP[P[i].MetID].PID, N_star);
	    fflush(stdout);
	    endrun(919190);
	  }
	if(MetP[P[i].MetID].iMass == 0)
	  {
	    printf("************\n  Warning [b] on Task %i at particle %i\n************\n", ThisTask, i);
	    fflush(stdout);
	    endrun(919191);
	  }
      }
#endif
}

/*! This function tries to find a split point in a range of cells in the
 *  domain-grid.  The range of cells starts at 'first', and ends at 'last'
 *  (inclusively). The number of cpus that holds the range is 'ncpu', with the
 *  first cpu given by 'cpustart'. If more than 2 cpus are to be split, the
 *  function calls itself recursively. The division tries to achieve a best
 *  particle-load balance under the constraint that 'maxload' and 'maxloadsph'
 *  may not be exceeded, and that each cpu holds at least one cell from the
 *  domaingrid. If such a decomposition cannot be achieved, a non-zero error
 *  code is returned.
 *
 *  After successful completion, DomainMyStart[] and DomainMyLast[] contain
 *  the first and last cell of the domaingrid assigned to the local task for
 *  the given type. Also, DomainTask[] contains for each cell the task it was
 *  assigned to.
 */
int domain_findSplit(int cpustart, int ncpu, int first, int last)
{
  int i, split, ok_left, ok_right;
  long long load, sphload, load_leftOfSplit, sphload_leftOfSplit;
  int ncpu_leftOfSplit;
  double maxAvgLoad_CurrentSplit, maxAvgLoad_NewSplit, dmax1, dmax2;

#ifdef LT_STELLAREVOLUTION
  long long starload, starload_leftOfSplit;
#endif

  ncpu_leftOfSplit = ncpu / 2;

#ifdef LT_STELLAREVOLUTION
  starload = 0;
#endif
  for(i = first, load = 0, sphload = 0; i <= last; i++)
    {
      load += DomainCount[i];
      sphload += DomainCountSph[i];
#ifdef LT_STELLAREVOLUTION
      starload += DomainCountStar[i];
#endif
    }

  split = first + ncpu_leftOfSplit;

#ifdef LT_STELLAREVOLUTION
  starload_leftOfSplit = 0;
#endif
  for(i = first, load_leftOfSplit = sphload_leftOfSplit = 0; i < split; i++)
    {
      load_leftOfSplit += DomainCount[i];
      sphload_leftOfSplit += DomainCountSph[i];
#ifdef LT_STELLAREVOLUTION
      starload_leftOfSplit += DomainCountStar[i];
#endif
    }

  /* find the best split point in terms of work-load balance */

  while(split < last - (ncpu - ncpu_leftOfSplit - 1) && split > 0)
    {
      maxAvgLoad_CurrentSplit =
	DMAX(load_leftOfSplit / ncpu_leftOfSplit, (load - load_leftOfSplit) / (ncpu - ncpu_leftOfSplit));

      maxAvgLoad_NewSplit =
	DMAX((load_leftOfSplit + DomainCount[split]) / ncpu_leftOfSplit,
	     (load - load_leftOfSplit - DomainCount[split]) / (ncpu - ncpu_leftOfSplit));

      if(maxAvgLoad_NewSplit <= maxAvgLoad_CurrentSplit)
	{
	  load_leftOfSplit += DomainCount[split];
	  sphload_leftOfSplit += DomainCountSph[split];
#ifdef LT_STELLAREVOLUTION
	  starload_leftOfSplit += DomainCountStar[split];
#endif
	  split++;
	}
      else
	break;
    }


  /* we will now have to check whether this solution is possible given the restrictions on the maximum load */

#ifdef LT_STELLAREVOLUTION
  starload_leftOfSplit = 0;
#endif
  for(i = first, load_leftOfSplit = 0, sphload_leftOfSplit = 0; i < split; i++)
    {
      load_leftOfSplit += DomainCount[i];
      sphload_leftOfSplit += DomainCountSph[i];
#ifdef LT_STELLAREVOLUTION
      starload_leftOfSplit += DomainCountStar[i];
#endif
    }

  if(load_leftOfSplit > maxload * ncpu_leftOfSplit ||
     (load - load_leftOfSplit) > maxload * (ncpu - ncpu_leftOfSplit))
    {
      /* we did not find a viable split */
      return -1;
    }

  if(sphload_leftOfSplit > maxloadsph * ncpu_leftOfSplit ||
     (sphload - sphload_leftOfSplit) > maxloadsph * (ncpu - ncpu_leftOfSplit))
    {
      /* we did not find a viable split */
      return -1;
    }

#ifdef LT_STELLAREVOLUTION
  if(((double) starload_leftOfSplit) / ncpu_leftOfSplit > maxloadstar ||
     ((double) (starload - starload_leftOfSplit)) / (ncpu - ncpu_leftOfSplit) > maxloadstar)
    {
      /* we did not find a viable split */
      if(ThisTask == 0)
	printf("*********************\n"
	       "i: %i split: %i\n"
	       "maxloadstar: %i rightload: %g leftload: %g",
	       i, split, maxloadstar,
	       ((double) starload_leftOfSplit) / ncpu_leftOfSplit,
	       ((double) (starload - starload_leftOfSplit)) / (ncpu - ncpu_leftOfSplit));

      return -1;
    }
#endif


  if(ncpu_leftOfSplit >= 2)
    ok_left = domain_findSplit(cpustart, ncpu_leftOfSplit, first, split - 1);
  else
    ok_left = 0;

  if((ncpu - ncpu_leftOfSplit) >= 2)
    ok_right = domain_findSplit(cpustart + ncpu_leftOfSplit, ncpu - ncpu_leftOfSplit, split, last);
  else
    ok_right = 0;

  if(ok_left == 0 && ok_right == 0)
    {
      /* found a viable split */

      if(ncpu_leftOfSplit == 1)
	{
	  for(i = first; i < split; i++)
	    DomainTask[i] = cpustart;

	  list_load[cpustart] = load_leftOfSplit;
	  list_loadsph[cpustart] = sphload_leftOfSplit;
	  DomainStartList[cpustart] = first;
	  DomainEndList[cpustart] = split - 1;
#ifdef LT_STELLAREVOLUTION
	  list_loadstar[cpustart] = starload_leftOfSplit;
#endif
	}

      if((ncpu - ncpu_leftOfSplit) == 1)
	{
	  for(i = split; i <= last; i++)
	    DomainTask[i] = cpustart + ncpu_leftOfSplit;

	  list_load[cpustart + ncpu_leftOfSplit] = load - load_leftOfSplit;
	  list_loadsph[cpustart + ncpu_leftOfSplit] = sphload - sphload_leftOfSplit;
	  DomainStartList[cpustart + ncpu_leftOfSplit] = split;
	  DomainEndList[cpustart + ncpu_leftOfSplit] = last;
#ifdef LT_STELLAREVOLUTION
	  list_loadstar[cpustart + ncpu_leftOfSplit] = starload - starload_leftOfSplit;
#endif
	}

      return 0;
    }

  /* we did not find a viable split */
  return -1;
}

int domain_findSplit_balanced(int cpustart, int ncpu, int first, int last)
{
  int i, split, ok_left, ok_right;
  long long load, sphload, load_leftOfSplit, sphload_leftOfSplit;
  int ncpu_leftOfSplit;
  double maxAvgLoad_CurrentSplit, maxAvgLoad_NewSplit, work, work_leftOfSplit;
  double speedfacleft, speedfacright, dmax1, dmax2;

#ifdef LT_STELLAREVOLUTION
  long long starload = 0, starload_leftOfSplit = 0;
#endif

  ncpu_leftOfSplit = ncpu / 2;

  for(i = first, load = 0, sphload = 0, work = 0; i <= last; i++)
    {
      load += DomainCount[i];
      sphload += DomainCountSph[i];
#ifdef LT_STELLAREVOLUTION
      starload += DomainCountStar[i];
#endif
      work += DomainWork[i];
    }

  split = first + ncpu_leftOfSplit;

  for(i = first, load_leftOfSplit = sphload_leftOfSplit = 0, work_leftOfSplit = 0; i < split; i++)
    {
      load_leftOfSplit += DomainCount[i];
      sphload_leftOfSplit += DomainCountSph[i];
#ifdef LT_STELLAREVOLUTION
      starload_leftOfSplit += DomainCountStar[i];
#endif
      work_leftOfSplit += DomainWork[i];
    }

  for(i = cpustart, speedfacleft = 0; i < cpustart + ncpu_leftOfSplit; i++)
    speedfacleft += 1 / list_speedfac[i];

  for(i = cpustart + ncpu_leftOfSplit, speedfacright = 0; i < cpustart + ncpu; i++)
    speedfacright += 1 / list_speedfac[i];

  /* find the best split point in terms of work-load balance */

  while(split < last - (ncpu - ncpu_leftOfSplit - 1) && split > 0)
    {
      maxAvgLoad_CurrentSplit =
	DMAX(work_leftOfSplit / speedfacleft, (work - work_leftOfSplit) / speedfacright);

      maxAvgLoad_NewSplit =
	DMAX((work_leftOfSplit + DomainWork[split]) / speedfacleft,
	     (work - work_leftOfSplit - DomainWork[split]) / speedfacright);

      if(maxAvgLoad_NewSplit <= maxAvgLoad_CurrentSplit)
	{
	  load_leftOfSplit += DomainCount[split];
	  sphload_leftOfSplit += DomainCountSph[split];
#ifdef LT_STELLAREVOLUTION
	  starload_leftOfSplit += DomainCountStar[split];
#endif
	  work_leftOfSplit += DomainWork[split];
	  split++;
	}
      else
	break;
    }


  /* we will now have to check whether this solution is possible given the restrictions on the maximum load */

#ifdef LT_STELLAREVOLUTION
  starload_leftOfSplit = 0;
#endif

  for(i = first, load_leftOfSplit = 0, sphload_leftOfSplit = 0; i < split; i++)
    {
      load_leftOfSplit += DomainCount[i];
      sphload_leftOfSplit += DomainCountSph[i];
#ifdef LT_STELLAREVOLUTION
      starload_leftOfSplit += DomainCountStar[i];
#endif
    }

  if(load_leftOfSplit > maxload * ncpu_leftOfSplit ||
     (load - load_leftOfSplit) > maxload * (ncpu - ncpu_leftOfSplit))
    {
      /* we did not find a viable split */
      return -1;
    }

  if(sphload_leftOfSplit > maxloadsph * ncpu_leftOfSplit ||
     (sphload - sphload_leftOfSplit) > maxloadsph * (ncpu - ncpu_leftOfSplit))
    {
      /* we did not find a viable split */
      return -1;
    }

#ifdef LT_STELLAREVOLUTION
  if(((double) starload_leftOfSplit) / ncpu_leftOfSplit > maxloadstar ||
     ((double) (starload - starload_leftOfSplit)) / (ncpu - ncpu_leftOfSplit) > maxloadstar)
    {
      /* we did not find a viable split */
      if(ThisTask == 0)
	printf("*********************\n"
	       "i: %i split: %i\n"
	       "maxloadstar: %i rightload: %g leftload: %g",
	       i, split, maxloadstar,
	       ((double) starload_leftOfSplit) / ncpu_leftOfSplit,
	       ((double) starload - starload_leftOfSplit) / (ncpu - ncpu_leftOfSplit));

      return -1;
    }
#endif

  if(ncpu_leftOfSplit >= 2)
    ok_left = domain_findSplit_balanced(cpustart, ncpu_leftOfSplit, first, split - 1);
  else
    ok_left = 0;

  if((ncpu - ncpu_leftOfSplit) >= 2)
    ok_right = domain_findSplit_balanced(cpustart + ncpu_leftOfSplit, ncpu - ncpu_leftOfSplit, split, last);
  else
    ok_right = 0;

  if(ok_left == 0 && ok_right == 0)
    {
      /* found a viable split */

      if(ncpu_leftOfSplit == 1)
	{
	  for(i = first; i < split; i++)
	    DomainTask[i] = cpustart;

	  list_load[cpustart] = load_leftOfSplit;
	  list_loadsph[cpustart] = sphload_leftOfSplit;
	  list_work[cpustart] = work_leftOfSplit;
	  DomainStartList[cpustart] = first;
	  DomainEndList[cpustart] = split - 1;
#ifdef LT_STELLAREVOLUTION
	  list_loadstar[cpustart] = starload_leftOfSplit;
#endif
	}

      if((ncpu - ncpu_leftOfSplit) == 1)
	{
	  for(i = split; i <= last; i++)
	    DomainTask[i] = cpustart + ncpu_leftOfSplit;

	  list_load[cpustart + ncpu_leftOfSplit] = load - load_leftOfSplit;
	  list_loadsph[cpustart + ncpu_leftOfSplit] = sphload - sphload_leftOfSplit;
	  list_work[cpustart + ncpu_leftOfSplit] = work - work_leftOfSplit;
#ifdef LT_STELLAREVOLUTION
	  list_loadstar[cpustart + ncpu_leftOfSplit] = starload - starload_leftOfSplit;
#endif
	  DomainStartList[cpustart + ncpu_leftOfSplit] = split;
	  DomainEndList[cpustart + ncpu_leftOfSplit] = last;
	}

      return 0;
    }

  /* we did not find a viable split */
  return -1;
}


/*! This function tries to improve the domain decomposition found by
 *  domain_findSplit() with respect to work-load balance.  To this end, the
 *  boundaries in the existing domain-split solution (which was found by
 *  trying to balance the particle load) are shifted as long as this leads
 *  to better work-load while still remaining within the allowed
 *  memory-imbalance constraints.
 */
void domain_shiftSplit(void)
{
  int i, task, iter = 0, moved;
  double maxw, newmaxw, dmax1, dmax2;

  for(task = 0; task < NTask; task++)
    list_work[task] = 0;

  for(i = 0; i < NTopleaves; i++)
    list_work[DomainTask[i]] += DomainWork[i];

  do
    {
      for(task = 0, moved = 0; task < NTask - 1; task++)
	{
	  maxw = DMAX(list_speedfac[task] * list_work[task], list_speedfac[task + 1] * list_work[task + 1]);

	  if(list_work[task] < list_work[task + 1])
	    {
	      newmaxw =
		DMAX(list_speedfac[task] * (list_work[task] + DomainWork[DomainStartList[task + 1]]),
		     list_speedfac[task + 1] * (list_work[task + 1] - DomainWork[DomainStartList[task + 1]]));
	      if(newmaxw <= maxw)
		{
		  if(DomainEndList[task + 1] > DomainStartList[task + 1])
		    {
		      if(list_load[task] + DomainCount[DomainStartList[task + 1]] <= maxload)
			{

			  if(list_loadsph[task] + DomainCountSph[DomainStartList[task + 1]] > maxloadsph)
			    continue;
#ifdef LT_STELLAREVOLUTION
			  if(list_loadstar[task] + DomainCountStar[DomainStartList[task + 1]] > maxloadstar)
			    continue;
#endif
			  /* ok, we can move one domain cell from right to left */
			  list_work[task] += DomainWork[DomainStartList[task + 1]];
			  list_load[task] += DomainCount[DomainStartList[task + 1]];
			  list_loadsph[task] += DomainCountSph[DomainStartList[task + 1]];
			  list_work[task + 1] -= DomainWork[DomainStartList[task + 1]];
			  list_load[task + 1] -= DomainCount[DomainStartList[task + 1]];
			  list_loadsph[task + 1] -= DomainCountSph[DomainStartList[task + 1]];
#ifdef LT_STELLAREVOLUTION
			  list_loadstar[task] += DomainCountStar[DomainStartList[task + 1]];
			  list_loadstar[task + 1] -= DomainCountStar[DomainStartList[task + 1]];
#endif

			  DomainTask[DomainStartList[task + 1]] = task;
			  DomainStartList[task + 1] += 1;
			  DomainEndList[task] += 1;

			  moved++;
			}
		    }
		}
	    }
	  else
	    {
	      newmaxw = DMAX(list_speedfac[task] * (list_work[task] - DomainWork[DomainEndList[task]]),
			     list_speedfac[task + 1] * (list_work[task + 1] +
							DomainWork[DomainEndList[task]]));
	      if(newmaxw <= maxw)
		{
		  if(DomainEndList[task] > DomainStartList[task])
		    {
		      if(list_load[task + 1] + DomainCount[DomainEndList[task]] <= maxload)
			{
			  if(list_loadsph[task + 1] + DomainCountSph[DomainEndList[task]] > maxloadsph)
			    continue;
#ifdef LT_STELLAREVOLUTION
			  if(list_loadstar[task + 1] + DomainCountStar[DomainEndList[task]] > maxloadstar)
			    continue;
#endif
			  /* ok, we can move one domain cell from left to right */
			  list_work[task] -= DomainWork[DomainEndList[task]];
			  list_load[task] -= DomainCount[DomainEndList[task]];
			  list_loadsph[task] -= DomainCountSph[DomainEndList[task]];
			  list_work[task + 1] += DomainWork[DomainEndList[task]];
			  list_load[task + 1] += DomainCount[DomainEndList[task]];
			  list_loadsph[task + 1] += DomainCountSph[DomainEndList[task]];
#ifdef LT_STELLAREVOLUTION
			  list_loadstar[task] -= DomainCountStar[DomainEndList[task]];
			  list_loadstar[task + 1] += DomainCountStar[DomainEndList[task]];
#endif

			  DomainTask[DomainEndList[task]] = task + 1;
			  DomainEndList[task] -= 1;
			  DomainStartList[task + 1] -= 1;

			  moved++;
			}
		    }
		}

	    }
	}

      iter++;
    }
  while(moved > 0 && iter < 10 * NTopleaves);
}


/*! This function counts how many particles have to be exchanged between two
 *  CPUs according to the domain split. If the CPUs are already quite full and
 *  hold data from other CPUs as well, not all the particles may be exchanged
 *  at once. In this case the communication phase has to be repeated, until
 *  enough of the third-party particles have moved off such that the
 *  decomposition can be completed.
 */
void domain_findExchangeNumbers(int task, int partner, int sphflag, int *send, int *recv)
{
  int numpartA, numpartsphA, ntobesentA, maxsendA, maxsendA_old;
  int numpartB, numpartsphB, ntobesentB, maxsendB, maxsendB_old;
  int imin1, imin2;

#ifdef LT_STELLAREVOLUTION
  int maxsendsA, maxsendsB, maxsendsA_old, maxsendsB_old;
  int numpartstarA, numpartstarB;
  double fracA, fracB;
  int delta_sendstar;
#endif

  numpartA = list_NumPart[task];
  numpartsphA = list_N_gas[task];

  numpartB = list_NumPart[partner];
  numpartsphB = list_N_gas[partner];

#ifdef LT_STELLAREVOLUTION
  numpartstarA = list_N_star[task];
  numpartstarB = list_N_star[partner];
#endif

  if(sphflag == 1)
    {
      ntobesentA = toGoSph[task * NTask + partner];
      ntobesentB = toGoSph[partner * NTask + task];
    }
  else
    {
      ntobesentA = toGo[task * NTask + partner] - toGoSph[task * NTask + partner];
      ntobesentB = toGo[partner * NTask + task] - toGoSph[partner * NTask + task];
    }

  maxsendA = IMIN(ntobesentA, All.BunchSizeDomain);
  maxsendB = IMIN(ntobesentB, All.BunchSizeDomain);

  do
    {
      maxsendA_old = maxsendA;
      maxsendB_old = maxsendB;

      maxsendA = IMIN(All.MaxPart - numpartB + maxsendB, maxsendA);
      maxsendB = IMIN(All.MaxPart - numpartA + maxsendA, maxsendB);
    }
  while((maxsendA != maxsendA_old) || (maxsendB != maxsendB_old));


#ifdef LT_STELLAREVOLUTION
  send_countstar = recv_countstar = 0;
  if(sphflag == 0)
    {
      maxsendsA = maxsendsB = 0;

      if(maxsendA > 0)
	{
	  fracA = (double) maxsendA / ntobesentA;
	  maxsendsA = toGoStar[task * NTask + partner] * fracA;
	  if(maxsendsA == 0 && toGoStar[task * NTask + partner] > 0)
	    maxsendsA++;
	}

      if(maxsendB > 0)
	{
	  fracB = (double) maxsendB / ntobesentB;
	  maxsendsB = toGoStar[partner * NTask + task] * fracB;
	  if(maxsendsB == 0 && toGoStar[partner * NTask + task] > 0)
	    maxsendsB++;
	}

      do
	{
	  maxsendsA_old = maxsendsA;
	  maxsendsB_old = maxsendsB;

	  maxsendsA = IMIN(All.MaxPartMet - numpartstarB + maxsendsB, maxsendsA);
	  maxsendsB = IMIN(All.MaxPartMet - numpartstarA + maxsendsA, maxsendsB);
	}
      while((maxsendsA != maxsendsA_old) || (maxsendsB != maxsendsB_old));

      send_countstar = maxsendsA;
      recv_countstar = maxsendsB;

      if((delta_sendstar = toGoStar[task * NTask + partner] - maxsendsA) > 0)
	{
	  if(delta_sendstar > (ntobesentA - maxsendA))
	    maxsendA = ntobesentA - delta_sendstar;
	}
      if((delta_sendstar = toGoStar[partner * NTask + task] - maxsendsB) > 0)
	{
	  if(delta_sendstar > (ntobesentB - maxsendB))
	    maxsendB = ntobesentB - delta_sendstar;
	}
    }
#endif

  /* now make also sure that there is enough space for SPH particeles */
  if(sphflag == 1)
    {
      do
	{
	  maxsendA_old = maxsendA;
	  maxsendB_old = maxsendB;

	  maxsendA = IMIN(All.MaxPartSph - numpartsphB + maxsendB, maxsendA);
	  maxsendB = IMIN(All.MaxPartSph - numpartsphA + maxsendA, maxsendB);
	}
      while((maxsendA != maxsendA_old) || (maxsendB != maxsendB_old));
    }

  *send = maxsendA;
  *recv = maxsendB;
}




/*! This function exchanges particles between two CPUs according to the domain
 *  split. In doing this, the memory boundaries which may restrict the exhange
 *  process are observed.
 */
void domain_exchangeParticles(int partner, int sphflag, int send_count, int recv_count)
{
  int i, no, n, count, rep;
  MPI_Status status;

#ifdef LT_STELLAREVOLUTION
  int countstar = 0;
  int j;
#endif

  for(n = 0, count = 0; count < send_count && n < NumPart; n++)
    {
      if(sphflag)
	{
	  if(P[n].Type != 0)
	    continue;
	}
      else
	{
	  if(P[n].Type == 0)
	    continue;
	}

      no = 0;

      while(TopNodes[no].Daughter >= 0)
	no = TopNodes[no].Daughter + (Key[n] - TopNodes[no].StartKey) / (TopNodes[no].Size / 8);

      no = TopNodes[no].Leaf;

      if(DomainTask[no] == partner)
	{
	  if(sphflag)		/* special reorder routine for SPH particles (need to stay at beginning) */
	    {
	      DomainPartBuf[count] = P[n];	/* copy particle and collect in contiguous memory */
	      DomainKeyBuf[count] = Key[n];
	      DomainSphBuf[count] = SphP[n];

	      P[n] = P[N_gas - 1];
	      P[N_gas - 1] = P[NumPart - 1];
#ifdef LT_STELLAREVOLUTION
	      if(P[NumPart - 1].Type == 4)
		MetP[P[NumPart - 1].MetID].PID = N_gas - 1;
#endif

	      Key[n] = Key[N_gas - 1];
	      Key[N_gas - 1] = Key[NumPart - 1];

	      SphP[n] = SphP[N_gas - 1];

	      N_gas--;
	    }
	  else
	    {
#ifdef LT_STELLAREVOLUTION
	      if(P[n].Type == 4)
		{
		  if(countstar == send_countstar)
		    /* maximum number of "sendable" stars already reached */
		    continue;
		  DomainMetBuf[countstar] = MetP[P[n].MetID];
		  DomainMetBuf[countstar].PID = count;
		  DomainPartBuf[count].MetID = countstar;
		  countstar++;
		  MetP[P[n].MetID] = MetP[N_star - 1];
		  P[MetP[N_star - 1].PID].MetID = P[n].MetID;

		  N_star--;
		}
	      else if((send_count - count) == (send_countstar - countstar))
		/* give the precedence to exchange stars */
		continue;
#endif
	      DomainPartBuf[count] = P[n];	/* copy particle and collect in contiguous memory */
	      DomainKeyBuf[count] = Key[n];

	      P[n] = P[NumPart - 1];
	      Key[n] = Key[NumPart - 1];
#ifdef LT_STELLAREVOLUTION
	      if((n < NumPart - 1) && P[NumPart - 1].Type == 4)
		MetP[P[NumPart - 1].MetID].PID = n;
#endif
	    }

	  count++;
	  NumPart--;
	  n--;
	}
    }

#ifdef LT_STELLAREVOLUTION
  if(countstar != send_countstar)
    {
      printf("Houston, we got a problem from stars...\n");
      printf("ThisTask=%d partner=%d countstar=%d send_countstar=%d\n[count=%d send_count=%d n=%d]\n",
	     ThisTask, partner, countstar, send_countstar, count, send_count, n);
    }
#endif
  if(count != send_count)
    {
      printf("Houston, we got a problem...\n");
      printf("ThisTask=%d count=%d send_count=%d\n", ThisTask, count, send_count);
      endrun(88);
    }
#ifdef LT_STELLAREVOLUTION
  if(countstar != send_countstar)
    endrun(89);
#endif

  /* transmit */

  for(rep = 0; rep < 2; rep++)
    {
      if((rep == 0 && ThisTask < partner) || (rep == 1 && ThisTask > partner))
	{
	  if(send_count > 0)
	    {
	      MPI_Ssend(&DomainPartBuf[0], send_count * sizeof(struct particle_data), MPI_BYTE, partner,
			TAG_PDATA, MPI_COMM_WORLD);

	      MPI_Ssend(&DomainKeyBuf[0], send_count * sizeof(peanokey), MPI_BYTE, partner, TAG_KEY,
			MPI_COMM_WORLD);

	      if(sphflag)
		MPI_Ssend(&DomainSphBuf[0], send_count * sizeof(struct sph_particle_data), MPI_BYTE, partner,
			  TAG_SPHDATA, MPI_COMM_WORLD);
#ifdef LT_STELLAREVOLUTION
	      else if(send_countstar > 0)
		MPI_Ssend(&DomainMetBuf[0], send_countstar * sizeof(struct met_particle_data), MPI_BYTE,
			  partner, TAG_METDATA, MPI_COMM_WORLD);
#endif
	    }
	}

      if((rep == 1 && ThisTask < partner) || (rep == 0 && ThisTask > partner))
	{
	  if(recv_count > 0)
	    {
	      if(sphflag)
		{
		  if((NumPart - N_gas) > recv_count)
		    {
		      for(i = 0; i < recv_count; i++)
			{
			  P[NumPart + i] = P[N_gas + i];
			  Key[NumPart + i] = Key[N_gas + i];
#ifdef LT_STELLAREVOLUTION
			  if(P[N_gas + i].Type == 4)
			    MetP[P[N_gas + i].MetID].PID = NumPart + i;
#endif
			}
		    }
		  else
		    {
		      for(i = NumPart - 1; i >= N_gas; i--)
			{
			  P[i + recv_count] = P[i];
			  Key[i + recv_count] = Key[i];
#ifdef LT_STELLAREVOLUTION
			  if(P[i].Type == 4)
			    MetP[P[i].MetID].PID = i + recv_count;
#endif
			}
		    }

		  MPI_Recv(&P[N_gas], recv_count * sizeof(struct particle_data), MPI_BYTE, partner, TAG_PDATA,
			   MPI_COMM_WORLD, &status);
		  MPI_Recv(&Key[N_gas], recv_count * sizeof(peanokey), MPI_BYTE, partner, TAG_KEY,
			   MPI_COMM_WORLD, &status);
		  MPI_Recv(&SphP[N_gas], recv_count * sizeof(struct sph_particle_data), MPI_BYTE, partner,
			   TAG_SPHDATA, MPI_COMM_WORLD, &status);

		  N_gas += recv_count;
		}
	      else
		{
		  MPI_Recv(&P[NumPart], recv_count * sizeof(struct particle_data), MPI_BYTE, partner,
			   TAG_PDATA, MPI_COMM_WORLD, &status);
		  MPI_Recv(&Key[NumPart], recv_count * sizeof(peanokey), MPI_BYTE, partner,
			   TAG_KEY, MPI_COMM_WORLD, &status);
#ifdef LT_STELLAREVOLUTION
		  if(recv_countstar > 0)
		    {
		      MPI_Recv(&MetP[N_star], recv_countstar * sizeof(struct met_particle_data), MPI_BYTE,
			       partner, TAG_METDATA, MPI_COMM_WORLD, &status);

		      j = NumPart;
		      for(i = 0; i < recv_countstar; i++)
			{

			  if(j == NumPart + recv_count)
			    {
			      printf("Houston, who knows what haps?\n");
			      endrun(90);
			    }
			  for(; j < NumPart + recv_count; j++)
			    if(P[j].Type == 4)
			      break;

			  MetP[N_star + i].PID = j;
			  P[j].MetID = N_star + i;
			  j++;
			  /*
			     P[(MetP[N_star + i].PID += NumPart)].MetID += N_star;
			   */
			}

		      N_star += recv_countstar;
		    }
#endif
		}

	      NumPart += recv_count;
	    }
	}
    }
}

/*! This function determines how many particles that are currently stored on
 *  the local CPU have to be moved off according to the domain decomposition.
 */
void domain_countToGo(void)
{
  int n, no;

  for(n = 0; n < NTask; n++)
    {
      local_toGo[n] = 0;
      local_toGoSph[n] = 0;
#ifdef LT_STELLAREVOLUTION
      local_toGoStar[n] = 0;
#endif
    }

  for(n = 0; n < NumPart; n++)
    {
      no = 0;

      while(TopNodes[no].Daughter >= 0)
	no = TopNodes[no].Daughter + (Key[n] - TopNodes[no].StartKey) / (TopNodes[no].Size / 8);

      no = TopNodes[no].Leaf;

      if(DomainTask[no] != ThisTask)
	{
	  local_toGo[DomainTask[no]] += 1;
	  if(P[n].Type == 0)
	    local_toGoSph[DomainTask[no]] += 1;
#ifdef LT_STELLAREVOLUTION
	  if(P[n].Type == 4)
	    local_toGoStar[DomainTask[no]] += 1;
#endif
	}
    }

  MPI_Allgather(local_toGo, NTask, MPI_INT, toGo, NTask, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgather(local_toGoSph, NTask, MPI_INT, toGoSph, NTask, MPI_INT, MPI_COMM_WORLD);
#ifdef LT_STELLAREVOLUTION
  MPI_Allgather(local_toGoStar, NTask, MPI_INT, toGoStar, NTask, MPI_INT, MPI_COMM_WORLD);
#endif
}


void domain_walktoptree(int no)
{
  int i;

  if(TopNodes[no].Daughter == -1)
    {
      TopNodes[no].Leaf = NTopleaves;
      NTopleaves++;
    }
  else
    {
      for(i = 0; i < 8; i++)
	domain_walktoptree(TopNodes[no].Daughter + i);
    }
}

/*! This routine bins the particles onto the domain-grid, i.e. it sums up the
 *  total number of particles and the total amount of work in each of the
 *  domain-cells. This information forms the basis for the actual decision on
 *  the adopted domain decomposition.
 */
void domain_sumCost(void)
{
  int i, n, no;
  double *local_DomainWork;
  int *local_DomainCount;
  int *local_DomainCountSph;

#ifdef LT_STELLAREVOLUTION
  int *local_DomainCountStar;
#endif

  local_DomainWork = (double *) mymalloc(NTopnodes * sizeof(double));
  local_DomainCount = (int *) mymalloc(NTopnodes * sizeof(int));
  local_DomainCountSph = (int *) mymalloc(NTopnodes * sizeof(int));
#ifdef LT_STELLAREVOLUTION
  local_DomainCountStar = (int *) mymalloc(NTopnodes * sizeof(int));
#endif


  NTopleaves = 0;

  domain_walktoptree(0);

  for(i = 0; i < NTopleaves; i++)
    {
      local_DomainWork[i] = 0;
      local_DomainCount[i] = 0;
      local_DomainCountSph[i] = 0;
#ifdef LT_STELLAREVOLUTION
      local_DomainCountStar[i] = 0;
#endif
    }

  if(ThisTask == 0)
    printf("NTopleaves= %d\n", NTopleaves);

  for(n = 0; n < NumPart; n++)
    {
      no = 0;

      while(TopNodes[no].Daughter >= 0)
	no = TopNodes[no].Daughter + (Key[n] - TopNodes[no].StartKey) / (TopNodes[no].Size / 8);

      no = TopNodes[no].Leaf;

      local_DomainWork[no] += domain_particle_costfactor(n);
      local_DomainCount[no] += 1;
      if(P[n].Type == 0)
	local_DomainCountSph[no] += 1;

#ifdef LT_STELLAREVOLUTION
      if(P[n].Type == 4)
	{
	  local_DomainCountStar[no] += 1;
#ifdef LT_COST_SE
	  local_DomainWork[no] += get_cost_SE(n);
#endif
	}
#endif
    }

  MPI_Allreduce(local_DomainWork, DomainWork, NTopleaves, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(local_DomainCount, DomainCount, NTopleaves, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(local_DomainCountSph, DomainCountSph, NTopleaves, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#ifdef LT_STELLAREVOLUTION
  MPI_Allreduce(local_DomainCountStar, DomainCountStar, NTopleaves, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  myfree(local_DomainCountStar);
#endif
  myfree(local_DomainCountSph);
  myfree(local_DomainCount);
  myfree(local_DomainWork);
}


/*! This routine finds the extent of the global domain grid. 
 */
void domain_findExtent(void)
{
  int i, j;
  double len, xmin[3], xmax[3], xmin_glob[3], xmax_glob[3];

  /* determine local extension */
  for(j = 0; j < 3; j++)
    {
      xmin[j] = MAX_REAL_NUMBER;
      xmax[j] = -MAX_REAL_NUMBER;
    }

  for(i = 0; i < NumPart; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  if(xmin[j] > P[i].Pos[j])
	    xmin[j] = P[i].Pos[j];

	  if(xmax[j] < P[i].Pos[j])
	    xmax[j] = P[i].Pos[j];
	}
    }

  MPI_Allreduce(xmin, xmin_glob, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(xmax, xmax_glob, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  len = 0;
  for(j = 0; j < 3; j++)
    if(xmax_glob[j] - xmin_glob[j] > len)
      len = xmax_glob[j] - xmin_glob[j];

  len *= 1.001;

  for(j = 0; j < 3; j++)
    {
      DomainCenter[j] = 0.5 * (xmin_glob[j] + xmax_glob[j]);
      DomainCorner[j] = 0.5 * (xmin_glob[j] + xmax_glob[j]) - 0.5 * len;
    }

  DomainLen = len;
  DomainFac = 1.0 / len * (((peanokey) 1) << (BITS_PER_DIMENSION));
}




void domain_determineTopTree(void)
{
  int i, ntop_local, ntop;
  int *ntopnodelist, *ntopoffset;

  for(i = 0; i < NumPart; i++)
    {
      KeySorted[i] = Key[i] = peano_hilbert_key((int) ((P[i].Pos[0] - DomainCorner[0]) * DomainFac),
						(int) ((P[i].Pos[1] - DomainCorner[1]) * DomainFac),
						(int) ((P[i].Pos[2] - DomainCorner[2]) * DomainFac),
						BITS_PER_DIMENSION);
    }

  qsort(KeySorted, NumPart, sizeof(peanokey), domain_compare_key);

  NTopnodes = 1;
  TopNodes[0].Daughter = -1;
  TopNodes[0].Size = PEANOCELLS;
  TopNodes[0].StartKey = 0;
  TopNodes[0].Count = NumPart;
  TopNodes[0].GravCost = gravcost;	/* local value */
  TopNodes[0].Pstart = 0;

  domain_topsplit_local(0, 0);

  toplist_local = (struct topnode_exchange *) mymalloc(NTopnodes * sizeof(struct topnode_exchange));

  for(i = 0, ntop_local = 0; i < NTopnodes; i++)
    {
      if(TopNodes[i].Daughter == -1)	/* only use leaves */
	{
	  toplist_local[ntop_local].Startkey = TopNodes[i].StartKey;
	  toplist_local[ntop_local].Count = TopNodes[i].Count;
	  toplist_local[ntop_local].GravCost = TopNodes[i].GravCost;
	  ntop_local++;
	}
    }

  ntopnodelist = (int *) mymalloc(sizeof(int) * NTask);
  ntopoffset = (int *) mymalloc(sizeof(int) * NTask);

  MPI_Allgather(&ntop_local, 1, MPI_INT, ntopnodelist, 1, MPI_INT, MPI_COMM_WORLD);

  for(i = 0, ntop = 0, ntopoffset[0] = 0; i < NTask; i++)
    {
      ntop += ntopnodelist[i];
      if(i > 0)
	ntopoffset[i] = ntopoffset[i - 1] + ntopnodelist[i - 1];
    }


  toplist = (struct topnode_exchange *) mymalloc(ntop * sizeof(struct topnode_exchange));

  for(i = 0; i < NTask; i++)
    {
      ntopnodelist[i] *= sizeof(struct topnode_exchange);
      ntopoffset[i] *= sizeof(struct topnode_exchange);
    }

  MPI_Allgatherv(toplist_local, ntop_local * sizeof(struct topnode_exchange), MPI_BYTE,
		 toplist, ntopnodelist, ntopoffset, MPI_BYTE, MPI_COMM_WORLD);

  qsort(toplist, ntop, sizeof(struct topnode_exchange), domain_compare_toplist);

  NTopnodes = 1;
  TopNodes[0].Daughter = -1;
  TopNodes[0].Size = PEANOCELLS;
  TopNodes[0].StartKey = 0;
  TopNodes[0].Count = All.TotNumPart;
  TopNodes[0].Pstart = 0;
  TopNodes[0].Blocks = ntop;

  domain_topsplit(0, 0);

  myfree(toplist);
  myfree(ntopoffset);
  myfree(ntopnodelist);
  myfree(toplist_local);

}




void domain_topsplit_local(int node, peanokey startkey)
{
  int i, p, sub, bin;

  if(TopNodes[node].Size >= 8)
    {
      TopNodes[node].Daughter = NTopnodes;

      for(i = 0; i < 8; i++)
	{
	  if(NTopnodes < MAXTOPNODES)
	    {
	      sub = TopNodes[node].Daughter + i;
	      TopNodes[sub].Size = TopNodes[node].Size / 8;
	      TopNodes[sub].Count = 0;
	      TopNodes[sub].GravCost = 0;
	      TopNodes[sub].Daughter = -1;
	      TopNodes[sub].StartKey = startkey + i * TopNodes[sub].Size;
	      TopNodes[sub].Pstart = TopNodes[node].Pstart;

	      NTopnodes++;
	    }
	  else
	    {
	      printf("task=%d: We are out of Topnodes. Increasing the constant MAXTOPNODES might help.\n",
		     ThisTask);
	      fflush(stdout);
	      endrun(13213);
	    }
	}

      for(p = TopNodes[node].Pstart; p < TopNodes[node].Pstart + TopNodes[node].Count; p++)
	{
	  bin = (KeySorted[p] - startkey) / (TopNodes[node].Size / 8);

	  if(bin < 0 || bin > 7)
	    {
	      printf("task=%d: something odd has happened here. bin=%d\n", ThisTask, bin);
	      fflush(stdout);
	      endrun(13123123);
	    }

	  sub = TopNodes[node].Daughter + bin;

	  if(TopNodes[sub].Count == 0)
	    TopNodes[sub].Pstart = p;

	  TopNodes[sub].Count++;
	  TopNodes[sub].GravCost += domain_particle_costfactor(p);
	}

      for(i = 0; i < 8; i++)
	{
	  sub = TopNodes[node].Daughter + i;

	  if(TopNodes[sub].GravCost > totgravcost / (8.0 * TOPNODEFACTOR * NTask))
	    domain_topsplit_local(sub, TopNodes[sub].StartKey);
	}
    }
}




void domain_topsplit(int node, peanokey startkey)
{
  int i, p, sub, bin;

  if(TopNodes[node].Size >= 8)
    {
      TopNodes[node].Daughter = NTopnodes;

      for(i = 0; i < 8; i++)
	{
	  if(NTopnodes < MAXTOPNODES)
	    {
	      sub = TopNodes[node].Daughter + i;
	      TopNodes[sub].Size = TopNodes[node].Size / 8;
	      TopNodes[sub].Count = 0;
	      TopNodes[sub].GravCost = 0;
	      TopNodes[sub].Blocks = 0;
	      TopNodes[sub].Daughter = -1;
	      TopNodes[sub].StartKey = startkey + i * TopNodes[sub].Size;
	      TopNodes[sub].Pstart = TopNodes[node].Pstart;
	      NTopnodes++;
	    }
	  else
	    {
	      printf("Task=%d: We are out of Topnodes. Increasing the constant MAXTOPNODES might help.\n",
		     ThisTask);
	      fflush(stdout);
	      endrun(137213);
	    }
	}

      for(p = TopNodes[node].Pstart; p < TopNodes[node].Pstart + TopNodes[node].Blocks; p++)
	{
	  bin = (toplist[p].Startkey - startkey) / (TopNodes[node].Size / 8);
	  sub = TopNodes[node].Daughter + bin;

	  if(bin < 0 || bin > 7)
	    endrun(77);

	  if(TopNodes[sub].Blocks == 0)
	    TopNodes[sub].Pstart = p;

	  TopNodes[sub].Count += toplist[p].Count;
	  TopNodes[sub].GravCost += toplist[p].GravCost;
	  TopNodes[sub].Blocks++;
	}

      for(i = 0; i < 8; i++)
	{
	  sub = TopNodes[node].Daughter + i;

	  if(TopNodes[sub].GravCost > totgravcost / (TOPNODEFACTOR * NTask))
	    domain_topsplit(sub, TopNodes[sub].StartKey);
	}
    }
}



int domain_compare_toplist(const void *a, const void *b)
{
  if(((struct topnode_exchange *) a)->Startkey < (((struct topnode_exchange *) b)->Startkey))
    return -1;

  if(((struct topnode_exchange *) a)->Startkey > (((struct topnode_exchange *) b)->Startkey))
    return +1;

  return 0;
}


int domain_compare_key(const void *a, const void *b)
{
  if(*(peanokey *) a < *(peanokey *) b)
    return -1;

  if(*(peanokey *) a > *(peanokey *) b)
    return +1;

  return 0;
}




/*

test()
{
  int i;
  peanokey key, keymin, keymax;
  
  keymax=0;
  keymin= (1<<30);

  for(i=0;i<NumPart;i++)
    {
      key=peano_hilbert_key((P[i].Pos[0] - DomainCorner[0]) * DomainFac,
			    (P[i].Pos[1] - DomainCorner[1]) * DomainFac,
			    (P[i].Pos[2] - DomainCorner[2]) * DomainFac, BITS_PER_DIMENSION);

      if(key<keymin)
	keymin=key;
      if(key>keymax)
	keymax=key;

    }

  printf("Task=%d keymin=%d keymax=%d\n", ThisTask, keymin, keymax);
  fflush(stdout);
}



dump_split()
{
  FILE *fd;
  int i;

  if(ThisTask==0)
    {
      fd= fopen("split.dat", "w");

      fwrite(&NTopleaves, sizeof(int), 1 , fd);
      fwrite(DomainCount, sizeof(int), NTopleaves , fd);
      fwrite(DomainWork, sizeof(double), NTopleaves, fd);
      fwrite(DomainTask, sizeof(int), NTopleaves, fd);
      fwrite(list_load, sizeof(int), NTask, fd);
      fwrite(list_loadsph, sizeof(int), NTask, fd);
      fwrite(list_work, sizeof(double), NTask, fd);
      fwrite(DomainStartList, sizeof(int), NTask, fd);
      fwrite(DomainEndList, sizeof(int), NTask, fd);

      for(i=0;i< NTopnodes;i++)
        {
          if(TopNodes[i].Daughter == -1)
            fwrite(&TopNodes[i].Size, sizeof(peanokey), 1, fd);
        }

      for(i=0;i< NTopnodes;i++)
        {
          if(TopNodes[i].Daughter == -1)
            fwrite(&TopNodes[i].Count, sizeof(long long), 1, fd);
        }

      for(i=0;i< NTopnodes;i++)
        {
          if(TopNodes[i].Daughter == -1)
            fwrite(&TopNodes[i].GravCost, sizeof(float), 1, fd);
        }

      fclose(fd);

      printf("maxload=%d\n", (int)maxload);
      printf("maxloadsph=%d\n", (int)maxloadsph);
    }

  endrun(0);
}


*/
