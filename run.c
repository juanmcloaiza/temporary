#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <unistd.h>
#include <ctype.h>

#include "allvars.h"
#include "proto.h"

/*! \file run.c
 *  \brief  iterates over timesteps, main loop
 */

/*! This routine contains the main simulation loop that iterates over
 * single timesteps. The loop terminates when the cpu-time limit is
 * reached, when a `stop' file is found in the output directory, or
 * when the simulation ends because we arrived at TimeMax.
 */
void run(void)
{
  FILE *fd;

#ifdef RADIATION
  int ifunc;
#endif
  int i, stopflag = 0;
  char buf[200], stopfname[200], contfname[200];
  double t0, t1;
  //Added by JM
  int * temp;
  //End of Added by JM


  sprintf(stopfname, "%sstop", All.OutputDir);
  sprintf(contfname, "%scont", All.OutputDir);
  unlink(contfname);

  do				/* main loop */
  {
    t0 = second();

    for(i = 0; i < CPU_PARTS; i++)
      CPU_Step[i] = 0;

    find_next_sync_point_and_drift();	/* find next synchronization point and drift particles to this time.
                       * If needed, this function will also write an output file
                       * at the desired time.
                       */

    every_timestep_stuff();	/* write some info to log-files */
    //Added by JM
    accrete_particles();
    //    All.NumForcesSinceLastDomainDecomp = 1 + All.TreeDomainUpdateFrequency * All.TotNumPart;
    //End of Added by JM

#ifdef RADIATION
    ifunc = init_rad(All.Time);
#endif

#ifdef COOLING
#ifndef BG_COOLING
    IonizeParams();		/* set UV background for the current time */
#endif
#endif

    DomainDecomposition();	/* do domain decomposition if needed */


    compute_accelerations(0);	/* compute accelerations for 
                   * the particles that are to be advanced  
                   */

    /* check whether we want a full energy statistics */
    if((All.Time - All.TimeLastStatistics) >= All.TimeBetStatistics)
    {
#ifdef COMPUTE_POTENTIAL_ENERGY
      compute_potential();
#endif
      energy_statistics();	/* compute and output energy statistics */
      All.TimeLastStatistics += All.TimeBetStatistics;
    }

    advance_and_find_timesteps();	/* 'kick' active particles in
                     * momentum space and compute new
                     * timesteps for them
                     */

    t1 = second();
    All.CPU_Total += timediff(t0, t1);
    CPU_Step[CPU_ALL] = timediff(t0, t1);

    write_cpu_log();		/* produce some CPU usage info */

    All.NumCurrentTiStep++;

    /* Check whether we need to interrupt the run */
    if(ThisTask == 0)
    {
      /* Is the stop-file present? If yes, interrupt the run. */
      if((fd = fopen(stopfname, "r")))
      {
        fclose(fd);
        stopflag = 1;
        unlink(stopfname);
      }

      /* are we running out of CPU-time ? If yes, interrupt run. */
      if(CPUThisRun > 0.85 * All.TimeLimitCPU)
      {
        printf("reaching time-limit. stopping.\n");
        stopflag = 2;
      }
    }

    MPI_Bcast(&stopflag, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(stopflag)
    {
      restart(0);		/* write restart file */
      MPI_Barrier(MPI_COMM_WORLD);

      if(stopflag == 2 && ThisTask == 0)
      {
        if((fd = fopen(contfname, "w")))
          fclose(fd);
      }

      if(stopflag == 2 && All.ResubmitOn && ThisTask == 0)
      {
        close_outputfiles();
        sprintf(buf, "%s", All.ResubmitCommand);
#ifndef NOCALLSOFSYSTEM
        system(buf);
#endif
      }
      return;
    }

    /* is it time to write a regular restart-file? (for security) */
    if(ThisTask == 0)
    {
      if((CPUThisRun - All.TimeLastRestartFile) >= All.CpuTimeBetRestartFile)
      {
        All.TimeLastRestartFile = CPUThisRun;
        stopflag = 3;
      }
      else
        stopflag = 0;
    }

    MPI_Bcast(&stopflag, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(stopflag == 3)
    {
      restart(0);		/* write an occasional restart file */
      stopflag = 0;
    }

    t1 = second();

    CPUThisRun += timediff(t0, t1);
  }
  while(All.Ti_Current < TIMEBASE && All.Time <= All.TimeMax);

  restart(0);

  savepositions(All.SnapshotFileCount++);	/* write a last snapshot
                       * file at final time (will
                       * be overwritten if
                       * All.TimeMax is increased
                       * and the run is continued)
                       */
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
#endif
}


/*! This function finds the next synchronization point of the system
 * (i.e. the earliest point of time any of the particles needs a force
 * computation), and drifts the system to this point of time.  If the
 * system dirfts over the desired time of a snapshot file, the
 * function will drift to this moment, generate an output, and then
 * resume the drift.
 */
void find_next_sync_point_and_drift(void)
{
  int n, min, min_glob, flag, *temp;
  double timeold;
  double t0, t1;

  t0 = second();

  timeold = All.Time;

  for(n = 1, min = P[0].Ti_endstep; n < NumPart; n++)
    if(min > P[n].Ti_endstep)
      min = P[n].Ti_endstep;

  MPI_Allreduce(&min, &min_glob, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  /* We check whether this is a full step where all particles are synchronized */
  flag = 1;
  for(n = 0; n < NumPart; n++)
    if(P[n].Ti_endstep > min_glob)
      flag = 0;

  MPI_Allreduce(&flag, &Flag_FullStep, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

#ifdef PMGRID
  if(min_glob >= All.PM_Ti_endstep)
  {
    min_glob = All.PM_Ti_endstep;
    Flag_FullStep = 1;
  }
#endif

  /* Determine 'NumForceUpdate', i.e. the number of particles on this processor that are going to be active */
  for(n = 0, NumForceUpdate = 0; n < NumPart; n++)
  {
    if(P[n].Ti_endstep == min_glob)
      NumForceUpdate++;
  }

  /* note: NumForcesSinceLastDomainDecomp has type "long long" */
  temp = (int *) mymalloc(NTask * sizeof(int));
  MPI_Allgather(&NumForceUpdate, 1, MPI_INT, temp, 1, MPI_INT, MPI_COMM_WORLD);
  for(n = 0; n < NTask; n++)
    All.NumForcesSinceLastDomainDecomp += temp[n];
  myfree(temp);

  t1 = second();
  All.CPU_Predict += timediff(t0, t1);
  CPU_Step[CPU_DRIFT] += timediff(t0, t1);


  while((min_glob >= All.Ti_nextoutput && All.Ti_nextoutput >= 0)
      || (min_glob >= All.Ti_nextlineofsight && All.Ti_nextlineofsight > 0))
  {
    if(All.Ti_nextlineofsight > 0 && All.Ti_nextlineofsight < All.Ti_nextoutput)
    {
      move_particles(All.Ti_Current, All.Ti_nextlineofsight);
      All.Ti_Current = All.Ti_nextlineofsight;
    }
    else
    {
      move_particles(All.Ti_Current, All.Ti_nextoutput);
      All.Ti_Current = All.Ti_nextoutput;
    }

    if(All.ComovingIntegrationOn)
      All.Time = All.TimeBegin * exp(All.Ti_Current * All.Timebase_interval);
    else
      All.Time = All.TimeBegin + All.Ti_Current * All.Timebase_interval;
#ifdef TIMEDEPGRAV
    All.G = All.Gini * dGfak(All.Time);
#endif


    if(All.Ti_Current == All.Ti_nextoutput)
    {
#ifdef OUTPUTPOTENTIAL
      All.NumForcesSinceLastDomainDecomp =
        (long long) (1 + All.TotNumPart * All.TreeDomainUpdateFrequency);
      DomainDecomposition();
      compute_potential();
#endif

#ifdef LT_STELLAREVOLUTION
      if(All.Ti_Current == TIMEBASE)
        evolve_SN(EVOLVE_SN);
#endif
#ifdef LT_SEvDbg
      get_metals_sumcheck(9);
#endif

      savepositions(All.SnapshotFileCount++);	/* write snapshot file */

      All.Ti_nextoutput = find_next_outputtime(All.Ti_nextoutput + 1);
    }

    if(All.Ti_Current == All.Ti_nextlineofsight)
    {
#ifdef OUTPUTLINEOFSIGHT
      lineofsight_output();
      All.Ti_nextlineofsight = find_next_lineofsighttime(All.Ti_nextlineofsight);
#endif
    }
  }

  move_particles(All.Ti_Current, min_glob);

  All.Ti_Current = min_glob;

  if(All.ComovingIntegrationOn)
    All.Time = All.TimeBegin * exp(All.Ti_Current * All.Timebase_interval);
  else
    All.Time = All.TimeBegin + All.Ti_Current * All.Timebase_interval;
#ifdef TIMEDEPGRAV
  All.G = All.Gini * dGfak(All.Time);
#endif

#ifdef LT_STELLAREVOLUTION
  All.Time_Age = get_age(All.Time);
#endif

  All.TimeStep = All.Time - timeold;
}



/*! this function returns the next output time that is equal or larger to
 *  ti_curr
 */
int find_next_outputtime(int ti_curr)
{
  int i, ti, ti_next, iter = 0;
  double next, time;

  ti_next = -1;

  if(All.OutputListOn)
  {
    for(i = 0; i < All.OutputListLength; i++)
    {
      time = All.OutputListTimes[i];

      if(time >= All.TimeBegin && time <= All.TimeMax)
      {
        if(All.ComovingIntegrationOn)
          ti = (int) (log(time / All.TimeBegin) / All.Timebase_interval);
        else
          ti = (int) ((time - All.TimeBegin) / All.Timebase_interval);

        if(ti >= ti_curr)
        {
          if(ti_next == -1)
            ti_next = ti;

          if(ti_next > ti)
            ti_next = ti;
        }
      }
    }
  }
  else
  {
    if(All.ComovingIntegrationOn)
    {
      if(All.TimeBetSnapshot <= 1.0)
      {
        printf("TimeBetSnapshot > 1.0 required for your simulation.\n");
        endrun(13123);
      }
    }
    else
    {
      if(All.TimeBetSnapshot <= 0.0)
      {
        printf("TimeBetSnapshot > 0.0 required for your simulation.\n");
        endrun(13123);
      }
    }

    time = All.TimeOfFirstSnapshot;

    iter = 0;

    while(time < All.TimeBegin)
    {
      if(All.ComovingIntegrationOn)
        time *= All.TimeBetSnapshot;
      else
        time += All.TimeBetSnapshot;

      iter++;

      if(iter > 1000000)
      {
        printf("Can't determine next output time.\n");
        endrun(110);
      }
    }

    while(time <= All.TimeMax)
    {
      if(All.ComovingIntegrationOn)
        ti = (int) (log(time / All.TimeBegin) / All.Timebase_interval);
      else
        ti = (int) ((time - All.TimeBegin) / All.Timebase_interval);

      if(ti >= ti_curr)
      {
        ti_next = ti;
        break;
      }

      if(All.ComovingIntegrationOn)
        time *= All.TimeBetSnapshot;
      else
        time += All.TimeBetSnapshot;

      iter++;

      if(iter > 1000000)
      {
        printf("Can't determine next output time.\n");
        endrun(111);
      }
    }
  }

  if(ti_next == -1)
  {
    ti_next = 2 * TIMEBASE;	/* this will prevent any further output */

    if(ThisTask == 0)
      printf("\nThere is no valid time for a further snapshot file.\n");
  }
  else
  {
    if(All.ComovingIntegrationOn)
      next = All.TimeBegin * exp(ti_next * All.Timebase_interval);
    else
      next = All.TimeBegin + ti_next * All.Timebase_interval;

    if(ThisTask == 0)
      printf("\nSetting next time for snapshot file to Time_next= %g\n\n", next);
  }

  return ti_next;
}




/*! This routine writes one line for every timestep to two log-files.
 * In FdInfo, we just list the timesteps that have been done, while in
 * FdCPU the cumulative cpu-time consumption in various parts of the
 * code is stored.
 */
void every_timestep_stuff(void)
{
  double z,hubble_a;

  if(ThisTask == 0)
  {
    if(All.ComovingIntegrationOn)
    {
      z = 1.0 / (All.Time) - 1;
      fprintf(FdInfo, "\nBegin Step %d, Time: %g, Redshift: %g, Systemstep: %g, Dloga: %g\n",
          All.NumCurrentTiStep, All.Time, z, All.TimeStep,
          log(All.Time) - log(All.Time - All.TimeStep));
      printf("\nBegin Step %d, Time: %g, Redshift: %g, Systemstep: %g, Dloga: %g\n", All.NumCurrentTiStep,
          All.Time, z, All.TimeStep, log(All.Time) - log(All.Time - All.TimeStep));

#ifdef CHEMISTRY
      printf("Abundances elec: %g, HM: %g, H2I: %g, H2II: %g\n",
          SphP[1].elec, SphP[1].HM, SphP[1].H2I, SphP[1].H2II);
#endif


      fflush(FdInfo);
    }
    else
    {
      fprintf(FdInfo, "\nBegin Step %d, Time: %g, Systemstep: %g\n", All.NumCurrentTiStep, All.Time,
          All.TimeStep);
      printf("\nBegin Step %d, Time: %g, Systemstep: %g\n", All.NumCurrentTiStep, All.Time, All.TimeStep);
      fflush(FdInfo);
    }

#ifdef XXLINFO
    if(Flag_FullStep == 1)
    {
      fprintf(FdXXL, "%d %g ", All.NumCurrentTiStep, All.Time);
#ifdef MAGNETIC
      fprintf(FdXXL, "%e ", MeanB);
#ifdef TRACEDIVB
      fprintf(FdXXL, "%e ", MaxDivB);
#endif
#endif
#ifdef TIME_DEP_ART_VISC
      fprintf(FdXXL, "%f ", MeanAlpha);
#endif
      fprintf(FdXXL, "\n");
      fflush(FdXXL);
    }
#endif

#ifdef DARKENERGY
    if(All.ComovingIntegrationOn == 1)
    {
      hubble_a = hubble_function(All.Time);
      fprintf(FdDE, "%d %g %e ", All.NumCurrentTiStep, All.Time, hubble_a);
#ifndef TIMEDEPDE
      fprintf(FdDE, "%e ", All.DarkEnergyParam);
#else
      fprintf(FdDE, "%e %e ", get_wa(All.Time), DarkEnergy_a(All.Time));
#endif
#ifdef TIMEDEPGRAV
      fprintf(FdDE, "%e %e", dHfak(All.Time), dGfak(All.Time));
#endif
      fprintf(FdDE, "\n");
      fflush(FdDE);
    }
#endif

  }

  set_random_numbers();
}



void write_cpu_log(void)
{
  double max_CPU_Step[CPU_PARTS], avg_CPU_Step[CPU_PARTS], t0, t1, tsum;
  int i;

  if(ThisTask == 0)
  {
    fprintf(FdCPU, "Step %d, Time: %g, CPUs: %d\n", All.NumCurrentTiStep, All.Time, NTask);
#ifdef SFR
#ifndef LT_STELLAREVOLUTION
    fprintf(FdCPU,
        "%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n",
        All.CPU_Total, All.CPU_Gravity, All.CPU_Hydro, All.CPU_Domain, All.CPU_Potential,
        All.CPU_Predict, All.CPU_TimeLine, All.CPU_Snapshot, All.CPU_TreeWalk, All.CPU_TreeConstruction,
        All.CPU_CommSum, All.CPU_Imbalance, All.CPU_HydCompWalk, All.CPU_HydCommSumm,
        All.CPU_HydImbalance, All.CPU_EnsureNgb, All.CPU_PM, All.CPU_Peano, All.CPU_SfrCool);
#else
    fprintf(FdCPU,
        "%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f ",
        All.CPU_Total, All.CPU_Gravity, All.CPU_Hydro, All.CPU_Domain, All.CPU_Potential,
        All.CPU_Predict, All.CPU_TimeLine, All.CPU_Snapshot, All.CPU_TreeWalk, All.CPU_TreeConstruction,
        All.CPU_CommSum, All.CPU_Imbalance, All.CPU_HydCompWalk, All.CPU_HydCommSumm,
        All.CPU_HydImbalance, All.CPU_EnsureNgb, All.CPU_PM, All.CPU_Peano, All.CPU_SfrCool);
#ifdef LT_SEv_INFO
    fprintf(FdCPU,
        "%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f ",
        All.CPU_SEv, All.CPU_SN_info, All.CPU_SN_find, All.CPU_SN_NeighFind, All.CPU_SN_NeighCheck,
        All.CPU_SN_Spread, All.CPU_SN_Calc, All.CPU_SN_Comm, All.CPU_SN_Imbalance);
#endif
#ifdef LT_EXTEGY_INFO
    fprintf(FdCPU, "%10.2f %10.2f", All.CPU_Eff_Iter, All.CPU_EE_info);
#endif
    fprintf(FdCPU, "\n");
#endif
#else
    fprintf(FdCPU,
        "%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n",
        All.CPU_Total, All.CPU_Gravity, All.CPU_Hydro, All.CPU_Domain, All.CPU_Potential,
        All.CPU_Predict, All.CPU_TimeLine, All.CPU_Snapshot, All.CPU_TreeWalk, All.CPU_TreeConstruction,
        All.CPU_CommSum, All.CPU_Imbalance, All.CPU_HydCompWalk, All.CPU_HydCommSumm,
        All.CPU_HydImbalance, All.CPU_EnsureNgb, All.CPU_PM, All.CPU_Peano);
#endif
    fflush(FdCPU);
  }

  MPI_Reduce(&CPU_Step, &max_CPU_Step, CPU_PARTS, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&CPU_Step, &avg_CPU_Step, CPU_PARTS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
  {
    for(i = 0; i < CPU_PARTS; i++)
      avg_CPU_Step[i] /= NTask;

    put_symbol(0.0, 1.0, '#');

    for(i = 1, tsum = 0.0; i < CPU_PARTS; i++)
    {
      t0 = tsum;
      t1 = tsum + avg_CPU_Step[i];
      put_symbol(t0 / max_CPU_Step[0], t1 / max_CPU_Step[0], CPU_Symbol[i]);
      tsum += t1 - t0;
    }

    for(i = 1; i < CPU_PARTS; i++)
    {
      t0 = tsum;
      t1 = tsum + max_CPU_Step[i] - avg_CPU_Step[i];
      put_symbol(t0 / max_CPU_Step[0], t1 / max_CPU_Step[0], CPU_SymbolImbalance[i]);
      tsum += t1 - t0;
    }

    put_symbol(tsum / max_CPU_Step[0], 1.0, '-');


    fprintf(FdBalance, "Step=%7d  sec=%10.3f  %s\n", All.NumCurrentTiStep, max_CPU_Step[0], CPU_String);
    fflush(FdBalance);
  }
}

void put_symbol(double t0, double t1, char c)
{
  int i, j;

  i = (int) (t0 * CPU_STRING_LEN + 0.5);
  j = (int) (t1 * CPU_STRING_LEN);

  if(i < 0)
    i = 0;
  if(j >= CPU_STRING_LEN)
    j = CPU_STRING_LEN;

  while(i <= j)
    CPU_String[i++] = c;

  CPU_String[CPU_STRING_LEN] = 0;
}


/*! This routine first calls a computation of various global
 * quantities of the particle distribution, and then writes some
 * statistics about the energies in the various particle components to
 * the file FdEnergy.
 */
void energy_statistics(void)
{
  compute_global_quantities_of_system();

  if(ThisTask == 0)
  {
    fprintf(FdEnergy,
        "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
        All.Time, SysState.EnergyInt, SysState.EnergyPot, SysState.EnergyKin, SysState.EnergyIntComp[0],
        SysState.EnergyPotComp[0], SysState.EnergyKinComp[0], SysState.EnergyIntComp[1],
        SysState.EnergyPotComp[1], SysState.EnergyKinComp[1], SysState.EnergyIntComp[2],
        SysState.EnergyPotComp[2], SysState.EnergyKinComp[2], SysState.EnergyIntComp[3],
        SysState.EnergyPotComp[3], SysState.EnergyKinComp[3], SysState.EnergyIntComp[4],
        SysState.EnergyPotComp[4], SysState.EnergyKinComp[4], SysState.EnergyIntComp[5],
        SysState.EnergyPotComp[5], SysState.EnergyKinComp[5], SysState.MassComp[0],
        SysState.MassComp[1], SysState.MassComp[2], SysState.MassComp[3], SysState.MassComp[4],
        SysState.MassComp[5]);

    fflush(FdEnergy);
  }
}
