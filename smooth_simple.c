#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"

#if ( defined(CONDUCTION) || defined(CR_DIFFUSION) || defined(SMOOTH_PHI) || defined(VOLUME_CORRECTION))

#ifdef VOLUME_CORRECTION
#define ANZ_VOL_FAC 500
static MyFloat rtab[ANZ_VOL_FAC];
static MyFloat ntab[ANZ_VOL_FAC];
#endif /* VOLUME_CORRECTION */


void compute_smoothed_values(void)
{
  int *noffset, *nbuffer, *nsend, *nsend_local, *numlist, *ndonelist;
  int i, j, n;
  int ndone;
  long long ntot, ntotleft;
  int maxfill, source;
  int level, ngrp, sendTask, recvTask;
  int place, nexport;
  MPI_Status status;
#ifdef VOLUME_CORRECTION
  double dV, hinv, hinv3, norm0, dr, rmax, u, wk;
  int ic;
#endif

  /* Display information message that this step is executed on Task 0 ... */
  if(ThisTask == 0)
    {
      printf("Updating SPH interpolants for:"
#ifdef CONDUCTION
	     " (temperature)"
#endif /* CONDUCTION */
#ifdef CR_DIFFUSION
	     " (CR diffusivity terms)"
#endif /* CR_DIFFUSION */
#ifdef SMOOTH_PHI
	     " (divB cleaning phi)"
#endif /* SMOOTH_PHI */
#ifdef VOLUME_CORRECTION
             " (Rho and GradRho)"
#endif /* VOLUME_CORRECTION */
	     "\n");
    }

  noffset = mymalloc(sizeof(int) * NTask);	/* offsets of bunches in common list */
  nbuffer = mymalloc(sizeof(int) * NTask);
  nsend_local = mymalloc(sizeof(int) * NTask);
  nsend = mymalloc(sizeof(int) * NTask * NTask);
  ndonelist = mymalloc(sizeof(int) * NTask);

  for(n = 0, NumSphUpdate = 0; n < N_gas; n++)
    {
      if(P[n].Type == 0)
	{
	  if(P[n].Ti_endstep == All.Ti_Current)
	    NumSphUpdate++;
	}
    }

  printf("Task %d: Updating %d particles\n",ThisTask,NumSphUpdate);

  numlist = mymalloc(NTask * sizeof(int) * NTask);
  MPI_Allgather(&NumSphUpdate, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 0, ntot = 0; i < NTask; i++)
    ntot += numlist[i];
  myfree(numlist);


  /* we will repeat the whole thing for those particles where we didn't find enough neighbours */

  i = 0;			/* beginn with this index */
  ntotleft = ntot;		/* particles left for all tasks together */

  while(ntotleft > 0)
    {
      for(j = 0; j < NTask; j++)
	nsend_local[j] = 0;

      /* do local particles and prepare export list */

      for(nexport = 0, ndone = 0; i < NumPart && nexport < All.BunchSizeDensity - NTask; i++)
	if(P[i].Type == 0)
	  if(P[i].Ti_endstep == All.Ti_Current)
	    {
	      ndone++;

	      for(j = 0; j < NTask; j++)
		Exportflag[j] = 0;

	      compute_smoothed_evaluate(i, 0);

	      for(j = 0; j < NTask; j++)
		{
		  if(Exportflag[j])
		    {
		      DensDataIn[nexport].Pos[0] = P[i].Pos[0];
		      DensDataIn[nexport].Pos[1] = P[i].Pos[1];
		      DensDataIn[nexport].Pos[2] = P[i].Pos[2];

		      DensDataIn[nexport].Hsml = PPP[i].Hsml;
#ifdef VOLUME_CORRECTION
		      DensDataIn[nexport].Density = SphP[i].a2.Density;
#endif
		      DensDataIn[nexport].Index = i;
		      DensDataIn[nexport].Task = j;
		      nexport++;
		      nsend_local[j]++;
		    }
		}
	    }


      qsort(DensDataIn, nexport, sizeof(struct densdata_in), dens_compare_key);

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
				   recvTask, TAG_CONDUCT_A,
				   &DensDataGet[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct densdata_in),
				   MPI_BYTE, recvTask, TAG_CONDUCT_A, MPI_COMM_WORLD, &status);
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }


	  for(j = 0; j < nbuffer[ThisTask]; j++)
	    compute_smoothed_evaluate(j, 1);


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
				   MPI_BYTE, recvTask, TAG_CONDUCT_B,
				   &DensDataPartialResult[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct densdata_out),
				   MPI_BYTE, recvTask, TAG_CONDUCT_B, MPI_COMM_WORLD, &status);

		      /* add the result to the particles */
		      for(j = 0; j < nsend_local[recvTask]; j++)
			{
			  source = j + noffset[recvTask];
			  place = DensDataIn[source].Index;
#ifdef CONDUCTION
			  SphP[place].SmoothedEntr += DensDataPartialResult[source].SmoothedEntr;
#ifdef CONDUCTION_SATURATION
			  SphP[place].GradEntr[0] += DensDataPartialResult[source].GradEntr[0];
			  SphP[place].GradEntr[1] += DensDataPartialResult[source].GradEntr[1];
			  SphP[place].GradEntr[2] += DensDataPartialResult[source].GradEntr[2];
#endif
#endif /* CONDUCTION */

#ifdef CR_DIFFUSION
			  SphP[place].CR_SmoothE0 += DensDataPartialResult[source].CR_SmoothE0;
			  SphP[place].CR_Smoothn0 += DensDataPartialResult[source].CR_Smoothn0;
#endif /* CR_DIFFUSION */

#ifdef SMOOTH_PHI
			  SphP[place].SmoothPhi += DensDataPartialResult[source].SmoothPhi;
#endif /* SMOOTH_PHI */
#ifdef VOLUME_CORRECTION
			  SphP[place].GradDensity[0] += DensDataPartialResult[source].GradDensity[0];
			  SphP[place].GradDensity[1] += DensDataPartialResult[source].GradDensity[1];
			  SphP[place].GradDensity[2] += DensDataPartialResult[source].GradDensity[2];
			  SphP[place].NormDensity += DensDataPartialResult[source].NormDensity;
#endif /* VOLUME_CORRECTION */
			}
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }
	  level = ngrp - 1;
	}

      MPI_Allgather(&ndone, 1, MPI_INT, ndonelist, 1, MPI_INT, MPI_COMM_WORLD);
      for(j = 0; j < NTask; j++)
	ntotleft -= ndonelist[j];
    }


  /* do final operations on results */
  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 0)
	if(P[i].Ti_endstep == All.Ti_Current)
	  {

#ifdef CONDUCTION
	    SphP[i].SmoothedEntr /= pow(SphP[i].a2.Density, GAMMA);
#ifdef CONDUCTION_SATURATION
	    SphP[i].GradEntr[0] /= pow(SphP[i].a2.Density, GAMMA);
	    SphP[i].GradEntr[1] /= pow(SphP[i].a2.Density, GAMMA);
	    SphP[i].GradEntr[2] /= pow(SphP[i].a2.Density, GAMMA);
#endif
#endif /* CONDUCTION */


#ifdef CR_DIFFUSION
	    SphP[i].CR_SmoothE0 /= SphP[i].a2.Density;
	    SphP[i].CR_Smoothn0 /= SphP[i].a2.Density;

	    /* Limit smoothed values so small-value particles
	     * will not be pulled into negative energy regimes
	     * by their neighbors due to far-too-high smoothed
	     * estimate */
	    if(SphP[i].CR_SmoothE0 > 2.0 * SphP[i].CR_E0)
	      SphP[i].CR_SmoothE0 = 2.0 * SphP[i].CR_E0;
	  
	    if(SphP[i].CR_Smoothn0 > 2.0 * SphP[i].CR_n0)
	      SphP[i].CR_Smoothn0 = 2.0 * SphP[i].CR_n0;

#endif /* CR_DIFFUSION */

#ifdef SMOOTH_PHI
	    SphP[i].SmoothPhi /= SphP[i].a2.Density;
#endif /* SMOOTH_PHI */

#ifdef VOLUME_CORRECTION
	    SphP[i].GradDensity[0] /= SphP[i].a2.Density;
	    SphP[i].GradDensity[1] /= SphP[i].a2.Density;
	    SphP[i].GradDensity[2] /= SphP[i].a2.Density;
	    /* 0th order density correction */

            ic = 0;
            dV = P[i].Mass / SphP[i].a2.Density;
            dr = pow(dV,1./3.);
            rmax = PPP[i].Hsml / dr;

            hinv = 1.0 / PPP[i].Hsml;
#ifndef  TWODIMS
            hinv3 = hinv * hinv * hinv;
#else
            hinv3 = hinv * hinv / boxSize_Z;
#endif
            norm0 = 0;
	    while(rtab[ic] < rmax)
	      {
		u = rtab[ic] * dr * hinv;
		if(u < 0.5)
		  wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
		else
		  wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

		norm0 += wk * ntab[ic] * dV;
                ic++;
	      }
	    /*
            printf("Task %d: %e %e\n",ThisTask,SphP[i].a2.Density,SphP[i].NormDensity);
	    */
#if defined(VOLUME_BIWEIGHT) || defined (VOLUME_QUADRIC) || defined (VOLUME_QUINTIC)
            SphP[i].a2.Density = SphP[i].NormDensity;
#else
            SphP[i].a2.Density = SphP[i].NormDensity;
#endif
#endif /* VOLUME_CORRECTION */
	  }
    }

  myfree(ndonelist);
  myfree(nsend);
  myfree(nsend_local);
  myfree(nbuffer);
  myfree(noffset);
}



/*! This function represents the core of the SPH density computation. The
*  target particle may either be local, or reside in the communication
*  buffer.
*/
void compute_smoothed_evaluate(int target, int mode)
{
  int j, n;
  int startnode, numngb_inbox;
  double h, h2, hinv, hinv3, hinv4;
  double wk, dwk;
  double dx, dy, dz, r, r2, u, mass_j;
  FLOAT *pos;

#ifdef CONDUCTION
  double smoothentr;

#ifdef CONDUCTION_SATURATION
  double gradentr[3];
#endif
#endif /* CONDUCTION */

#ifdef CR_DIFFUSION
  double rCR_SmoothE0 = 0.0;
  double rCR_Smoothn0 = 0.0;
#endif /* CR_DIFFUSION */

#ifdef SMOOTH_PHI
  double SmoothPhi = 0.0;
#endif /* SMOOTH_PHI */

#ifdef VOLUME_CORRECTION
  double graddensity[3],density,normdensity,dV;
#endif /* VOLUME_CORRECTION */

#ifdef CONDUCTION
  smoothentr = 0;

#ifdef CONDUCTION_SATURATION
  gradentr[0] = gradentr[1] = gradentr[2] = 0;
#endif
#endif /* CONDUCTION */

#ifdef VOLUME_CORRECTION
  graddensity[0] = graddensity[1] = graddensity[2] = normdensity = 0;
#endif /* VOLUME_CORRECTION */

  if(mode == 0)
    {
      pos = P[target].Pos;
      h = PPP[target].Hsml;
#ifdef VOLUME_CORRECTION
      density = SphP[target].a2.Density;
#endif
    }
  else
    {
      pos = DensDataGet[target].Pos;
      h = DensDataGet[target].Hsml;
#ifdef VOLUME_CORRECTION
      density = DensDataGet[target].Density;
#endif
    }

#ifdef VOLUME_CORRECTION
    dV = P[target].Mass / density;
#endif

  h2 = h * h;
  hinv = 1.0 / h;
#ifndef  TWODIMS
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv * hinv / boxSize_Z;
#endif
  hinv4 = hinv3 * hinv;



  startnode = All.MaxPart;

  do
    {
      numngb_inbox = ngb_treefind_variable(&pos[0], h, &startnode);

      for(n = 0; n < numngb_inbox; n++)
	{
	  j = Ngblist[n];

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

#ifdef CONDUCTION
	      smoothentr += mass_j * wk * pow(SphP[j].a2.Density, GAMMA_MINUS1) * SphP[j].Entropy;

#ifdef CONDUCTION_SATURATION
	      if(r > 0)
		{
		  gradentr[0] +=
		    mass_j * dwk * dx / r * SphP[j].Entropy * pow(SphP[j].a2.Density, GAMMA_MINUS1);
		  gradentr[1] +=
		    mass_j * dwk * dy / r * SphP[j].Entropy * pow(SphP[j].a2.Density, GAMMA_MINUS1);
		  gradentr[2] +=
		    mass_j * dwk * dz / r * SphP[j].Entropy * pow(SphP[j].a2.Density, GAMMA_MINUS1);
		}
#endif
#endif /* CONDUCTION */

#ifdef CR_DIFFUSION
	      rCR_SmoothE0 += mass_j * wk * SphP[j].CR_E0;
	      rCR_Smoothn0 += mass_j * wk * SphP[j].CR_n0;
#endif /* CR_DIFFUSION */

#ifdef SMOOTH_PHI
	      SmoothPhi += mass_j * wk * SphP[j].PhiPred;
#endif /* SMOOTH_PHI */

#ifdef VOLUME_CORRECTION
	      if(r > 0)
		{
		  graddensity[0] += mass_j * dwk * dx / r * (density - SphP[j].a2.Density);
		  graddensity[1] += mass_j * dwk * dy / r * (density - SphP[j].a2.Density);
		  graddensity[2] += mass_j * dwk * dz / r * (density - SphP[j].a2.Density);
		}
#ifdef VOLUME_BIWEIGHT
	      wk =  NORM_COEFF_2 * (1 - u * u) * (1 - u * u);
#endif
#ifdef VOLUME_QUADRIC
	      if(u < 0.2)
		wk = NORM_COEFF_4 * ( pow(5 * u + 5,4) - 5 * pow(5 * u + 3,4) + 10 * pow(5 * u + 1,4));
	      else
		if(u < 0.6)
		  wk = NORM_COEFF_4 * ( pow(5 - 5 * u,4) - 5 * pow(3 - 5 * u,4));
		else
		  wk = NORM_COEFF_4 * ( pow(5 - 5 * u,4));
#endif
#ifdef VOLUME_QUINTIC
	      if(u < 1./3.)
		wk = NORM_COEFF_5 * ( pow(3 - 3 * u,5) - 6 * pow(2 - 3 * u,5) + 15 * pow(1 - 3 * u,5));
	      else
		if(u < 2./3)
		  wk = NORM_COEFF_5 * ( pow(3 - 3 * u,5) - 6 * pow(2 - 3 * u,5));
		else
		  wk = NORM_COEFF_5 * ( pow(3 - 3 * u,5));
#endif
#if defined(VOLUME_BIWEIGHT) || defined(VOLUME_QUADRIC) || defined(VOLUME_QUINTIC)
	      normdensity += mass_j * hinv3 * pow(density/SphP[j].a2.Density,1.75) * wk;
#else
	      normdensity += mass_j * pow(density/SphP[j].a2.Density,1.75) * wk;
#endif
#endif /* VOLUME_CORRECTION */
	    }
	}
    }
  while(startnode >= 0);


  if(mode == 0)
    {

#ifdef CONDUCTION
      SphP[target].SmoothedEntr = smoothentr;
#ifdef CONDUCTION_SATURATION
      SphP[target].GradEntr[0] = gradentr[0];
      SphP[target].GradEntr[1] = gradentr[1];
      SphP[target].GradEntr[2] = gradentr[2];
#endif
#endif /* CONDUCTION */

#ifdef CR_DIFFUSION
      SphP[target].CR_SmoothE0 = rCR_SmoothE0;
      SphP[target].CR_Smoothn0 = rCR_Smoothn0;
#endif /* CR_DIFFUSION */

#ifdef SMOOTH_PHI
      SphP[target].SmoothPhi = SmoothPhi;
#endif /* SMOOTH_PHI */

#ifdef VOLUME_CORRECTION
      SphP[target].GradDensity[0] = graddensity[0];
      SphP[target].GradDensity[1] = graddensity[1];
      SphP[target].GradDensity[2] = graddensity[2];
      SphP[target].NormDensity = normdensity;
#endif /* VOLUME_CORRECTION */
    }
  else
    {

#ifdef CONDUCTION
      DensDataResult[target].SmoothedEntr = smoothentr;
#ifdef CONDUCTION_SATURATION
      DensDataResult[target].GradEntr[0] = gradentr[0];
      DensDataResult[target].GradEntr[1] = gradentr[1];
      DensDataResult[target].GradEntr[2] = gradentr[2];
#endif
#endif /* CONDUCTION */

#ifdef CR_DIFFUSION
      DensDataResult[target].CR_SmoothE0 = rCR_SmoothE0;
      DensDataResult[target].CR_Smoothn0 = rCR_Smoothn0;
#endif /* CR_DIFFUSION */

#ifdef SMOOTH_PHI
      DensDataResult[target].SmoothPhi = SmoothPhi;
#endif /* SMOOTH_PHI */

#ifdef VOLUME_CORRECTION
      DensDataResult[target].GradDensity[0] = graddensity[0];
      DensDataResult[target].GradDensity[1] = graddensity[1];
      DensDataResult[target].GradDensity[2] = graddensity[2];
      DensDataResult[target].NormDensity = normdensity;
#endif /* VOLUME_CORRECTION */
    }
}

#ifdef VOLUME_CORRECTION

/* set-up table with radial distribution of weights
 */

void vol_weights_init(void)
{
  int count = 0;
  char buf[200], buf1[200], buf2[200];
  FILE *fd;

  if((fd = fopen(All.VolCorrectFile, "r")))
    {
      if(ThisTask == 0)
	{
	  printf("\nreading volume corection from file `%s'\n", All.VolCorrectFile);
	  fflush(stdout);
	}
      count = 0;
      while(!feof(fd) && count < ANZ_VOL_FAC)
	{
	  if(fgets(buf, 200, fd))
	    {
	      if(sscanf(buf, "%s%s", buf1, buf2) < 2)
		{
		  if(ThisTask == 0)
		    {
		      printf("Wrong syntax in file `%s', line %d\n", All.VolCorrectFile, count);
		      fflush(stdout);
		    }
		  endrun(0);
		}
	      rtab[count] = atof(buf1);
	      ntab[count] = atof(buf2);
	      count++;
	    }
	}
      fclose(fd);
      if(count >= ANZ_VOL_FAC - 1)
	{
	  if(ThisTask == 0)
	    {
	      printf("File `%s' contains to many datapoints, increase ANZ_VOL_FAC !\n", All.VolCorrectFile);
	      fflush(stdout);
	    }
	  endrun(0);
	}
      else
	{
	  while(count < ANZ_VOL_FAC)
	    {
	      rtab[count] = 1e20;
	      ntab[count] = 0.0;
	      count++;
	    }
	}
    }
  else
    {
      if(ThisTask == 0)
	{
	  printf("\nFile `%s' not found !\n", All.VolCorrectFile);
	  fflush(stdout);
	}
      endrun(0);
    }
}

#endif /* VOLUME_CORRECTION */


#endif
