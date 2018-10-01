#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <mpi.h>

#include "allvars.h"
#include "proto.h"



void read_ic_cluster_gas(char *fname)
{
#define BLOCKSIZE 50000
  FILE *fd = 0;
  int i, j, k, n;
  double sqr_a;
  float *fp, dummy[BLOCKSIZE][3];
  int pc, id, pc_here;
  int type, pr, left, groupsize;
  MPI_Status status;
  int npart;
  float a0, pmhr;
  char buf[100];
  int n_for_this_task, n_in_file;
  double massfac, posfac, velfac, TempTemp, BaryonFrac;
  int blocks, subblock;
  int nhr, nlr;
  int counttype3, counttype2, counttype0, N_init_gas;
  int nhr_blocks, nlr_blocks;
  double rho, d;
  int sumnumpart;

#ifdef T3E
  short int dummy4byte;		/* Note: int has 8 Bytes on the T3E ! */
#else
  int dummy4byte;
#endif

#define SKIP fread(&dummy4byte,sizeof(int),1,fd);


  BaryonFrac = All.OmegaBaryon / All.Omega0;


  massfac = 0.3 * 3 * 0.1 * 0.1 / (8 * M_PI * All.G) * pow(141300.0 / 760, 3);
  posfac = 141300.0;
  velfac = 14130.0;

  /*
     massfac= 0.3 * 3*0.1*0.1/ (8*M_PI*All.G) * pow(479000.0/512, 3);
     posfac=  479000.0;
     velfac=  47900.0;
   */



  /* Below, Bepi's new format is assumed !!!!!! */
  /* for the old one, the HR particle mass 'pmhr' has to be set
     by hand. SEE BELOW.
   */

  if(ThisTask == 0)
    {
      for(i = 0; i < 5; i++)
	{
	  All.MassTable[i] = 0;
	}


      if(!(fd = fopen(fname, "r")))
	{
	  printf("can't open file `%s'.\n", fname);
	  endrun(123);
	}

      printf("\nREADING FILE '%s'....\n\n", fname);
      fflush(stdout);

      SKIP;
      fread(&nhr, sizeof(int), 1, fd);
      fread(&nlr, sizeof(int), 1, fd);
      fread(&a0, sizeof(float), 1, fd);
      if(dummy4byte == 16)
	fread(&pmhr, sizeof(float), 1, fd);
      else
	{
	  pmhr = 12.8440488;	/* S1 *//* here set by hand, if necessary */
	  /* pmhr=  2.55189071; *//* S2 */
	}
      SKIP;


      All.MassTable[1] = (1 - BaryonFrac) * pmhr * massfac;	/* high-res particles */
      All.MassTable[0] = BaryonFrac * pmhr * massfac;

      printf("All.MassTable[1]=%g\n", All.MassTable[1]);

      printf("file contains %d HR and %d LR particles.\n", nhr, nlr);

      All.TotN_gas = nhr;	/* gas parts are put on top of the high-res particles */


      printf("\nN_sph: %d\nN_halo: %d\nN_disk: %d\n\n", (int) All.TotN_gas, nhr, nlr);



      All.TotNumPart = 2 * nhr + nlr;

      fclose(fd);
    }

  MPI_Bcast(&All, sizeof(struct global_data_all_processes), MPI_BYTE, 0, MPI_COMM_WORLD);

  All.MaxPart = (int) (All.PartAllocFactor * (All.TotNumPart / NTask));	/* sets the maximum number of particles that may 
									   reside on a processor */
  All.MaxPartSph = (int) (All.PartAllocFactor * (All.TotN_gas / NTask));

  allocate_memory();


  if(ThisTask == 0)		/*re-open and re-read the header */
    {
      sprintf(buf, "%s", fname);

      if(!(fd = fopen(buf, "r")))
	{
	  printf("can't open file `%s'.\n", buf);
	  endrun(123);
	}

      SKIP;
      fread(&nhr, sizeof(int), 1, fd);
      fread(&nlr, sizeof(int), 1, fd);
      fread(&a0, sizeof(float), 1, fd);
      if(dummy4byte == 16)
	fread(&pmhr, sizeof(float), 1, fd);
      SKIP;


      nhr_blocks = nhr / 1000000 + 1;	/* for lcdmtest, blocks=1 */
      nlr_blocks = nlr / 1000000 + 1;
    }

  MPI_Bcast(&nhr, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nhr_blocks, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nlr_blocks, 1, MPI_INT, 0, MPI_COMM_WORLD);


  for(blocks = 0, N_init_gas = 0; blocks < nhr_blocks; blocks++)
    {
      if(blocks < (nhr_blocks - 1))
	n_in_file = 1000000;
      else
	n_in_file = nhr % 1000000;

      n_for_this_task = n_in_file / NTask;
      if(ThisTask < (n_in_file % NTask))	/* simple division to get the correct number */
	n_for_this_task++;

      N_init_gas += n_for_this_task;	/* N_init_gas  will then hold the total(!) number of gas particles 
					   for this processor */
    }


  pc = 0 + N_init_gas;
  id = 1 + nhr;
  NumPart = N_init_gas;

  for(blocks = 0; blocks < (nhr_blocks + nlr_blocks); blocks++)
    {
      if(ThisTask == 0)
	{
	  SKIP;
	  fread(&npart, sizeof(int), 1, fd);
	  SKIP;
	  n_in_file = npart;

	  if(blocks < nhr_blocks)
	    type = 1;
	  else
	    type = 2;

	  /* printf("block=%d n_in_file=%d type=%d\n",blocks,n_in_file,type); */
	}

      MPI_Bcast(&n_in_file, 1, MPI_INT, 0, MPI_COMM_WORLD);	/* receive type of particle and total number */
      MPI_Bcast(&type, 1, MPI_INT, 0, MPI_COMM_WORLD);

      n_for_this_task = n_in_file / NTask;
      if(ThisTask < (n_in_file % NTask))	/* simple division to get the correct number */
	n_for_this_task++;


      for(subblock = 0; subblock < 3; subblock++)
	{
	  /* subblock 0 => pos
	     subblock 1 => vel
	     subblock 2 => mass */

	  if(type == 1 && subblock == 2)
	    continue;		/*  HR part's have no mass array  */

	  if(ThisTask == 0)
	    {
	      SKIP;
	    }

	  pc_here = pc;

	  for(pr = 0; pr < NTask; pr++)	/* go through all processes, note: pr is the receiving process */
	    {
	      if(ThisTask == 0 || ThisTask == pr)
		{
		  n = n_for_this_task;	/* number of particles for this process */

		  if(ThisTask == 0 && pr > 0)	/* pre-communication */
		    MPI_Recv(&n, 1, MPI_INT, pr, TAG_N, MPI_COMM_WORLD, &status);

		  if(ThisTask == pr && pr > 0)
		    MPI_Ssend(&n, 1, MPI_INT, 0, TAG_N, MPI_COMM_WORLD);

		  left = n;

		  while(left > 0)
		    {
		      if(left > BLOCKSIZE)	/* restrict transmission size to buffer length */
			groupsize = BLOCKSIZE;
		      else
			groupsize = left;

		      if(ThisTask == 0)
			{
			  if(subblock < 2)
			    {
			      fread(&dummy[0][0], sizeof(float), 3 * groupsize, fd);
			    }
			  else
			    {
			      if(blocks > 0)
				{
				  fread(&dummy[0][0], sizeof(float), groupsize, fd);
				}
			    }
			}

		      if(ThisTask == 0 && pr != 0)
			{
			  if(subblock < 2)
			    {
			      MPI_Ssend(&dummy[0][0], 3 * groupsize, MPI_FLOAT, pr, TAG_PDATA,
					MPI_COMM_WORLD);
			    }
			  else
			    {
			      if(blocks > 0)
				{
				  MPI_Ssend(&dummy[0][0], groupsize, MPI_FLOAT, pr, TAG_PDATA,
					    MPI_COMM_WORLD);
				}
			    }
			}

		      if(ThisTask != 0 && pr != 0)
			{
			  if(subblock < 2)
			    {
			      MPI_Recv(&dummy[0][0], 3 * groupsize, MPI_FLOAT, 0, TAG_PDATA, MPI_COMM_WORLD,
				       &status);
			    }
			  else
			    {
			      if(blocks > 0)
				{
				  MPI_Recv(&dummy[0][0], groupsize, MPI_FLOAT, 0, TAG_PDATA, MPI_COMM_WORLD,
					   &status);
				}
			    }
			}

		      if(ThisTask == pr)
			{
			  for(i = 0, fp = &dummy[0][0]; i < groupsize; i++)
			    {
			      if(subblock == 0)
				{
				  P[pc_here].Type = type;
				  P[pc_here].ID = id;	/* now set ID */
				  for(k = 0; k < 3; k++)
				    P[pc_here].Pos[k] = dummy[i][k];
				  pc_here++;
				  id++;
				}

			      if(subblock == 1)
				{
				  for(k = 0; k < 3; k++)
				    P[pc_here].Vel[k] = dummy[i][k];
				  pc_here++;
				}
			      if(subblock == 2)
				{
				  P[pc_here].Mass = fp[i];
				  pc_here++;
				}
			    }

			}

		      left -= groupsize;
		    }
		}

	      MPI_Barrier(MPI_COMM_WORLD);
	      MPI_Bcast(&id, 1, MPI_INT, pr, MPI_COMM_WORLD);
	    }

	  if(ThisTask == 0)
	    {
	      SKIP;
	    }
	}

      pc += n_for_this_task;
      NumPart += n_for_this_task;
    }

  if(ThisTask == 0)
    {
      fclose(fd);
    }

  MPI_Barrier(MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      printf("\nREADING DONE.\n\n");
      fflush(stdout);
    }

  /* now convert the units */

  sqr_a = sqrt(All.Time);

  counttype0 = counttype2 = counttype3 = 0;

  for(i = N_init_gas; i < NumPart; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  P[i].Pos[j] = P[i].Pos[j] * posfac;	/* here in units of kpc/h */

	  P[i].Vel[j] = P[i].Vel[j] * velfac;	/* comoving velocity xdot on km/sec */

	  P[i].Vel[j] *= sqr_a;	/* transform to velocity variable u */
	}

      if(P[i].Type == 1)
	{
	  P[i].Mass = All.MassTable[1];
	}

      if(P[i].Type == 2)
	{
	  P[i].Mass *= massfac;

	  counttype2++;

	  /*
	     r2=P[i].Pos[0]*P[i].Pos[0] + P[i].Pos[1]*P[i].Pos[1] + P[i].Pos[2]*P[i].Pos[2];
	     if(sqrt(r2) > 24000.0) 
	     {                    
	     P[i].Type=3;
	     counttype3++;
	     counttype2--;
	     }
	   */
	}
    }


  rho = 0.3 * 3 * 0.1 * 0.1 / (8 * M_PI * 43007.1);
  d = pow((All.MassTable[1] + All.MassTable[0]) / rho, 1.0 / 3);
  if(ThisTask == 0)
    printf("d= %g\n", d);

  for(i = 0; i < N_init_gas; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  P[i].Pos[j] = P[N_init_gas + i].Pos[j];
	  P[i].Vel[j] = P[N_init_gas + i].Vel[j];

	  P[i].Pos[j] += All.MassTable[1] / (All.MassTable[1] + All.MassTable[0]) * 0.1 * d;
	  P[N_init_gas + i].Pos[j] -= All.MassTable[0] / (All.MassTable[1] + All.MassTable[0]) * 0.1 * d;
	}
      P[i].Mass = All.MassTable[0];
      P[i].Type = 0;
      P[i].ID = P[N_init_gas + i].ID - nhr;
      counttype0++;
    }

  /* assign the initial internal energy to SPH particles */

  TempTemp = (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.InitGasTemp;

  /* unit conversion */
  TempTemp *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

  for(i = 0; i < N_init_gas; i++)
    {
      SphP[i].Entropy = TempTemp;
      /* Note: the coversion to real entropy will be done in the function init(),
         after the densities have been computed */
    }

  N_gas = N_init_gas;

  MPI_Allreduce(&counttype0, &sumnumpart, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  All.TotN_gas = sumnumpart;


  MPI_Allreduce(&NumPart, &sumnumpart, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  All.TotNumPart = sumnumpart;


  fflush(stdout);

  if(ThisTask == 0)
    {
      printf("u (internal)=%g\n", TempTemp);

      printf("\n\n");
      printf("particles loaded: %d \n\n", (int) All.TotNumPart);
      printf("particles of type 0: %d   (gas)\n", (int) All.TotN_gas);

      printf("\n");
      printf("Collisionless particles   :  %d\n", (int) (All.TotNumPart - All.TotN_gas));
      printf("Baryonic particles        :  %d\n", (int) All.TotN_gas);
      printf("                             ---------\n");
      printf("Total number of particles :  %d\n\n", (int) All.TotNumPart);
    }

#undef BLOCKSIZE
}
