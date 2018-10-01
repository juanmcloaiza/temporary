#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"


/*
 * This function read in the cooling tables
 *
 * !!! to be rewritten and generalized !!!
 */

#ifdef LT_METAL_COOLING

void read_metalcool_table(void)
{
  FILE *file;
  char *fnames[8] = { "suth8", "suth7", "suth6", "suth5",
    "suth4", "suth3", "suth2", "suth1"
  };
  int i, j;
  double m, t, Ln, L;

  Zvalue[0] = -4.0;
  Zvalue[1] = -3.0;
  Zvalue[2] = -2.0;
  Zvalue[3] = -1.5;
  Zvalue[4] = -1.0;
  Zvalue[5] = -0.5;
  Zvalue[6] = 0.0;
  Zvalue[7] = 0.5;

  if(ThisTask == 0)
    {
      printf("\nreading Sutherland & Dopita tables");
      fflush(stdout);

      for(i = 0; i < 8; i++)
	{
	  if((file = fopen(fnames[i], "r")) == NULL)
	    {
	      printf("\ni can't open Suth. & Dop input file:" "%s\n", fnames[i]);
	      MPI_Finalize();
	      exit(404040);
	    }
	  else
	    {
	      printf(".. %d", i);
	      for(j = 0; j < Zlength; j++)
		{
		  fscanf(file, "%lf %lf %lf %lf\n", &m, &t, &Ln, &L);
		  Lsuthdop[i][j] = Ln;
		  if(i == 0)
		    ZTemp[j] = t;
		}
	    }
	  fclose(file);
	}
      printf("\n");
    }


  MPI_Bcast((void *) &Lsuthdop[0][0], ZBins * Zlength * sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast((void *) &ZTemp[0], Zlength * sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);


  if(ThisTask == 0)
    printf("\n\n");
}

#endif
