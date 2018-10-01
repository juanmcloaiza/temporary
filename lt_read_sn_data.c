#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"


/*

  This function read in the table used by metal pollution

*/

#ifdef LT_STELLAREVOLUTION

int *Zbins_dim, *Mbins_dim;
int *account_just_produced;
double *Zbins, *Mbins;
double *DataSpace, *Yields[LT_NMet];

#define LL 1500

void read_Yields(FILE * file)
{
  char s[LL], name[20];
  int i, j, k;

  /*
   * read in Sn Ia yields
   */

  /* read the number of metal bins (1 for metal independent yields) */
  do
    fgets(s, LL, file);
  while(strchr("%# \n", s[0]) != 0x0);
  /* allocate space for metal bins and read/store them */
  *Zbins_dim = atoi(s);
  Zbins = (double *) calloc(*Zbins_dim, sizeof(double));
  if(*Zbins_dim > 1)
    {
      do
	fgets(s, LL, file);
      while(strchr("%# \n", s[0]) != 0x0);
      for(j = 0; j < *Zbins_dim; j++)
	sscanf(&s[0], "%lg%[^n]s", &Zbins[j], &s[0]);
    }
  else
    Zbins[0] = 0;
  /* allocate space for mass bins and read/store them */
  do
    fgets(s, LL, file);
  while(strchr("%# \n", s[0]) != 0x0);
  *Mbins_dim = atoi(s);
  Mbins = (double *) calloc(*Mbins_dim, sizeof(double));
  if(*Mbins_dim > 1)
    {
      do
	fgets(s, LL, file);
      while(strchr("%# \n", s[0]) != 0x0);
      for(j = 0; j < *Mbins_dim; j++)
	sscanf(&s[0], "%lg%[^n]s", &Mbins[j], &s[0]);
    }
  else
    Mbins[0] = 0;

  /*
     do
     fgets(s, LL, file);
     while(strchr("%# \n", s[0]) != 0x0);
     *account_just_produced = atoi(s);
   */

  DataSpace = (double *) calloc(LT_NMet * *Zbins_dim * *Mbins_dim, sizeof(double));

  Yields[0] = &DataSpace[0];
  for(j = 1; j < LT_NMet; j++)
    Yields[j] = Yields[j - 1] + *Zbins_dim * *Mbins_dim;

  /* actually read yields. they are organized in subsequent blocks, one for each
     metal bin. each block is a table, whose rows refer to a single element and
     columns to the mass array */
  j = 0;
  while(j < *Zbins_dim)
    {
      do
	fgets(s, LL, file);
      while(!feof(file) && (strchr("%# \n", s[0]) != 0x0));

      while(!feof(file) && (strchr("%# \n", s[0]) == 0x0))
	{
	  /* find the element name */
	  sscanf(s, "%s %[^\n]s", name, &s[0]);

	  for(k = 0; k < LT_NMet; k++)
	    if(strcmp(name, MetNames[k]) == 0)
	      break;

	  /* this is an element you want to use (as specified in metals.dat) */
	  if(k < LT_NMet)
	    for(i = 0; i < *Mbins_dim; i++)
	      sscanf(s, "%lg%[^\n]s", &Yields[k][*Mbins_dim * j + i], &s[0]);

	  fgets(s, LL, file);
	}
      j++;
    }
  return;
}


#ifdef LT_SNIa
void read_SnIa_yields(void)
{
  char buff[300];
  int set, i, j;
  FILE *file;

  /*
   * read in Sn Ia yields
   */

  IaZbins_dim = (int *) mymalloc(All.Ia_Nset_ofYields * sizeof(int));
  IaZbins = (double **) mymalloc(All.Ia_Nset_ofYields * sizeof(double *));
  IaMbins_dim = (int *) mymalloc(All.Ia_Nset_ofYields * sizeof(int));
  IaMbins = (double **) mymalloc(All.Ia_Nset_ofYields * sizeof(double *));
  SnIaYields = (double ***) mymalloc(All.Ia_Nset_ofYields * sizeof(double **));
  /*  SnIaY_Give_Produced = (int*)mymalloc(All.Ia_Nset_ofYields * sizeof(int)); */
  for(set = 0; set < All.Ia_Nset_ofYields; set++)
    SnIaYields[set] = (double **) mymalloc(LT_NMet * sizeof(double *));

  for(set = 0; set < All.Ia_Nset_ofYields; set++)
    {
      if(ThisTask == 0)
	{
	  if(All.Ia_Nset_ofYields > 1)
	    sprintf(buff, "%s.%03d", All.SnIaDataFile, set);
	  else
	    strcpy(buff, All.SnIaDataFile);
	  if((file = fopen(buff, "r")) == NULL)
	    {
	      printf("I can't open SnIa data input file: <%s>\n", buff);
	      MPI_Finalize();
	      exit(0);
	    }
	  else
	    {
	      Zbins_dim = &IaZbins_dim[set];
	      Mbins_dim = &IaMbins_dim[set];
	      /*account_just_produced = &SnIaY_Give_Produced[set]; */
	      read_Yields(file);
	      fclose(file);

	      IaZbins[set] = Zbins;
	      IaMbins[set] = Mbins;
	      for(i = 0; i < LT_NMet; i++)
		SnIaYields[set][i] = Yields[i];
	    }
	}

      MPI_Bcast(&IaZbins_dim[set], 1, MPI_INT, 0, MPI_COMM_WORLD);
      if(ThisTask != 0)
	IaZbins[set] = (double *) calloc(IaZbins_dim[set], sizeof(double));
      MPI_Bcast(IaZbins[set], IaZbins_dim[set] * sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);

      MPI_Bcast(&IaMbins_dim[set], 1, MPI_INT, 0, MPI_COMM_WORLD);
      if(ThisTask != 0)
	IaMbins[set] = (double *) calloc(IaMbins_dim[set], sizeof(double));
      MPI_Bcast(IaMbins[set], IaMbins_dim[set] * sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);

      if(ThisTask != 0)
	{
	  DataSpace = (double *) calloc(LT_NMet * IaZbins_dim[set] * IaMbins_dim[set], sizeof(double));

	  SnIaYields[set][0] = DataSpace;
	  for(j = 1; j < LT_NMet; j++)
	    SnIaYields[set][j] = SnIaYields[set][j - 1] + IaZbins_dim[set] * IaMbins_dim[set];
	}
      MPI_Bcast(&DataSpace[0], LT_NMet * IaZbins_dim[set] * IaMbins_dim[set] * sizeof(double),
		MPI_BYTE, 0, MPI_COMM_WORLD);

    }

  return;
}
#endif


#ifdef LT_SNII
void read_SnII_yields(void)
{
  char buff[300];
  int set, i, j, k;
  FILE *file;

  /*
   * read in Sn II yields
   */

  IIZbins_dim = (int *) mymalloc(All.II_Nset_ofYields * sizeof(int));
  IIZbins = (double **) mymalloc(All.II_Nset_ofYields * sizeof(double *));
  IIMbins_dim = (int *) mymalloc(All.II_Nset_ofYields * sizeof(int));
  IIMbins = (double **) mymalloc(All.II_Nset_ofYields * sizeof(double *));
  SnIIYields = (double ***) mymalloc(All.II_Nset_ofYields * sizeof(double **));
  /* SnIIY_Give_Produced = (int*)mymalloc(All.Ia_Nset_ofYields * sizeof(int)); */
  for(set = 0; set < All.II_Nset_ofYields; set++)
    SnIIYields[set] = (double **) mymalloc(LT_NMet * sizeof(double *));
  SnIIEj = (double **) calloc(All.II_Nset_ofYields, sizeof(double *));


  for(set = 0; set < All.II_Nset_ofYields; set++)
    {
      if(ThisTask == 0)
	{
	  if(All.II_Nset_ofYields > 1)
	    sprintf(buff, "%s.%03d", All.SnIIDataFile, set);
	  else
	    strcpy(buff, All.SnIIDataFile);
	  if((file = fopen(buff, "r")) == NULL)
	    {
	      printf("I can't open SnII data input file: <%s>\n", buff);
	      MPI_Finalize();
	      exit(0);
	    }
	  else
	    {
	      Zbins_dim = &IIZbins_dim[set];
	      Mbins_dim = &IIMbins_dim[set];
	      /*account_just_produced = &SnIIY_Give_Produced[set]; */
	      read_Yields(file);
	      fclose(file);

	      IIZbins[set] = Zbins;
	      IIMbins[set] = Mbins;
	      for(i = 0; i < LT_NMet; i++)
		SnIIYields[set][i] = Yields[i];
	    }
	}

      MPI_Bcast(&IIZbins_dim[set], 1, MPI_INT, 0, MPI_COMM_WORLD);
      if(ThisTask != 0)
	IIZbins[set] = (double *) calloc(IIZbins_dim[set], sizeof(double));
      MPI_Bcast(&IIZbins[set][0], IIZbins_dim[set] * sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);

      MPI_Bcast(&IIMbins_dim[set], 1, MPI_INT, 0, MPI_COMM_WORLD);
      if(ThisTask != 0)
	IIMbins[set] = (double *) calloc(IIMbins_dim[set], sizeof(double));
      MPI_Bcast(&IIMbins[set][0], IIMbins_dim[set] * sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);

      if(ThisTask != 0)
	{
	  DataSpace = (double *) calloc(LT_NMet * IIZbins_dim[set] * IIMbins_dim[set], sizeof(double));

	  SnIIYields[set][0] = DataSpace;
	  for(j = 1; j < LT_NMet; j++)
	    SnIIYields[set][j] = SnIIYields[set][j - 1] + IIZbins_dim[set] * IIMbins_dim[set];
	}
      MPI_Bcast(&DataSpace[0], LT_NMet * IIZbins_dim[set] * IIMbins_dim[set] * sizeof(double),
		MPI_BYTE, 0, MPI_COMM_WORLD);

      MPI_Barrier(MPI_COMM_WORLD);

      /* SnIIEj[set] will contain the total ejected mass in all element present in file for each couple
         (Zbin,Mbin) */
      SnIIEj[set] = (double *) calloc(IIZbins_dim[set] * IIMbins_dim[set], sizeof(double));
      memcpy((void *) SnIIEj[set], (void *) SnIIYields[set][FillEl],
	     (size_t) (IIZbins_dim[set] * IIMbins_dim[set] * sizeof(double)));

      /* the last element FillEl will contain the difference between the total ejected mass and the sum
         of ejecta from the used elements (which can be less than those present in the file) */
      for(i = 0; i < IIZbins_dim[set]; i++)
	for(j = 0; j < IIMbins_dim[set]; j++)
	  {
	    SnIIYields[set][FillEl][i * IIMbins_dim[set] + j] = 0;
	    for(k = 0; k < LT_NMet; k++)
	      if(k != FillEl)
		SnIIYields[set][FillEl][i * IIMbins_dim[set] + j] +=
		  SnIIYields[set][k][i * IIMbins_dim[set] + j];
	    if(fabs
	       ((SnIIEj[set][i * IIMbins_dim[set] + j] -
		 SnIIYields[set][FillEl][i * IIMbins_dim[set] + j]) / SnIIEj[set][i * IIMbins_dim[set] + j]) <
	       1e-2)
	      SnIIYields[set][FillEl][i * IIMbins_dim[set] + j] = 0;
	    else
	      SnIIYields[set][FillEl][i * IIMbins_dim[set] + j] =
		SnIIEj[set][i * IIMbins_dim[set] + j] - SnIIYields[set][FillEl][i * IIMbins_dim[set] + j];
	  }
    }
  return;
}
#endif


#ifdef LT_AGB
void read_AGB_yields(void)
{
  char buff[300];
  int set, i, j, k;
  FILE *file;


  /*
   * read in Sn AGB yields
   */

  AGBZbins_dim = (int *) mymalloc(All.AGB_Nset_ofYields * sizeof(int));
  AGBZbins = (double **) mymalloc(All.AGB_Nset_ofYields * sizeof(double *));
  AGBMbins_dim = (int *) mymalloc(All.AGB_Nset_ofYields * sizeof(int));
  AGBMbins = (double **) mymalloc(All.AGB_Nset_ofYields * sizeof(double *));
  AGBYields = (double ***) mymalloc(All.AGB_Nset_ofYields * sizeof(double **));
  /* AGBY_Give_Produced = (int*)mymalloc(All.Ia_Nset_ofYields * sizeof(int)); */
  for(set = 0; set < All.AGB_Nset_ofYields; set++)
    AGBYields[set] = (double **) mymalloc(LT_NMet * sizeof(double *));
  AGBEj = (double **) calloc(All.AGB_Nset_ofYields, sizeof(double *));

  for(set = 0; set < All.AGB_Nset_ofYields; set++)
    {
      if(ThisTask == 0)
	{
	  if(All.AGB_Nset_ofYields > 1)
	    sprintf(buff, "%s.%03d", All.AGBDataFile, set);
	  else
	    strcpy(buff, All.AGBDataFile);
	  if((file = fopen(buff, "r")) == NULL)
	    {
	      printf("I can't open AGB data input file: <%s>\n", buff);
	      MPI_Finalize();
	      exit(1);
	    }
	  else
	    {
	      Zbins_dim = &AGBZbins_dim[set];
	      Mbins_dim = &AGBMbins_dim[set];
	      /*account_just_produced = &AGBY_Give_Produced[set]; */
	      read_Yields(file);
	      fclose(file);

	      AGBZbins[set] = Zbins;
	      AGBMbins[set] = Mbins;
	      for(i = 0; i < LT_NMet; i++)
		AGBYields[set][i] = Yields[i];
	    }
	}

      MPI_Bcast(&AGBZbins_dim[set], 1, MPI_INT, 0, MPI_COMM_WORLD);
      if(ThisTask != 0)
	AGBZbins[set] = (double *) calloc(AGBZbins_dim[set], sizeof(double));
      MPI_Bcast(&AGBZbins[set][0], AGBZbins_dim[set] * sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);

      MPI_Bcast(&AGBMbins_dim[set], 1, MPI_INT, 0, MPI_COMM_WORLD);
      if(ThisTask != 0)
	AGBMbins[set] = (double *) calloc(AGBMbins_dim[set], sizeof(double));
      MPI_Bcast(&AGBMbins[set][0], AGBMbins_dim[set] * sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);

      if(ThisTask != 0)
	{
	  DataSpace = (double *) calloc(LT_NMet * AGBZbins_dim[set] * AGBMbins_dim[set], sizeof(double));
	  AGBYields[set][0] = DataSpace;
	  for(j = 1; j < LT_NMet; j++)
	    AGBYields[set][j] = AGBYields[set][j - 1] + AGBZbins_dim[set] * AGBMbins_dim[set];
	}
      MPI_Bcast(&DataSpace[0], LT_NMet * AGBZbins_dim[set] * AGBMbins_dim[set] * sizeof(double),
		MPI_BYTE, 0, MPI_COMM_WORLD);

      /* AGBEj[set] will contain the total ejected mass in all element present in file for each couple
         (Zbin,Mbin) */
      AGBEj[set] = (double *) calloc(AGBZbins_dim[set] * AGBMbins_dim[set], sizeof(double));
      memcpy((void *) AGBEj[set], (void *) AGBYields[set][FillEl],
	     (size_t) (AGBZbins_dim[set] * AGBMbins_dim[set] * sizeof(double)));

      /* the last element FillEl will contain the difference between the total ejected mass and the sum
         of ejecta from the used elements (which can be less than those present in the file) */
      for(i = 0; i < AGBZbins_dim[set]; i++)
	for(j = 0; j < AGBMbins_dim[set]; j++)
	  {
	    AGBYields[set][FillEl][i * AGBMbins_dim[set] + j] = 0;
	    for(k = 0; k < LT_NMet; k++)
	      if(k != FillEl)
		AGBYields[set][FillEl][i * AGBMbins_dim[set] + j] +=
		  AGBYields[set][k][i * AGBMbins_dim[set] + j];
	    if(fabs
	       ((AGBEj[set][i * AGBMbins_dim[set] + j] -
		 AGBYields[set][FillEl][i * AGBMbins_dim[set] + j]) / AGBEj[set][i * AGBMbins_dim[set] + j]) <
	       1e-2)
	      AGBYields[set][FillEl][i * AGBMbins_dim[set] + j] = 0;
	    else
	      AGBYields[set][FillEl][i * AGBMbins_dim[set] + j] =
		AGBEj[set][i * AGBMbins_dim[set] + j] - AGBYields[set][FillEl][i * AGBMbins_dim[set] + j];

	  }
    }
  return;
}
#endif


#endif
