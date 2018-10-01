#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"


/*
 * This function read in the names of metals
 */

#ifdef LT_STELLAREVOLUTION

void read_metals(void)
{
  char s[200], name[200];
  int j;
  FILE *file;

  Hel = -1;
  Iron = -1;
  Oxygen = -1;
  FillEl = -1;

  if(ThisTask == 0)
    {
      if((file = fopen("metals.dat", "r")) == NULL)
	{
	  printf("I can't open metals data input file:" "%s\n", "metals.dat");
	  endrun(88888);
	}

      for(j = 0; j < LT_NMet; j++)
	{
	  if(feof(file))
	    {
	      printf("something wrong with <metals.dat>\n");
	      endrun(88889);
	    }
	  fgets(s, 200, file);
	  if((sscanf(s, "%[a-zA-Z]s %lg", &name[0], &MetSolarValues[j])) == 1)
	    printf("it seems that no solar abundance is present for element %s\n", name);

	  MetNames[j] = (char *) mymalloc(sizeof(char) * strlen(name));
	  strcpy(MetNames[j], name);
	  if(strcmp(name, "He") == 0)
	    Hel = j;
	  if(strcmp(name, "Fe") == 0)
	    Iron = j;
	  if(strcmp(name, "O") == 0)
	    Oxygen = j;
	  if(strcmp(name, "Ej") == 0)
	    FillEl = j;
	}

      Hyd = LT_NMet - 1;
      MetNames[Hyd] = (char *) mymalloc(sizeof(char));
      strcpy(MetNames[Hyd], "H");

      if(FillEl == -1)
	{
	  printf("you don't have specified the FillEl position.. better to do\n");
	  endrun(10001000);
	}

      fclose(file);

#ifdef LT_METAL_COOLING
      if(Iron == -1)
	{
	  if(Oxygen >= 0)
	    printf("you don't trace IRON. metal cooling will be calculated inferring X_Fe from X_O\n");
	  else
	    {
	      printf("you don't trace neither IRON nor OXYGEN. So far, metal cooling cannot be used\n");
	      endrun(993399);
	    }
	}
#endif

    }

  MPI_Bcast(&Iron, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Oxygen, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Hel, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&FillEl, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);

#ifdef LT_METAL_COOLING

/*   if(ThisTask == 0) */
/*     { */
/*       if((file = fopen("metals_for_cooling.dat", "r")) == NULL) */
/*	{ */
/*	  printf("I can't open metals data input file:" "%s\n", "metals_for_cooling.dat"); */
/*           endrun(88888); */
/*	} */

/*       fgets(s, 5, file); */
/*       NMet4C = atoi(s); */

/*       if(NMet4C > LT_NMet) */
/*         { */
/*           printf("oops.. to calculate the cooling function you need more elements than\n" */
/*                  "those present in metals.dat\n"); */
/*           endrun(88890); */
/*         } */

/*       Met4cNames = (char**)mymalloc(NMet4C * sizeof(char*)); */
/*       Z4Cool = (int*)mymalloc(NMet4C * sizeof(char*)); */
/*       for(j = 0; j < NMet4C; j++) */
/*	{ */
/*           if(feof(file)) */
/*             { */
/*               printf("something wrong with <metals_for_cooling.dat>\n"); */
/*               endrun(88891); */
/*             } */
/*	  fgets(s, 5, file); */
/*           sscanf(s, "%[^ \n]s ", &s[0]); */

/*	  Met4cNames[j] = (char *)mymalloc(sizeof(char) * strlen(s)); */
/*	  strcpy(Met4cNames[j], s); */

/*           for(i = 0; i < LT_NMet; i++) */
/*             if(strcmp(MetNames[i], Met4cNames[j]) == 0) */
/*               break; */
/*           if(i < LT_NMet) */
/*             Z4Cool[j] = i; */
/*           else */
/*             printf("oops, it seems that the element %3s needed by cooling is not present!\n", */
/*                    Met4cNames[j]); */
/*	} */
/*     } */

/*   MPI_Bcast(&NMet4C, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD); */
/*   if(ThisTask != 0) */
/*     Z4Cool = (int*)mymalloc(NMet4C * sizeof(char*)); */
/*   MPI_Bcast(&Z4Cool[0], NMet4C * sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD); */

#endif
  return;
}
#endif
