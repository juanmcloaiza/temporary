#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"
//Added by JM:

void accrete_particles(void){
  float my_r, my_x, my_y, my_z;
  float my_v, my_vrad, my_vtan, my_vx, my_vy, my_vz;
  float my_angmom, my_energy;
  float my_angmomacc = sqrt(All.G * my_Mbh * my_racc);
  int   i, my_numpartate , n ;
  struct particle_data * Ptemp;
  struct sph_particle_data * SPHtemp;
  int * temp;
  my_numpartate = 0;
  flag_accretion = 0;
  for(i = 0; i < NumPart; i++){
    if(P[i].Ti_endstep != All.Ti_Current) continue;
      /*positions*/
    my_x = P[i].Pos[0];
    my_y = P[i].Pos[1];
    my_z = P[i].Pos[2];
    my_r = sqrt( my_x*my_x + my_y*my_y + my_z*my_z ) ;
     /*velocities*/
    my_vx = P[i].Vel[0];
    my_vy = P[i].Vel[1];
    my_vz = P[i].Vel[2];
    my_v = sqrt( my_vx*my_vx + my_vy*my_vy + my_vz*my_vz ) ;
    my_vrad = ( my_x*my_vx + my_y*my_vy + my_z*my_vz )/my_r ;
    my_vtan = sqrt( my_v * my_v - my_vrad * my_vrad );
     /*angular momentum*/
    my_angmom = my_r * my_vtan;
    my_energy = 0.5 * my_v*my_v - All.G * my_Mbh / my_r; 
    /*Flag particle for accretion*/
    if( (my_r < 0.1*my_racc) || ( (my_r < my_racc) && (my_angmom < my_angmomacc) && (my_energy < 0.0) ) ){
      my_numpartate++;
      P[i].Mass = 0;
    }  
  }
  /*Call this function even if in this process there was no accretion*/
  rearrange_particle_sequence();

  if(my_numpartate>0){
    printf("%d particles accreted in process %d\n", my_numpartate, ThisTask);
    printf("NumPart in process %d now: %d \n",ThisTask,NumPart );
    printf("N_gas in process %d now: %d \n",ThisTask,N_gas );
    printf("TotNumPart now: %u \n",(int)All.TotNumPart );
    printf("TotN_gas now: %u \n",(int)All.TotN_gas );
    fflush(stdout);
  }
  
    /*Update number of particles in all processes*/
/*
    temp = malloc(NTask * sizeof(int)); 
    All.TotNumPart = 0; 
    All.TotN_gas = 0; 
  
    MPI_Allgather(&NumPart, 1, MPI_INT, temp, 1, MPI_INT, MPI_COMM_WORLD); 
    for(i = 0; i< NTask; i++) { 
      All.TotNumPart += temp[i]; 
    } 
  
    for(i=0; i<NTask; i++) temp[i] = 0; 
    MPI_Allgather(&N_gas, 1, MPI_INT, temp, 1, MPI_INT, MPI_COMM_WORLD); 
  
    for(i = 0; i< NTask; i++) { 
      All.TotN_gas += temp[i]; 
    } 
    free(temp);
  
    //Print the updated value
    if(my_numpartate > 0){
      printf("TotNumPart after accretion: %u \n",(int)All.TotNumPart );
      printf("TotN_gas after accretion: %u \n",(int)All.TotN_gas );
    }
*/

  if(ThisTask==0){
    fprintf(FdAccretion, "%.4E \t %u \t %.4E \n",All.Time, All.TotNumPart, All.TotNumPart * All.MassTable[0]);
    fflush(FdAccretion);
  }
} 
//END of Added by JM

void rearrange_particle_sequence(void)
{
  int i, j;
  struct particle_data psave;

#if defined(BLACK_HOLES) || defined(MYSWITCH)
  int count_elim, count_gaselim, tot_elim, tot_gaselim;
  int * temp;
#endif

/*#ifdef SFR
  if(Stars_converted)
    {
      N_gas -= Stars_converted;
      Stars_converted = 0;

      for(i = 0; i < N_gas; i++)
    if(P[i].Type != 0)
      {
        for(j = N_gas; j < NumPart; j++)
          if(P[j].Type == 0)
        break;

        if(j >= NumPart)
          endrun(181170);

        psave = P[i];
        P[i] = P[j];
        SphP[i] = SphP[j];
        P[j] = psave;
      }
    }
#endif*/

#if defined(BLACK_HOLES) || defined(MYSWITCH)
  count_elim = 0;
  count_gaselim = 0;

  for(i = 0; i < NumPart; i++)
    if(P[i].Mass == 0){
      if(P[i].Type == 0){
        P[i] = P[N_gas - 1];
        SphP[i] = SphP[N_gas - 1];
        P[N_gas - 1] = P[NumPart - 1];
        N_gas--;
        count_gaselim++;
      }
      else{
        P[i] = P[NumPart - 1];
      }
      NumPart--;
      i--;
      count_elim++;
    }

  MPI_Allreduce(&count_elim, &tot_elim, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&count_gaselim, &tot_gaselim, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(ThisTask == 0){
      printf("Accrete-particles: Eliminated %d gas particles and merged away %d black holes.\n",
         tot_gaselim, tot_elim - tot_gaselim);
      fflush(stdout);
  }

    temp = malloc(NTask * sizeof(int)); 
    All.TotNumPart = 0; 
    All.TotN_gas = 0; 
  
    MPI_Allgather(&NumPart, 1, MPI_INT, temp, 1, MPI_INT, MPI_COMM_WORLD); 
    for(i = 0; i< NTask; i++) { 
      All.TotNumPart += temp[i]; 
    } 
  
    for(i=0; i<NTask; i++) temp[i] = 0; 
    MPI_Allgather(&N_gas, 1, MPI_INT, temp, 1, MPI_INT, MPI_COMM_WORLD); 
  
    for(i = 0; i< NTask; i++) { 
      All.TotN_gas += temp[i]; 
    } 
    free(temp);


#endif

}
