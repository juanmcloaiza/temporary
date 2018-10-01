#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"

#define r_acc 0.01
static int Nsink_here, Nsink_all;
struct my_sink
{
  float Pos[3];
  float Vel[3];
  float Mass;
}
*Sink_here,
*Sink_all;

void get_all_sink_particles(void);
void get_local_sink_particles(void);

void accrete_particles(void){
  float r_sink, x_sink, y_sink, z_sink, m_sink;
  float v_sink, v_rad_sink, v_tan_sink, vx_sink, vy_sink, vz_sink;
  float r_rel, x_rel, y_rel, z_rel;
  float v_rel, v_rad_rel, v_tan_rel, vx_rel, vy_rel, vz_rel;
  float my_angmom, my_energy;
  float my_angmomacc;
  float accretion_threshold = N0;
  int   i, my_numpartate, n, k;
  struct particle_data * Ptemp;
  struct sph_particle_data * SPHtemp;
  int * temp;

  /*First of all, get info about all the sink particles*/
  Nsink_here = NtypeLocal[5];
  Nsink_all = Ntype[5];
  get_local_sink_particles();
  get_all_sink_particles();
  /*Set counters to zero*/
  my_numpartate = 0;
  flag_accretion_threshold = 0;
  //printf("Task %d: \n N_gas = %d, NumPart = %d, NtypeLocal[0] = %d, NtypeLocal[3] = %d \n",ThisTask,N_gas,NumPart,NtypeLocal[0],NtypeLocal[4]);
  //printf("Task %d: \n All.TotN_gas = %d, All.TotNumPart = %d, Ntype[0] = %d, Ntype[3] = %d \n",ThisTask,All.TotN_gas,All.TotNumPart,Ntype[0],Ntype[4]);
  //Sink_particle stuff:
  for( k = 0; k < Nsink_all; k++ ){
    m_sink = Sink_all[k].Mass;
    x_sink = Sink_all[k].Pos[0];
    y_sink = Sink_all[k].Pos[1];
    z_sink = Sink_all[k].Pos[2];
     
    vx_sink = Sink_all[k].Vel[0];
    vy_sink = Sink_all[k].Vel[1];
    vz_sink = Sink_all[k].Vel[2];
    my_angmomacc = sqrt(All.G * m_sink * r_acc);
    /******************************/
    for(i = 0; i < N_gas; i++){
      if( P[i].ID > 0 ){
        /*positions*/
        x_rel = P[i].Pos[0] - x_sink;
        y_rel = P[i].Pos[1] - y_sink;
        z_rel = P[i].Pos[2] - z_sink;
        r_rel = sqrt( x_rel*x_rel + y_rel*y_rel + z_rel*z_rel ) ;
        /*velocities*/
        vx_rel = P[i].Vel[0] - vx_sink;
        vy_rel = P[i].Vel[1] - vy_sink;;
        vz_rel = P[i].Vel[2] - vz_sink;;
        v_rel = sqrt( vx_rel*vx_rel + vy_rel*vy_rel + vz_rel*vz_rel ) ;
        v_rad_rel = ( x_rel*vx_rel + y_rel*vy_rel + z_rel*vz_rel )/r_rel ;
        v_tan_rel = sqrt( v_rel * v_rel - v_rad_rel * v_rad_rel );
        /*angular momentum*/
        my_angmom = r_rel * v_tan_rel;
        my_energy = 0.5 * v_rel*v_rel - All.G * m_sink / r_rel; 
        /*Flag particle for accretion*/
        if( (r_rel < 0.1*r_acc) || ( (r_rel < r_acc) && (my_angmom < my_angmomacc) && (my_energy < 0.0) ) ){
          if( P[i].ID == 2 ){
              printf("P[i].ID = 2\n");
              exit(EXIT_FAILURE);
          }
          my_numpartate++;
          P[i].ID = 0;
        }
      }
    }
  }
  /*Update the particle list in the current processor*/
  if(my_numpartate > 0){
    Ptemp = malloc(All.MaxPart * sizeof(struct particle_data));
    if( Ptemp == 0){
      printf("Error: Out of memory for malloc Ptemp\n");
      exit(EXIT_FAILURE);
    }
    SPHtemp = malloc(All.MaxPartSph*sizeof(struct sph_particle_data));
    if( SPHtemp == 0){
      printf("Error: Out of memory for malloc SPHtemp\n");
      exit(EXIT_FAILURE);
    }
    for(n = 0 ,i = 0; i < NumPart; i++){
      if(P[i].ID != 0){
        Ptemp[n] = P[i];
        if(i < N_gas) SPHtemp[n] = SphP[i];
        n+=1;
      } 
    }
    free(P);
    P = Ptemp;
    free(SphP);
    SphP = SPHtemp;
    NumPart = NumPart - my_numpartate;
    N_gas = N_gas - my_numpartate;
    
    /*If the Number of particles is less than 90% of the original particles present,
    a domain decomposition is done*/
    if( (float)NumPart < accretion_threshold) flag_accretion_threshold = 1;
    printf("%d particles accreted in process %d\n", my_numpartate, ThisTask);
    printf("Numpart/N0 = %f\n", (float)NumPart/(float)N0 );
    fflush(stdout);
  }

  /*Update number of particles in all processes*/
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

  if(ThisTask==0){
    fprintf(FdAccretion, "%.4E \t %u \t %.4E \n",All.Time, All.TotNumPart, All.TotNumPart * All.MassTable[0]);
    fflush(FdAccretion);
  }
} 

void get_local_sink_particles(void){
  int i,j;
  Sink_here = malloc(Nsink_here*sizeof(struct my_sink));
  for(i = 0 ; i < Nsink_here ; i++){
      for( j = 0 ; j < 3 ; j ++)
        Sink_here[i].Pos[j] = P[N_gas + i].Pos[j];
      for( j = 0 ; j < 3 ; j ++)
        Sink_here[i].Vel[j] = P[N_gas + i].Vel[j];
      Sink_here[i].Mass = P[N_gas + i].Mass;
  }
}

void get_all_sink_particles(void){
  int *N_in_task, *offsets;
  int *Nb_in_task, *offsetsb;
  int i, s; 
  /*Send-receive number of particles from each processors to all processors*/
  N_in_task = (int*)malloc(NTask*sizeof(int));
  Nb_in_task = (int*)malloc(NTask*sizeof(int));
  offsets = (int*)malloc(NTask*sizeof(int));
  offsetsb = (int*)malloc(NTask*sizeof(int));
  /*N_in_task[I] is the number of sink particles in processor I*/
  MPI_Allgather(&Nsink_here,1,MPI_INT,N_in_task,1,MPI_INT,MPI_COMM_WORLD);
  s=0;
  for (i=0; i < NTask; i++){
      offsets[i] = s; 
      s += N_in_task[i];
  }
  if(s != Nsink_all){
	  printf("\nFAILURE in getting tottal number of sink particles\n");
	  exit(EXIT_FAILURE);
  }

  for (i = 0; i < NTask; i++){
      Nb_in_task[i] = N_in_task[i] * sizeof(struct my_sink);
      offsetsb[i] = offsets[i] * sizeof(struct my_sink);
  }
  Sink_all = malloc(Nsink_all*sizeof(struct my_sink));
  MPI_Allgatherv
      ((void*)Sink_here,sizeof(struct my_sink)*Nsink_here, MPI_BYTE,(void*)Sink_all,
       Nb_in_task,offsetsb,MPI_BYTE,MPI_COMM_WORLD);
  free(N_in_task); 
  free(Nb_in_task); 
  free(offsets);
  free(offsetsb);
}
