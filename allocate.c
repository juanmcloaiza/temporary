#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"



/* Allocates communication buffers, plus the list of the interior
 * domains (for neighbour search and hydro). 
 */
void allocate_commbuffers(void)
{
  size_t bytes;

  Exportflag = (char *) mymalloc(NTask * sizeof(char));
  DomainStartList = (int *) mymalloc(NTask * sizeof(int));
  DomainEndList = (int *) mymalloc(NTask * sizeof(int));

  TopNodes = (struct topnode_data *) mymalloc(MAXTOPNODES * sizeof(struct topnode_data));

  DomainWork = (double *) mymalloc(MAXTOPNODES * sizeof(double));
  DomainCount = (int *) mymalloc(MAXTOPNODES * sizeof(int));
  DomainCountSph = (int *) mymalloc(MAXTOPNODES * sizeof(int));
  DomainTask = (int *) mymalloc(MAXTOPNODES * sizeof(int));
  DomainNodeIndex = (int *) mymalloc(MAXTOPNODES * sizeof(int));
  DomainHmax = (FLOAT *) mymalloc(MAXTOPNODES * sizeof(FLOAT));
  DomainMoment = (struct DomainNODE *) mymalloc(MAXTOPNODES * sizeof(struct DomainNODE));
#ifdef LT_STELLAREVOLUTION
  DomainCountStar = (int *) mymalloc(MAXTOPNODES * sizeof(int));
#endif



  if(!(CommBuffer = mymalloc(bytes = All.BufferSize * 1024 * 1024)))
    {
      printf("failed to allocate memory for `CommBuffer' (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(2);
    }

  All.BunchSizeForce =
    (All.BufferSize * 1024 * 1024) / (sizeof(struct gravdata_index) + 2 * sizeof(struct gravdata_in));

  if(All.BunchSizeForce & 1)
    All.BunchSizeForce -= 1;	/* make sure that All.BunchSizeForce is even 
				   --> 8-byte alignment for 64bit processors */

  GravDataIndexTable = (struct gravdata_index *) CommBuffer;
  GravDataIn = (struct gravdata_in *) (GravDataIndexTable + +All.BunchSizeForce);
  GravDataGet = GravDataIn + All.BunchSizeForce;
  GravDataOut = GravDataIn;	/* this will overwrite the GravDataIn-Table */
  GravDataResult = GravDataGet;	/* this will overwrite the GravDataGet-Table */


  All.BunchSizeDensity =
    (All.BufferSize * 1024 * 1024) / (2 * sizeof(struct densdata_in) + 2 * sizeof(struct densdata_out));

  DensDataIn = (struct densdata_in *) CommBuffer;
  DensDataGet = DensDataIn + All.BunchSizeDensity;
  DensDataResult = (struct densdata_out *) (DensDataGet + All.BunchSizeDensity);
  DensDataPartialResult = DensDataResult + All.BunchSizeDensity;

#if defined(SFR_METALS) || defined (BG_SFR)
  All.BunchSizeMetal = (All.BufferSize * 1024 * 1024) / (2 * sizeof(struct metaldata_in));
  MetalDataIn = (struct metaldata_in *) CommBuffer;
  MetalDataGet = MetalDataIn + All.BunchSizeMetal;
#ifdef SFR_PROMOTION
  All.BunchSizeHotNgbs =
    (All.BufferSize * 1024 * 1024) / (2 * sizeof(struct hotngbs_in) + 2 * sizeof(struct hotngbs_out));

  HotNgbsIn = (struct hotngbs_in *) CommBuffer;
  HotNgbsGet = HotNgbsIn + All.BunchSizeHotNgbs;
  HotNgbsResult = (struct hotngbs_out *) (HotNgbsGet + All.BunchSizeHotNgbs);
  HotNgbsPartialResult = HotNgbsResult + All.BunchSizeHotNgbs;
#endif
#endif

#ifdef BLACK_HOLES
  All.BunchSizeBlackhole =
    (All.BufferSize * 1024 * 1024) / (2 * sizeof(struct blackholedata_in) +
				      2 * sizeof(struct blackholedata_out));
  BlackholeDataIn = (struct blackholedata_in *) CommBuffer;
  BlackholeDataGet = BlackholeDataIn + All.BunchSizeBlackhole;
  BlackholeDataResult = (struct blackholedata_out *) (BlackholeDataGet + All.BunchSizeBlackhole);
  BlackholeDataPartialResult = BlackholeDataResult + All.BunchSizeBlackhole;
#endif

#ifdef LT_STELLAREVOLUTION
  All.BunchSizeMetal_Ngb = (All.BufferSize * 1024 * 1024) /
    (2 * sizeof(struct metal_ngbfindingdata_in) + 2 * sizeof(struct metal_ngbfindingdata_out));

  MetalNgbIn = (struct metal_ngbfindingdata_in *) CommBuffer;
  MetalNgbGet = MetalNgbIn + All.BunchSizeMetal_Ngb;
  MetalNgbResultIn = (struct metal_ngbfindingdata_out *) (MetalNgbGet + All.BunchSizeMetal_Ngb);
  MetalNgbResultGet = MetalNgbResultIn + All.BunchSizeMetal_Ngb;

  All.BunchSizeMetal_Spread = (All.BufferSize * 1024 * 1024) / (2 * sizeof(struct metaldata_in));

  MetalDataIn = (struct metaldata_in *) CommBuffer;
  MetalDataGet = MetalDataIn + All.BunchSizeMetal_Spread;

#endif

  All.BunchSizeHydro =
    (All.BufferSize * 1024 * 1024) / (2 * sizeof(struct hydrodata_in) + 2 * sizeof(struct hydrodata_out));

  HydroDataIn = (struct hydrodata_in *) CommBuffer;
  HydroDataGet = HydroDataIn + All.BunchSizeHydro;
  HydroDataResult = (struct hydrodata_out *) (HydroDataGet + All.BunchSizeHydro);
  HydroDataPartialResult = HydroDataResult + All.BunchSizeHydro;

#ifndef LT_STELLAREVOLUTION
  All.BunchSizeDomain =
    (All.BufferSize * 1024 * 1024) / (sizeof(struct particle_data) + sizeof(struct sph_particle_data) +
				      sizeof(peanokey));
#else
  All.BunchSizeDomain =
    (All.BufferSize * 1024 * 1024) / (sizeof(struct particle_data) + sizeof(struct sph_particle_data) +
				      sizeof(peanokey) + sizeof(struct met_particle_data));
#endif
  if(All.BunchSizeDomain & 1)
    All.BunchSizeDomain -= 1;	/* make sure that All.BunchSizeDomain is even 
				   --> 8-byte alignment of DomainKeyBuf for 64bit processors */

  DomainPartBuf = (struct particle_data *) CommBuffer;
  DomainSphBuf = (struct sph_particle_data *) (DomainPartBuf + All.BunchSizeDomain);
  DomainKeyBuf = (peanokey *) (DomainSphBuf + All.BunchSizeDomain);
#ifdef LT_STELLAREVOLUTION
  DomainMetBuf = (struct met_particle_data *) (DomainKeyBuf + All.BunchSizeDomain);
#endif


#ifdef FOF
  All.BunchSizeFoF =
    (All.BufferSize * 1024 * 1024) / (2 * sizeof(struct fofdata_in) + 2 * sizeof(struct fofdata_out));

  FoFDataIn = (struct fofdata_in *) CommBuffer;
  FoFDataGet = FoFDataIn + All.BunchSizeFoF;
  FoFDataResult = (struct fofdata_out *) (FoFDataGet + All.BunchSizeFoF);
  FoFDataPartialResult = FoFDataResult + All.BunchSizeFoF;
#endif

#ifdef MHM
  All.BunchSizeKinetic = (All.BufferSize * 1024 * 1024) / (2 * sizeof(struct kindata_in));
  KinDataIn = (struct kindata_in *) CommBuffer;
  KinDataGet = KinDataIn + All.BunchSizeKinetic;
#endif


#ifdef BG_SFR
  All.BunchSizeSolidAngle =
    (All.BufferSize * 1024 * 1024) / (2 * sizeof(struct solidangledata_in) +
				      2 * sizeof(struct solidangledata_out));

  SolidAngleDataIn = (struct solidangledata_in *) CommBuffer;
  SolidAngleDataGet = SolidAngleDataIn + All.BunchSizeSolidAngle;
  SolidAngleDataResult = (struct solidangledata_out *) (SolidAngleDataGet + All.BunchSizeSolidAngle);
  SolidAngleDataPartialResult = SolidAngleDataResult + All.BunchSizeSolidAngle;
#endif



  if(ThisTask == 0)
    {
      printf("\nAllocated %d MByte communication buffer per processor.\n\n", All.BufferSize);
      printf("Communication buffer has room for %d particles in gravity computation\n", All.BunchSizeForce);
      printf("Communication buffer has room for %d particles in density computation\n", All.BunchSizeDensity);
      printf("Communication buffer has room for %d particles in hydro computation\n", All.BunchSizeHydro);
      printf("Communication buffer has room for %d particles in domain decomposition\n", All.BunchSizeDomain);
#ifdef BG_SFR
      printf("Communication buffer has room for %d particles in metal spreading\n", All.BunchSizeMetal);
#endif

#ifdef LT_STELLAREVOLUTION
      printf("Communication buffer has room for %d particles in metal ngb searching\n",
	     All.BunchSizeMetal_Ngb);
      printf("Communication buffer has room for %d particles in metal spreading\n",
	     All.BunchSizeMetal_Spread);
#endif
      printf("\n");
    }
}





/* This routine allocates memory for 
 * particle storage, both the collisionless and the SPH particles.
 * The memory for the ordered binary tree of the timeline
 * is also allocated.
 */
void allocate_memory(void)
{
  size_t bytes;
  double bytes_tot = 0;

  if(All.MaxPart > 0)
    {
      if(!(P = (struct particle_data *) mymalloc(bytes = All.MaxPart * sizeof(struct particle_data))))
	{
	  printf("failed to allocate memory for `P' (%g MB).\n", bytes / (1024.0 * 1024.0));
	  endrun(1);
	}
      bytes_tot += bytes;

      if(ThisTask == 0)
	printf("\nAllocated %g MByte for particle storage.\n\n", bytes_tot / (1024.0 * 1024.0));
    }

  if(All.MaxPartSph > 0)
    {
      bytes_tot = 0;

      if(!
	 (SphP =
	  (struct sph_particle_data *) mymalloc(bytes = All.MaxPartSph * sizeof(struct sph_particle_data))))
	{
	  printf("failed to allocate memory for `SphP' (%g MB).\n", bytes / (1024.0 * 1024.0));
	  endrun(1);
	}
      bytes_tot += bytes;

      if(ThisTask == 0)
	printf("Allocated %g MByte for storage of SPH data.\n\n", bytes_tot / (1024.0 * 1024.0));

#ifdef LT_STELLAREVOLUTION
      bytes_tot = 0;

      if(!
	 (MetP =
	  (struct met_particle_data *) mymalloc(bytes = All.MaxPartMet * sizeof(struct met_particle_data))))
	{
	  printf("failed to allocate memory for `MetP' (%g MB).\n", bytes / (1024.0 * 1024.0));
	  endrun(1);
	}
      bytes_tot += bytes;

      if(ThisTask == 0)
	printf("Allocated %g MByte for storage of metal data.\n\n", bytes_tot / (1024.0 * 1024.0));
#endif
    }
}




/* This routine frees the memory for the particle storage
 * and for the communication buffers,
 * but we don't actually call it in the code. 
 * When the program terminats, the memory will be automatically
 * freed by the operating system.
 */
void free_memory(void)
{
  if(All.MaxPart > 0)
    myfree(P);

  if(All.MaxPartSph > 0)
    myfree(SphP);

  myfree(CommBuffer);
}
