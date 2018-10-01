#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <mpi.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"

/*! \file fof.c
 *  \brief parallel FoF group finder
 */

#ifdef FOF


static struct id_list
{
  long long GrIndex;
  int ID;
  int Type;
  int Task;
  int Index;
  FLOAT Density;
  FLOAT Pos[3];
  FLOAT Mass;
  FLOAT Sfr;
#ifdef BLACK_HOLES
  FLOAT BH_Mass;
  FLOAT BH_Mdot;
#endif
#ifdef WINDS
  FLOAT DelayTime;
#endif
}
 *local_ids, *export_ids;

#ifdef BUBBLES

static struct biggest_group
{
  int Len;
  float CM[3];
  double Mass[6];
}
 *BiggestGroup;

#endif



static double LinkL;
static int Ngroups, TotNgroups, Nids;
static long long TotNids;
static int Nlocal;

static int *local_toGo, *toGo, *export_offset, *import_offset;
static int *Head, *Len, *Next, *Tail, *MinID, *MinIDTask;
static int *GroupLen, *GroupOffset, *GroupLenType;
static double *GroupMassType;
static float *GroupCM, *GroupSfr;
static int *GroupIDs;

#ifdef BLACK_HOLES
static float *GroupMbh;
static float *GroupMdot;
static int num_bh;
static int *NBHs;
static float **GroupMbh_ind;
static float **GroupMdot_ind;
static float ***GroupPosbh_ind;
#endif

static float *fof_nearest_distance;
static float *fof_nearest_hsml;



void fof_fof(int num)
{
  int i, ndm, ndmtot, largestgroup;
  double mass, masstot, rhodm;
#ifdef BLACK_HOLES
  int k;
#endif

  if(ThisTask == 0)
    {
      printf("\nBegin to compute FoF group catalogues...");
      fflush(stdout);
    }


  if(All.NumForcesSinceLastDomainDecomp > All.TotNumPart * All.TreeDomainUpdateFrequency)
    {
      DomainDecomposition();
      if(ThisTask == 0)
	printf("Tree construction.\n");
      force_treebuild(NumPart);
      All.NumForcesSinceLastDomainDecomp = 0;
      TreeReconstructFlag = 0;
    }

  for(i = 0, ndm = 0, mass = 0; i < NumPart; i++)
    if(P[i].Type == 1)
      {
	ndm++;
	mass += P[i].Mass;
      }

  MPI_Allreduce(&ndm, &ndmtot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  rhodm = (All.Omega0 - All.OmegaBaryon) * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

  LinkL = LINKLENGTH * pow(masstot / ndmtot / rhodm, 1.0 / 3);

  if(ThisTask == 0)
    {
      printf("\nComoving linking length: %g\n", LinkL);
      fflush(stdout);
    }

  Head = mymalloc(NumPart * sizeof(int));
  Len = mymalloc(NumPart * sizeof(int));
  Next = mymalloc(NumPart * sizeof(int));
  Tail = mymalloc(NumPart * sizeof(int));
  MinID = mymalloc(NumPart * sizeof(int));
  MinIDTask = mymalloc(NumPart * sizeof(int));

  local_toGo = mymalloc(sizeof(int) * NTask);
  toGo = mymalloc(sizeof(int) * NTask * NTask);
  export_offset = mymalloc(sizeof(int) * NTask);
  import_offset = mymalloc(sizeof(int) * NTask);


  for(i = 0; i < NumPart; i++)
    {
      Head[i] = Tail[i] = i;
      Len[i] = 1;
      Next[i] = -1;
      MinID[i] = P[i].ID;
      MinIDTask[i] = ThisTask;
    }

  fof_find_groups();

  fof_find_nearest_dmparticle();

  fof_exchange_id_lists();

  fof_compile_catalogue();



  MPI_Allreduce(&Ngroups, &TotNgroups, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allgather(&Nids, 1, MPI_INT, local_toGo, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 0, TotNids = 0; i < NTask; i++)
    TotNids += local_toGo[i];

  if(TotNgroups > 0)
    MPI_Allreduce(&GroupLen[0], &largestgroup, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  else
    largestgroup = 0;


  if(ThisTask == 0)
    {
      printf("\nTotal number of groups with at least %d particles: %d\n", GROUP_MIN_LEN, TotNgroups);
      if(TotNgroups > 0)
	{
	  printf("Largest group has %d particles.\n", largestgroup);
	  printf("Total number of particles in groups: %d%09d\n\n",
		 (int) (TotNids / 1000000000), (int) (TotNids % 1000000000));
	}
    }

#ifdef BLACK_HOLES
  if(num < 0)
    fof_make_black_holes();
#endif

#ifdef BUBBLES
  if(num < 0)
    find_CM_of_biggest_group();
#endif

#ifdef MULTI_BUBBLES
  if(num < 0)
    multi_bubbles();
#endif

  if(num >= 0)
    {
      fof_save_groups(num);

      if(ThisTask == 0)
	{
	  printf("Group catalogues saved.\n\n");
	  fflush(stdout);
	}
    }

#ifdef BLACK_HOLES
  for(k = 2; k >= 0; k--)
    {
      for(i = Ngroups - 1; i >= 0; i--) 
	myfree(GroupPosbh_ind[k][i]);
      myfree(GroupPosbh_ind[k]);
    } 
  myfree(GroupPosbh_ind);
  
  for(i = Ngroups - 1; i >= 0; i--)
    myfree(GroupMdot_ind[i]);
  myfree(GroupMdot_ind);
  
  for(i = Ngroups - 1; i >= 0; i--)
    myfree(GroupMbh_ind[i]);
  myfree(GroupMbh_ind);
#endif
  myfree(GroupIDs);
  myfree(GroupOffset);
  myfree(GroupSfr);
  myfree(GroupCM);
  myfree(GroupMassType);
  myfree(GroupLenType);
  myfree(GroupLen);
#ifdef BLACK_HOLES
  myfree(GroupMdot);
  myfree(GroupMbh);
  myfree(NBHs);
#endif
  
#ifdef BUBBLES
  myfree(BiggestGroup);
#endif


  myfree(local_ids);

  myfree(import_offset);
  myfree(export_offset);
  myfree(toGo);
  myfree(local_toGo);


  myfree(MinIDTask);
  myfree(MinID);
  myfree(Tail);
  myfree(Next);
  myfree(Len);
  myfree(Head);


}



void fof_find_groups(void)
{
  int i, j, n, count, ntot, ntotleft, ndonetot, level, maxfill, link_inner, link_inner_tot, link_count;
  int nexport, ndone, ngrp, sendTask, recvTask, link_across, link_across_tot;
  int *noffset, *nbuffer, *nsend_local, *nsend;
  MPI_Status status;

  noffset = mymalloc(sizeof(int) * NTask);	/* offsets of bunches in common list */
  nbuffer = mymalloc(sizeof(int) * NTask);
  nsend_local = mymalloc(sizeof(int) * NTask);
  nsend = mymalloc(sizeof(int) * NTask * NTask);

  if(ThisTask == 0)
    {
      printf("\nStart linking particles\n");
      fflush(stdout);
    }

  for(n = 0, count = 0; n < NumPart; n++)
    {
      if(P[n].Type == 1)
	count++;
    }

  MPI_Allreduce(&count, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);


  do
    {
      i = 0;			/* beginn with this index */
      ntotleft = ntot;		/* particles left for all tasks together */
      link_across = 0;

      while(ntotleft > 0)
	{
	  for(j = 0; j < NTask; j++)
	    nsend_local[j] = 0;

	  /* do local particles and prepare export list */

	  for(nexport = 0, ndone = 0; i < NumPart && nexport < All.BunchSizeFoF - NTask; i++)
	    if(P[i].Type == 1)
	      {
		ndone++;

		for(j = 0; j < NTask; j++)
		  Exportflag[j] = 0;

		fof_find_dmparticles_evaluate(i, 0);

		for(j = 0; j < NTask; j++)
		  {
		    if(Exportflag[j])
		      {
			FoFDataIn[nexport].Pos[0] = P[i].Pos[0];
			FoFDataIn[nexport].Pos[1] = P[i].Pos[1];
			FoFDataIn[nexport].Pos[2] = P[i].Pos[2];
			FoFDataIn[nexport].Index = i;
			FoFDataIn[nexport].Task = j;
			nexport++;
			nsend_local[j]++;
		      }
		  }
	      }

	  qsort(FoFDataIn, nexport, sizeof(struct fofdata_in), fof_compare_key);

	  for(j = 1, noffset[0] = 0; j < NTask; j++)
	    noffset[j] = noffset[j - 1] + nsend_local[j - 1];

	  MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);

	  /* now do the particles that need to be exported */

	  do
	    {
	      for(j = 0; j < nexport; j++)
		{
		  FoFDataIn[j].MinID = MinID[Head[FoFDataIn[j].Index]];
		  FoFDataIn[j].MinIDTask = MinIDTask[Head[FoFDataIn[j].Index]];
		}

	      link_inner = 0;

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
		      if(maxfill >= All.BunchSizeFoF)
			break;

		      sendTask = ThisTask;
		      recvTask = ThisTask ^ ngrp;

		      if(recvTask < NTask)
			{
			  if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
			    {
			      /* get the particles */
			      MPI_Sendrecv(&FoFDataIn[noffset[recvTask]],
					   nsend_local[recvTask] * sizeof(struct fofdata_in), MPI_BYTE,
					   recvTask, TAG_FOF_F,
					   &FoFDataGet[nbuffer[ThisTask]],
					   nsend[recvTask * NTask + ThisTask] * sizeof(struct fofdata_in),
					   MPI_BYTE, recvTask, TAG_FOF_F, MPI_COMM_WORLD, &status);
			    }
			}

		      for(j = 0; j < NTask; j++)
			if((j ^ ngrp) < NTask)
			  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
		    }


		  for(j = 0; j < nbuffer[ThisTask]; j++)
		    {
		      link_count = fof_find_dmparticles_evaluate(j, 1);
		      link_inner += link_count;
		      link_across += link_count;
		    }

		  level = ngrp - 1;
		}

	      MPI_Allreduce(&link_inner, &link_inner_tot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	      if(ThisTask == 0)
		{
		  printf("  link inner %d\n", link_inner_tot);
		  fflush(stdout);
		}
	    }
	  while(link_inner_tot > 0);

	  MPI_Allreduce(&ndone, &ndonetot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	  ntotleft -= ndonetot;
	}

      MPI_Allreduce(&link_across, &link_across_tot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      if(ThisTask == 0)
	{
	  printf("have done %d cross links\n", link_across_tot);
	  fflush(stdout);
	}
    }
  while(link_across_tot > 0);


  myfree(nsend);
  myfree(nsend_local);
  myfree(nbuffer);
  myfree(noffset);

  if(ThisTask == 0)
    {
      printf("Local groups found.\n");
      fflush(stdout);
    }
}


int fof_find_dmparticles_evaluate(int target, int mode)
{
  int j, n, links, p, s, ss;
  int startnode, numngb_inbox;
  double h2, dx, dy, dz, r2;
  FLOAT pos[3];

  links = 0;

  if(mode == 0)
    {
      pos[0] = P[target].Pos[0];
      pos[1] = P[target].Pos[1];
      pos[2] = P[target].Pos[2];
    }
  else
    {
      pos[0] = FoFDataGet[target].Pos[0];
      pos[1] = FoFDataGet[target].Pos[1];
      pos[2] = FoFDataGet[target].Pos[2];
    }

  startnode = All.MaxPart;
  h2 = LinkL * LinkL;

  do
    {
      numngb_inbox = ngb_treefind_darkmatter(&pos[0], LinkL, &startnode);

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

	  if(r2 <= h2)		/* this is a pair */
	    {
	      if(mode == 0)
		{
		  if(Head[target] != Head[j])	/* only if not yet linked */
		    {
		      if(Len[Head[target]] > Len[Head[j]])	/* p group is longer */
			{
			  p = target;
			  s = j;
			}
		      else
			{
			  p = j;
			  s = target;
			}
		      Next[Tail[Head[p]]] = Head[s];

		      Tail[Head[p]] = Tail[Head[s]];

		      Len[Head[p]] += Len[Head[s]];

		      ss = Head[s];
		      do
			Head[ss] = Head[p];
		      while((ss = Next[ss]) >= 0);

		      if(MinID[Head[s]] < MinID[Head[p]])
			{
			  MinID[Head[p]] = MinID[Head[s]];
			  MinIDTask[Head[p]] = MinIDTask[Head[s]];
			}
		    }
		}
	      else		/* mode is 1 */
		{
		  if(MinID[Head[j]] > FoFDataGet[target].MinID)
		    {
		      MinID[Head[j]] = FoFDataGet[target].MinID;
		      MinIDTask[Head[j]] = FoFDataGet[target].MinIDTask;
		      links++;
		    }
		}

	    }
	}
    }
  while(startnode >= 0);

  return links;
}


void fof_exchange_id_lists(void)
{
  int n, task, nsend, nget, level, sendTask, recvTask, j, k;

  MPI_Status status;

  for(n = 0; n < NTask; n++)
    local_toGo[n] = 0;

  for(n = 0; n < NumPart; n++)
    {
      task = MinIDTask[Head[n]];

      if(task != ThisTask)
	local_toGo[task]++;
    }

  MPI_Allgather(local_toGo, NTask, MPI_INT, toGo, NTask, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nsend = nget = 0; j < NTask; j++)
    {
      nsend += toGo[ThisTask * NTask + j];
      nget += toGo[j * NTask + ThisTask];
    }

  printf("Task=%d needs to send %d particles for FoF-grouplist and gets %d particles.\n", ThisTask, nsend,
	 nget);
  fflush(stdout);

  Nlocal = NumPart - nsend + nget;

  local_ids = mymalloc(1 + Nlocal * sizeof(struct id_list));
  export_ids = mymalloc(1 + nsend * sizeof(struct id_list));

  if(!(export_ids) || !(local_ids))
    {
      printf("failed to allocate memory\n");
      endrun(152311);
    }

  import_offset[0] = NumPart - nsend;
  export_offset[0] = 0;

  for(n = 1; n < NTask; n++)
    {
      export_offset[n] = export_offset[n - 1] + toGo[ThisTask * NTask + (n - 1)];
      import_offset[n] = import_offset[n - 1] + toGo[(n - 1) * NTask + ThisTask];
    }

  for(n = 0; n < NTask; n++)
    local_toGo[n] = 0;

  for(n = 0; n < NumPart; n++)
    {
      task = MinIDTask[Head[n]];

      if(task != ThisTask)
	{
	  export_ids[export_offset[task] + local_toGo[task]].GrIndex = MinID[Head[n]];
	  export_ids[export_offset[task] + local_toGo[task]].ID = P[n].ID;
	  export_ids[export_offset[task] + local_toGo[task]].Type = P[n].Type;
	  export_ids[export_offset[task] + local_toGo[task]].Mass = P[n].Mass;
	  for(k = 0; k < 3; k++)
	    export_ids[export_offset[task] + local_toGo[task]].Pos[k] = P[n].Pos[k];
	  export_ids[export_offset[task] + local_toGo[task]].Task = ThisTask;
	  export_ids[export_offset[task] + local_toGo[task]].Index = n;
	  export_ids[export_offset[task] + local_toGo[task]].Sfr = 0;
#ifdef BLACK_HOLES
	  if(P[n].Type == 5)
	    {
	      export_ids[export_offset[task] + local_toGo[task]].BH_Mass = P[n].BH_Mass;
	      export_ids[export_offset[task] + local_toGo[task]].BH_Mdot = P[n].BH_Mdot;
	    }
	  else
	    {
	      export_ids[export_offset[task] + local_toGo[task]].BH_Mass = 0;
	      export_ids[export_offset[task] + local_toGo[task]].BH_Mdot = 0;
	    }
#endif
	  if(P[n].Type == 0)
	    {
	      export_ids[export_offset[task] + local_toGo[task]].Density = SphP[n].a2.Density;
#ifdef SFR
	      export_ids[export_offset[task] + local_toGo[task]].Sfr = SphP[n].Sfr;
#endif
#ifdef WINDS
	      export_ids[export_offset[task] + local_toGo[task]].DelayTime = SphP[n].DelayTime;
#endif
	    }
	}
      else
	{
	  local_ids[local_toGo[task]].GrIndex = MinID[Head[n]];
	  local_ids[local_toGo[task]].ID = P[n].ID;
	  local_ids[local_toGo[task]].Type = P[n].Type;
	  local_ids[local_toGo[task]].Mass = P[n].Mass;
	  for(k = 0; k < 3; k++)
	    local_ids[local_toGo[task]].Pos[k] = P[n].Pos[k];
	  local_ids[local_toGo[task]].Task = ThisTask;
	  local_ids[local_toGo[task]].Index = n;
	  local_ids[local_toGo[task]].Sfr = 0;

	  if(P[n].Type == 0)
	    {
	      local_ids[local_toGo[task]].Density = SphP[n].a2.Density;
#ifdef SFR
	      local_ids[local_toGo[task]].Sfr = SphP[n].Sfr;
#endif
#ifdef WINDS
	      local_ids[local_toGo[task]].DelayTime = SphP[n].DelayTime;
#endif
	    }

#ifdef BLACK_HOLES
	  if(P[n].Type == 5)
	    {
	      local_ids[local_toGo[task]].BH_Mass = P[n].BH_Mass;
	      local_ids[local_toGo[task]].BH_Mdot = P[n].BH_Mdot;
	    }
	  else
	    {
	      local_ids[local_toGo[task]].BH_Mass = 0;
	      local_ids[local_toGo[task]].BH_Mdot = 0;
	    }
#endif
	}

      local_toGo[task]++;
    }

  for(level = 1; level < (1 << PTask); level++)
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
	MPI_Sendrecv(&export_ids[export_offset[recvTask]],
		     toGo[ThisTask * NTask + recvTask] * sizeof(struct id_list),
		     MPI_BYTE, recvTask, TAG_FOF_D,
		     &local_ids[import_offset[recvTask]],
		     toGo[recvTask * NTask + ThisTask] * sizeof(struct id_list),
		     MPI_BYTE, recvTask, TAG_FOF_D, MPI_COMM_WORLD, &status);
    }

  myfree(export_ids);
}



int fof_grid_compare(const void *a, const void *b)
{
  if(((struct id_list *) a)->GrIndex > (((struct id_list *) b)->GrIndex))
    return -1;

  if(((struct id_list *) a)->GrIndex < (((struct id_list *) b)->GrIndex))
    return +1;

  if(((struct id_list *) a)->Type < (((struct id_list *) b)->Type))
    return -1;

  if(((struct id_list *) a)->Type > (((struct id_list *) b)->Type))
    return +1;


  if(((struct id_list *) a)->ID < (((struct id_list *) b)->ID))
    return -1;

  if(((struct id_list *) a)->ID > (((struct id_list *) b)->ID))
    return +1;

  return 0;
}




void fof_compile_catalogue(void)
{
  int i, j, k, start, groups;

#ifdef BLACK_HOLES
  int l;
#endif
  
  qsort(local_ids, Nlocal, sizeof(struct id_list), fof_grid_compare);

  Ngroups = 0;
  Nids = 0;
  groups = 0;

  for(i = 0, start = 0; i < Nlocal; i++)
    {
      if(local_ids[i].GrIndex != local_ids[start].GrIndex)
	{
	  for(j = start; j < i; j++)
	    local_ids[j].GrIndex = groups + (((long long) (i - start)) << 32);

	  if((i - start) >= GROUP_MIN_LEN)
	    {
	      Ngroups++;
	      Nids += (i - start);
	    }

	  groups++;
	  start = i;
	}
    }

  /* finish last group */
  for(j = start; j < Nlocal; j++)
    local_ids[j].GrIndex = groups + (((long long) (Nlocal - start)) << 32);

  if((Nlocal - start) >= GROUP_MIN_LEN)
    {
      Ngroups++;
      Nids += (Nlocal - start);
    }


  /* now sort the groups by size */
  qsort(local_ids, Nlocal, sizeof(struct id_list), fof_grid_compare);

#ifdef BLACK_HOLES
  NBHs = mymalloc(1 + Ngroups * sizeof(int));
  GroupMbh = mymalloc((1 + Ngroups * sizeof(float)));
  GroupMdot = mymalloc((1 + Ngroups * sizeof(float)));
#endif
  GroupLen = mymalloc(1 + Ngroups * sizeof(int));
  GroupLenType = mymalloc(6 * (1 + Ngroups * sizeof(int)));
  GroupMassType = mymalloc(6 * (1 + Ngroups * sizeof(double)));
  GroupCM = mymalloc(3 * (1 + Ngroups * sizeof(float)));
  GroupSfr = mymalloc((1 + Ngroups * sizeof(float)));
  GroupOffset = mymalloc(1 + Ngroups * sizeof(int));
  GroupIDs = mymalloc(1 + Nids * sizeof(int));


#ifdef BUBBLES
  BiggestGroup = mymalloc(NTask * sizeof(struct biggest_group));
#endif


  if(!(GroupLen) || !(GroupOffset) || !(GroupIDs))
    {
      printf("failed to allocate memory\n");
      endrun(1523711);
    }

  GroupLen[0] = 0;

#ifdef BLACK_HOLES
  num_bh = 0;
#endif

  for(i = 0, j = 0, GroupOffset[0] = 0; i < Ngroups; i++)
    {
      GroupLen[i] = (local_ids[j].GrIndex >> 32);
      GroupSfr[i] = 0;
#ifdef BLACK_HOLES
      GroupMdot[i] = 0;
      GroupMbh[i] = 0;
#endif
#ifdef BLACK_HOLES
      NBHs[i] = 0;
#endif
      for(k = 0; k < 3; k++)
	{
	  GroupCM[i * 3 + k] = 0;
	}

      for(k = 0; k < 6; k++)
	{
	  GroupLenType[i * 6 + k] = 0;
	  GroupMassType[i * 6 + k] = 0;
	}

      for(k = 0; k < GroupLen[i]; k++)
	{
	  GroupLenType[i * 6 + local_ids[j + k].Type]++;
	  GroupMassType[i * 6 + local_ids[j + k].Type] += local_ids[j + k].Mass;

#ifdef SFR
	  if(local_ids[j + k].Type == 0)
	    GroupSfr[i] += local_ids[j + k].Sfr;
#endif
#ifdef BLACK_HOLES
	  if(local_ids[j + k].Type == 5)
	    {
	      GroupMdot[i] += local_ids[j + k].BH_Mdot;
	      GroupMbh[i] += local_ids[j + k].BH_Mass;
	      NBHs[i]++;
	      num_bh++;
	    }
#endif
	}

      fof_get_group_center(&GroupCM[3 * i], j, GroupLen[i]);

      j += GroupLen[i];
      if(i > 0)
	GroupOffset[i] = GroupOffset[i - 1] + GroupLen[i - 1];
    }

  for(i = 0; i < Nids; i++)
    GroupIDs[i] = local_ids[i].ID;

#ifdef BLACK_HOLES
  GroupMbh_ind = mymalloc(1 + Ngroups * sizeof(float *));
  for(i = 0; i < Ngroups; i++)
    GroupMbh_ind[i] = mymalloc(1 + NBHs[i] * sizeof(float));
  
  GroupMdot_ind = mymalloc(1 + Ngroups * sizeof(float *));
  for(i = 0; i < Ngroups; i++)
    GroupMdot_ind[i] = mymalloc(1 + NBHs[i] * sizeof(float));

  GroupPosbh_ind = mymalloc(3 * sizeof(float **));
  for(j = 0; j < 3; j++)
    {
      GroupPosbh_ind[j] = mymalloc(1 + Ngroups * sizeof(float *));
      for(i = 0; i < Ngroups; i++)
	GroupPosbh_ind[j][i] = mymalloc(1 + NBHs[i] * sizeof(float));
    }


  GroupLen[0] = 0;
  
  for(i = 0, j = 0, GroupOffset[0] = 0; i < Ngroups; i++)
    {
      GroupLen[i] = (local_ids[j].GrIndex >> 32);
      
      for(l = 0; l < NBHs[i]; l++)
	{
	  GroupMbh_ind[i][l] = 0;
	  GroupMdot_ind[i][l] = 0;
	  GroupPosbh_ind[0][i][l] = 0;
	  GroupPosbh_ind[1][i][l] = 0;
	  GroupPosbh_ind[2][i][l] = 0;
	}
	  
      for(k = 0, l = 0; k < GroupLen[i]; k++)
	{
	  if(local_ids[j + k].Type == 5)
	    { 
	      GroupMbh_ind[i][l] = local_ids[j + k].BH_Mass;
	      GroupMdot_ind[i][l] = local_ids[j + k].BH_Mdot;
	      GroupPosbh_ind[0][i][l] = local_ids[j + k].Pos[0];
	      GroupPosbh_ind[1][i][l] = local_ids[j + k].Pos[1];
	      GroupPosbh_ind[2][i][l] = local_ids[j + k].Pos[2];
	      l++;
	    }
	}	
      
            
      j += GroupLen[i];
      if(i > 0)
	GroupOffset[i] = GroupOffset[i - 1] + GroupLen[i - 1];
    }
  
#endif
}

void fof_get_group_center(float *cm, int first, int len)
{
  double masstot, x, y, z, xx, yy, zz, cx, cy, cz, r, rmin, rmax;
  int i, count;


#ifdef PERIODIC
  cx = cy = cz = 0;
#else
  cx = local_ids[first].Pos[0];
  cy = local_ids[first].Pos[1];
  cz = local_ids[first].Pos[2];
#endif
  rmax = 1.0e30;

  do
    {
      rmin = 0;
      count = 0;
      xx = yy = zz = masstot = 0;

      for(i = 0; i < len; i++)
	{
	  x = local_ids[first + i].Pos[0];
	  y = local_ids[first + i].Pos[1];
	  z = local_ids[first + i].Pos[2];
#ifdef PERIODIC
	  x = fof_periodic(x - local_ids[first].Pos[0]);
	  y = fof_periodic(y - local_ids[first].Pos[1]);
	  z = fof_periodic(z - local_ids[first].Pos[2]);
#endif
	  r = sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy) + (z - cz) * (z - cz));
	  if(r < rmax)
	    {
	      xx += local_ids[first + i].Mass * x;
	      yy += local_ids[first + i].Mass * y;
	      zz += local_ids[first + i].Mass * z;
	      masstot += local_ids[first + i].Mass;
	      count++;

	      if(r > rmin)
		rmin = r;
	    }
	}

      if(count)
	{
	  cx = xx / masstot;
	  cy = yy / masstot;
	  cz = zz / masstot;
	  rmax = 0.8 * rmin;

	  if(len > 100000)
	    {
	      printf("%g %g %g  rmax=%g  count=%d\n", cx, cy, cz, rmax, count);
	      fflush(stdout);
	    }
	}
    }
  while(count > 100);

#ifdef PERIODIC
  cx = fof_periodic_wrap(cx + local_ids[first].Pos[0]);
  cy = fof_periodic_wrap(cy + local_ids[first].Pos[1]);
  cz = fof_periodic_wrap(cz + local_ids[first].Pos[2]);
#endif

  cm[0] = cx;
  cm[1] = cy;
  cm[2] = cz;
}


#ifdef BLACK_HOLES

void fof_make_black_holes(void)
{
  int i, j, k, n, index, nblackhole, ntot, task;
  int nsend, nget, nlocal, sendTask, recvTask, level;
  int *local_indices, *export_indices;
  double maxdens, massDMpart;
  MPI_Status status;

  for(n = 0; n < NTask; n++)
    local_toGo[n] = 0;

  nblackhole = 0;

  if(All.MassTable[1] > 0)
    massDMpart = All.MassTable[1];
  else
    massDMpart = All.massDMpart;

  for(i = 0; i < Ngroups; i++)
    {

      if(GroupLenType[i * 6 + 1] * massDMpart >=
	 (All.Omega0 - All.OmegaBaryon) / All.Omega0 * All.MinFoFMassForNewSeed)
	if(GroupLenType[i * 6 + 5] == 0)
	  {
	    maxdens = 0;
	    index = -1;
	    for(k = 0; k < GroupLen[i]; k++)
	      if(local_ids[GroupOffset[i] + k].Type == 0)
		{
#ifdef WINDS
		  if(local_ids[GroupOffset[i] + k].DelayTime == 0)
#endif
		    if(local_ids[GroupOffset[i] + k].Density > maxdens)
		      {
			maxdens = local_ids[GroupOffset[i] + k].Density;
			index = GroupOffset[i] + k;
		      }

		}

	    if(index >= 0)
	      {
		nblackhole++;
		if(local_ids[index].Task != ThisTask)
		  local_toGo[local_ids[index].Task]++;
	      }
	  }
    }


  MPI_Allgather(local_toGo, NTask, MPI_INT, toGo, NTask, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nsend = nget = 0; j < NTask; j++)
    {
      nsend += toGo[ThisTask * NTask + j];
      nget += toGo[j * NTask + ThisTask];
    }

  nlocal = nblackhole - nsend + nget;

  local_indices = mymalloc(nlocal * sizeof(int));
  export_indices = mymalloc(nsend * sizeof(int));

  import_offset[0] = nblackhole - nsend;
  export_offset[0] = 0;

  for(n = 1; n < NTask; n++)
    {
      export_offset[n] = export_offset[n - 1] + toGo[ThisTask * NTask + (n - 1)];
      import_offset[n] = import_offset[n - 1] + toGo[(n - 1) * NTask + ThisTask];
    }

  for(n = 0; n < NTask; n++)
    local_toGo[n] = 0;


  for(i = 0; i < Ngroups; i++)
    {
      if(GroupLenType[i * 6 + 1] * massDMpart >=
	 (All.Omega0 - All.OmegaBaryon) / All.Omega0 * All.MinFoFMassForNewSeed)
	if(GroupLenType[i * 6 + 5] == 0)
	  {
	    maxdens = 0;
	    index = -1;
	    for(k = 0; k < GroupLen[i]; k++)
	      if(local_ids[GroupOffset[i] + k].Type == 0)
		{
#ifdef WINDS
		  if(local_ids[GroupOffset[i] + k].DelayTime == 0)
#endif
		    if(local_ids[GroupOffset[i] + k].Density > maxdens)
		      {
			maxdens = local_ids[GroupOffset[i] + k].Density;
			index = GroupOffset[i] + k;
		      }

		}
	    if(index >= 0)
	      {
		if(local_ids[index].Task != ThisTask)
		  {
		    task = local_ids[index].Task;
		    export_indices[export_offset[task] + local_toGo[task]] = local_ids[index].Index;
		  }
		else
		  {
		    local_indices[local_toGo[ThisTask]] = local_ids[index].Index;
		  }
		local_toGo[local_ids[index].Task]++;


		GroupLenType[i * 6 + 1]--;
		GroupLenType[i * 6 + 5]++;
	      }
	  }
    }


  for(level = 1; level < (1 << PTask); level++)
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
	MPI_Sendrecv(&export_indices[export_offset[recvTask]],
		     toGo[ThisTask * NTask + recvTask] * sizeof(int),
		     MPI_BYTE, recvTask, TAG_FOF_E,
		     &local_indices[import_offset[recvTask]],
		     toGo[recvTask * NTask + ThisTask] * sizeof(int),
		     MPI_BYTE, recvTask, TAG_FOF_E, MPI_COMM_WORLD, &status);
    }

  MPI_Allreduce(&nlocal, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      printf("\nMaking %d new black hole particles\n\n", ntot);
      fflush(stdout);
    }

  for(n = 0; n < nlocal; n++)
    {
      if(P[local_indices[n]].Type != 0)
	endrun(7772);

      P[local_indices[n]].Type = 5;	/* make it a black hole particle */
      P[local_indices[n]].StellarAge = All.Time;
      P[local_indices[n]].BH_Mass = All.SeedBlackHoleMass;
      P[local_indices[n]].BH_Mdot = 0;
#ifdef BH_BUBBLES
      P[local_indices[n]].BH_Mass_ini = All.SeedBlackHoleMass;
#endif
      Stars_converted++;
    }

  All.TotN_gas -= ntot;
  myfree(export_indices);
  myfree(local_indices);
}

#endif






void fof_save_groups(int num)
{
  int nprocgroup, masterTask, groupTask;
  char buf[500];

  if(ThisTask == 0)
    {
      sprintf(buf, "%s/groups_%03d", All.OutputDir, num);
      mkdir(buf, 02755);
    }
  MPI_Barrier(MPI_COMM_WORLD);


  if(NTask < All.NumFilesWrittenInParallel)

    {
      printf
	("Fatal error.\nNumber of processors must be a smaller or equal than `NumFilesWrittenInParallel'.\n");
      endrun(241931);
    }

  nprocgroup = NTask / All.NumFilesWrittenInParallel;
  if((NTask % All.NumFilesWrittenInParallel))
    nprocgroup++;
  masterTask = (ThisTask / nprocgroup) * nprocgroup;
  for(groupTask = 0; groupTask < nprocgroup; groupTask++)
    {
      if(ThisTask == (masterTask + groupTask))	/* ok, it's this processor's turn */
	fof_save_local_catalogue(num);
      MPI_Barrier(MPI_COMM_WORLD);	/* wait inside the group */
    }
}




void fof_save_local_catalogue(int num)
{
  FILE *fd;
  char buf[500];
#ifdef BLACK_HOLES
  int i,k;
#endif
  
  sprintf(buf, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "group_tab", num, ThisTask);
  if(!(fd = fopen(buf, "w")))
    {
      printf("can't open file `%s`\n", buf);
      endrun(1183);
    }

  fwrite(&Ngroups, sizeof(int), 1, fd);
  fwrite(&Nids, sizeof(int), 1, fd);
  fwrite(&TotNgroups, sizeof(int), 1, fd);
  fwrite(&NTask, sizeof(int), 1, fd);
  fwrite(GroupLen, sizeof(int), Ngroups, fd);
  fwrite(GroupOffset, sizeof(int), Ngroups, fd);
  fwrite(GroupLenType, 6 * sizeof(int), Ngroups, fd);
  fwrite(GroupMassType, 6 * sizeof(double), Ngroups, fd);
  fwrite(GroupCM, 3 * sizeof(float), Ngroups, fd);
  fwrite(GroupSfr, sizeof(float), Ngroups, fd);
#ifdef BLACK_HOLES
  fwrite(GroupMbh, sizeof(float), Ngroups, fd);
  fwrite(GroupMdot, sizeof(float), Ngroups, fd);
  fwrite(NBHs, sizeof(int), Ngroups, fd);
  
  for(i = 0; i < Ngroups; i++)
    fwrite(GroupMbh_ind[i], sizeof(float), NBHs[i], fd);
  
  for(i = 0; i < Ngroups; i++)
    fwrite(GroupMdot_ind[i], sizeof(float), NBHs[i], fd);
  
  for(k = 0; k < 3; k++)
    for(i = 0; i < Ngroups; i++)
      fwrite(GroupPosbh_ind[k][i], sizeof(float), NBHs[i], fd);
#endif
  fclose(fd);
  sprintf(buf, "%s/groups_%03d/%s_%03d.%d", All.OutputDir, num, "group_ids", num, ThisTask);
  if(!(fd = fopen(buf, "w")))
    {
      printf("can't open file `%s`\n", buf);
      endrun(1184);
    }

  fwrite(&Ngroups, sizeof(int), 1, fd);
  fwrite(&Nids, sizeof(int), 1, fd);
  fwrite(&TotNgroups, sizeof(int), 1, fd);
  fwrite(&NTask, sizeof(int), 1, fd);
  fwrite(GroupIDs, sizeof(int), Nids, fd);
  fclose(fd);
}



void fof_find_nearest_dmparticle(void)
{
  int i, j, n, count, ntot, ntotleft, ndonetot, level, maxfill;
  int nexport, ndone, ngrp, sendTask, recvTask, place, source, npleft, iter;
  int *noffset, *nbuffer, *nsend_local, *nsend;
  MPI_Status status;
  fof_nearest_distance = mymalloc(sizeof(float) * NumPart);
  fof_nearest_hsml = mymalloc(sizeof(float) * NumPart);
  noffset = mymalloc(sizeof(int) * NTask);	/* offsets of bunches in common list */
  nbuffer = mymalloc(sizeof(int) * NTask);
  nsend_local = mymalloc(sizeof(int) * NTask);
  nsend = mymalloc(sizeof(int) * NTask * NTask);
  if(ThisTask == 0)
    {
      printf("\nStart finding nearest dm-particle\n");
      fflush(stdout);
    }

  for(n = 0, count = 0; n < NumPart; n++)
    {
      if(P[n].Type == 0 || P[n].Type >= 4)
	{
	  fof_nearest_distance[n] = 1.0e30;
	  if(P[n].Type == 0)
	    fof_nearest_hsml[n] = PPP[n].Hsml;
	  else
	    fof_nearest_hsml[n] = 0.1 * LinkL;
	  count++;
	}
    }

  MPI_Allreduce(&count, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  iter = 0;
  /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
  do
    {
      i = 0;			/* beginn with this index */
      ntotleft = ntot;		/* particles left for all tasks together */
      while(ntotleft > 0)
	{
	  for(j = 0; j < NTask; j++)
	    nsend_local[j] = 0;
	  /* do local particles and prepare export list */
	  for(nexport = 0, ndone = 0; i < NumPart && nexport < All.BunchSizeFoF - NTask; i++)
	    if(P[i].Type == 0 || P[i].Type >= 4)
	      {
		ndone++;
		for(j = 0; j < NTask; j++)
		  Exportflag[j] = 0;
		fof_find_nearest_dmparticle_evaluate(i, 0);
		for(j = 0; j < NTask; j++)
		  {
		    if(Exportflag[j])
		      {
			FoFDataIn[nexport].Pos[0] = P[i].Pos[0];
			FoFDataIn[nexport].Pos[1] = P[i].Pos[1];
			FoFDataIn[nexport].Pos[2] = P[i].Pos[2];
			FoFDataIn[nexport].Hsml = fof_nearest_hsml[i];
			FoFDataIn[nexport].Index = i;
			FoFDataIn[nexport].Task = j;
			nexport++;
			nsend_local[j]++;
		      }
		  }
	      }

	  qsort(FoFDataIn, nexport, sizeof(struct fofdata_in), fof_compare_key);
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
		  if(maxfill >= All.BunchSizeFoF)
		    break;
		  sendTask = ThisTask;
		  recvTask = ThisTask ^ ngrp;
		  if(recvTask < NTask)
		    {
		      if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
			{
			  /* get the particles */
			  MPI_Sendrecv(&FoFDataIn[noffset[recvTask]],
				       nsend_local[recvTask] * sizeof(struct fofdata_in), MPI_BYTE,
				       recvTask, TAG_FOF_F,
				       &FoFDataGet[nbuffer[ThisTask]],
				       nsend[recvTask * NTask + ThisTask] * sizeof(struct fofdata_in),
				       MPI_BYTE, recvTask, TAG_FOF_F, MPI_COMM_WORLD, &status);
			}
		    }

		  for(j = 0; j < NTask; j++)
		    if((j ^ ngrp) < NTask)
		      nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
		}


	      for(j = 0; j < nbuffer[ThisTask]; j++)
		fof_find_nearest_dmparticle_evaluate(j, 1);
	      /* get the result */
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
		  if(maxfill >= All.BunchSizeFoF)
		    break;
		  sendTask = ThisTask;
		  recvTask = ThisTask ^ ngrp;
		  if(recvTask < NTask)
		    {
		      if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
			{
			  /* send the results */
			  MPI_Sendrecv(&FoFDataResult[nbuffer[ThisTask]],
				       nsend[recvTask * NTask + ThisTask] * sizeof(struct fofdata_out),
				       MPI_BYTE, recvTask, TAG_FOF_G,
				       &FoFDataPartialResult[noffset[recvTask]],
				       nsend_local[recvTask] * sizeof(struct fofdata_out),
				       MPI_BYTE, recvTask, TAG_FOF_G, MPI_COMM_WORLD, &status);
			  /* add the result to the particles */
			  for(j = 0; j < nsend_local[recvTask]; j++)
			    {
			      source = j + noffset[recvTask];
			      place = FoFDataIn[source].Index;
			      if(FoFDataPartialResult[source].Distance < fof_nearest_distance[place])
				{
				  fof_nearest_distance[place] = FoFDataPartialResult[source].Distance;
				  MinID[place] = FoFDataPartialResult[source].MinID;
				  MinIDTask[place] = FoFDataPartialResult[source].MinIDTask;
				}
			    }
			}
		    }

		  for(j = 0; j < NTask; j++)
		    if((j ^ ngrp) < NTask)
		      nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
		}

	      level = ngrp - 1;
	    }

	  MPI_Allreduce(&ndone, &ndonetot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	  ntotleft -= ndonetot;
	}

      /* do final operations on results */
      for(i = 0, npleft = 0; i < NumPart; i++)
	{
	  if(P[i].Type == 0 || P[i].Type >= 4)
	    {
	      if(fof_nearest_distance[i] > 1.0e29)
		{
		  /* need to redo this particle */
		  npleft++;
		  fof_nearest_hsml[i] *= 2.0;
		  if(iter >= MAXITER - 10)
		    {
		      printf("i=%d task=%d ID=%d Hsml=%g  pos=(%g|%g|%g)\n",
			     i, ThisTask, P[i].ID, fof_nearest_hsml[i],
			     P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
		      fflush(stdout);
		    }
		}
	    }
	}

      MPI_Allreduce(&npleft, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      if(ntot > 0)
	{
	  iter++;
	  if(iter > 0 && ThisTask == 0)
	    {
	      printf("fof-nearest iteration %d: need to repeat for %d particles.\n", iter, ntot);
	      fflush(stdout);
	    }

	  if(iter > MAXITER)
	    {
	      printf("failed to converge in fof-nearest\n");
	      fflush(stdout);
	      endrun(1159);
	    }
	}
    }
  while(ntot > 0);
  myfree(nsend);
  myfree(nsend_local);
  myfree(nbuffer);
  myfree(noffset);
  myfree(fof_nearest_hsml);
  myfree(fof_nearest_distance);
  if(ThisTask == 0)
    {
      printf("\ndone finding nearest dm-particle\n");
      fflush(stdout);
    }
}


void fof_find_nearest_dmparticle_evaluate(int target, int mode)
{
  int j, n, index;
  int startnode, numngb_inbox;
  double h, r2max;
  double dx, dy, dz, r2;
  FLOAT pos[3];

  if(mode == 0)
    {
      pos[0] = P[target].Pos[0];
      pos[1] = P[target].Pos[1];
      pos[2] = P[target].Pos[2];
      h = fof_nearest_hsml[target];
    }
  else
    {
      pos[0] = FoFDataGet[target].Pos[0];
      pos[1] = FoFDataGet[target].Pos[1];
      pos[2] = FoFDataGet[target].Pos[2];
      h = FoFDataGet[target].Hsml;
    }

  startnode = All.MaxPart;
  index = -1;
  r2max = 1.0e30;
  do
    {
      numngb_inbox = ngb_treefind_darkmatter(&pos[0], h, &startnode);
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
	  if(r2 < r2max && r2 < h * h)
	    {
	      index = j;
	      r2max = r2;
	    }
	}
    }
  while(startnode >= 0);
  if(mode == 0)
    {
      if(index >= 0)
	{
	  fof_nearest_distance[target] = sqrt(r2max);
	  MinID[target] = MinID[Head[index]];
	  MinIDTask[target] = MinIDTask[Head[index]];
	}
    }
  else
    {
      if(index >= 0)
	{
	  FoFDataResult[target].Distance = sqrt(r2max);
	  FoFDataResult[target].MinID = MinID[Head[index]];
	  FoFDataResult[target].MinIDTask = MinIDTask[Head[index]];
	}
      else
	FoFDataResult[target].Distance = 2.0e30;
    }
}



/*! This routine is a comparison kernel used in a sort routine to group
 *  particles that are exported to the same processor.
 */
int fof_compare_key(const void *a, const void *b)
{
  if(((struct fofdata_in *) a)->Task < (((struct fofdata_in *) b)->Task))
    return -1;
  if(((struct fofdata_in *) a)->Task > (((struct fofdata_in *) b)->Task))
    return +1;
  return 0;
}

double fof_periodic(double x)
{
  if(x >= 0.5 * All.BoxSize)
    x -= All.BoxSize;
  if(x < -0.5 * All.BoxSize)
    x += All.BoxSize;
  return x;
}


double fof_periodic_wrap(double x)
{
  while(x >= All.BoxSize)
    x -= All.BoxSize;
  while(x < 0)
    x += All.BoxSize;
  return x;
}

#ifdef BUBBLES

int compare_length_values(const void *a, const void *b)
{
  if(((struct biggest_group *) a)->Len < (((struct biggest_group *) b)->Len))
    return -1;

  if(((struct biggest_group *) a)->Len > (((struct biggest_group *) b)->Len))
    return +1;

  return 0;
}


void find_CM_of_biggest_group(void)
{
  int i, j;
  int *BGroupLen, *BGroupLen_common;
  float *BGroupCM, *BGroupCM_common;
  double *BGroupMass, *BGroupMass_common;

  BGroupLen = mymalloc(sizeof(int));
  BGroupLen_common = mymalloc(NTask * sizeof(int));
  BGroupCM = mymalloc(3 * sizeof(float));
  BGroupCM_common = mymalloc(3 * NTask * sizeof(float));
  BGroupMass = mymalloc(6 * sizeof(double));
  BGroupMass_common = mymalloc(6 * NTask * sizeof(double));

  BGroupLen[0] = GroupLen[0];

  for(i = 0; i < 3; i++)
    BGroupCM[i] = GroupCM[0 * 3 + i];

  for(i = 0; i < 6; i++)
    BGroupMass[i] = GroupMassType[0 * 6 + i];

  MPI_Allgather(BGroupLen, 1, MPI_INT, BGroupLen_common, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgather(BGroupCM, 3, MPI_FLOAT, BGroupCM_common, 3, MPI_FLOAT, MPI_COMM_WORLD);
  MPI_Allgather(BGroupMass, 6, MPI_DOUBLE, BGroupMass_common, 6, MPI_DOUBLE, MPI_COMM_WORLD);


  for(i = 0; i < NTask; i++)
    {
      BiggestGroup[i].Len = BGroupLen_common[i];
      BiggestGroup[i].CM[0] = BGroupCM_common[i * 3 + 0];
      BiggestGroup[i].CM[1] = BGroupCM_common[i * 3 + 1];
      BiggestGroup[i].CM[2] = BGroupCM_common[i * 3 + 2];
      for(j = 0; j < 6; j++)
	BiggestGroup[i].Mass[j] = BGroupMass_common[i * 6 + j];
    }

  qsort(BiggestGroup, NTask, sizeof(struct biggest_group), compare_length_values);

  All.BiggestGroupLen = BiggestGroup[NTask - 1].Len;

  for(i = 0; i < 3; i++)
    All.BiggestGroupCM[i] = BiggestGroup[NTask - 1].CM[i];

  All.BiggestGroupMass = BiggestGroup[NTask - 1].Mass[0] + BiggestGroup[NTask - 1].Mass[1] + BiggestGroup[NTask - 1].Mass[4];	/* added stars mass */

  if(ThisTask == 0)
    {
      printf("Biggest group length has %d particles.\n", All.BiggestGroupLen);
      printf("CM of biggest group is: (%g|%g|%g)\n", All.BiggestGroupCM[0], All.BiggestGroupCM[1],
	     All.BiggestGroupCM[2]);
      printf("Mass of biggest group is: %g\n", All.BiggestGroupMass);
    }

  myfree(BGroupLen);
  myfree(BGroupLen_common);
  myfree(BGroupCM);
  myfree(BGroupCM_common);
  myfree(BGroupMass);
  myfree(BGroupMass_common);

}

#endif



#ifdef MULTI_BUBBLES

void multi_bubbles(void)
{
  double phi, theta;
  double dx, dy, dz, rr, r2, dE;
  double E_bubble, totE_bubble, hubble_a = 0.0;
  double BubbleDistance, BubbleRadius, BubbleEnergy;
  FLOAT Mass_bubble, totMass_bubble;
  FLOAT pos[3];
  int numngb, tot_numngb, startnode, numngb_inbox;
  int n, i, j, k, l;
  int nheat, tot_nheat;
  int eff_nheat, tot_eff_nheat;
  double *GroupMassType_common, *GroupMassType_dum;
  float *GroupCM_common_x, *GroupCM_dum_x;
  float *GroupCM_common_y, *GroupCM_dum_y;
  float *GroupCM_common_z, *GroupCM_dum_z;
  int logical;
  int *nn_heat, *disp;
  double massDMpart;

  if(All.MassTable[1] > 0)
    massDMpart = All.MassTable[1];
  else
    massDMpart = All.massDMpart;


  if(All.ComovingIntegrationOn)
    {
      hubble_a = hubble_function(All.Time) / All.Hubble;
    }

  nheat = tot_nheat = 0;
  eff_nheat = tot_eff_nheat = 0;

  logical = 0;

  for(k = 0; k < Ngroups; k++)
    {
      if(massDMpart > 0)
	{
	  if(GroupLenType[k * 6 + 1] * massDMpart >=
	     (All.Omega0 - All.OmegaBaryon) / All.Omega0 * All.MinFoFMassForNewSeed)
	    nheat++;
	}
      else
	{
	  printf("The DM particles mass is zero! I will stop.\n");
	  endrun(0);
	}

    }

  MPI_Allreduce(&nheat, &tot_nheat, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(ThisTask == 0)
    printf("The total number of clusters to heat is: %d\n", tot_nheat);


  nn_heat = mymalloc(NTask * sizeof(int));
  disp = mymalloc(NTask * sizeof(int));

  MPI_Allgather(&nheat, 1, MPI_INT, nn_heat, 1, MPI_INT, MPI_COMM_WORLD);

  for(k = 1, disp[0] = 0; k < NTask; k++)
    disp[k] = disp[k - 1] + nn_heat[k - 1];


  if(tot_nheat > 0)
    {
      GroupMassType_common = mymalloc(tot_nheat * sizeof(double));
      GroupMassType_dum = mymalloc(nheat * sizeof(double));

      GroupCM_common_x = mymalloc(tot_nheat * sizeof(float));
      GroupCM_dum_x = mymalloc(nheat * sizeof(float));

      GroupCM_common_y = mymalloc(tot_nheat * sizeof(float));
      GroupCM_dum_y = mymalloc(nheat * sizeof(float));

      GroupCM_common_z = mymalloc(tot_nheat * sizeof(float));
      GroupCM_dum_z = mymalloc(nheat * sizeof(float));


      for(k = 0, i = 0; k < Ngroups; k++)
	{
	  if(massDMpart > 0)
	    {
	      if(GroupLenType[k * 6 + 1] * massDMpart >=
		 (All.Omega0 - All.OmegaBaryon) / All.Omega0 * All.MinFoFMassForNewSeed)
		{
		  GroupCM_dum_x[i] = GroupCM[k * 3 + 0];
		  GroupCM_dum_y[i] = GroupCM[k * 3 + 1];
		  GroupCM_dum_z[i] = GroupCM[k * 3 + 2];

		  GroupMassType_dum[i] = GroupMassType[k * 6 + 0] + GroupMassType[k * 6 + 1] + GroupMassType[k * 6 + 4];	/* gas+dm+stars */

		  i++;
		}
	    }
	  else
	    {
	      printf("The DM particles mass is zero! I will stop.\n");
	      endrun(0);
	    }
	}

      MPI_Allgatherv(GroupMassType_dum, nheat, MPI_DOUBLE, GroupMassType_common, nn_heat, disp, MPI_DOUBLE,
		     MPI_COMM_WORLD);

      MPI_Allgatherv(GroupCM_dum_x, nheat, MPI_FLOAT, GroupCM_common_x, nn_heat, disp, MPI_FLOAT,
		     MPI_COMM_WORLD);

      MPI_Allgatherv(GroupCM_dum_y, nheat, MPI_FLOAT, GroupCM_common_y, nn_heat, disp, MPI_FLOAT,
		     MPI_COMM_WORLD);

      MPI_Allgatherv(GroupCM_dum_z, nheat, MPI_FLOAT, GroupCM_common_z, nn_heat, disp, MPI_FLOAT,
		     MPI_COMM_WORLD);


      for(l = 0; l < tot_nheat; l++)
	{
	  if(All.ComovingIntegrationOn > 0)
	    {
	      BubbleDistance =
		All.BubbleDistance * 1. / All.Time * pow(GroupMassType_common[l] / All.ClusterMass200,
							 1. / 3.) / pow(hubble_a, 2. / 3.);

	      BubbleRadius =
		All.BubbleRadius * 1. / All.Time * pow(GroupMassType_common[l] / All.ClusterMass200,
						       1. / 3.) / pow(hubble_a, 2. / 3.);

	      BubbleEnergy =
		All.BubbleEnergy * pow(GroupMassType_common[l] / All.ClusterMass200, 5. / 3.) * pow(hubble_a,
												    2. / 3.);

	      phi = theta = rr = 0.0;

	      phi = 2 * M_PI * get_random_number(0);
	      theta = acos(2 * get_random_number(0) - 1);
	      rr = pow(get_random_number(0), 1. / 3.) * BubbleDistance;

	      pos[0] = pos[1] = pos[2] = 0.0;

	      pos[0] = sin(theta) * cos(phi);
	      pos[1] = sin(theta) * sin(phi);
	      pos[2] = cos(theta);

	      for(k = 0; k < 3; k++)
		pos[k] *= rr;

	      pos[0] += GroupCM_common_x[l];
	      pos[1] += GroupCM_common_y[l];
	      pos[2] += GroupCM_common_z[l];


	      /* First, let's see how many particles are in the bubble */

	      numngb = 0;
	      E_bubble = 0.;
	      Mass_bubble = 0.;

	      startnode = All.MaxPart;
	      do
		{
		  numngb_inbox = ngb_treefind_variable(pos, BubbleRadius, &startnode);

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

		      if(r2 < BubbleRadius * BubbleRadius)
			{
			  numngb++;

			  if(All.ComovingIntegrationOn)
			    E_bubble +=
			      SphP[j].Entropy * P[j].Mass * pow(SphP[j].a2.Density / pow(All.Time, 3),
								GAMMA_MINUS1) / GAMMA_MINUS1;
			  else
			    E_bubble +=
			      SphP[j].Entropy * P[j].Mass * pow(SphP[j].a2.Density,
								GAMMA_MINUS1) / GAMMA_MINUS1;

			  Mass_bubble += P[j].Mass;

			}
		    }
		}
	      while(startnode >= 0);


	      tot_numngb = totE_bubble = totMass_bubble = 0.0;

	      MPI_Allreduce(&numngb, &tot_numngb, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	      MPI_Allreduce(&E_bubble, &totE_bubble, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	      MPI_Allreduce(&Mass_bubble, &totMass_bubble, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);


	      if(tot_numngb == 0)
		{
		  tot_numngb = 1;
		  totMass_bubble = 1;
		  totE_bubble = 1;
		  logical = 0;
		}
	      else
		{
		  logical = 1;
		  eff_nheat++;
		}

	      totE_bubble *= All.UnitEnergy_in_cgs;


	      if(ThisTask == 0)
		{
		  if(logical == 1)
		    {
		      printf("%g, %g, %g, %g, %d, %g, %g, %g\n", GroupMassType_common[l], GroupCM_common_x[l],
			     GroupCM_common_y[l], GroupCM_common_z[l], tot_numngb, BubbleRadius, BubbleEnergy,
			     (BubbleEnergy + totE_bubble) / totE_bubble);

		    }
		  fflush(stdout);
		}

	      /* now find particles in Bubble again, and inject energy */

	      startnode = All.MaxPart;

	      do
		{
		  numngb_inbox = ngb_treefind_variable(pos, BubbleRadius, &startnode);

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

		      if(r2 < BubbleRadius * BubbleRadius)
			{
			  /* with sf on gas particles have mass that is not fixed */
			  /* energy we want to inject in this particle */

			  if(logical == 1)
			    dE = ((BubbleEnergy / All.UnitEnergy_in_cgs) / totMass_bubble) * P[j].Mass;
			  else
			    dE = 0;


			  if(All.ComovingIntegrationOn)
			    SphP[j].Entropy +=
			      GAMMA_MINUS1 * dE / P[j].Mass / pow(SphP[j].a2.Density / pow(All.Time, 3),
								  GAMMA_MINUS1);
			  else
			    SphP[j].Entropy +=
			      GAMMA_MINUS1 * dE / P[j].Mass / pow(SphP[j].a2.Density, GAMMA_MINUS1);
			  if(dE > 0 && P[j].ID > 0)
			    P[j].ID = -P[j].ID;

			}
		    }
		}
	      while(startnode >= 0);
	    }
	}

    }
  else
    {
      printf("There are no clusters to heat. I will stop.\n");
      endrun(0);
    }

  if(tot_nheat > 0)
    {
      myfree(GroupCM_dum_z);
      myfree(GroupCM_common_z);
      myfree(GroupCM_dum_y); 
      myfree(GroupCM_common_y);
      myfree(GroupCM_dum_x);
      myfree(GroupCM_common_x);
      myfree(GroupMassType_dum);
      myfree(GroupMassType_common);

      MPI_Allreduce(&eff_nheat, &tot_eff_nheat, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      if(ThisTask == 0)
	printf("The total effective! number of clusters to heat is: %d\n", eff_nheat);
    }

  myfree(disp);
  myfree(nn_heat);

}

#endif




#endif /* of FOF */
