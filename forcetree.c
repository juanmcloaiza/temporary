#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"


/*! \file forcetree.c
 *  \brief gravitational tree and code for Ewald correction
 *
 *  This file contains the computation of the gravitational force by means
 *  of a tree. The type of tree implemented is a geometrical oct-tree,
 *  starting from a cube encompassing all particles. This cube is
 *  automatically found in the domain decomposition, which also splits up
 *  the global "top-level" tree along node boundaries, moving the particles
 *  of different parts of the tree to separate processors. Tree nodes can
 *  be dynamically updated in drift/kick operations to avoid having to
 *  reconstruct the tree every timestep.
 */

/*! auxialiary variable used to set-up non-recursive walk */
static int last;



/*! length of lock-up table for short-range force kernel in TreePM algorithm */
#define NTAB 1000
/*! variables for short-range lookup table */
static float tabfac, shortrange_table[NTAB], shortrange_table_potential[NTAB];

/*! toggles after first tree-memory allocation, has only influence on log-files */
static int first_flag = 0;




#ifdef PERIODIC
/*! Macro that maps a distance to the nearest periodic neighbour */
#define NEAREST(x) (((x)>boxhalf)?((x)-boxsize):(((x)<-boxhalf)?((x)+boxsize):(x)))
/*! Size of 3D lock-up table for Ewald correction force */
#define EN  64
/*! 3D lock-up table for Ewald correction to force and potential. Only one
 *  octant is stored, the rest constructed by using the symmetry
 */
static FLOAT fcorrx[EN + 1][EN + 1][EN + 1];
static FLOAT fcorry[EN + 1][EN + 1][EN + 1];
static FLOAT fcorrz[EN + 1][EN + 1][EN + 1];
static FLOAT potcorr[EN + 1][EN + 1][EN + 1];
static double fac_intp;
#endif



/*! This function is a driver routine for constructing the gravitational
 *  oct-tree, which is done by calling a small number of other functions.
 */
int force_treebuild(int npart)
{
  double t0, t1;

  t0 = second();

  Numnodestree = force_treebuild_single(npart);
  force_flag_localnodes();

  t1 = second();
  CPU_Step[CPU_TREEBUILD] += timediff(t0, t1);
  All.Cadj_Cpu += timediff(t0, t1);


  force_exchange_pseudodata();


  t0 = second();

  force_treeupdate_pseudos(All.MaxPart);

  TimeOfLastTreeConstruction = All.Time;

  t1 = second();
  CPU_Step[CPU_TREEBUILD] += timediff(t0, t1);
  All.Cadj_Cpu += timediff(t0, t1);

  return Numnodestree;
}



/*! Constructs the gravitational oct-tree.
 *
 *  The index convention for accessing tree nodes is the following: the
 *  indices 0...NumPart-1 reference single particles, the indices
 *  All.MaxPart.... All.MaxPart+nodes-1 reference tree nodes. `Nodes_base'
 *  points to the first tree node, while `nodes' is shifted such that
 *  nodes[All.MaxPart] gives the first tree node. Finally, node indices
 *  with values 'All.MaxPart + MaxNodes' and larger indicate "pseudo
 *  particles", i.e. multipole moments of top-level nodes that lie on
 *  different CPUs. If such a node needs to be opened, the corresponding
 *  particle must be exported to that CPU. The 'Extnodes' structure
 *  parallels that of 'Nodes'. Its information is only needed for the SPH
 *  part of the computation. (The data is split onto these two structures
 *  as a tuning measure.  If it is merged into 'Nodes' a somewhat bigger
 *  size of the nodes also for gravity would result, which would reduce
 *  cache utilization slightly.
 */
int force_treebuild_single(int npart)
{
  int i, j, subnode = 0, shift, parent, numnodes, rep;
  int nfree, th, nn, no;
  struct NODE *nfreep;
  FLOAT lenhalf;
  FLOAT minbound[3], maxbound[3];
  peanokey key, morton, th_key, *morton_list;


  /* create an empty root node  */
  nfree = All.MaxPart;		/* index of first free node */
  nfreep = &Nodes[nfree];	/* select first node */

  nfreep->len = DomainLen;
  for(j = 0; j < 3; j++)
    nfreep->center[j] = DomainCenter[j];
  for(j = 0; j < 8; j++)
    nfreep->u.suns[j] = -1;


  numnodes = 1;
  nfreep++;
  nfree++;

  /* create a set of empty nodes corresponding to the top-level domain
   * grid. We need to generate these nodes first to make sure that we have a
   * complete top-level tree which allows the easy insertion of the
   * pseudo-particles at the right place
   */

  force_create_empty_nodes(All.MaxPart, 0, 1, 0, 0, 0, &numnodes, &nfree);


  /* if a high-resolution region in a global tree is used, we need to generate
   * an additional set empty nodes to make sure that we have a complete
   * top-level tree for the high-resolution inset
   */

  nfreep = &Nodes[nfree];
  parent = -1;			/* note: will not be used below before it is changed */

  morton_list = (peanokey *) mymalloc(npart * sizeof(peanokey));

  /* now we insert all particles */
  for(i = 0; i < npart; i++)
    {
      rep = 0;

      key = peano_and_morton_key((int) ((P[i].Pos[0] - DomainCorner[0]) * DomainFac),
				 (int) ((P[i].Pos[1] - DomainCorner[1]) * DomainFac),
				 (int) ((P[i].Pos[2] - DomainCorner[2]) * DomainFac), BITS_PER_DIMENSION,
				 &morton);
      morton_list[i] = morton;

      shift = 3 * (BITS_PER_DIMENSION - 1);

      no = 0;
      while(TopNodes[no].Daughter >= 0)
	{
	  no = TopNodes[no].Daughter + (key - TopNodes[no].StartKey) / (TopNodes[no].Size / 8);
	  shift -= 3;
	  rep++;
	}

      no = TopNodes[no].Leaf;
      th = DomainNodeIndex[no];

      while(1)
	{
	  if(th >= All.MaxPart)	/* we are dealing with an internal node */
	    {
	      if(shift >= 0)
		{
		  subnode = ((morton >> shift) & 7);
		}
	      else
		{
		  subnode = 0;
		  if(P[i].Pos[0] > Nodes[th].center[0])
		    subnode += 1;
		  if(P[i].Pos[1] > Nodes[th].center[1])
		    subnode += 2;
		  if(P[i].Pos[2] > Nodes[th].center[2])
		    subnode += 4;
		}

#ifndef NOTREERND
	      if(Nodes[th].len < 1.0e-3 * All.ForceSoftening[P[i].Type])
		{
		  /* seems like we're dealing with particles at identical (or extremely close)
		   * locations. Randomize subnode index to allow tree construction. Note: Multipole moments
		   * of tree are still correct, but this will only happen well below gravitational softening
		   * length-scale anyway.
		   */
		  subnode = (int) (8.0 * get_random_number((P[i].ID + rep) % (RNDTABLE + (rep & 3))));

		  if(subnode >= 8)
		    subnode = 7;
		}
#endif

	      nn = Nodes[th].u.suns[subnode];

	      shift -= 3;

	      if(nn >= 0)	/* ok, something is in the daughter slot already, need to continue */
		{
		  parent = th;
		  th = nn;
		  rep++;
		}
	      else
		{
		  /* here we have found an empty slot where we can attach
		   * the new particle as a leaf.
		   */
		  Nodes[th].u.suns[subnode] = i;
		  break;	/* done for this particle */
		}
	    }
	  else
	    {
	      /* We try to insert into a leaf with a single particle.  Need
	       * to generate a new internal node at this point.
	       */
	      Nodes[parent].u.suns[subnode] = nfree;

	      nfreep->len = 0.5 * Nodes[parent].len;
	      lenhalf = 0.25 * Nodes[parent].len;

	      if(subnode & 1)
		nfreep->center[0] = Nodes[parent].center[0] + lenhalf;
	      else
		nfreep->center[0] = Nodes[parent].center[0] - lenhalf;

	      if(subnode & 2)
		nfreep->center[1] = Nodes[parent].center[1] + lenhalf;
	      else
		nfreep->center[1] = Nodes[parent].center[1] - lenhalf;

	      if(subnode & 4)
		nfreep->center[2] = Nodes[parent].center[2] + lenhalf;
	      else
		nfreep->center[2] = Nodes[parent].center[2] - lenhalf;

	      nfreep->u.suns[0] = -1;
	      nfreep->u.suns[1] = -1;
	      nfreep->u.suns[2] = -1;
	      nfreep->u.suns[3] = -1;
	      nfreep->u.suns[4] = -1;
	      nfreep->u.suns[5] = -1;
	      nfreep->u.suns[6] = -1;
	      nfreep->u.suns[7] = -1;

	      if(shift >= 0)
		{
		  th_key = morton_list[th];
		  subnode = ((th_key >> shift) & 7);
		}
	      else
		{
		  subnode = 0;
		  if(P[th].Pos[0] > nfreep->center[0])
		    subnode += 1;
		  if(P[th].Pos[1] > nfreep->center[1])
		    subnode += 2;
		  if(P[th].Pos[2] > nfreep->center[2])
		    subnode += 4;
		}

#ifndef NOTREERND
	      if(nfreep->len < 1.0e-3 * All.ForceSoftening[P[th].Type])
		{
		  /* seems like we're dealing with particles at identical (or extremely close)
		   * locations. Randomize subnode index to allow tree construction. Note: Multipole moments
		   * of tree are still correct, but this will only happen well below gravitational softening
		   * length-scale anyway.
		   */
		  subnode = (int) (8.0 * get_random_number((P[th].ID + rep) % (RNDTABLE + (rep & 3))));

		  if(subnode >= 8)
		    subnode = 7;
		}
#endif
	      nfreep->u.suns[subnode] = th;

	      th = nfree;	/* resume trying to insert the new particle at
				 * the newly created internal node
				 */

	      numnodes++;
	      nfree++;
	      nfreep++;

	      if((numnodes) >= MaxNodes)
		{
		  printf("task %d: maximum number %d of tree-nodes reached.\n", ThisTask, MaxNodes);
		  printf("for particle %d\n", i);
		  dump_particles();
		  endrun(1);
		}
	    }
	}
    }

  myfree(morton_list);


  /* insert the pseudo particles that represent the mass distribution of other domains */
  force_insert_pseudo_particles();


  /* now compute the multipole moments recursively */
  last = -1;

  force_update_node_recursive(All.MaxPart, -1, -1, minbound, maxbound);

  if(last >= All.MaxPart)
    {
      if(last >= All.MaxPart + MaxNodes)	/* a pseudo-particle */
	Nextnode[last - MaxNodes] = -1;
      else
	Nodes[last].u.d.nextnode = -1;
    }
  else
    Nextnode[last] = -1;

  return numnodes;
}



/*! This function recursively creates a set of empty tree nodes which
 *  corresponds to the top-level tree for the domain grid. This is done to
 *  ensure that this top-level tree is always "complete" so that we can easily
 *  associate the pseudo-particles of other CPUs with tree-nodes at a given
 *  level in the tree, even when the particle population is so sparse that
 *  some of these nodes are actually empty.
*/
void force_create_empty_nodes(int no, int topnode, int bits, int x, int y, int z, int *nodecount,
			      int *nextfree)
{
  int i, j, k, n, sub, count;
  FLOAT lenhalf;

  if(TopNodes[topnode].Daughter >= 0)
    {
      for(i = 0; i < 2; i++)
	for(j = 0; j < 2; j++)
	  for(k = 0; k < 2; k++)
	    {
	      sub = 7 & peano_hilbert_key((x << 1) + i, (y << 1) + j, (z << 1) + k, bits);

	      count = i + 2 * j + 4 * k;

	      Nodes[no].u.suns[count] = *nextfree;

	      lenhalf = 0.25 * Nodes[no].len;
	      Nodes[*nextfree].len = 0.5 * Nodes[no].len;
	      Nodes[*nextfree].center[0] = Nodes[no].center[0] + (2 * i - 1) * lenhalf;
	      Nodes[*nextfree].center[1] = Nodes[no].center[1] + (2 * j - 1) * lenhalf;
	      Nodes[*nextfree].center[2] = Nodes[no].center[2] + (2 * k - 1) * lenhalf;

	      for(n = 0; n < 8; n++)
		Nodes[*nextfree].u.suns[n] = -1;

	      if(TopNodes[TopNodes[topnode].Daughter + sub].Daughter == -1)
		DomainNodeIndex[TopNodes[TopNodes[topnode].Daughter + sub].Leaf] = *nextfree;

	      *nextfree = *nextfree + 1;
	      *nodecount = *nodecount + 1;

	      if((*nodecount) >= MaxNodes || (*nodecount) >= MAXTOPNODES)
		{
		  printf("task %d: maximum number %d/%d of tree-nodes reached.\n", ThisTask, MaxNodes,
			 MAXTOPNODES);
		  printf("in create empty nodes\n");
		  dump_particles();
		  endrun(11);
		}

	      force_create_empty_nodes(*nextfree - 1, TopNodes[topnode].Daughter + sub,
				       bits + 1, 2 * x + i, 2 * y + j, 2 * z + k, nodecount, nextfree);
	    }
    }
}



/*! this function inserts pseudo-particles which will represent the mass
 *  distribution of the other CPUs. Initially, the mass of the
 *  pseudo-particles is set to zero, and their coordinate is set to the center
 *  of the domain-cell they correspond to. These quantities will be updated
 *  later on.
 */
void force_insert_pseudo_particles(void)
{
  int i, index, subnode, nn, th;

  for(i = 0; i < NTopleaves; i++)
    {
      index = DomainNodeIndex[i];

      DomainMoment[i].mass = 0;
      DomainMoment[i].s[0] = Nodes[index].center[0];
      DomainMoment[i].s[1] = Nodes[index].center[1];
      DomainMoment[i].s[2] = Nodes[index].center[2];
    }

  for(i = 0; i < NTopleaves; i++)
    {
      if(i < DomainMyStart || i > DomainMyLast)
	{
	  th = All.MaxPart;	/* select index of first node in tree */

	  while(1)
	    {
	      if(th >= All.MaxPart)	/* we are dealing with an internal node */
		{
		  if(th >= All.MaxPart + MaxNodes)
		    endrun(888);	/* this can't be */

		  subnode = 0;
		  if(DomainMoment[i].s[0] > Nodes[th].center[0])
		    subnode += 1;
		  if(DomainMoment[i].s[1] > Nodes[th].center[1])
		    subnode += 2;
		  if(DomainMoment[i].s[2] > Nodes[th].center[2])
		    subnode += 4;

		  nn = Nodes[th].u.suns[subnode];

		  if(nn >= 0)	/* ok, something is in the daughter slot already, need to continue */
		    {
		      th = nn;
		    }
		  else
		    {
		      /* here we have found an empty slot where we can
		       * attach the pseudo particle as a leaf
		       */
		      Nodes[th].u.suns[subnode] = All.MaxPart + MaxNodes + i;

		      break;	/* done for this pseudo particle */
		    }
		}
	      else
		{
		  endrun(889);	/* this can't be */
		}
	    }
	}
    }
}


/*! this routine determines the multipole moments for a given internal node
 *  and all its subnodes using a recursive computation.  The result is
 *  stored in the Nodes[] structure in the sequence of this tree-walk.
 *
 *  Note that the bitflags-variable for each node is used to store in the
 *  lowest bits some special information: Bit 0 flags whether the node
 *  belongs to the top-level tree corresponding to the domain
 *  decomposition, while Bit 1 signals whether the top-level node is
 *  dependent on local mass.
 *
 *  If UNEQUALSOFTENINGS is set, bits 2-4 give the particle type with
 *  the maximum softening among the particles in the node, and bit 5
 *  flags whether the node contains any particles with lower softening
 *  than that.
 */
void force_update_node_recursive(int no, int sib, int father, FLOAT * minbound, FLOAT * maxbound)
{
  int j, jj, k, p, pp, nextsib, suns[8], count_daughters, count_particles, multiple_flag, flag_single;
  FLOAT hmax, lenmax;
  FLOAT single_s[3], single_mass = 0;
  FLOAT s[3], mass;
  FLOAT mindaughter[3], maxdaughter[3];

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
  int maxsofttype, diffsoftflag;
#else
  FLOAT maxsoft;
#endif
#endif
  struct particle_data *pa;


  if(no >= All.MaxPart && no < All.MaxPart + MaxNodes)	/* internal node */
    {
      for(k = 0; k < 3; k++)
	{
	  minbound[k] = MAX_REAL_NUMBER;
	  maxbound[k] = -MAX_REAL_NUMBER;
	}

      for(j = 0; j < 8; j++)
	suns[j] = Nodes[no].u.suns[j];	/* this "backup" is necessary because the nextnode entry will
					   overwrite one element (union!) */
      if(last >= 0)
	{
	  if(last >= All.MaxPart)
	    {
	      if(last >= All.MaxPart + MaxNodes)	/* a pseudo-particle */
		Nextnode[last - MaxNodes] = no;
	      else
		Nodes[last].u.d.nextnode = no;
	    }
	  else
	    Nextnode[last] = no;
	}

      last = no;

      mass = 0;
      s[0] = 0;
      s[1] = 0;
      s[2] = 0;
      hmax = 0;
      count_daughters = 0;
      count_particles = 0;
      flag_single = 0;

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
      maxsofttype = 7;
      diffsoftflag = 0;
#else
      maxsoft = 0;
#endif
#endif

      for(j = 0; j < 8; j++)
	{
	  if((p = suns[j]) >= 0)
	    {
	      /* check if we have a sibling on the same level */
	      for(jj = j + 1; jj < 8; jj++)
		if((pp = suns[jj]) >= 0)
		  break;

	      if(jj < 8)	/* yes, we do */
		nextsib = pp;
	      else
		nextsib = sib;

	      force_update_node_recursive(p, nextsib, no, mindaughter, maxdaughter);

	      if(p >= All.MaxPart)	/* an internal node or pseudo particle */
		{
		  if(p >= All.MaxPart + MaxNodes)	/* a pseudo particle */
		    {
		      /* nothing to be done here because the mass of the
		       * pseudo-particle is still zero. This will be changed
		       * later.
		       */
		    }
		  else
		    {
		      count_daughters++;

		      if(flag_single == 0)
			{
			  single_s[0] = Nodes[p].u.d.s[0];
			  single_s[1] = Nodes[p].u.d.s[1];
			  single_s[2] = Nodes[p].u.d.s[2];
			  single_mass = Nodes[p].u.d.mass;
			}

		      mass += (Nodes[p].u.d.mass);
		      s[0] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[0]);
		      s[1] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[1]);
		      s[2] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[2]);

		      if(Nodes[p].u.d.mass > 0)
			{
			  flag_single++;

			  if((Nodes[p].u.d.bitflags & 128))
			    count_particles += 2;
			  else
			    count_particles++;
			}


		      if(Extnodes[p].hmax > hmax)
			hmax = Extnodes[p].hmax;

		      for(k = 0; k < 3; k++)
			{
			  if(minbound[k] > mindaughter[k])
			    minbound[k] = mindaughter[k];

			  if(maxbound[k] < maxdaughter[k])
			    maxbound[k] = maxdaughter[k];
			}

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
		      diffsoftflag |= (Nodes[p].u.d.bitflags >> 5) & 1;

		      if(maxsofttype == 7)
			{
			  maxsofttype = (Nodes[p].u.d.bitflags >> 2) & 7;
			}
		      else
			{
			  if(((Nodes[p].u.d.bitflags >> 2) & 7) != 7)
			    {
			      if(All.ForceSoftening[((Nodes[p].u.d.bitflags >> 2) & 7)] >
				 All.ForceSoftening[maxsofttype])
				{
				  maxsofttype = ((Nodes[p].u.d.bitflags >> 2) & 7);
				  diffsoftflag = 1;
				}
			      else
				{
				  if(All.ForceSoftening[((Nodes[p].u.d.bitflags >> 2) & 7)] <
				     All.ForceSoftening[maxsofttype])
				    diffsoftflag = 1;
				}
			    }
			}
#else
		      if(Nodes[p].maxsoft > maxsoft)
			maxsoft = Nodes[p].maxsoft;
#endif
#endif
		    }
		}
	      else		/* a particle */
		{
		  count_daughters++;

		  pa = &P[p];
		  if(flag_single == 0)
		    {
		      single_s[0] = pa->Pos[0];
		      single_s[1] = pa->Pos[1];
		      single_s[2] = pa->Pos[2];
		      single_mass = pa->Mass;
		    }
		  flag_single++;

		  mass += (pa->Mass);
		  s[0] += (pa->Mass * pa->Pos[0]);
		  s[1] += (pa->Mass * pa->Pos[1]);
		  s[2] += (pa->Mass * pa->Pos[2]);

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
		  if(maxsofttype == 7)
		    {
		      maxsofttype = pa->Type;
		    }
		  else
		    {
		      if(All.ForceSoftening[pa->Type] > All.ForceSoftening[maxsofttype])
			{
			  maxsofttype = pa->Type;
			  diffsoftflag = 1;
			}
		      else
			{
			  if(All.ForceSoftening[pa->Type] < All.ForceSoftening[maxsofttype])
			    diffsoftflag = 1;
			}
		    }
#else
		  if(pa->Type == 0)
		    {
#ifdef ADAPTIVE_GRAVSOFT_FORGAS_HSML
		      if(dmin(All.SofteningTable[P[i].Type],ADAPTIVE_GRAVSOFT_FORGAS_HSML * PPP[i].Hsml) > maxsoft)
			maxsoft = dmin(All.SofteningTable[P[i].Type],ADAPTIVE_GRAVSOFT_FORGAS_HSML * PPP[i].Hsml);
#else
		      if(All.ForceSoftening[0] * pow(pa->Mass / All.ReferenceGasMass, 1.0 / 3) > maxsoft)
			maxsoft = All.ForceSoftening[0] * pow(pa->Mass / All.ReferenceGasMass, 1.0 / 3);
#endif
		    }
		  else
		    {
		      if(All.ForceSoftening[pa->Type] > maxsoft)
			maxsoft = All.ForceSoftening[pa->Type];
		    }
#endif
#endif
		  if(pa->Type == 0)
		    if(PPP[p].Hsml > hmax)
		      hmax = PPP[p].Hsml;

		  for(k = 0; k < 3; k++)
		    {
		      if(minbound[k] > pa->Pos[k])
			minbound[k] = pa->Pos[k];

		      if(maxbound[k] < pa->Pos[k])
			maxbound[k] = pa->Pos[k];
		    }

		  count_particles++;
		}
	    }
	}


      if(mass)
	{
	  s[0] /= mass;
	  s[1] /= mass;
	  s[2] /= mass;
	}
      else
	{
	  s[0] = Nodes[no].center[0];
	  s[1] = Nodes[no].center[1];
	  s[2] = Nodes[no].center[2];
	}

      if(flag_single == 1)
	{
	  Nodes[no].u.d.s[0] = single_s[0];
	  Nodes[no].u.d.s[1] = single_s[1];
	  Nodes[no].u.d.s[2] = single_s[2];
	  Nodes[no].u.d.mass = single_mass;
	}
      else
	{
	  Nodes[no].u.d.s[0] = s[0];
	  Nodes[no].u.d.s[1] = s[1];
	  Nodes[no].u.d.s[2] = s[2];
	  Nodes[no].u.d.mass = mass;
	}

      for(k = 0, lenmax = Nodes[no].len; k < 3; k++)
	{
	  if(lenmax < 2.0 * (Nodes[no].center[k] - minbound[k]))
	    lenmax = 2.0 * (Nodes[no].center[k] - minbound[k]);
	  if(lenmax < 2.0 * (maxbound[k] - Nodes[no].center[k]))
	    lenmax = 2.0 * (maxbound[k] - Nodes[no].center[k]);
	}

      Nodes[no].len = lenmax;

      Extnodes[no].hmax = hmax;

      if(count_particles > 1)	/* this flags that the node represents more than one particle */
	multiple_flag = 1;
      else
	multiple_flag = 0;

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
      Nodes[no].u.d.bitflags =
	4 * maxsofttype + 32 * diffsoftflag + 128 * multiple_flag + 256 * count_daughters;
#else
      Nodes[no].u.d.bitflags = 128 * multiple_flag + 256 * count_daughters;
      Nodes[no].maxsoft = maxsoft;
#endif
#else
      Nodes[no].u.d.bitflags = 128 * multiple_flag + 256 * count_daughters;
#endif

      Nodes[no].u.d.sibling = sib;
      Nodes[no].u.d.father = father;
    }
  else				/* single particle or pseudo particle */
    {
      if(last >= 0)
	{
	  if(last >= All.MaxPart)
	    {
	      if(last >= All.MaxPart + MaxNodes)	/* a pseudo-particle */
		Nextnode[last - MaxNodes] = no;
	      else
		Nodes[last].u.d.nextnode = no;
	    }
	  else
	    Nextnode[last] = no;
	}

      last = no;

      if(no < All.MaxPart)	/* only set it for single particles */
	Father[no] = father;
    }

}





/*! This function communicates the values of the multipole moments of the
 *  top-level tree-nodes of the domain grid.  This data can then be used to
 *  update the pseudo-particles on each CPU accordingly.
 */
void force_exchange_pseudodata(void)
{
  int i, no, recvTask;
  int *recvcounts, *recvoffset;
  double t0, t1;

  t0 = second();

  for(i = DomainMyStart; i <= DomainMyLast; i++)
    {
      no = DomainNodeIndex[i];

      /* read out the multipole moments from the local base cells */
      DomainMoment[i].s[0] = Nodes[no].u.d.s[0];
      DomainMoment[i].s[1] = Nodes[no].u.d.s[1];
      DomainMoment[i].s[2] = Nodes[no].u.d.s[2];
      DomainMoment[i].mass = Nodes[no].u.d.mass;
      DomainMoment[i].len = Nodes[no].len;
      DomainMoment[i].bitflags = Nodes[no].u.d.bitflags;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
      DomainMoment[i].maxsoft = Nodes[no].maxsoft;
#endif
    }

  /* share the pseudo-particle data accross CPUs */

    recvcounts = (int *)mymalloc(sizeof(int) * NTask);
    recvoffset = (int *)mymalloc(sizeof(int) * NTask);

  for(recvTask = 0; recvTask < NTask; recvTask++)
    {
      recvcounts[recvTask] =
	(DomainEndList[recvTask] - DomainStartList[recvTask] + 1) * sizeof(struct DomainNODE);
      recvoffset[recvTask] = DomainStartList[recvTask] * sizeof(struct DomainNODE);
    }

  t1 = second();
  CPU_Step[CPU_TREEBUILD] += timediff(t0, t1);
  All.Cadj_Cpu += timediff(t0, t1);
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = second();


  MPI_Allgatherv(&DomainMoment[DomainStartList[ThisTask]], recvcounts[ThisTask], MPI_BYTE,
		 &DomainMoment[0], recvcounts, recvoffset, MPI_BYTE, MPI_COMM_WORLD);


  t1 = second();
  CPU_Step[CPU_TREECOMM] += timediff(t0, t1);
  t0 = second();


  myfree(recvoffset);
  myfree(recvcounts);

  for(i = 0; i <= DomainEndList[NTask - 1]; i++)
    {
      if(i < DomainMyStart || i > DomainMyLast)
	{
	  no = DomainNodeIndex[i];

	  Nodes[no].u.d.s[0] = DomainMoment[i].s[0];
	  Nodes[no].u.d.s[1] = DomainMoment[i].s[1];
	  Nodes[no].u.d.s[2] = DomainMoment[i].s[2];
	  Nodes[no].u.d.mass = DomainMoment[i].mass;
	  Nodes[no].len = DomainMoment[i].len;
	  Nodes[no].u.d.bitflags =
	    (Nodes[no].u.d.bitflags & (~((15 << 2) + 128))) | (DomainMoment[i].bitflags & ((15 << 2) + 128));
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
	  Nodes[no].maxsoft = DomainMoment[i].maxsoft;
#endif
	}
    }

  t1 = second();
  CPU_Step[CPU_TREEBUILD] += timediff(t0, t1);
  All.Cadj_Cpu += timediff(t0, t1);
}

/*! This function updates the top-level tree after the multipole moments of
 *  the pseudo-particles have been updated.
 */
void force_treeupdate_pseudos(int no)
{
  int j, p, flag_single, count_daughters, count_particles, multiple_flag;
  FLOAT hmax;
  FLOAT s[3], mass;
  FLOAT single_s[3], single_mass = 0;

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
  int maxsofttype, diffsoftflag;
#else
  FLOAT maxsoft;
#endif
#endif

  mass = 0;
  s[0] = 0;
  s[1] = 0;
  s[2] = 0;
  hmax = 0;
  flag_single = 0;
  count_particles = 0;
#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
  maxsofttype = 7;
  diffsoftflag = 0;
#else
  maxsoft = 0;
#endif
#endif

  count_daughters = ((Nodes[no].u.d.bitflags >> 8) & 15);

  p = Nodes[no].u.d.nextnode;

  for(j = 0; j < count_daughters; j++)
    {
      if(p >= All.MaxPart && p < All.MaxPart + MaxNodes)	/* internal node */
	{
	  if(Nodes[p].u.d.bitflags & 64)
	    force_treeupdate_pseudos(p);

	  if(flag_single == 0)
	    {
	      single_s[0] = Nodes[p].u.d.s[0];
	      single_s[1] = Nodes[p].u.d.s[1];
	      single_s[2] = Nodes[p].u.d.s[2];
	      single_mass = Nodes[p].u.d.mass;
	    }

	  mass += (Nodes[p].u.d.mass);
	  s[0] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[0]);
	  s[1] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[1]);
	  s[2] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[2]);
	  if(Extnodes[p].hmax > hmax)
	    hmax = Extnodes[p].hmax;

	  if(Nodes[p].u.d.mass > 0)
	    {
	      flag_single++;

	      if((Nodes[p].u.d.bitflags & 128))
		count_particles += 2;
	      else
		count_particles++;
	    }

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
	  diffsoftflag |= (Nodes[p].u.d.bitflags >> 5) & 1;

	  if(maxsofttype == 7)
	    {
	      maxsofttype = (Nodes[p].u.d.bitflags >> 2) & 7;
	    }
	  else
	    {
	      if(((Nodes[p].u.d.bitflags >> 2) & 7) != 7)
		{
		  if(All.ForceSoftening[((Nodes[p].u.d.bitflags >> 2) & 7)] > All.ForceSoftening[maxsofttype])
		    {
		      maxsofttype = ((Nodes[p].u.d.bitflags >> 2) & 7);
		      diffsoftflag = 1;
		    }
		  else
		    {
		      if(All.ForceSoftening[((Nodes[p].u.d.bitflags >> 2) & 7)] <
			 All.ForceSoftening[maxsofttype])
			diffsoftflag = 1;
		    }
		}
	    }
#else
	  if(Nodes[p].maxsoft > maxsoft)
	    maxsoft = Nodes[p].maxsoft;
#endif
#endif
	}
      else
	endrun(6767);		/* may not happen */

      p = Nodes[p].u.d.sibling;
    }

  if(flag_single == 1)
    {
      Nodes[no].u.d.s[0] = single_s[0];
      Nodes[no].u.d.s[1] = single_s[1];
      Nodes[no].u.d.s[2] = single_s[2];
      Nodes[no].u.d.mass = single_mass;
    }
  else
    {
      if(mass)
	{
	  s[0] /= mass;
	  s[1] /= mass;
	  s[2] /= mass;
	}
      else
	{
	  s[0] = Nodes[no].center[0];
	  s[1] = Nodes[no].center[1];
	  s[2] = Nodes[no].center[2];
	}

      Nodes[no].u.d.s[0] = s[0];
      Nodes[no].u.d.s[1] = s[1];
      Nodes[no].u.d.s[2] = s[2];
      Nodes[no].u.d.mass = mass;
    }

  Extnodes[no].hmax = hmax;

  if(count_particles > 1)
    multiple_flag = 1;
  else
    multiple_flag = 0;


#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
  Nodes[no].u.d.bitflags =
    (Nodes[no].u.d.bitflags & (~(4 * 7 + 32 + 128))) + 4 * maxsofttype + 32 * diffsoftflag +
    128 * multiple_flag;
#else
  Nodes[no].u.d.bitflags = (Nodes[no].u.d.bitflags & (~128)) + 128 * multiple_flag;
  Nodes[no].maxsoft = maxsoft;
#endif
#else
  Nodes[no].u.d.bitflags = (Nodes[no].u.d.bitflags & (~128)) + 128 * multiple_flag;
#endif
}



/*! This function flags nodes in the top-level tree that are dependent on
 *  local particle data.
 */
void force_flag_localnodes(void)
{
  int no, i;

  /* mark all top-level nodes */

  for(i = 0; i < NTopleaves; i++)
    {
      no = DomainNodeIndex[i];

      Nodes[no].u.d.bitflags |= (i << 12);

      while(no >= 0)
	{
	  if((Nodes[no].u.d.bitflags & 1))
	    break;

	  Nodes[no].u.d.bitflags |= 1;

	  no = Nodes[no].u.d.father;
	}

      /* mark also internal top level nodes */

      no = DomainNodeIndex[i];
      no = Nodes[no].u.d.father;

      while(no >= 0)
	{
	  if((Nodes[no].u.d.bitflags & 64))
	    break;

	  Nodes[no].u.d.bitflags |= 64;

	  no = Nodes[no].u.d.father;
	}
    }

  /* mark top-level nodes that contain local particles */

  for(i = DomainMyStart; i <= DomainMyLast; i++)
    {
      no = DomainNodeIndex[i];

      while(no >= 0)
	{
	  if((Nodes[no].u.d.bitflags & 2))
	    break;

	  Nodes[no].u.d.bitflags |= 2;

	  no = Nodes[no].u.d.father;
	}
    }
}


/*! When a new additional star particle is created, we can put it into the
 *  tree at the position of the spawning gas particle. This is possible
 *  because the Nextnode[] array essentially describes the full tree walk as a
 *  link list. Multipole moments of tree nodes need not be changed.
 */
void force_add_star_to_tree(int igas, int istar)
{
  int no;

  no = Nextnode[igas];
  Nextnode[igas] = istar;
  Nextnode[istar] = no;
  Father[istar] = Father[igas];
  Nodes[Father[istar]].u.d.bitflags |= 128;
  Nodes[Father[istar]].u.d.bitflags += 256;
}


void force_dynamic_update(void)
{
  int i, k, no;
  FLOAT minbound[3], maxbound[3];
  double t0, t1;

  t0 = second();

  for(i = DomainMyStart; i <= DomainMyLast; i++)
    {
      no = DomainNodeIndex[i];

      force_dynamic_update_node(no, 0, minbound, maxbound);

      for(k = 0; k < 3; k++)
	{
	  DomainMoment[i].minbound[k] = minbound[k];
	  DomainMoment[i].maxbound[k] = maxbound[k];
	}
    }

  t1 = second();
  CPU_Step[CPU_TREEUPDATE] += timediff(t0, t1);
  All.Cadj_Cpu += timediff(t0, t1);


  force_exchange_pseudodata();


  t0 = second();

  force_dynamic_update_node(All.MaxPart, 1, minbound, maxbound);

  t1 = second();
  CPU_Step[CPU_TREEUPDATE] += timediff(t0, t1);
  All.Cadj_Cpu += timediff(t0, t1);
}



void force_dynamic_update_node(int no, int mode, FLOAT * minbound, FLOAT * maxbound)
{
  int i, j, k, p, single_count, count_daughters;
  FLOAT hmax, lenmax;
  FLOAT single_s[3], single_mass = 0;
  FLOAT s[3], mass;
  FLOAT mindaughter[3], maxdaughter[3];
  struct particle_data *pa;


  if(no >= All.MaxPart && no < All.MaxPart + MaxNodes)	/* internal node */
    {
      mass = 0;
      s[0] = 0;
      s[1] = 0;
      s[2] = 0;
      hmax = 0;
      single_count = 0;

      for(k = 0; k < 3; k++)
	{
	  minbound[k] = MAX_REAL_NUMBER;
	  maxbound[k] = -MAX_REAL_NUMBER;
	}

      count_daughters = ((Nodes[no].u.d.bitflags >> 8) & 15);

      p = Nodes[no].u.d.nextnode;

      for(j = 0; j < count_daughters; j++)
	{
	  if(p >= All.MaxPart)	/* an internal node or pseudo particle */
	    {
	      if(p >= All.MaxPart + MaxNodes)	/* a pseudo particle */
		{
		  /* nothing to be done here because the mass of the
		   * pseudo-particle is still zero. This will be changed
		   * later.
		   */
		  endrun(2000);	/* should not occur */
		}
	      else
		{
		  if(mode == 0 || (mode == 1 && (Nodes[p].u.d.bitflags & 64)))
		    force_dynamic_update_node(p, mode, mindaughter, maxdaughter);

		  if(mode == 1 && (Nodes[p].u.d.bitflags & 64) == 0)
		    {
		      i = Nodes[p].u.d.bitflags >> 12;

		      if(DomainNodeIndex[i] != p)
			endrun(7777);

		      for(k = 0; k < 3; k++)
			{
			  mindaughter[k] = DomainMoment[i].minbound[k];
			  maxdaughter[k] = DomainMoment[i].maxbound[k];
			}
		    }

		  if(single_count == 0)
		    {
		      single_s[0] = Nodes[p].u.d.s[0];
		      single_s[1] = Nodes[p].u.d.s[1];
		      single_s[2] = Nodes[p].u.d.s[2];
		      single_mass = Nodes[p].u.d.mass;
		    }
		  if(Nodes[p].u.d.mass > 0)
		    single_count++;

		  mass += (Nodes[p].u.d.mass);
		  s[0] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[0]);
		  s[1] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[1]);
		  s[2] += (Nodes[p].u.d.mass * Nodes[p].u.d.s[2]);

		  if(Extnodes[p].hmax > hmax)
		    hmax = Extnodes[p].hmax;

		  for(k = 0; k < 3; k++)
		    {
		      if(minbound[k] > mindaughter[k])
			minbound[k] = mindaughter[k];

		      if(maxbound[k] < maxdaughter[k])
			maxbound[k] = maxdaughter[k];
		    }

		  p = Nodes[p].u.d.sibling;
		}
	    }
	  else			/* a particle */
	    {
	      pa = &P[p];
	      if(single_count == 0)
		{
		  single_s[0] = pa->Pos[0];
		  single_s[1] = pa->Pos[1];
		  single_s[2] = pa->Pos[2];
		  single_mass = pa->Mass;
		}
	      single_count++;
	      mass += (pa->Mass);
	      s[0] += (pa->Mass * pa->Pos[0]);
	      s[1] += (pa->Mass * pa->Pos[1]);
	      s[2] += (pa->Mass * pa->Pos[2]);

	      if(pa->Type == 0)
		if(PPP[p].Hsml > hmax)
		  hmax = PPP[p].Hsml;

	      for(k = 0; k < 3; k++)
		{
		  if(minbound[k] > pa->Pos[k])
		    minbound[k] = pa->Pos[k];

		  if(maxbound[k] < pa->Pos[k])
		    maxbound[k] = pa->Pos[k];
		}

	      p = Nextnode[p];
	    }
	}

      if(single_count == 1)
	{
	  Nodes[no].u.d.s[0] = single_s[0];
	  Nodes[no].u.d.s[1] = single_s[1];
	  Nodes[no].u.d.s[2] = single_s[2];
	  Nodes[no].u.d.mass = single_mass;
	}
      else
	{
	  if(mass)
	    {
	      s[0] /= mass;
	      s[1] /= mass;
	      s[2] /= mass;
	    }
	  else
	    {
	      s[0] = Nodes[no].center[0];
	      s[1] = Nodes[no].center[1];
	      s[2] = Nodes[no].center[2];
	    }

	  Nodes[no].u.d.s[0] = s[0];
	  Nodes[no].u.d.s[1] = s[1];
	  Nodes[no].u.d.s[2] = s[2];
	  Nodes[no].u.d.mass = mass;
	}

      for(k = 0, lenmax = Nodes[no].len; k < 3; k++)
	{
	  if(lenmax < 2.0 * (Nodes[no].center[k] - minbound[k]))
	    lenmax = 2.0 * (Nodes[no].center[k] - minbound[k]);
	  if(lenmax < 2.0 * (maxbound[k] - Nodes[no].center[k]))
	    lenmax = 2.0 * (maxbound[k] - Nodes[no].center[k]);
	}

      Nodes[no].len = lenmax;

      Extnodes[no].hmax = hmax;
    }
  else
    {
      endrun(2001);		/* should not occur */
    }
}





/*! This function updates the hmax-values in tree nodes that hold SPH
 *  particles. These values are needed to find all neighbors in the
 *  hydro-force computation.  Since the Hsml-values are potentially changed in
 *  the SPH-denity computation, force_update_hmax() should be carried out just
 *  before the hydrodynamical SPH forces are computed, i.e. after density().
 */
void force_update_hmax(void)
{
  int i, no, recvTask;
  int *recvcounts, *recvoffset;
  double t0, t1;

  t0 = second();

  for(i = DomainMyStart; i <= DomainMyLast; i++)
    {
      no = DomainNodeIndex[i];

      force_update_hmax_of_node(no, 0);

      DomainHmax[i] = Extnodes[no].hmax;
    }


  /* share the hmax-data of the pseudo-particles accross CPUs */

    recvcounts = (int *)mymalloc(sizeof(int) * NTask);
    recvoffset = (int *)mymalloc(sizeof(int) * NTask);

  for(recvTask = 0; recvTask < NTask; recvTask++)
    {
      recvcounts[recvTask] = (DomainEndList[recvTask] - DomainStartList[recvTask] + 1) * sizeof(FLOAT);
      recvoffset[recvTask] = DomainStartList[recvTask] * sizeof(FLOAT);
    }


  t1 = second();
  CPU_Step[CPU_TREEHMAXUPDATE] += timediff(t0, t1);
  All.Cadj_Cpu += timediff(t0, t1);
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = second();


  MPI_Allgatherv(&DomainHmax[DomainStartList[ThisTask]], recvcounts[ThisTask], MPI_BYTE,
		 &DomainHmax[0], recvcounts, recvoffset, MPI_BYTE, MPI_COMM_WORLD);


  t1 = second();
  CPU_Step[CPU_TREECOMM] += timediff(t0, t1);
  t0 = second();


  myfree(recvoffset);
  myfree(recvcounts);

  for(i = 0; i <= DomainEndList[NTask - 1]; i++)
    if(i < DomainMyStart || i > DomainMyLast)
      {
	no = DomainNodeIndex[i];

	Extnodes[no].hmax = DomainHmax[i];
      }

  force_update_hmax_of_node(All.MaxPart, 1);


  t1 = second();
  CPU_Step[CPU_TREEHMAXUPDATE] += timediff(t0, t1);
  All.Cadj_Cpu += timediff(t0, t1);
}



/*! This routine updates the hmax-value of a node
 */
void force_update_hmax_of_node(int no, int mode)
{
  int j, p, count_daughters;
  FLOAT hmax;

  if(no >= All.MaxPart && no < All.MaxPart + MaxNodes)	/* internal node */
    {
      hmax = 0;

      count_daughters = ((Nodes[no].u.d.bitflags >> 8) & 15);

      p = Nodes[no].u.d.nextnode;

      for(j = 0; j < count_daughters; j++)
	{
	  if(p >= All.MaxPart)	/* an internal node or pseudo particle */
	    {
	      if(p >= All.MaxPart + MaxNodes)	/* a pseudo particle */
		{
		  /* nothing to be done here because the mass of the
		   * pseudo-particle is still zero. This will be changed
		   * later.
		   */
		  endrun(1000);	/* should not occur */
		}
	      else
		{
		  if(mode == 0 || (mode == 1 && (Nodes[p].u.d.bitflags & 64)))
		    force_update_hmax_of_node(p, mode);

		  if(Extnodes[p].hmax > hmax)
		    hmax = Extnodes[p].hmax;

		  p = Nodes[p].u.d.sibling;
		}
	    }
	  else			/* a particle */
	    {
	      if(P[p].Type == 0)
		if(PPP[p].Hsml > hmax)
		  hmax = PPP[p].Hsml;

	      p = Nextnode[p];
	    }
	}
      Extnodes[no].hmax = hmax;
    }
  else
    {
      endrun(1001);		/* should not occur */
    }
}





/*! This routine computes the gravitational force for a given local
 *  particle, or for a particle in the communication buffer. Depending on
 *  the value of TypeOfOpeningCriterion, either the geometrical BH
 *  cell-opening criterion, or the `relative' opening criterion is used.
 */
int force_treeevaluate(int target, int mode, double *ewaldcountsum)
{
  struct NODE *nop = 0;
  int no, ninteractions, ptype;
  double r2, dx, dy, dz, mass, r, fac, u, h, h_inv, h3_inv;
  double pos_x, pos_y, pos_z, aold;
  DOUBLE acc_x, acc_y, acc_z;

#ifdef EVALPOTENTIAL
  double wp;
  DOUBLE pot;

  pot = 0.0;
#endif
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
  double soft = 0;
#endif
#ifdef PERIODIC
  double boxsize, boxhalf;

  boxsize = All.BoxSize;
  boxhalf = 0.5 * All.BoxSize;
#endif

  acc_x = 0;
  acc_y = 0;
  acc_z = 0;
  ninteractions = 0;

  if(mode == 0)
    {
      pos_x = P[target].Pos[0];
      pos_y = P[target].Pos[1];
      pos_z = P[target].Pos[2];
      ptype = P[target].Type;
      aold = All.ErrTolForceAcc * P[target].OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
      if(ptype == 0)
	{
#ifdef ADAPTIVE_GRAVSOFT_FORGAS_HSML
	  soft = dmin(All.SofteningTable[P[target].Type],ADAPTIVE_GRAVSOFT_FORGAS_HSML * PPP[target].Hsml);
#else
	  soft = All.ForceSoftening[0] * pow(P[target].Mass / All.ReferenceGasMass, 1.0 / 3);
#endif
	}
      else
	soft = All.ForceSoftening[ptype];
#endif
    }
  else
    {
      pos_x = GravDataGet[target].u.Pos[0];
      pos_y = GravDataGet[target].u.Pos[1];
      pos_z = GravDataGet[target].u.Pos[2];
#ifdef UNEQUALSOFTENINGS
      ptype = GravDataGet[target].v.Type;
#else
      ptype = P[0].Type;
#endif
      aold = All.ErrTolForceAcc * GravDataGet[target].w.OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
      soft = GravDataGet[target].Soft;
#endif
    }



#ifndef UNEQUALSOFTENINGS
  h = All.ForceSoftening[ptype];
  h_inv = 1.0 / h;
  h3_inv = h_inv * h_inv * h_inv;
#endif
  no = All.MaxPart;		/* root node */

  while(no >= 0)
    {
      if(no < All.MaxPart)	/* single particle */
	{
	  /* the index of the node is the index of the particle */
	  /* observe the sign */

	  dx = P[no].Pos[0] - pos_x;
	  dy = P[no].Pos[1] - pos_y;
	  dz = P[no].Pos[2] - pos_z;

	  mass = P[no].Mass;
	}
      else
	{
	  if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
	    {
	      if(mode == 0)
		{
		  Exportflag[DomainTask[no - (All.MaxPart + MaxNodes)]] = 1;
		}
	      no = Nextnode[no - MaxNodes];
	      continue;
	    }
	  nop = &Nodes[no];
	  dx = nop->u.d.s[0] - pos_x;
	  dy = nop->u.d.s[1] - pos_y;
	  dz = nop->u.d.s[2] - pos_z;

	  mass = nop->u.d.mass;
	}
#ifdef PERIODIC
      dx = NEAREST(dx);
      dy = NEAREST(dy);
      dz = NEAREST(dz);
#endif
      r2 = dx * dx + dy * dy + dz * dz;

      if(no < All.MaxPart)
	{
#ifdef UNEQUALSOFTENINGS
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
	  h = soft;

	  if(P[no].Type == 0)
	    {
#ifdef ADAPTIVE_GRAVSOFT_FORGAS_HSML
	      if(h < dmin(All.SofteningTable[P[no].Type],ADAPTIVE_GRAVSOFT_FORGAS_HSML * PPP[no].Hsml);
		h = dmin(All.SofteningTable[P[no].Type],ADAPTIVE_GRAVSOFT_FORGAS_HSML * PPP[no].Hsml);
#else
	      if(h < All.ForceSoftening[0] * pow(P[no].Mass / All.ReferenceGasMass, 1.0 / 3))
		h = All.ForceSoftening[0] * pow(P[no].Mass / All.ReferenceGasMass, 1.0 / 3);
#endif
	    }
	  else
	    {
	      if(h < All.ForceSoftening[P[no].Type])
		h = All.ForceSoftening[P[no].Type];
	    }
#else
	  h = All.ForceSoftening[ptype];
	  if(h < All.ForceSoftening[P[no].Type])
	    h = All.ForceSoftening[P[no].Type];
#endif
#endif
	  no = Nextnode[no];
	}
      else			/* we have an  internal node. Need to check opening criterion */
	{
	  if(mode == 1)
	    {
	      if((nop->u.d.bitflags & 3) == 1)	/* if it's a top-level node
						 * which does not contain
						 * local particles we can
						 * continue to do a short-cut */
		{
		  no = nop->u.d.sibling;
		  continue;
		}
	    }


	  if(All.ErrTolTheta)	/* check Barnes-Hut opening criterion */
	    {
	      if(nop->len * nop->len > r2 * All.ErrTolTheta * All.ErrTolTheta)
		{
		  /* open cell */
		  no = nop->u.d.nextnode;
		  continue;
		}
	    }
	  else			/* check relative opening criterion */
	    {
	      if(mass * nop->len * nop->len > r2 * r2 * aold)
		{
		  /* open cell */
		  no = nop->u.d.nextnode;
		  continue;
		}

	      /* check in addition whether we lie inside the cell */

	      if(fabs(nop->center[0] - pos_x) < 0.60 * nop->len)
		{
		  if(fabs(nop->center[1] - pos_y) < 0.60 * nop->len)
		    {
		      if(fabs(nop->center[2] - pos_z) < 0.60 * nop->len)
			{
			  no = nop->u.d.nextnode;
			  continue;
			}
		    }
		}
	    }


#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
	  h = All.ForceSoftening[ptype];
	  if(h < All.ForceSoftening[(nop->u.d.bitflags >> 2) & 7])
	    {
	      h = All.ForceSoftening[(nop->u.d.bitflags >> 2) & 7];
	      if(r2 < h * h)
		{
		  if(((nop->u.d.bitflags >> 5) & 1))	/* bit-5 signals that there are particles of different softening in the node */
		    {
		      no = nop->u.d.nextnode;
		      continue;
		    }
		}
	    }
#else
	  h = soft;

	  if(h < nop->maxsoft)
	    {
	      h = nop->maxsoft;
	      if(r2 < h * h)
		{
		  no = nop->u.d.nextnode;
		  continue;
		}
	    }
#endif
#endif

	  no = nop->u.d.sibling;	/* ok, node can be used */

	  if(mode == 1)
	    {
	      if(((nop->u.d.bitflags) & 1))	/* Bit 0 signals that this node belongs to top-level tree */
		continue;
	    }
	}

      r = sqrt(r2);

      if(r >= h)
	{
	  fac = mass / (r2 * r);
#ifdef EVALPOTENTIAL
	  pot += FLT(-mass / r);
#endif
	}
      else
	{
#ifdef UNEQUALSOFTENINGS
	  h_inv = 1.0 / h;
	  h3_inv = h_inv * h_inv * h_inv;
#endif
	  u = r * h_inv;
	  if(u < 0.5)
	    fac = mass * h3_inv * (10.666666666667 + u * u * (32.0 * u - 38.4));
	  else
	    fac =
	      mass * h3_inv * (21.333333333333 - 48.0 * u +
			       38.4 * u * u - 10.666666666667 * u * u * u - 0.066666666667 / (u * u * u));

#ifdef EVALPOTENTIAL
	  if(u < 0.5)
	    wp = -2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6));
	  else
	    wp =
	      -3.2 + 0.066666666667 / u + u * u * (10.666666666667 +
						   u * (-16.0 + u * (9.6 - 2.133333333333 * u)));
	  pot += FLT(mass * h_inv * wp);
#endif
	}

#ifdef EVALPOTENTIAL
#ifdef PERIODIC
      pot += FLT(mass * ewald_pot_corr(dx, dy, dz));
#endif
#endif

      acc_x += FLT(dx * fac);
      acc_y += FLT(dy * fac);
      acc_z += FLT(dz * fac);

      ninteractions++;
    }


  /* store result at the proper place */
  if(mode == 0)
    {
      P[target].g.dGravAccel[0] = acc_x;
      P[target].g.dGravAccel[1] = acc_y;
      P[target].g.dGravAccel[2] = acc_z;
      P[target].GravCost = ninteractions;
#ifdef EVALPOTENTIAL
      P[target].p.dPotential = pot;
#endif
    }
  else
    {
      GravDataResult[target].u.Acc[0] = acc_x;
      GravDataResult[target].u.Acc[1] = acc_y;
      GravDataResult[target].u.Acc[2] = acc_z;
      GravDataResult[target].w.Ninteractions = ninteractions;
#ifdef EVALPOTENTIAL
      GravDataResult[target].v.Potential = pot;
#endif
    }

#ifdef PERIODIC
  *ewaldcountsum += force_treeevaluate_ewald_correction(target, mode, pos_x, pos_y, pos_z, aold);
#endif

  return ninteractions;
}






#ifdef PMGRID
/*! In the TreePM algorithm, the tree is walked only locally around the target
 *  coordinate.  Tree nodes that fall outside a box of half side-length Rcut=
 *  RCUT*ASMTH*MeshSize can be discarded. The short-range potential is
 *  modified by a complementary error function compared to the Newtonian
 *  form. The resulting short-range suppression compared to the Newtonian
 *  force is tabulated, because looking up from this table is faster than
 *  recomputing the corresponding factor, despite the memory-access panelty
 *  (which is bad for the cache performance) incurred by the table.
 */
int force_treeevaluate_shortrange(int target, int mode)
{
  struct NODE *nop = 0;
  int no, ptype, ninteractions, tabindex;
  double r2, dx, dy, dz, mass, r, fac, u, h, h_inv, h3_inv;
  double pos_x, pos_y, pos_z, aold;
  double eff_dist;
  double rcut, asmth, asmthfac, rcut2, dist;
  DOUBLE acc_x, acc_y, acc_z;

#ifdef ADAPTIVE_GRAVSOFT_FORGAS
  double soft = 0;
#endif
#ifdef EVALPOTENTIAL
  double wp, facpot;
  DOUBLE pot = 0;
#endif

#ifdef PERIODIC
  double boxsize, boxhalf;

  boxsize = All.BoxSize;
  boxhalf = 0.5 * All.BoxSize;
#endif


  acc_x = 0;
  acc_y = 0;
  acc_z = 0;
  ninteractions = 0;

  if(mode == 0)
    {
      pos_x = P[target].Pos[0];
      pos_y = P[target].Pos[1];
      pos_z = P[target].Pos[2];
      ptype = P[target].Type;
      aold = All.ErrTolForceAcc * P[target].OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
      if(ptype == 0)
	{
#ifdef ADAPTIVE_GRAVSOFT_FORGAS_HSML
	  soft = dmin(All.SofteningTable[P[target].Type],ADAPTIVE_GRAVSOFT_FORGAS_HSML * PPP[target].Hsml);
#else
	  soft = All.ForceSoftening[0] * pow(P[target].Mass / All.ReferenceGasMass, 1.0 / 3);
#endif
	}
      else
	soft = All.ForceSoftening[ptype];
#endif
    }
  else
    {
      pos_x = GravDataGet[target].u.Pos[0];
      pos_y = GravDataGet[target].u.Pos[1];
      pos_z = GravDataGet[target].u.Pos[2];
#ifdef UNEQUALSOFTENINGS
      ptype = GravDataGet[target].v.Type;
#else
      ptype = P[0].Type;
#endif
      aold = All.ErrTolForceAcc * GravDataGet[target].w.OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
      soft = GravDataGet[target].Soft;
#endif
    }

  rcut = All.Rcut[0];
  asmth = All.Asmth[0];
#ifdef PLACEHIGHRESREGION
  if(((1 << ptype) & (PLACEHIGHRESREGION)))
    {
      rcut = All.Rcut[1];
      asmth = All.Asmth[1];
    }
#endif
  rcut2 = rcut * rcut;

  asmthfac = 0.5 / asmth * (NTAB / 3.0);

#ifndef UNEQUALSOFTENINGS
  h = All.ForceSoftening[ptype];
  h_inv = 1.0 / h;
  h3_inv = h_inv * h_inv * h_inv;
#endif
  no = All.MaxPart;		/* root node */

  while(no >= 0)
    {
      if(no < All.MaxPart)
	{
	  /* the index of the node is the index of the particle */
	  dx = P[no].Pos[0] - pos_x;
	  dy = P[no].Pos[1] - pos_y;
	  dz = P[no].Pos[2] - pos_z;
#ifdef PERIODIC
	  dx = NEAREST(dx);
	  dy = NEAREST(dy);
	  dz = NEAREST(dz);
#endif
	  r2 = dx * dx + dy * dy + dz * dz;

	  mass = P[no].Mass;
#ifdef UNEQUALSOFTENINGS
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
	  h = soft;

	  if(P[no].Type == 0)
	    {
#ifdef ADAPTIVE_GRAVSOFT_FORGAS_HSML
	      if(h < dmin(All.SofteningTable[P[no].Type],ADAPTIVE_GRAVSOFT_FORGAS_HSML * PPP[no].Hsml);
		h = dmin(All.SofteningTable[P[no].Type],ADAPTIVE_GRAVSOFT_FORGAS_HSML * PPP[no].Hsml);
#else
	      if(h < All.ForceSoftening[0] * pow(P[no].Mass / All.ReferenceGasMass, 1.0 / 3))
		h = All.ForceSoftening[0] * pow(P[no].Mass / All.ReferenceGasMass, 1.0 / 3);
#endif
	    }
	  else
	    {
	      if(h < All.ForceSoftening[P[no].Type])
		h = All.ForceSoftening[P[no].Type];
	    }
#else
	  h = All.ForceSoftening[ptype];
	  if(h < All.ForceSoftening[P[no].Type])
	    h = All.ForceSoftening[P[no].Type];
#endif
#endif
	  no = Nextnode[no];
	}
      else			/* we have an  internal node */
	{
	  if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
	    {
	      if(mode == 0)
		{
		  Exportflag[DomainTask[no - (All.MaxPart + MaxNodes)]] = 1;
		}
	      no = Nextnode[no - MaxNodes];
	      continue;
	    }

	  nop = &Nodes[no];

	  if(mode == 1)
	    {
	      if((nop->u.d.bitflags & 3) == 1)	/* if it's a top-level node
						 * which does not contain
						 * local particles we can
						 * continue at this point
						 */
		{
		  no = nop->u.d.sibling;
		  continue;
		}
	    }

	  mass = nop->u.d.mass;

	  dx = nop->u.d.s[0] - pos_x;
	  dy = nop->u.d.s[1] - pos_y;
	  dz = nop->u.d.s[2] - pos_z;
#ifdef PERIODIC
	  dx = NEAREST(dx);
	  dy = NEAREST(dy);
	  dz = NEAREST(dz);
#endif
	  r2 = dx * dx + dy * dy + dz * dz;

	  if(r2 > rcut2)
	    {
	      if(nop->u.d.bitflags & 128)
		{
		  /* check whether we can stop walking along this branch */
		  eff_dist = rcut + 0.5 * nop->len;
#ifdef PERIODIC
		  dist = NEAREST(nop->center[0] - pos_x);
#else
		  dist = nop->center[0] - pos_x;
#endif
		  if(dist < -eff_dist || dist > eff_dist)
		    {
		      no = nop->u.d.sibling;
		      continue;
		    }
#ifdef PERIODIC
		  dist = NEAREST(nop->center[1] - pos_y);
#else
		  dist = nop->center[1] - pos_y;
#endif
		  if(dist < -eff_dist || dist > eff_dist)
		    {
		      no = nop->u.d.sibling;
		      continue;
		    }
#ifdef PERIODIC
		  dist = NEAREST(nop->center[2] - pos_z);
#else
		  dist = nop->center[2] - pos_z;
#endif
		  if(dist < -eff_dist || dist > eff_dist)
		    {
		      no = nop->u.d.sibling;
		      continue;
		    }
		}
	    }


	  if(All.ErrTolTheta)	/* check Barnes-Hut opening criterion */
	    {
	      if(nop->len * nop->len > r2 * All.ErrTolTheta * All.ErrTolTheta)
		{
		  /* open cell */
		  no = nop->u.d.nextnode;
		  continue;
		}
	    }
	  else			/* check relative opening criterion */
	    {
	      if(mass * nop->len * nop->len > r2 * r2 * aold)
		{
		  /* open cell */
		  no = nop->u.d.nextnode;
		  continue;
		}

	      /* check in addition whether we lie inside the cell */

	      if(fabs(nop->center[0] - pos_x) < 0.60 * nop->len)
		{
		  if(fabs(nop->center[1] - pos_y) < 0.60 * nop->len)
		    {
		      if(fabs(nop->center[2] - pos_z) < 0.60 * nop->len)
			{
			  no = nop->u.d.nextnode;
			  continue;
			}
		    }
		}
	    }

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
	  h = All.ForceSoftening[ptype];
	  if(h < All.ForceSoftening[(nop->u.d.bitflags >> 2) & 7])
	    {
	      h = All.ForceSoftening[(nop->u.d.bitflags >> 2) & 7];
	      if(r2 < h * h)
		{
		  if(((nop->u.d.bitflags >> 5) & 1))	/* bit-5 signals that there are particles of different softening in the node */
		    {
		      no = nop->u.d.nextnode;

		      continue;
		    }
		}
	    }
#else
	  h = soft;

	  if(h < nop->maxsoft)
	    {
	      h = nop->maxsoft;
	      if(r2 < h * h)
		{
		  no = nop->u.d.nextnode;
		  continue;
		}
	    }
#endif
#endif
	  no = nop->u.d.sibling;	/* ok, node can be used */

	  if(mode == 1)
	    {
	      if(((nop->u.d.bitflags) & 1))	/* Bit 0 signals that this node belongs to top-level tree */
		continue;
	    }
	}

      r = sqrt(r2);

      if(r >= h)
	{
	  fac = mass / (r2 * r);
#ifdef EVALPOTENTIAL
	  facpot = -mass / r;
#endif
	}
      else
	{
#ifdef UNEQUALSOFTENINGS
	  h_inv = 1.0 / h;
	  h3_inv = h_inv * h_inv * h_inv;
#endif
	  u = r * h_inv;
	  if(u < 0.5)
	    fac = mass * h3_inv * (10.666666666667 + u * u * (32.0 * u - 38.4));
	  else
	    fac =
	      mass * h3_inv * (21.333333333333 - 48.0 * u +
			       38.4 * u * u - 10.666666666667 * u * u * u - 0.066666666667 / (u * u * u));
#ifdef EVALPOTENTIAL
	  if(u < 0.5)
	    wp = -2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6));
	  else
	    wp =
	      -3.2 + 0.066666666667 / u + u * u * (10.666666666667 +
						   u * (-16.0 + u * (9.6 - 2.133333333333 * u)));

	  facpot = mass * h_inv * wp;
#endif
	}

      tabindex = (int) (asmthfac * r);

      if(tabindex < NTAB)
	{
	  fac *= shortrange_table[tabindex];

	  acc_x += FLT(dx * fac);
	  acc_y += FLT(dy * fac);
	  acc_z += FLT(dz * fac);

#ifdef EVALPOTENTIAL
	  pot += FLT(facpot * shortrange_table_potential[tabindex]);
#endif
	  ninteractions++;
	}
    }


  /* store result at the proper place */
  if(mode == 0)
    {
      P[target].g.dGravAccel[0] = acc_x;
      P[target].g.dGravAccel[1] = acc_y;
      P[target].g.dGravAccel[2] = acc_z;
      P[target].GravCost = ninteractions;
#ifdef EVALPOTENTIAL
      P[target].p.dPotential = pot;
#endif
    }
  else
    {
      GravDataResult[target].u.Acc[0] = acc_x;
      GravDataResult[target].u.Acc[1] = acc_y;
      GravDataResult[target].u.Acc[2] = acc_z;
      GravDataResult[target].w.Ninteractions = ninteractions;
#ifdef EVALPOTENTIAL
      GravDataResult[target].v.Potential = pot;
#endif
    }

  return ninteractions;
}

#endif



#ifdef PERIODIC
/*! This function computes the Ewald correction, and is needed if periodic
 *  boundary conditions together with a pure tree algorithm are used. Note
 *  that the ordinary tree walk does not carry out this correction directly as
 *  it was done in Gadget-1.1. Instead, the tree is walked a second time. This
 *  is actually faster because the "Ewald-Treewalk" can use a different
 *  opening criterion than the normal tree walk. In particular, the Ewald
 *  correction is negligible for particles that are very close, but it is
 *  large for particles that are far away (this is quite different for the
 *  normal direct force). So we can here use a different opening
 *  criterion. Sufficient accuracy is usually obtained if the node length has
 *  dropped to a certain fraction ~< 0.25 of the BoxLength. However, we may
 *  only short-cut the interaction list of the normal full Ewald tree walk if
 *  we are sure that the whole node and all daughter nodes "lie on the same
 *  side" of the periodic boundary, i.e. that the real tree walk would not
 *  find a daughter node or particle that was mapped to a different nearest
 *  neighbour position when the tree walk would be further refined.
 */
int force_treeevaluate_ewald_correction(int target, int mode, double pos_x, double pos_y, double pos_z,
					double aold)
{
  struct NODE *nop = 0;
  int no, cost;
  double dx, dy, dz, mass, r2;
  int signx, signy, signz;
  int i, j, k, openflag;
  double u, v, w;
  double f1, f2, f3, f4, f5, f6, f7, f8;
  DOUBLE acc_x, acc_y, acc_z;
  double boxsize, boxhalf;

  boxsize = All.BoxSize;
  boxhalf = 0.5 * All.BoxSize;

  acc_x = 0;
  acc_y = 0;
  acc_z = 0;
  cost = 0;

  no = All.MaxPart;

  while(no >= 0)
    {
      if(no < All.MaxPart)	/* single particle */
	{
	  /* the index of the node is the index of the particle */
	  /* observe the sign */

	  dx = P[no].Pos[0] - pos_x;
	  dy = P[no].Pos[1] - pos_y;
	  dz = P[no].Pos[2] - pos_z;
	  mass = P[no].Mass;
	}
      else
	{
	  if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
	    {
	      if(mode == 0)
		{
		  Exportflag[DomainTask[no - (All.MaxPart + MaxNodes)]] = 1;
		}

	      no = Nextnode[no - MaxNodes];
	      continue;
	    }

	  nop = &Nodes[no];
	  dx = nop->u.d.s[0] - pos_x;
	  dy = nop->u.d.s[1] - pos_y;
	  dz = nop->u.d.s[2] - pos_z;
	  mass = nop->u.d.mass;
	}

      dx = NEAREST(dx);
      dy = NEAREST(dy);
      dz = NEAREST(dz);

      if(no < All.MaxPart)
	no = Nextnode[no];
      else			/* we have an  internal node. Need to check opening criterion */
	{
	  openflag = 0;

	  r2 = dx * dx + dy * dy + dz * dz;

	  if(All.ErrTolTheta)	/* check Barnes-Hut opening criterion */
	    {
	      if(nop->len * nop->len > r2 * All.ErrTolTheta * All.ErrTolTheta)
		{
		  openflag = 1;
		}
	    }
	  else			/* check relative opening criterion */
	    {
	      if(mass * nop->len * nop->len > r2 * r2 * aold)
		{
		  openflag = 1;
		}
	      else
		{
		  if(fabs(nop->center[0] - pos_x) < 0.60 * nop->len)
		    {
		      if(fabs(nop->center[1] - pos_y) < 0.60 * nop->len)
			{
			  if(fabs(nop->center[2] - pos_z) < 0.60 * nop->len)
			    {
			      openflag = 1;
			    }
			}
		    }
		}
	    }

	  if(openflag)
	    {
	      /* now we check if we can avoid opening the cell */

	      u = nop->center[0] - pos_x;
	      if(u > boxhalf)
		u -= boxsize;
	      if(u < -boxhalf)
		u += boxsize;

	      if(fabs(u) > 0.5 * (boxsize - nop->len))
		{
		  no = nop->u.d.nextnode;
		  continue;
		}

	      u = nop->center[1] - pos_y;
	      if(u > boxhalf)
		u -= boxsize;
	      if(u < -boxhalf)
		u += boxsize;

	      if(fabs(u) > 0.5 * (boxsize - nop->len))
		{
		  no = nop->u.d.nextnode;
		  continue;
		}

	      u = nop->center[2] - pos_z;
	      if(u > boxhalf)
		u -= boxsize;
	      if(u < -boxhalf)
		u += boxsize;

	      if(fabs(u) > 0.5 * (boxsize - nop->len))
		{
		  no = nop->u.d.nextnode;
		  continue;
		}

	      /* if the cell is too large, we need to refine
	       * it further
	       */
	      if(nop->len > 0.20 * boxsize)
		{
		  /* cell is too large */
		  no = nop->u.d.nextnode;
		  continue;
		}
	    }

	  no = nop->u.d.sibling;	/* ok, node can be used */

	  if(mode == 1)
	    {
	      if((nop->u.d.bitflags & 1))	/* Bit 0 signals that this node belongs to top-level tree */
		continue;
	    }
	}

      /* compute the Ewald correction force */

      if(dx < 0)
	{
	  dx = -dx;
	  signx = +1;
	}
      else
	signx = -1;

      if(dy < 0)
	{
	  dy = -dy;
	  signy = +1;
	}
      else
	signy = -1;

      if(dz < 0)
	{
	  dz = -dz;
	  signz = +1;
	}
      else
	signz = -1;

      u = dx * fac_intp;
      i = (int) u;
      if(i >= EN)
	i = EN - 1;
      u -= i;
      v = dy * fac_intp;
      j = (int) v;
      if(j >= EN)
	j = EN - 1;
      v -= j;
      w = dz * fac_intp;
      k = (int) w;
      if(k >= EN)
	k = EN - 1;
      w -= k;

      /* compute factors for trilinear interpolation */

      f1 = (1 - u) * (1 - v) * (1 - w);
      f2 = (1 - u) * (1 - v) * (w);
      f3 = (1 - u) * (v) * (1 - w);
      f4 = (1 - u) * (v) * (w);
      f5 = (u) * (1 - v) * (1 - w);
      f6 = (u) * (1 - v) * (w);
      f7 = (u) * (v) * (1 - w);
      f8 = (u) * (v) * (w);

      acc_x += FLT(mass * signx * (fcorrx[i][j][k] * f1 +
				   fcorrx[i][j][k + 1] * f2 +
				   fcorrx[i][j + 1][k] * f3 +
				   fcorrx[i][j + 1][k + 1] * f4 +
				   fcorrx[i + 1][j][k] * f5 +
				   fcorrx[i + 1][j][k + 1] * f6 +
				   fcorrx[i + 1][j + 1][k] * f7 + fcorrx[i + 1][j + 1][k + 1] * f8));

      acc_y += FLT(mass * signy * (fcorry[i][j][k] * f1 +
				   fcorry[i][j][k + 1] * f2 +
				   fcorry[i][j + 1][k] * f3 +
				   fcorry[i][j + 1][k + 1] * f4 +
				   fcorry[i + 1][j][k] * f5 +
				   fcorry[i + 1][j][k + 1] * f6 +
				   fcorry[i + 1][j + 1][k] * f7 + fcorry[i + 1][j + 1][k + 1] * f8));

      acc_z += FLT(mass * signz * (fcorrz[i][j][k] * f1 +
				   fcorrz[i][j][k + 1] * f2 +
				   fcorrz[i][j + 1][k] * f3 +
				   fcorrz[i][j + 1][k + 1] * f4 +
				   fcorrz[i + 1][j][k] * f5 +
				   fcorrz[i + 1][j][k + 1] * f6 +
				   fcorrz[i + 1][j + 1][k] * f7 + fcorrz[i + 1][j + 1][k + 1] * f8));
      cost++;
    }


  /* add the result at the proper place */

  if(mode == 0)
    {
      P[target].g.dGravAccel[0] += acc_x;
      P[target].g.dGravAccel[1] += acc_y;
      P[target].g.dGravAccel[2] += acc_z;
      P[target].GravCost += cost;
    }
  else
    {
      GravDataResult[target].u.Acc[0] += acc_x;
      GravDataResult[target].u.Acc[1] += acc_y;
      GravDataResult[target].u.Acc[2] += acc_z;
      GravDataResult[target].w.Ninteractions += cost;
    }

  return cost;
}

#endif





#if defined(COMPUTE_POTENTIAL_ENERGY) || defined(OUTPUTPOTENTIAL)

/*! This routine computes the gravitational potential by walking the
 *  tree. The same opening criteria is used as for the gravitational force
 *  walk.
 */
void force_treeevaluate_potential(int target, int mode)
{
  struct NODE *nop = 0;
  DOUBLE pot;
  int no, ptype;
  double r2, dx, dy, dz, mass, r, u, h, h_inv, wp;
  double pos_x, pos_y, pos_z, aold;

#ifdef ADAPTIVE_GRAVSOFT_FORGAS
  double soft = 0;
#endif
#ifdef PERIODIC
  double boxsize, boxhalf;

  boxsize = All.BoxSize;
  boxhalf = 0.5 * All.BoxSize;
#endif

  pot = 0;

  if(mode == 0)
    {
      pos_x = P[target].Pos[0];
      pos_y = P[target].Pos[1];
      pos_z = P[target].Pos[2];
      ptype = P[target].Type;
      aold = All.ErrTolForceAcc * P[target].OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
      if(ptype == 0)
	{
#ifdef ADAPTIVE_GRAVSOFT_FORGAS_HSML
	  soft = dmin(All.SofteningTable[P[target].Type],ADAPTIVE_GRAVSOFT_FORGAS_HSML * PPP[target].Hsml);
#else
	  soft = All.ForceSoftening[0] * pow(P[target].Mass / All.ReferenceGasMass, 1.0 / 3);
#endif
	}
      else
	soft = All.ForceSoftening[ptype];
#endif
    }
  else
    {
      pos_x = GravDataGet[target].u.Pos[0];
      pos_y = GravDataGet[target].u.Pos[1];
      pos_z = GravDataGet[target].u.Pos[2];
#ifdef UNEQUALSOFTENINGS
      ptype = GravDataGet[target].v.Type;
#else
      ptype = P[0].Type;
#endif
      aold = All.ErrTolForceAcc * GravDataGet[target].w.OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
      soft = GravDataGet[target].Soft;
#endif
    }


#ifndef UNEQUALSOFTENINGS
  h = All.ForceSoftening[ptype];
  h_inv = 1.0 / h;
#endif
  no = All.MaxPart;

  while(no >= 0)
    {
      if(no < All.MaxPart)	/* single particle */
	{
	  /* the index of the node is the index of the particle */
	  /* observe the sign */

	  dx = P[no].Pos[0] - pos_x;
	  dy = P[no].Pos[1] - pos_y;
	  dz = P[no].Pos[2] - pos_z;
	  mass = P[no].Mass;
	}
      else
	{
	  if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
	    {
	      if(mode == 0)
		{
		  Exportflag[DomainTask[no - (All.MaxPart + MaxNodes)]] = 1;
		}
	      no = Nextnode[no - MaxNodes];
	      continue;
	    }

	  nop = &Nodes[no];
	  dx = nop->u.d.s[0] - pos_x;
	  dy = nop->u.d.s[1] - pos_y;
	  dz = nop->u.d.s[2] - pos_z;
	  mass = nop->u.d.mass;
	}

#ifdef PERIODIC
      dx = NEAREST(dx);
      dy = NEAREST(dy);
      dz = NEAREST(dz);
#endif
      r2 = dx * dx + dy * dy + dz * dz;

      if(no < All.MaxPart)
	{
#ifdef UNEQUALSOFTENINGS
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
	  h = soft;

	  if(P[no].Type == 0)
	    {
#ifdef ADAPTIVE_GRAVSOFT_FORGAS_HSML
	      if(h < dmin(All.SofteningTable[P[no].Type],ADAPTIVE_GRAVSOFT_FORGAS_HSML * PPP[no].Hsml);
		h = dmin(All.SofteningTable[P[no].Type],ADAPTIVE_GRAVSOFT_FORGAS_HSML * PPP[no].Hsml);
#else
	      if(h < All.ForceSoftening[0] * pow(P[no].Mass / All.ReferenceGasMass, 1.0 / 3))
		h = All.ForceSoftening[0] * pow(P[no].Mass / All.ReferenceGasMass, 1.0 / 3);
#endif

	    }
	  else
	    {
	      if(h < All.ForceSoftening[P[no].Type])
		h = All.ForceSoftening[P[no].Type];
	    }
#else
	  h = All.ForceSoftening[ptype];
	  if(h < All.ForceSoftening[P[no].Type])
	    h = All.ForceSoftening[P[no].Type];
#endif
#endif
	  no = Nextnode[no];
	}
      else			/* we have an internal node. Need to check opening criterion */
	{
	  if(mode == 1)
	    {
	      if((nop->u.d.bitflags & 3) == 1)	/* if it's a top-level node
						 * which does not contain
						 * local particles we can make
						 * a short-cut
						 */
		{
		  no = nop->u.d.sibling;
		  continue;
		}
	    }

	  if(All.ErrTolTheta)	/* check Barnes-Hut opening criterion */
	    {
	      if(nop->len * nop->len > r2 * All.ErrTolTheta * All.ErrTolTheta)
		{
		  /* open cell */
		  no = nop->u.d.nextnode;
		  continue;
		}
	    }
	  else			/* check relative opening criterion */
	    {
	      if(mass * nop->len * nop->len > r2 * r2 * aold)
		{
		  /* open cell */
		  no = nop->u.d.nextnode;
		  continue;
		}

	      if(fabs(nop->center[0] - pos_x) < 0.60 * nop->len)
		{
		  if(fabs(nop->center[1] - pos_y) < 0.60 * nop->len)
		    {
		      if(fabs(nop->center[2] - pos_z) < 0.60 * nop->len)
			{
			  no = nop->u.d.nextnode;
			  continue;
			}
		    }
		}
	    }
#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
	  h = All.ForceSoftening[ptype];
	  if(h < All.ForceSoftening[(nop->u.d.bitflags >> 2) & 7])
	    {
	      h = All.ForceSoftening[(nop->u.d.bitflags >> 2) & 7];
	      if(r2 < h * h)
		{
		  if(((nop->u.d.bitflags >> 5) & 1))	/* bit-5 signals that there are particles of different softening in the node */
		    {
		      no = nop->u.d.nextnode;
		      continue;
		    }
		}
	    }
#else
	  h = soft;

	  if(h < nop->maxsoft)
	    {
	      h = nop->maxsoft;
	      if(r2 < h * h)
		{
		  no = nop->u.d.nextnode;
		  continue;
		}
	    }
#endif
#endif

	  no = nop->u.d.sibling;	/* node can be used */

	  if(mode == 1)
	    {
	      if(((nop->u.d.bitflags) & 1))	/* Bit 0 signals that this node belongs to top-level tree */
		continue;
	    }
	}

      r = sqrt(r2);

      if(r >= h)
	pot += FLT(-mass / r);
      else
	{
#ifdef UNEQUALSOFTENINGS
	  h_inv = 1.0 / h;
#endif
	  u = r * h_inv;

	  if(u < 0.5)
	    wp = -2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6));
	  else
	    wp =
	      -3.2 + 0.066666666667 / u + u * u * (10.666666666667 +
						   u * (-16.0 + u * (9.6 - 2.133333333333 * u)));

	  pot += FLT(mass * h_inv * wp);
	}
#ifdef PERIODIC
      pot += FLT(mass * ewald_pot_corr(dx, dy, dz));
#endif
    }

  /* store result at the proper place */

  if(mode == 0)
    P[target].p.dPotential = pot;
  else
    GravDataResult[target].v.Potential = pot;
}




#ifdef PMGRID
/*! This function computes the short-range potential when the TreePM algorithm
 *  is used. This potential is the Newtonian potential, modified by a
 *  complementary error function.
 */
void force_treeevaluate_potential_shortrange(int target, int mode)
{
  struct NODE *nop = 0;
  DOUBLE pot;
  int no, ptype, tabindex;
  double r2, dx, dy, dz, mass, r, u, h, h_inv, wp;
  double pos_x, pos_y, pos_z, aold;
  double eff_dist, fac, rcut, asmth, asmthfac;
  double dxx, dyy, dzz;

#ifdef ADAPTIVE_GRAVSOFT_FORGAS
  double soft = 0;
#endif


#ifdef PERIODIC
  double boxsize, boxhalf;

  boxsize = All.BoxSize;
  boxhalf = 0.5 * All.BoxSize;
#endif

  pot = 0;

  if(mode == 0)
    {
      pos_x = P[target].Pos[0];
      pos_y = P[target].Pos[1];
      pos_z = P[target].Pos[2];
      ptype = P[target].Type;
      aold = All.ErrTolForceAcc * P[target].OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
      if(ptype == 0)
	{
#ifdef ADAPTIVE_GRAVSOFT_FORGAS_HSML
	  soft = dmin(All.SofteningTable[P[target].Type],ADAPTIVE_GRAVSOFT_FORGAS_HSML * PPP[target].Hsml);
#else
	  soft = All.ForceSoftening[0] * pow(P[target].Mass / All.ReferenceGasMass, 1.0 / 3);
	}
      else
	soft = All.ForceSoftening[ptype];
#endif
#endif
    }
  else
    {
      pos_x = GravDataGet[target].u.Pos[0];
      pos_y = GravDataGet[target].u.Pos[1];
      pos_z = GravDataGet[target].u.Pos[2];
#ifdef UNEQUALSOFTENINGS
      ptype = GravDataGet[target].v.Type;
#else
      ptype = P[0].Type;
#endif
      aold = All.ErrTolForceAcc * GravDataGet[target].w.OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
      soft = GravDataGet[target].Soft;
#endif
    }


  rcut = All.Rcut[0];
  asmth = All.Asmth[0];
#ifdef PLACEHIGHRESREGION
  if(((1 << ptype) & (PLACEHIGHRESREGION)))
    {
      rcut = All.Rcut[1];
      asmth = All.Asmth[1];
    }
#endif
  asmthfac = 0.5 / asmth * (NTAB / 3.0);

#ifndef UNEQUALSOFTENINGS
  h = All.ForceSoftening[ptype];
  h_inv = 1.0 / h;
#endif

  no = All.MaxPart;

  while(no >= 0)
    {
      if(no < All.MaxPart)	/* single particle */
	{
	  /* the index of the node is the index of the particle */
	  /* observe the sign  */

	  dx = P[no].Pos[0] - pos_x;
	  dy = P[no].Pos[1] - pos_y;
	  dz = P[no].Pos[2] - pos_z;
	  mass = P[no].Mass;
	}
      else
	{
	  if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
	    {
	      if(mode == 0)
		{
		  Exportflag[DomainTask[no - (All.MaxPart + MaxNodes)]] = 1;
		}
	      no = Nextnode[no - MaxNodes];
	      continue;
	    }

	  nop = &Nodes[no];
	  dx = nop->u.d.s[0] - pos_x;
	  dy = nop->u.d.s[1] - pos_y;
	  dz = nop->u.d.s[2] - pos_z;
	  mass = nop->u.d.mass;
	}

#ifdef PERIODIC
      dx = NEAREST(dx);
      dy = NEAREST(dy);
      dz = NEAREST(dz);
#endif
      r2 = dx * dx + dy * dy + dz * dz;

      if(no < All.MaxPart)
	{
#ifdef UNEQUALSOFTENINGS
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
	  h = soft;

	  if(P[no].Type == 0)
	    {
#ifdef ADAPTIVE_GRAVSOFT_FORGAS_HSML
	      if(h < dmin(All.SofteningTable[P[no].Type],ADAPTIVE_GRAVSOFT_FORGAS_HSML * PPP[no].Hsml);
		h = dmin(All.SofteningTable[P[no].Type],ADAPTIVE_GRAVSOFT_FORGAS_HSML * PPP[no].Hsml);
#else
	      if(h < All.ForceSoftening[0] * pow(P[no].Mass / All.ReferenceGasMass, 1.0 / 3))
		h = All.ForceSoftening[0] * pow(P[no].Mass / All.ReferenceGasMass, 1.0 / 3);
#endif
	    }
	  else
	    {
	      if(h < All.ForceSoftening[P[no].Type])
		h = All.ForceSoftening[P[no].Type];
	    }
#else
	  h = All.ForceSoftening[ptype];
	  if(h < All.ForceSoftening[P[no].Type])
	    h = All.ForceSoftening[P[no].Type];
#endif
#endif
	  no = Nextnode[no];
	}
      else			/* we have an  internal node. Need to check opening criterion */
	{
	  /* check whether we can stop walking along this branch */
	  if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
	    {
	      if(mode == 0)
		{
		  Exportflag[DomainTask[no - (All.MaxPart + MaxNodes)]] = 1;
		}
	      no = Nextnode[no - MaxNodes];
	      continue;
	    }

	  if(mode == 1)
	    {
	      if((nop->u.d.bitflags & 3) == 1)	/* if it's a top-level node which does not contain local particles */
		{
		  no = nop->u.d.sibling;
		  continue;
		}
	    }

	  eff_dist = rcut + 0.5 * nop->len;

	  dxx = nop->center[0] - pos_x;	/* observe the sign ! */
	  dyy = nop->center[1] - pos_y;	/* this vector is -y in my thesis notation */
	  dzz = nop->center[2] - pos_z;
#ifdef PERIODIC
	  dxx = NEAREST(dxx);
	  dyy = NEAREST(dyy);
	  dzz = NEAREST(dzz);
#endif
	  if(dxx < -eff_dist || dxx > eff_dist)
	    {
	      no = nop->u.d.sibling;
	      continue;
	    }

	  if(dyy < -eff_dist || dyy > eff_dist)
	    {
	      no = nop->u.d.sibling;
	      continue;
	    }

	  if(dzz < -eff_dist || dzz > eff_dist)
	    {
	      no = nop->u.d.sibling;
	      continue;
	    }

	  if(All.ErrTolTheta)	/* check Barnes-Hut opening criterion */
	    {
	      if(nop->len * nop->len > r2 * All.ErrTolTheta * All.ErrTolTheta)
		{
		  /* open cell */
		  no = nop->u.d.nextnode;
		  continue;
		}
	    }
	  else			/* check relative opening criterion */
	    {
	      if(mass * nop->len * nop->len > r2 * r2 * aold)
		{
		  /* open cell */
		  no = nop->u.d.nextnode;
		  continue;
		}

	      if(fabs(nop->center[0] - pos_x) < 0.60 * nop->len)
		{
		  if(fabs(nop->center[1] - pos_y) < 0.60 * nop->len)
		    {
		      if(fabs(nop->center[2] - pos_z) < 0.60 * nop->len)
			{
			  no = nop->u.d.nextnode;
			  continue;
			}
		    }
		}
	    }

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
	  h = All.ForceSoftening[ptype];
	  if(h < All.ForceSoftening[(nop->u.d.bitflags >> 2) & 7])
	    {
	      h = All.ForceSoftening[(nop->u.d.bitflags >> 2) & 7];
	      if(r2 < h * h)
		{
		  /* bit-5 signals that there are particles of
		   * different softening in the node
		   */
		  if(((nop->u.d.bitflags >> 5) & 1))
		    {
		      no = nop->u.d.nextnode;
		      continue;
		    }
		}
	    }
#else
	  h = soft;

	  if(h < nop->maxsoft)
	    {
	      h = nop->maxsoft;
	      if(r2 < h * h)
		{
		  no = nop->u.d.nextnode;
		  continue;
		}
	    }
#endif
#endif
	  no = nop->u.d.sibling;	/* node can be used */

	  if(mode == 1)
	    {
	      if(((nop->u.d.bitflags) & 1))	/* Bit 0 signals that this node belongs to top-level tree */
		continue;
	    }
	}

      r = sqrt(r2);

      tabindex = (int) (r * asmthfac);

      if(tabindex < NTAB)
	{
	  fac = shortrange_table_potential[tabindex];

	  if(r >= h)
	    pot += FLT(-fac * mass / r);
	  else
	    {
#ifdef UNEQUALSOFTENINGS
	      h_inv = 1.0 / h;
#endif
	      u = r * h_inv;

	      if(u < 0.5)
		wp = -2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6));
	      else
		wp =
		  -3.2 + 0.066666666667 / u + u * u * (10.666666666667 +
						       u * (-16.0 + u * (9.6 - 2.133333333333 * u)));
	      pot += FLT(fac * mass * h_inv * wp);
	    }
	}
    }


  /* store result at the proper place */
  if(mode == 0)
    P[target].p.dPotential = pot;
  else
    GravDataResult[target].v.Potential = pot;
}

#endif

#endif /*  end of COMPUTE_POTENTIAL_ENERGY/OUTPUTPOTENTIAL block */



/*! This function allocates the memory used for storage of the tree and of
 *  auxiliary arrays needed for tree-walk and link-lists.  Usually, maxnodes
 *  approximately equal to 0.7*maxpart is sufficient to store the tree for up
 *  to maxpart particles.
 */
void force_treeallocate(int maxnodes, int maxpart)
{
  int i;
  size_t bytes;
  double allbytes = 0;
  double u;

  MaxNodes = maxnodes;

  if(!(Nodes_base = (struct NODE *) mymalloc(bytes = (MaxNodes + 1) * sizeof(struct NODE))))
    {
      printf("failed to allocate memory for %d tree-nodes (%g MB).\n", MaxNodes, bytes / (1024.0 * 1024.0));
      endrun(3);
    }
  allbytes += bytes;

  if(!(Extnodes_base = (struct extNODE *) mymalloc(bytes = (MaxNodes + 1) * sizeof(struct extNODE))))
    {
      printf("failed to allocate memory for %d tree-extnodes (%g MB).\n", MaxNodes,
	     bytes / (1024.0 * 1024.0));
      endrun(3);
    }
  allbytes += bytes;

  Nodes = Nodes_base - All.MaxPart;
  Extnodes = Extnodes_base - All.MaxPart;

  if(!(Nextnode = (int *) mymalloc(bytes = (maxpart + MAXTOPNODES) * sizeof(int))))
    {
      printf("Failed to allocate %d spaces for 'Nextnode' array (%g MB)\n", maxpart + MAXTOPNODES,
	     bytes / (1024.0 * 1024.0));
      exit(0);
    }
  allbytes += bytes;

  if(!(Father = (int *) mymalloc(bytes = (maxpart) * sizeof(int))))
    {
      printf("Failed to allocate %d spaces for 'Father' array (%g MB)\n", maxpart, bytes / (1024.0 * 1024.0));
      exit(0);
    }
  allbytes += bytes;

  if(first_flag == 0)
    {
      first_flag = 1;

      if(ThisTask == 0)
	printf("\nAllocated %g MByte for BH-tree. %d\n\n", allbytes / (1024.0 * 1024.0),
	       (int) (sizeof(struct NODE) + sizeof(struct extNODE)));

      tabfac = NTAB / 3.0;

      for(i = 0; i < NTAB; i++)
	{
	  u = 3.0 / NTAB * (i + 0.5);
	  shortrange_table[i] = erfc(u) + 2.0 * u / sqrt(M_PI) * exp(-u * u);
	  shortrange_table_potential[i] = erfc(u);
	}
    }
}


/*! This function frees the memory allocated for the tree, i.e. it frees the
 *  space allocated by the function force_treeallocate().
 */
void force_treefree(void)
{
  myfree(Father);
  myfree(Nextnode);
  myfree(Extnodes_base);
  myfree(Nodes_base);
}




/*! This function does the force computation with direct summation for the
 *  specified particle in the communication buffer. This can be useful for
 *  debugging purposes, in particular for explicit checks of the force
 *  accuracy.
 */
#ifdef FORCETEST
int force_treeevaluate_direct(int target, int mode)
{
  double epsilon;
  double h, h_inv, dx, dy, dz, r, r2, u, r_inv, fac, dmax1, dmax2;
  int i, ptype;
  double pos_x, pos_y, pos_z;
  double acc_x, acc_y, acc_z;

#ifdef PERIODIC
  double fcorr[3];
#endif
#ifdef PERIODIC
  double boxsize, boxhalf;

  boxsize = All.BoxSize;
  boxhalf = 0.5 * All.BoxSize;
#endif

  acc_x = 0;
  acc_y = 0;
  acc_z = 0;

  if(mode == 0)
    {
      pos_x = P[target].Pos[0];
      pos_y = P[target].Pos[1];
      pos_z = P[target].Pos[2];
      ptype = P[target].Type;
    }
  else
    {
      pos_x = GravDataGet[target].u.Pos[0];
      pos_y = GravDataGet[target].u.Pos[1];
      pos_z = GravDataGet[target].u.Pos[2];
      ptype = GravDataGet[target].v.Type;
    }

  for(i = 0; i < NumPart; i++)
    {
      epsilon = DMAX(All.ForceSoftening[P[i].Type], All.ForceSoftening[ptype]);

      h = epsilon;
      h_inv = 1 / h;

      dx = P[i].Pos[0] - pos_x;
      dy = P[i].Pos[1] - pos_y;
      dz = P[i].Pos[2] - pos_z;

#ifdef PERIODIC
      while(dx > boxhalf)
	dx -= boxsize;
      while(dy > boxhalf)
	dy -= boxsize;
      while(dz > boxhalf)
	dz -= boxsize;
      while(dx < -boxhalf)
	dx += boxsize;
      while(dy < -boxhalf)
	dy += boxsize;
      while(dz < -boxhalf)
	dz += boxsize;
#endif
      r2 = dx * dx + dy * dy + dz * dz;

      r = sqrt(r2);

      u = r * h_inv;

      if(u >= 1)
	{
	  r_inv = 1 / r;

	  fac = P[i].Mass * r_inv * r_inv * r_inv;
	}
      else
	{
	  if(u < 0.5)
	    fac = P[i].Mass * h_inv * h_inv * h_inv * (10.666666666667 + u * u * (32.0 * u - 38.4));
	  else
	    fac =
	      P[i].Mass * h_inv * h_inv * h_inv * (21.333333333333 -
						   48.0 * u + 38.4 * u * u -
						   10.666666666667 * u * u *
						   u - 0.066666666667 / (u * u * u));
	}

      acc_x += dx * fac;
      acc_y += dy * fac;
      acc_z += dz * fac;

#ifdef PERIODIC
      if(u > 1.0e-5)
	{
	  ewald_corr(dx, dy, dz, fcorr);

	  acc_x += P[i].Mass * fcorr[0];
	  acc_y += P[i].Mass * fcorr[1];
	  acc_z += P[i].Mass * fcorr[2];
	}
#endif
    }


  if(mode == 0)
    {
      P[target].GravAccelDirect[0] = acc_x;
      P[target].GravAccelDirect[1] = acc_y;
      P[target].GravAccelDirect[2] = acc_z;
    }
  else
    {
      GravDataResult[target].u.Acc[0] = acc_x;
      GravDataResult[target].u.Acc[1] = acc_y;
      GravDataResult[target].u.Acc[2] = acc_z;
    }


  return NumPart;
}
#endif


/*! This function dumps some of the basic particle data to a file. If the tree
 * construction fails, it is called just before the run terminates with an
 * error message. Examination of the generated file may then give clues to
 * what caused the problem.
 */
void dump_particles(void)
{
  FILE *fd;
  char buffer[200];
  int i;

  sprintf(buffer, "particles%d.dat", ThisTask);
  fd = fopen(buffer, "w");
  my_fwrite(&NumPart, 1, sizeof(int), fd);

  for(i = 0; i < NumPart; i++)
    my_fwrite(&P[i].Pos[0], 3, sizeof(FLOAT), fd);

  for(i = 0; i < NumPart; i++)
    my_fwrite(&P[i].Vel[0], 3, sizeof(FLOAT), fd);

  for(i = 0; i < NumPart; i++)
    my_fwrite(&P[i].ID, 1, sizeof(int), fd);

  fclose(fd);
}



#ifdef PERIODIC

/*! This function initializes tables with the correction force and the
 *  correction potential due to the periodic images of a point mass located at
 * the origin. These corrections are obtained by Ewald summation. (See
 * Hernquist, Bouchet, Suto, ApJS, 1991, 75, 231) The correction fields are
 * used to obtain the full periodic force if periodic boundaries combined with
 * the pure tree algorithm are used. For the TreePM algorithm, the Ewald
 * correction is not used.
 *
 * The correction fields are stored on disk once they are computed. If a
 * corresponding file is found, they are loaded from disk to speed up the
 * initialization.  The Ewald summation is done in parallel, i.e. the
 * processors share the work to compute the tables if needed.
 */
void ewald_init(void)
{
  int i, j, k, beg, len, size, n, task, count;
  double x[3], force[3];
  char buf[200];
  FILE *fd;

  if(ThisTask == 0)
    {
      printf("initialize Ewald correction...\n");
      fflush(stdout);
    }

#ifdef DOUBLEPRECISION
  sprintf(buf, "ewald_spc_table_%d_dbl.dat", EN);
#else
  sprintf(buf, "ewald_spc_table_%d.dat", EN);
#endif

  if((fd = fopen(buf, "r")))
    {
      if(ThisTask == 0)
	{
	  printf("\nreading Ewald tables from file `%s'\n", buf);
	  fflush(stdout);
	}

      my_fread(&fcorrx[0][0][0], sizeof(FLOAT), (EN + 1) * (EN + 1) * (EN + 1), fd);
      my_fread(&fcorry[0][0][0], sizeof(FLOAT), (EN + 1) * (EN + 1) * (EN + 1), fd);
      my_fread(&fcorrz[0][0][0], sizeof(FLOAT), (EN + 1) * (EN + 1) * (EN + 1), fd);
      my_fread(&potcorr[0][0][0], sizeof(FLOAT), (EN + 1) * (EN + 1) * (EN + 1), fd);
      fclose(fd);
    }
  else
    {
      if(ThisTask == 0)
	{
	  printf("\nNo Ewald tables in file `%s' found.\nRecomputing them...\n", buf);
	  fflush(stdout);
	}

      /* ok, let's recompute things. Actually, we do that in parallel. */

      size = (EN + 1) * (EN + 1) * (EN + 1) / NTask;


      beg = ThisTask * size;
      len = size;
      if(ThisTask == (NTask - 1))
	len = (EN + 1) * (EN + 1) * (EN + 1) - beg;

      for(i = 0, count = 0; i <= EN; i++)
	for(j = 0; j <= EN; j++)
	  for(k = 0; k <= EN; k++)
	    {
	      n = (i * (EN + 1) + j) * (EN + 1) + k;
	      if(n >= beg && n < (beg + len))
		{
		  if(ThisTask == 0)
		    {
		      if((count % (len / 20)) == 0)
			{
			  printf("%4.1f percent done\n", count / (len / 100.0));
			  fflush(stdout);
			}
		    }

		  x[0] = 0.5 * ((double) i) / EN;
		  x[1] = 0.5 * ((double) j) / EN;
		  x[2] = 0.5 * ((double) k) / EN;

		  ewald_force(i, j, k, x, force);

		  fcorrx[i][j][k] = force[0];
		  fcorry[i][j][k] = force[1];
		  fcorrz[i][j][k] = force[2];

		  if(i + j + k == 0)
		    potcorr[i][j][k] = 2.8372975;
		  else
		    potcorr[i][j][k] = ewald_psi(x);

		  count++;
		}
	    }

      for(task = 0; task < NTask; task++)
	{
	  beg = task * size;
	  len = size;
	  if(task == (NTask - 1))
	    len = (EN + 1) * (EN + 1) * (EN + 1) - beg;

#ifdef DOUBLEPRECISION
	  MPI_Bcast(&fcorrx[0][0][beg], len, MPI_DOUBLE, task, MPI_COMM_WORLD);
	  MPI_Bcast(&fcorry[0][0][beg], len, MPI_DOUBLE, task, MPI_COMM_WORLD);
	  MPI_Bcast(&fcorrz[0][0][beg], len, MPI_DOUBLE, task, MPI_COMM_WORLD);
	  MPI_Bcast(&potcorr[0][0][beg], len, MPI_DOUBLE, task, MPI_COMM_WORLD);
#else
	  MPI_Bcast(&fcorrx[0][0][beg], len, MPI_FLOAT, task, MPI_COMM_WORLD);
	  MPI_Bcast(&fcorry[0][0][beg], len, MPI_FLOAT, task, MPI_COMM_WORLD);
	  MPI_Bcast(&fcorrz[0][0][beg], len, MPI_FLOAT, task, MPI_COMM_WORLD);
	  MPI_Bcast(&potcorr[0][0][beg], len, MPI_FLOAT, task, MPI_COMM_WORLD);
#endif
	}

      if(ThisTask == 0)
	{
	  printf("\nwriting Ewald tables to file `%s'\n", buf);
	  fflush(stdout);

	  if((fd = fopen(buf, "w")))
	    {
	      my_fwrite(&fcorrx[0][0][0], sizeof(FLOAT), (EN + 1) * (EN + 1) * (EN + 1), fd);
	      my_fwrite(&fcorry[0][0][0], sizeof(FLOAT), (EN + 1) * (EN + 1) * (EN + 1), fd);
	      my_fwrite(&fcorrz[0][0][0], sizeof(FLOAT), (EN + 1) * (EN + 1) * (EN + 1), fd);
	      my_fwrite(&potcorr[0][0][0], sizeof(FLOAT), (EN + 1) * (EN + 1) * (EN + 1), fd);
	      fclose(fd);
	    }
	}
    }

  fac_intp = 2 * EN / All.BoxSize;

  for(i = 0; i <= EN; i++)
    for(j = 0; j <= EN; j++)
      for(k = 0; k <= EN; k++)
	{
	  potcorr[i][j][k] /= All.BoxSize;
	  fcorrx[i][j][k] /= All.BoxSize * All.BoxSize;
	  fcorry[i][j][k] /= All.BoxSize * All.BoxSize;
	  fcorrz[i][j][k] /= All.BoxSize * All.BoxSize;
	}

  if(ThisTask == 0)
    {
      printf("initialization of periodic boundaries finished.\n");
      fflush(stdout);
    }
}


/*! This function looks up the correction force due to the infinite number of
 * periodic particle/node images. We here use trilinear interpolation to get
 * it from the precomputed tables, which contain one octant around the target
 * particle at the origin. The other octants are obtained from it by
 * exploiting the symmetry properties.
 */
#ifdef FORCETEST
void ewald_corr(double dx, double dy, double dz, double *fper)
{
  int signx, signy, signz;
  int i, j, k;
  double u, v, w;
  double f1, f2, f3, f4, f5, f6, f7, f8;

  if(dx < 0)
    {
      dx = -dx;
      signx = +1;
    }
  else
    signx = -1;

  if(dy < 0)
    {
      dy = -dy;
      signy = +1;
    }
  else
    signy = -1;

  if(dz < 0)
    {
      dz = -dz;
      signz = +1;
    }
  else
    signz = -1;

  u = dx * fac_intp;
  i = (int) u;
  if(i >= EN)
    i = EN - 1;
  u -= i;
  v = dy * fac_intp;
  j = (int) v;
  if(j >= EN)
    j = EN - 1;
  v -= j;
  w = dz * fac_intp;
  k = (int) w;
  if(k >= EN)
    k = EN - 1;
  w -= k;

  f1 = (1 - u) * (1 - v) * (1 - w);
  f2 = (1 - u) * (1 - v) * (w);
  f3 = (1 - u) * (v) * (1 - w);
  f4 = (1 - u) * (v) * (w);
  f5 = (u) * (1 - v) * (1 - w);
  f6 = (u) * (1 - v) * (w);
  f7 = (u) * (v) * (1 - w);
  f8 = (u) * (v) * (w);

  fper[0] = signx * (fcorrx[i][j][k] * f1 +
		     fcorrx[i][j][k + 1] * f2 +
		     fcorrx[i][j + 1][k] * f3 +
		     fcorrx[i][j + 1][k + 1] * f4 +
		     fcorrx[i + 1][j][k] * f5 +
		     fcorrx[i + 1][j][k + 1] * f6 +
		     fcorrx[i + 1][j + 1][k] * f7 + fcorrx[i + 1][j + 1][k + 1] * f8);

  fper[1] = signy * (fcorry[i][j][k] * f1 +
		     fcorry[i][j][k + 1] * f2 +
		     fcorry[i][j + 1][k] * f3 +
		     fcorry[i][j + 1][k + 1] * f4 +
		     fcorry[i + 1][j][k] * f5 +
		     fcorry[i + 1][j][k + 1] * f6 +
		     fcorry[i + 1][j + 1][k] * f7 + fcorry[i + 1][j + 1][k + 1] * f8);

  fper[2] = signz * (fcorrz[i][j][k] * f1 +
		     fcorrz[i][j][k + 1] * f2 +
		     fcorrz[i][j + 1][k] * f3 +
		     fcorrz[i][j + 1][k + 1] * f4 +
		     fcorrz[i + 1][j][k] * f5 +
		     fcorrz[i + 1][j][k + 1] * f6 +
		     fcorrz[i + 1][j + 1][k] * f7 + fcorrz[i + 1][j + 1][k + 1] * f8);
}
#endif


/*! This function looks up the correction potential due to the infinite number
 * of periodic particle/node images. We here use tri-linear interpolation to
 * get it from the precomputed table, which contains one octant around the
 * target particle at the origin. The other octants are obtained from it by
 * exploiting symmetry properties.
 */
double ewald_pot_corr(double dx, double dy, double dz)
{
  int i, j, k;
  double u, v, w;
  double f1, f2, f3, f4, f5, f6, f7, f8;

  if(dx < 0)
    dx = -dx;

  if(dy < 0)
    dy = -dy;

  if(dz < 0)
    dz = -dz;

  u = dx * fac_intp;
  i = (int) u;
  if(i >= EN)
    i = EN - 1;
  u -= i;
  v = dy * fac_intp;
  j = (int) v;
  if(j >= EN)
    j = EN - 1;
  v -= j;
  w = dz * fac_intp;
  k = (int) w;
  if(k >= EN)
    k = EN - 1;
  w -= k;

  f1 = (1 - u) * (1 - v) * (1 - w);
  f2 = (1 - u) * (1 - v) * (w);
  f3 = (1 - u) * (v) * (1 - w);
  f4 = (1 - u) * (v) * (w);
  f5 = (u) * (1 - v) * (1 - w);
  f6 = (u) * (1 - v) * (w);
  f7 = (u) * (v) * (1 - w);
  f8 = (u) * (v) * (w);

  return potcorr[i][j][k] * f1 +
    potcorr[i][j][k + 1] * f2 +
    potcorr[i][j + 1][k] * f3 +
    potcorr[i][j + 1][k + 1] * f4 +
    potcorr[i + 1][j][k] * f5 +
    potcorr[i + 1][j][k + 1] * f6 + potcorr[i + 1][j + 1][k] * f7 + potcorr[i + 1][j + 1][k + 1] * f8;
}



/*! This function computes the potential correction term by means of Ewald
 * summation.
 */
double ewald_psi(double x[3])
{
  double alpha, psi;
  double r, sum1, sum2, hdotx;
  double dx[3];
  int i, n[3], h[3], h2;

  alpha = 2.0;

  for(n[0] = -4, sum1 = 0; n[0] <= 4; n[0]++)
    for(n[1] = -4; n[1] <= 4; n[1]++)
      for(n[2] = -4; n[2] <= 4; n[2]++)
	{
	  for(i = 0; i < 3; i++)
	    dx[i] = x[i] - n[i];

	  r = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
	  sum1 += erfc(alpha * r) / r;
	}

  for(h[0] = -4, sum2 = 0; h[0] <= 4; h[0]++)
    for(h[1] = -4; h[1] <= 4; h[1]++)
      for(h[2] = -4; h[2] <= 4; h[2]++)
	{
	  hdotx = x[0] * h[0] + x[1] * h[1] + x[2] * h[2];
	  h2 = h[0] * h[0] + h[1] * h[1] + h[2] * h[2];
	  if(h2 > 0)
	    sum2 += 1 / (M_PI * h2) * exp(-M_PI * M_PI * h2 / (alpha * alpha)) * cos(2 * M_PI * hdotx);
	}

  r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);

  psi = M_PI / (alpha * alpha) - sum1 - sum2 + 1 / r;

  return psi;
}


/*! This function computes the force correction term (difference between full
 *  force of infinite lattice and nearest image) by Ewald summation.
 */
void ewald_force(int iii, int jjj, int kkk, double x[3], double force[3])
{
  double alpha, r2;
  double r, val, hdotx, dx[3];
  int i, h[3], n[3], h2;

  alpha = 2.0;

  for(i = 0; i < 3; i++)
    force[i] = 0;

  if(iii == 0 && jjj == 0 && kkk == 0)
    return;

  r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];

  for(i = 0; i < 3; i++)
    force[i] += x[i] / (r2 * sqrt(r2));

  for(n[0] = -4; n[0] <= 4; n[0]++)
    for(n[1] = -4; n[1] <= 4; n[1]++)
      for(n[2] = -4; n[2] <= 4; n[2]++)
	{
	  for(i = 0; i < 3; i++)
	    dx[i] = x[i] - n[i];

	  r = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);

	  val = erfc(alpha * r) + 2 * alpha * r / sqrt(M_PI) * exp(-alpha * alpha * r * r);

	  for(i = 0; i < 3; i++)
	    force[i] -= dx[i] / (r * r * r) * val;
	}

  for(h[0] = -4; h[0] <= 4; h[0]++)
    for(h[1] = -4; h[1] <= 4; h[1]++)
      for(h[2] = -4; h[2] <= 4; h[2]++)
	{
	  hdotx = x[0] * h[0] + x[1] * h[1] + x[2] * h[2];
	  h2 = h[0] * h[0] + h[1] * h[1] + h[2] * h[2];

	  if(h2 > 0)
	    {
	      val = 2.0 / ((double) h2) * exp(-M_PI * M_PI * h2 / (alpha * alpha)) * sin(2 * M_PI * hdotx);

	      for(i = 0; i < 3; i++)
		force[i] -= h[i] * val;
	    }
	}
}

#endif
