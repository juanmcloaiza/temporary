#ifndef INLINE_FUNC
#ifdef INLINE
#define INLINE_FUNC inline
#else
#define INLINE_FUNC
#endif
#endif


void force_flag_localnodes(void);

int  force_treeevaluate(int target, int mode, double *ewaldcountsum);
int  force_treeevaluate_shortrange(int target, int mode);
int  force_treeevaluate_ewald_correction(int target, int mode, double pos_x, double pos_y, double pos_z, double aold);
void force_treeevaluate_potential_shortrange(int target, int mode);
void force_treeevaluate_potential(int target, int type);


void force_tree_discardpartials(void);
void force_treeupdate_pseudos(int);
void force_update_pseudoparticles(void);


void force_dynamic_update(void);
void force_dynamic_update_node(int no, int mode, FLOAT *minbound, FLOAT *maxbound);

void force_update_hmax(void);
void force_update_hmax_of_node(int no, int mode);


void force_create_empty_nodes(int no, int topnode, int bits, int x, int y, int z, int *nodecount, int *nextfree);

void force_exchange_pseudodata(void);

void force_insert_pseudo_particles(void);

void force_add_star_to_tree(int igas, int istar);

void   force_costevaluate(void);
int    force_getcost_single(void);
int    force_getcost_quadru(void);
void   force_resetcost(void);
void   force_setupnonrecursive(int no);
void   force_treeallocate(int maxnodes, int maxpart);  
int    force_treebuild(int npart);
int    force_treebuild_single(int npart);

int    force_treeevaluate_direct(int target, int mode);


void   force_treefree(void);
void   force_update_node(int no, int flag);

void force_update_node_recursive(int no, int sib, int father, FLOAT *minbound, FLOAT *maxbound);

void   force_update_size_of_parent_node(int no);

void   dump_particles(void);

FLOAT  INLINE_FUNC ngb_periodic(FLOAT x);
FLOAT  INLINE_FUNC ngb_periodic_longbox(FLOAT x);
FLOAT  ngb_select_closest(int k, int n, FLOAT *arr, int *ind);
void   ngb_treeallocate(int npart);
void   ngb_treebuild(void);


void   ngb_treefree(void);
void   ngb_treesearch(int);
void   ngb_treesearch_pairs(int);
void   ngb_update_nodes(void);
void   ngb_treesearch_notsee(int no);


#ifdef SFR_DECOUPLING
int ngb_treefind_variable(FLOAT searchcenter[3], FLOAT hguess, int *startnode, FLOAT densityold, FLOAT entropy, FLOAT *vel);
#else
int ngb_treefind_variable(FLOAT searchcenter[3], FLOAT hguess, int *startnode);
#endif

#ifdef SFR_DECOUPLING
int ngb_treefind_pairs(FLOAT searchcenter[3], FLOAT hsml, int *startnode, FLOAT densityold, FLOAT entropy, FLOAT *vel);
#else
int ngb_treefind_pairs(FLOAT searchcenter[3], FLOAT hsml, int *startnode);
#endif





