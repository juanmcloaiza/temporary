#ifndef ALLVARS_H
#include "allvars.h"
#endif

void   DomainDecomposition(void); 

void domain_shiftSplit(void);
int  domain_compare_key(const void *a, const void *b);
void domain_decompose(void);
void domain_countToGo(void);
int  domain_findSplit(int cpustart, int ncpu, int first, int last);
void domain_findExchangeNumbers(int task, int partner, int sphflag, int *send, int *recv);
void domain_exchangeParticles(int partner, int sphflag, int send_count, int recv_count);
void domain_sumCost(void);
void domain_findExtent(void);

int domain_findSplit_balanced(int cpustart, int ncpu, int first, int last);
double domain_particle_costfactor(int i);

void domain_topsplit(int node, peanokey startkey);
void domain_topsplit_local(int node, peanokey startkey);
void domain_determineTopTree(void);

int domain_compare_key(const void *a, const void *b);
int domain_compare_toplist(const void *a, const void *b);





