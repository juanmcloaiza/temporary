/*! \file allvars.c
 *  \brief creates global variables.
 *
 *  This file creates all the global variables that are declared in allvars.h
 *  
 *  To produce 'allvars.c' from 'allvars.h', do the following:
 *
 *     - Erase all #define's
 *     - add #include "allvars.h" 
 *     - delete all keywords 'extern'
 *     - delete all struct definitions enclosed in {...}, e.g.
 *        "extern struct global_data_all_processes {....} All;"
 *        becomes "struct global_data_all_processes All;"
 */

#include "allvars.h"


int ThisTask;			/* the local processors  */
int NTask, PTask;		/* note: NTask = 2^PTask */

double CPUThisRun;

int NumForceUpdate;
int NumSphUpdate;
int TimeTreeRoot;
int RestartFlag;
int Flag_FullStep;
int TreeReconstructFlag;

int NumPart;			/* Note: this is the LOCAL processor value */
int N_gas;			/* Note: this is the LOCAL processor value */
//Added by JM
int flag_accretion; /*!< LOCAL flag. If set, domain decomposition will take place */
//Added by JM
long long Ntype[6];
int NtypeLocal[6];

#ifdef LT_STELLAREVOLUTION
int N_star;
#endif

#ifdef SFR
int Stars_converted;

#ifdef SFR_METALS
double TotalEnergy;
double DEnergy_spawned, DEnergy_converted;
double DEnergy_radiation, DEnergy_promotion;
double DEnergy_feedback, TotalReservoir;

#ifdef SFR_FEEDBACK
double ESN;
int nhot, ncold;
#endif
#endif
#endif

gsl_rng *random_generator;

#ifdef SFR_METALS
int Flag_phase;
int Flag_promotion;
#endif

double CPU_Step[CPU_PARTS];
char CPU_Symbol[CPU_PARTS] =
  { '-', '*', '=', ';', '<', '[', '^', ':', '.', '~', '|', '+', '"', '/', '`', ',', '>', '@', '#', '&', '<' };
char CPU_SymbolImbalance[CPU_PARTS] =
  { 'a', 't', 'u', 'v', 'b', 'w', 'd', 'r', 'h', 'm', 'n', 'l', 'o', 'p', 's', 'f', 'i', 'g', 'c', 'e', 'z' };
char CPU_String[CPU_STRING_LEN + 1];

struct topnode_data *TopNodes;

int NTopnodes, NTopleaves;


double TimeOfLastTreeConstruction;

int *Ngblist;

#ifdef PERIODIC
FLOAT boxSize_X, boxHalf_X;
FLOAT boxSize_Y, boxHalf_Y;
FLOAT boxSize_Z, boxHalf_Z;
#endif

peanokey *DomainKeyBuf;

peanokey *Key, *KeySorted;

double DomainCorner[3], DomainCenter[3], DomainLen, DomainFac;
int DomainMyStart, DomainMyLast;
int *DomainStartList, *DomainEndList;

double *DomainWork;
int *DomainCount;
int *DomainCountSph;
int *DomainTask;

#ifdef LT_STELLAREVOLUTION
int *DomainCountStar;
#endif


int *DomainNodeIndex;

struct DomainNODE *DomainMoment;


FLOAT *DomainHmax;
FLOAT *DomainTreeNodeLen;



double RndTable[RNDTABLE];



/* variables for input/output ,  usually only used on process 0 
 */
char ParameterFile[100];
FILE *FdInfo, *FdEnergy, *FdTimings, *FdCPU, *FdBalance;
//Added by JM
FILE *FdAccretion;
//End of Added by JM

#ifdef SFR
FILE *FdSfr;
#endif

#ifdef BG_STELLAR_EVOLUTION
FILE *FdMetGas;			/*!< file handle for metals_gas.txt log-file. */
FILE *FdMetStars;		/*!< file handle for metals_star.txt log-file. */
FILE *FdMetSF;			/*!< file handle for metals_sf.txt log-file. */
FILE *FdMetTot;			/*!< file handle for metals_tot.txt log-file. */
#endif

#ifdef BLACK_HOLES
FILE *FdBlackHoles;
FILE *FdBlackHolesDetails;
#endif



#ifdef SFR_METALS
FILE *FdMphase;
FILE *FdSNE;

#if defined(SFR_SNI) || defined(SFR_SNII)
FILE *FdSN;
#endif
#ifdef SFR_PROMOTION
FILE *FdPromotion;
#endif
#endif


#ifdef FORCETEST
FILE *FdForceTest;
#endif

#ifdef XXLINFO
FILE *FdXXL;

#ifdef MAGNETIC
double MeanB;

#ifdef TRACEDIVB
double MaxDivB;
#endif
#endif
#ifdef TIME_DEP_ART_VISC
double MeanAlpha;
#endif
#endif

#ifdef DARKENERGY
FILE *FdDE;  /*!< file handle for darkenergy.txt log-file. */
#endif


#ifdef LT_CharT_INFO
FILE *FdCharT;
#endif

#ifdef LT_STELLAREVOLUTION

gsl_error_handler_t *old_error_handler;
int my_gslstatus;

FILE *FdSnInit, *FdWarn;

gsl_function F;
int gsl_status;
gsl_integration_workspace *w;

gsl_interp_accel *accel;
gsl_spline *spline;

/* ----------------- *
 *   metal species   *
 * ----------------- */

char *MetNames[LT_NMet];
double MetSolarValues[LT_NMet];
int Hyd, Hel, Iron, Oxygen, FillEl;

float xclouds;
int search_for_metalspread;

char SnIaDataFile[300], SnIIDataFile[300], AGBDataFile[300];

#ifdef LT_TRACK_CONTRIBUTES
unsigned int   Packing_Factor;
unsigned int   *Power10_Factors, Max_Power10;
unsigned int   TrackMask;
float Max_Packed_Int, UnPacking_Factor;
#endif

/* -------------- *
 *   SF related   *
 * -------------- */

double *PhysDensTh, *FEVP;
int sfrrate_filenum;

/* -------------- *
 *   Sn related   *
 * -------------- */

struct SDtype SD;

#ifdef LT_SNIa
/*
 * SnIaData defines the yields for SnIa.
 * we define it as a pointer so that in the future we can easily extend the
 * number of parameters on which they depend (e.g. mass of the system, time).
 * currently they are fixed
 */
double ***SnIaYields;

/*
 * it is common that yields tables use the same mass array and Z array for
 * all elements. in this case, we simply use IaZbins, IaMbins to store these 
 * common data. otherwise (different mass/Z array for each element), you can
 * just use the already existent Yield structure, just allocating as twice as
 * the amount of memory needed for the yields' value and using half of the
 * memory to store also the array's value. You can also use SnIaY_dim to store
 * the dimensions.
 */
double **IaZbins, **IaMbins;
int *IaZbins_dim, *IaMbins_dim, *SnIaY_dim[LT_NMet];
#endif

#ifdef LT_SNII
/*
 * the same as for SnIa.
 */
double ***SnIIYields, **SnIIEj;
double ***SnII_ShortLiv_Yields;
double **IIZbins, **IIMbins;
int *IIZbins_dim, *IIMbins_dim, *SnIIY_dim[LT_NMet];
double ***SnII_steps;
int *SnII_Nsteps;
#ifdef LT_TRACK_CONTRIBUTES
float ***SnII_ShortLiv_Yield_Contrib;
#endif
#endif

#ifdef LT_AGB
double ***AGBYields, **AGBEj;
double **AGBZbins, **AGBMbins;
int *AGBZbins_dim, *AGBMbins_dim, *AGB_dim[LT_NMet];
#endif

#if defined(LT_AGB) || defined(LT_SnIa)
double ***LLv_steps;
int *LLv_Nsteps;
#endif

char *IMF_Spec_Labels[2]= {"power law", "general"};
IMF_SPEC IMF_Spec;
int IsThere_TimeDep_IMF;
int IsThere_ZDep_IMF;

IMF_Type *IMFs, *IMFp, IMFu;
int IMFs_dim;

double IMF_func_arg;
double *IMF_CommBuff;

FILE *FdIMFin, *FdIMF;


#ifdef LT_SEv_INFO

/*   temporary   */

/* double delta_Trun, *Tpop_num, *Tpop_mass; */
/* double *tot_Tpop_num, *tot_Tpop_mass; */
/* double delta_TCrun, *TCpop_num, *TCpop_mass; */
/* double *tot_TCpop_num, *tot_TCpop_mass; */
/* int Trun_size, TCrun_size; */
/* FILE *FdTrun, *FdTCrun; */

/* ............. */


FILE *FdSn, *FdMetals, *FdSnLost, *FdSnDetails;

#ifdef LT_SEvDbg
FILE *FdMetSumCheck;
unsigned int FirstID;
int checkFirstID;
#endif

#ifdef WINDS
FILE *FdWinds;
#endif

#endif

#ifdef LT_SEv_INFO_DETAILS
double DetailsW[LT_NMet], DetailsWo[LT_NMet];
#endif

#ifdef LT_EXTEGY_INFO
FILE *FdExtEgy;
#endif
#endif

double DriftTable[DRIFT_TABLE_LENGTH], GravKickTable[DRIFT_TABLE_LENGTH], HydroKickTable[DRIFT_TABLE_LENGTH];


void *CommBuffer;		/* communication buffer, used at a number of places */

char *Exportflag;

/* ---==={{( ( METAL COOLING SECTION ) )}}===--- */
double sphpZ, R;
double ZTemp[Zlength], Zvalue[8], Lsuthdop[8][Zlength];

#ifdef LT_METAL_COOLING
/*double ThInst_onset[ZBins] = { 4.9, 4.9, 4.95, 5.35, 5.35, 5.35, 5.35, 5.35 };*/
double ThInst_onset[ZBins] = { 5, 5, 5, 5.3, 5.3, 5.3, 5.3, 5.3 };

#else

/*double ThInst_onset[ZBins] = { 4.9, 4.9, 4.95, 5.35, 5.35, 5.35, 5.35, 5.35 };*/
double ThInst_onset[ZBins] = { 5 };

#endif
/* ---===--- */

/* ---==={{( ( CHEMICAL ENRICHMENT SECTION ) )}}===--- */
#ifdef LT_STELLAREVOLUTION

float xclouds;

#ifdef LT_SEvDbg
unsigned int FirstID;
int checkFirstID;
#endif

#endif
/* ---===--- */


/* this structure contains data which is the SAME for all 
 * tasks (mostly code parameters read from the parameter file). 
 * Holding this data in a structure is convenient for writing/reading
 * the restart file, and it allows the introduction of new global
 * variables in a simple way. The only thing to do is to introduce them
 * into this structure.
 */
struct global_data_all_processes All;



/* The following structure holds all the information that is
 * stored for each particle of the simulation.
 */
struct particle_data *P, *DomainPartBuf;

#ifdef LT_STELLAREVOLUTION
/* the following struture holds data that is stored for each star particle
 * in addition to the collisionless variables.
 */
struct met_particle_data *MetP, *DomainMetBuf;
#endif


/* the following struture holds data that is stored for each SPH particle
 * in addition to the collisionless variables.
 */
struct sph_particle_data *SphP, *DomainSphBuf;




/* global state of system 
*/
struct state_of_system SysState, SysStateAtStart, SysStateAtEnd;






/* Various structure for communication during the gravity 
 * computation.
 */
struct gravdata_in *GravDataIn, *GravDataGet, *GravDataResult, *GravDataOut;

struct gravdata_index *GravDataIndexTable;

#ifdef SFR_METALS
struct metaldata_in *MetalDataIn, *MetalDataGet;
#endif

#ifdef BG_SFR
struct metaldata_in *MetalDataIn, *MetalDataGet;
#endif

/* Various structure for communication during the density
 * computation.
 */
struct densdata_in *DensDataIn, *DensDataGet;

#ifdef BG_SFR
struct solidangledata_in *SolidAngleDataIn, *SolidAngleDataGet;
struct solidangledata_out *SolidAngleDataResult, *SolidAngleDataPartialResult;
#endif


struct densdata_out *DensDataResult, *DensDataPartialResult;


#ifdef BLACK_HOLES
struct blackholedata_in *BlackholeDataIn, *BlackholeDataGet;
struct blackholedata_out *BlackholeDataResult, *BlackholeDataPartialResult;
#endif


#ifdef FOF
struct fofdata_in *FoFDataIn, *FoFDataGet;
struct fofdata_out *FoFDataResult, *FoFDataPartialResult;
#endif




/* Various structures for communication during the 
 * computation of hydrodynamical forces.
 */
struct hydrodata_in *HydroDataIn, *HydroDataGet;

struct hydrodata_out *HydroDataResult, *HydroDataPartialResult;


#ifdef SFR_PROMOTION
struct hotngbs_in *HotNgbsIn, *HotNgbsGet;
struct hotngbs_out *HotNgbsResult, *HotNgbsPartialResult;
#endif

#ifdef MHM
struct kindata_in *KinDataIn, *KinDataGet;
#endif

#ifdef LT_STELLAREVOLUTION
struct metal_ngbfindingdata_in *MetalNgbIn, *MetalNgbGet;
struct metal_ngbfindingdata_out *MetalNgbResultIn, *MetalNgbResultGet;
struct metaldata_in *MetalDataIn, *MetalDataGet;
#endif

/* Header for the standard file format.
 */
struct io_header header;



/*******************
 ** Variables for Tree
 ********************
 */


struct NODE *Nodes, *Nodes_base;

struct extNODE *Extnodes, *Extnodes_base;


int MaxNodes;			/* maximum allowed number of internal nodes */
int Numnodestree;		/* number of (internal) nodes in each tree */


int *Nextnode;			/* gives next node in tree walk  (nodes array) */
int *Father;			/* gives parent node in tree (Prenodes array) */


#ifdef STATICNFW
double Rs, R200;
double Dc;
double RhoCrit, V200;
double fac;
#endif



#ifdef CHEMISTRY
/* ----- Tables ------- */
double T[N_T], J0_nu[N_nu], J_nu[N_nu], nu[N_nu];
double k1a[N_T], k2a[N_T], k3a[N_T], k4a[N_T], k5a[N_T], k6a[N_T], k7a[N_T], k8a[N_T], k9a[N_T], k10a[N_T],
  k11a[N_T];
double k12a[N_T], k13a[N_T], k14a[N_T], k15a[N_T], k16a[N_T], k17a[N_T], k18a[N_T], k19a[N_T], k20a[N_T],
  k21a[N_T];
double ciHIa[N_T], ciHeIa[N_T], ciHeIIa[N_T], ciHeISa[N_T], reHIIa[N_T], brema[N_T];
double ceHIa[N_T], ceHeIa[N_T], ceHeIIa[N_T], reHeII1a[N_T], reHeII2a[N_T], reHeIIIa[N_T];

/* cross-sections */
#ifdef RADIATION
double sigma24[N_nu], sigma25[N_nu], sigma26[N_nu], sigma27[N_nu], sigma28[N_nu], sigma29[N_nu],
  sigma30[N_nu], sigma31[N_nu];
#endif

#endif /* chemisitry */
