/* Definitions and declarations for Metallicity Model  */


/* constants */
#define Nab 8      /* Number of Sutherland & Dopita tables  */
#define delta 0.05 /* Interval of cooling tables  */


/* variables */
extern int *Nlines, cont_sm;
extern  double **logLambda_i;
extern  double **logT_i;
extern double **yield1, **yield2, **yield3, **yield4, **yield5;
extern double metal, sm2, FeHgas;

extern float **nsimfww;
extern double XH, yhelium;
extern int numenrich;


/* functions */
void find_low_density_tail(void);
int  compare_density_values(const void *a, const void *b);
void promote_particles(void);
void decide_on_promotion(void);
void add_in_energies(void);
void read_coolrate_table(void);
double CoolingRate_SD(double logT, double rho, double *nelec);
double get_Lambda_SD(double logT, double abund);
void read_yield_table(void);
double SNII_yields(int indice);
void imf(void);
void flag_SN_starparticles(void);
void update_weights(void);
double SNI_yields(int indice);
void enrichment(void);
float integrated_time(int indice, double time_hubble_a);
void phase_mass(void);
void  copy_densities(void);
int hotngbs_compare_key(const void *a, const void *b);
void hotngbs_evaluate(int target, int mode);
void find_hot_neighbours(void);
int ngb_treefind_hotngbs(FLOAT searchcenter[3], FLOAT hguess, int *startnode, FLOAT entropy);

/* input/output */
 extern FILE *fpcool;




