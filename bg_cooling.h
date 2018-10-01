#ifdef BG_COOLING

#define EL_NAME_LENGTH 20 // maximum length of character string

// cooling tables
// table nmber 1:
float ****net_heating; /// to be allocated as [redshifts(2)][number of species+number helium fractions][number of density bins][number of temperature bins]
// table number 2:
float ****u_to_t; /// to be allocated as [redshifts(2)][number helium fractions][number of density bins][number of temperature bins]
// electron abundance as function of Helium abundance and thermal energy (collisional tables)
float **collisional_electron_abundance; /// to be allocated as [number helium fractions][number of temperature bins]
//

float *tempgrid, *n_hgrid, *ugrid, *Hefracs, *cooling_redshifts, *solar_abundances; 
char **Henames, **CoolHeat_Element_Names, **Solar_Abundance_Names;
int nHes, ntemps, nhdens, nCoolHeats, nSolar_Abundances, nZs, cooling_redshift_index;

int  *SPH_Name_Pointers, *Solar_Name_Pointers;  // arrays that will map the cooling tables to the abundances, etc.

char SPH_Element_Name[BG_NELEMENTS][EL_NAME_LENGTH];

#define STEFAN       7.57e-15 //???????
#define THOMSON      6.6524587e-25 //???????
#define COMPTON_RATE 5.406e-36 //???????
#define T_CMB        2.728 //???????
#define COMP_COEFF   (4.0*STEFAN*THOMSON*BOLTZMANN/(ELECTRONMASS*C))

#define SULPSILI 0.505273
#define CALCSILI 0.0794201

///shared functions ///
void InitCool(void); //???????
void InitChemistry(void);
void get_redshift_index(float z, int z_index_old, float *dz);
void get_redshift_index(float z, int z_index_old, float *dz);
void GetCollisTable(void);
void GetCoolingTables(int z_index_old);
void BroadcastCoolingTables();
double CoolingRateFromU(double u, double inn_h, double inZ, float d_z, double particle_Z[]);
double DoCooling(double u_old, double rho, double dt, double inz, float d_z, double particle_Z[]);
void allocate_header_arrays(void);

#endif /* BG_COOLING */
