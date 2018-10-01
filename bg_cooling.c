//Implementing version
#ifdef BG_COOLING
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string.h>
#include <hdf5.h>
#include "allvars.h"
#include "bg_cooling.h"

float Hefrac;

/// Used to interpolate the cooling given cooling grids at two redshifts ///
double threewayinterp(int spec_id, int temp_i, int n_h_i, float d_temp, float d_n_h, float d_z)
{
  double X;

  X = d_z   *   d_n_h   *   d_temp   *   net_heating[0][spec_id][n_h_i][temp_i]
    + d_z   *   d_n_h   *   (1.0-d_temp)*net_heating[0][spec_id][n_h_i][temp_i+1]
    + d_z   *   (1.0-d_n_h)*d_temp   *   net_heating[0][spec_id][n_h_i+1][temp_i]
    + d_z   *   (1.0-d_n_h)*(1.0-d_temp)*net_heating[0][spec_id][n_h_i+1][temp_i+1] 
    + (1.0-d_z)*d_n_h   *   d_temp   *   net_heating[1][spec_id][n_h_i][temp_i] 
    + (1.0-d_z)*d_n_h   *   (1.0-d_temp)*net_heating[1][spec_id][n_h_i][temp_i+1]
    + (1.0-d_z)*(1.0-d_n_h)*d_temp   *   net_heating[1][spec_id][n_h_i+1][temp_i]
    + (1.0-d_z)*(1.0-d_n_h)*(1.0-d_temp)*net_heating[1][spec_id][n_h_i+1][temp_i+1];

  return X;
}

/// Interpolated from cloudy tables that use an mmw based on the Helium fraction ///
double convert_u_to_temp(int *n_h_i, int *He_i, float *d_n_h, float d_z, float *d_He, double u, double inn_h, double inHe)
{

  /// This involves a crazy amount of interpolation (in 4 variables!!), but is worth it. ///
  double T, d_u;
  int u_i ;

  /// We use the indices we find later on for the Cooling ///
  if (inHe <Hefracs[0])
  {
    *d_He = 1.0;
    *He_i = 0;
    puts("Helium fraction out of range, using lowest value");
  }
  else if (inHe > Hefracs[nHes-1])
  {
    *d_He = 0.0;
    *He_i = nHes-2;
    puts("Helium fraction out of range, using highest value");
  }
  else
  {
    *He_i = (int)floor((inHe - Hefracs[0])/(Hefracs[1] - Hefracs[0]));
    *d_He = (Hefracs[*He_i+1]- inHe)/(Hefracs[*He_i+1]-Hefracs[*He_i]);
  }    
    
  if (cooling_redshift_index != -1)
  {  
    if (inn_h <n_hgrid[0])
    {
      *d_n_h = 1.0;
      *n_h_i = 0;
      puts("Density out of range, using lowest value");
    }
    else if (inn_h > n_hgrid[nhdens-1])
    {
      *d_n_h = 0.0;
      *n_h_i = nhdens-2;
      puts("Density out of range, using highest value");
    }
    else
    {
      *n_h_i = (int)floor((log10(inn_h) - log10(n_hgrid[0]))/(log10(n_hgrid[1]) - log10(n_hgrid[0])));
      *d_n_h = (log10(n_hgrid[*n_h_i+1]) - log10(inn_h))/(log10(n_hgrid[*n_h_i+1])-log10(n_hgrid[*n_h_i]));
    }    
  }
  else /// No density dependence for the collisional cooling ///
  { 
    *n_h_i = 0;
    *d_n_h = 0.0;
  }
 
  if (u < ugrid[0])
  {
    d_u = 1.0;
    u_i = 0;
    puts("Energy Density out of range, using lowest value");
  }
  else if (u > ugrid[ntemps-1])
  {
    d_u = 0.0;
    u_i = ntemps-2;
    puts("Energy Density out of range, using highest value");
  }
  else
  {
    u_i = (int)floor((log10(u) - log10(ugrid[0]))/(log10(ugrid[1]) - log10(ugrid[0])));
    d_u = (log10(ugrid[u_i+1]) - log10(u))/(log10(ugrid[u_i+1])-log10(ugrid[u_i]));
  }    

  /// returns log_10 of T!! ///
  if (cooling_redshift_index != -1)
  {
    T =       (*d_He)   *   d_z   *   (*d_n_h)   *   d_u      *u_to_t[0][*He_i][*n_h_i][u_i]
            + (*d_He)   *   d_z   *   (*d_n_h)   *   (1.0-d_u)*u_to_t[0][*He_i][*n_h_i][u_i+1]
            + (*d_He)   *   d_z   *   (1.0-(*d_n_h))*d_u   *   u_to_t[0][*He_i][(*n_h_i)+1][u_i]
            + (*d_He)   *   d_z   *   (1.0-(*d_n_h))*(1.0-d_u)*u_to_t[0][*He_i][(*n_h_i)+1][u_i+1]
            + (*d_He)   *   (1.0-d_z)*(*d_n_h)   *   d_u   *   u_to_t[1][*He_i][*n_h_i][u_i]
 	    + (*d_He)   *   (1.0-d_z)*(*d_n_h)   *   (1.0-d_u)*u_to_t[1][*He_i][*n_h_i][u_i+1]
	    + (*d_He)   *   (1.0-d_z)*(1.0-(*d_n_h))*d_u   *   u_to_t[1][*He_i][(*n_h_i)+1][u_i]
  	    + (*d_He)   *   (1.0-d_z)*(1.0-(*d_n_h))*(1.0-d_u)*u_to_t[1][*He_i][(*n_h_i)+1][u_i+1]
	    + (1.0-(*d_He))*d_z   *   (*d_n_h)   *   d_u   *   u_to_t[0][(*He_i)+1][*n_h_i][u_i]
	    + (1.0-(*d_He))*d_z   *   (*d_n_h)   *   (1.0-d_u)*u_to_t[0][(*He_i)+1][*n_h_i][u_i+1]
	    + (1.0-(*d_He))*d_z   *   (1.0-(*d_n_h))*d_u   *   u_to_t[0][(*He_i)+1][(*n_h_i)+1][u_i]
	    + (1.0-(*d_He))*d_z   *   (1.0-(*d_n_h))*(1.0-d_u)*u_to_t[0][(*He_i)+1][(*n_h_i)+1][u_i+1]
	    + (1.0-(*d_He))*(1.0-d_z)*(*d_n_h)   *   d_u   *   u_to_t[1][(*He_i)+1][*n_h_i][u_i]
	    + (1.0-(*d_He))*(1.0-d_z)*(*d_n_h)   *   (1.0-d_u)*u_to_t[1][(*He_i)+1][*n_h_i][u_i+1]
	    + (1.0-(*d_He))*(1.0-d_z)*(1.0-(*d_n_h))*d_u   *   u_to_t[1][(*He_i)+1][(*n_h_i)+1][u_i]
	    + (1.0-(*d_He))*(1.0-d_z)*(1.0-(*d_n_h))*(1.0-d_u)*u_to_t[1][(*He_i)+1][(*n_h_i)+1][u_i+1];
  }
  else
    T =       (*d_He)   *   d_u      *u_to_t[0][*He_i][0][u_i]
            + (*d_He)   *   (1.0-d_u)*u_to_t[0][*He_i][0][u_i+1]
            + (1.0-(*d_He))*d_u   *   u_to_t[0][(*He_i)+1][0][u_i]
            + (1.0-(*d_He))*(1.0-d_u)*u_to_t[0][(*He_i)+1][0][u_i+1];

  return T;
}

///  Calculates (heating rate-cooling rate)/n_h^2 in cgs units  ///
double CoolingRate(int temp_i, int n_h_i, int He_i, float d_temp, float d_n_h, float d_He, float d_z, float T, float n_h, double inZ, double particle_Z[])
{
  double Lamb_he1, Lamb_he2, T_gam, LambdaNet = 0.0, n_e, n_e_he1, n_e_he2;
  float fact;
  int i;

  // Heating - Cooling is computed by adding up all the contributions from the different elements
  //Metal-free cooling//
  if (cooling_redshift_index != -1) ///Radiative///
  {
    Lamb_he1 = threewayinterp(nCoolHeats+He_i, temp_i, n_h_i, d_temp, d_n_h, d_z); 
    Lamb_he2 = threewayinterp(nCoolHeats+He_i+1, temp_i, n_h_i, d_temp, d_n_h, d_z);
  }
  else ///Collisional and compton cooling///
  {
    n_e_he1 = d_temp*collisional_electron_abundance[He_i][temp_i]+(1.0-d_temp)*collisional_electron_abundance[He_i][temp_i+1];
    n_e_he2 = d_temp*collisional_electron_abundance[He_i+1][temp_i]+(1.0-d_temp)*collisional_electron_abundance[He_i+1][temp_i+1];
    n_e = (d_He*n_e_he1 +(1.0 - d_He)*n_e_he2)*n_h;
    T_gam = T_CMB*(1+inZ);
    LambdaNet += COMP_COEFF*pow(T_gam,4.0)*n_e*(T_gam-T)/(pow(n_h,2.0)); /// Inverse Compton Scattering ///
    Lamb_he1 = d_temp*net_heating[0][nCoolHeats+He_i][0][temp_i]+(1.0-d_temp)*net_heating[0][nCoolHeats+He_i][0][temp_i+1];
    Lamb_he2 = d_temp*net_heating[0][nCoolHeats+He_i+1][0][temp_i]+(1.0-d_temp)*net_heating[0][nCoolHeats+He_i+1][0][temp_i+1];
  }
  LambdaNet += d_He*Lamb_he1 +(1.0 - d_He)*Lamb_he2;

  /// For each element, find the abundance and multiply it by the interpolated heating-cooling ///
  for (i = 0; i < nCoolHeats; i++)
  {
    fact = 0.0;
    if (SPH_Name_Pointers[i] >= 0)
      fact = particle_Z[SPH_Name_Pointers[i]]/solar_abundances[Solar_Name_Pointers[i]];
    else if (SPH_Name_Pointers[i] > -999)
    {
      if(strcmp(CoolHeat_Element_Names[i], "Sulphur") == 0)
	fact = SULPSILI*particle_Z[abs(SPH_Name_Pointers[i])]/solar_abundances[Solar_Name_Pointers[i]];
      else if(strcmp(CoolHeat_Element_Names[i], "Calcium") == 0)
	fact = CALCSILI*particle_Z[abs(SPH_Name_Pointers[i])]/solar_abundances[Solar_Name_Pointers[i]];
    }	  
    
    if (cooling_redshift_index != -1)
      LambdaNet += threewayinterp(i, temp_i, n_h_i, d_temp, d_n_h, d_z)*fact;	  
    else
      LambdaNet += (d_temp*net_heating[0][i][0][temp_i]+(1.0-d_temp)*net_heating[0][i][0][temp_i+1])*fact;
  }
	  
  return (LambdaNet);
}

/// General routine to calculate the cooling from the energy density ///
double CoolingRateFromU(double u, double inn_h, double inZ, float d_z, double particle_Z[])
{
  double temp;
  int temp_i, n_h_i, He_i;
  float d_n_h, d_He, d_temp;

  temp = convert_u_to_temp(&n_h_i, &He_i, &d_n_h, d_z, &d_He, u, inn_h, Hefrac);

  if (temp <log10(tempgrid[0])) /// Get temperature index ///
  {
    d_temp = 1.0;
    temp_i = 0;
  }
  else if (temp > log10(tempgrid[ntemps-1]))
  {
    d_temp = 0.0;
    temp_i = ntemps-2;
  }
  else
  {
    temp_i = (int)floor((temp - log10(tempgrid[0]))/(log10(tempgrid[1]) - log10(tempgrid[0])));
    d_temp = (log10(tempgrid[temp_i+1])- temp)/(log10(tempgrid[temp_i+1])-log10(tempgrid[temp_i]));
  }    

  return CoolingRate(temp_i, n_h_i, He_i, d_temp, d_n_h, d_He, d_z, pow(10.0, temp), inn_h, inZ, particle_Z);;
}


double DoCooling(double u_old, double rho, double dt, double inz, float d_z, double particle_Z[])
{
  float inn_h, XH;
  double du, ratefact, u, u_upper, u_lower, LambdaNet;
  int iter;
  

  /// Find the Hydrogen and Helium fractions ///
  for (iter = 0; iter <BG_NELEMENTS; iter++)
    {
      if (strcmp(SPH_Element_Name[iter], "Hydrogen") == 0)
	XH = particle_Z[iter];

      if (strcmp(SPH_Element_Name[iter], "Helium") == 0)
	Hefrac = particle_Z[iter];
    }

  /// initialize some variables ///
  inn_h = rho*XH/PROTONMASS;
  ratefact = inn_h * inn_h / rho;
  
  u = u_old;
  u_lower = u;
  u_upper = u;
  
  
  /// Iterative, implicit cooling ///
  LambdaNet = CoolingRateFromU(u, inn_h, inz, d_z, particle_Z);

  // first check whether cooling rate is worth taking into account.
  du = ratefact * LambdaNet * dt;
  if(fabs(du / u) < 1.0e-6)
  {
    return u;
  }

  // bracketing 
  if(u - u_old - ratefact * LambdaNet * dt < 0)	// heating 
  {
    u_upper *= sqrt(1.1);
    u_lower /= sqrt(1.1);
    while(u_upper - u_old - ratefact * CoolingRateFromU(u_upper, inn_h, inz,  d_z, particle_Z) * dt < 0)
    {
      u_upper *= 1.1;
      u_lower *= 1.1;
    }    
  }
  
  if(u - u_old - ratefact * LambdaNet * dt > 0) // cooling
  {
    u_lower /= sqrt(1.1);
    u_upper *= sqrt(1.1);
    while(u_lower - u_old - ratefact * CoolingRateFromU(u_lower, inn_h, inz, d_z, particle_Z) * dt > 0)
    {
      u_upper /= 1.1;
      u_lower /= 1.1;
    }
  }

  iter = 0;

  do // iterate to convergence
  {
    u = 0.5 * (u_lower + u_upper);
    
    LambdaNet = CoolingRateFromU(u, inn_h, inz, d_z, particle_Z);
    
    if(u - u_old - ratefact * LambdaNet * dt > 0)
    {
      u_upper = u;
    }
    else
    {
      u_lower = u;
    }
    
    du = u_upper - u_lower;
    
    iter++;
    
    if(iter >= (MAXITER - 10))
      printf("u= %g\n", u);
  }
  while(fabs(du / u) > 1.0e-6 && iter < MAXITER);
  
  if(iter >= MAXITER)
    printf("failed to converge in DoCooling()\n");
	
  return u;
}

 
#endif  /* closes #ifdef BG_COOLING at the beginning of the file */
