/// This contains functions pertaining to reading in cooling tables, and finding redshifts ///
#ifdef BG_COOLING

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string.h>
#include <hdf5.h>
#include <mpi.h>
#include "allvars.h"
#include "bg_cooling.h"


void get_redshift_index(float z, int z_index_old, float *dz)
{
  // assumes cooling_redshifts table is in increasing order.

  if(z > cooling_redshifts[nZs-1]) /// before the earliest redshift, flag for collisional cooling
  {
    cooling_redshift_index = -1;
    *dz = 0.0;
  }
  else if (z < cooling_redshifts[0]) /// at the end, just use the last value
  {
    *dz = 1.0;
    cooling_redshift_index = 0;
  }
  else
  {
    if (z_index_old < 0)
      z_index_old = nZs - 1;

    *dz = 0.0;

    for (cooling_redshift_index = z_index_old; cooling_redshift_index > 0; cooling_redshift_index--) /// start at the previous index and search
    {
      if (z > cooling_redshifts[cooling_redshift_index])
      {
	*dz = (cooling_redshifts[cooling_redshift_index+1] - z)/(cooling_redshifts[cooling_redshift_index+1] - cooling_redshifts[cooling_redshift_index]);
	break;
      }
    }
  }
}


/// Get the cooling table for collisional cooling (before redshift ~8) ///
void GetCollisTable()
{
  hid_t tempfile_id, dataset;
  herr_t status;
  char fname[100], set_name[100];
  int Hes, specs, j;
  float hgrid[ntemps], cgrid[ntemps], temp[ntemps], n_e[ntemps];
  
  sprintf(fname, "%sz_collis.hdf5", All.CoolTablePath);

  tempfile_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  // For normal elements ///
  for (specs = 0; specs<nCoolHeats; specs++)
  {
    sprintf(set_name, "/%s/Heating", CoolHeat_Element_Names[specs]);
    dataset = H5Dopen(tempfile_id, set_name);
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, hgrid);
    status = H5Dclose(dataset);
 
    sprintf(set_name, "/%s/Cooling", CoolHeat_Element_Names[specs]);
    dataset = H5Dopen(tempfile_id, set_name);
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cgrid);
    status = H5Dclose(dataset);
    for (j = 0; j< ntemps; j++)
    {
      net_heating[0][specs][0][j] = hgrid[j] -cgrid[j];
    }
  }

  /// Helium ///
  for (Hes =0; Hes < nHes; Hes++) 
  {
    sprintf(set_name, "/Metal_free/%s/Heating", Henames[Hes]);
    dataset = H5Dopen(tempfile_id, set_name);
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, hgrid);
    status = H5Dclose(dataset);
 
    sprintf(set_name, "/Metal_free/%s/Cooling", Henames[Hes]);
    dataset = H5Dopen(tempfile_id, set_name);
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cgrid);
    status = H5Dclose(dataset);

    sprintf(set_name, "/Metal_free/%s/Temperature", Henames[Hes]);
    dataset = H5Dopen(tempfile_id, set_name);
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp);
    status = H5Dclose(dataset);

    sprintf(set_name, "/Metal_free/%s/Electron_density_over_n_h", Henames[Hes]);
    dataset = H5Dopen(tempfile_id, set_name);
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_e);
    status = H5Dclose(dataset);

    for (j = 0; j< ntemps; j++)
    {
      net_heating[0][nCoolHeats+Hes][0][j] = hgrid[j] -cgrid[j];
      u_to_t[0][Hes][0][j] = log10(temp[j]);
      collisional_electron_abundance[Hes][j] = n_e[j];

    }
  }
  status = H5Fclose(tempfile_id);
    
}

/// Get the cooling tables that bound the given redshift ///
void GetCoolingTables(int z_index_old)
{
  hid_t tempfile_id, dataset;
  herr_t status;
  char fname[100], set_name[100];
  int Hes, specs, i, j, skipme = 0;
  float hgrid[nhdens][ntemps], cgrid[nhdens][ntemps], temp[nhdens][ntemps], heattune_factor;

  /// check to see if we can use the last one ///
  if(cooling_redshift_index == z_index_old -1) 
  {
    skipme = 1;
    for (specs = 0; specs<nCoolHeats+nHes; specs++)
    {
      for (i =0; i <nhdens; i++)
      {
	for (j = 0; j< ntemps; j++)
	{
	  net_heating[1][specs][i][j] = net_heating[0][specs][i][j];
	  if (specs<nHes)
	    u_to_t[1][specs][i][j] = u_to_t[0][specs][i][j];
	}
      }
    }
  }
 
  
  sprintf(fname, "%sz_%1.3f.hdf5", All.CoolTablePath, cooling_redshifts[cooling_redshift_index]);

  tempfile_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  dataset = H5Dopen(tempfile_id, "/Header/Heat_tuning_factor");
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &heattune_factor);
  status = H5Dclose(dataset);

  // Redshift number one
  // For normal elements
  for (specs = 0; specs<nCoolHeats; specs++)
  {
    sprintf(set_name, "/%s/Heating", CoolHeat_Element_Names[specs]);
    dataset = H5Dopen(tempfile_id, set_name);
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, hgrid);
    status = H5Dclose(dataset);
 
    sprintf(set_name, "/%s/Cooling", CoolHeat_Element_Names[specs]);
    dataset = H5Dopen(tempfile_id, set_name);
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cgrid);
    status = H5Dclose(dataset);
    for (i =0; i <nhdens; i++)
    {
      for (j = 0; j< ntemps; j++)
      {
	net_heating[0][specs][i][j] = hgrid[i][j] -cgrid[i][j];
      }
    }
  }

  // Helium    
  for (Hes =0; Hes < nHes; Hes++) 
  {
    sprintf(set_name, "/Metal_free/%s/Heating", Henames[Hes]);
    dataset = H5Dopen(tempfile_id, set_name);
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, hgrid);
    status = H5Dclose(dataset);
 
    sprintf(set_name, "/Metal_free/%s/Cooling", Henames[Hes]);
    dataset = H5Dopen(tempfile_id, set_name);
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cgrid);
    status = H5Dclose(dataset);

    sprintf(set_name, "/Metal_free/%s/Temperature", Henames[Hes]);
    dataset = H5Dopen(tempfile_id, set_name);
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp);
    status = H5Dclose(dataset);

    for (i =0; i <nhdens; i++)
    {
      for (j = 0; j< ntemps; j++)
      {
	/// multiplied by a tuning factor to account for reionization ///
	net_heating[0][nCoolHeats+Hes][i][j] = hgrid[i][j] *heattune_factor -cgrid[i][j];
	u_to_t[0][Hes][i][j] = log10(temp[i][j]);
      }
    }
  }
  status = H5Fclose(tempfile_id);

  ///redshift two///
  if(skipme == 0)
  {
    sprintf(fname, "%sz_%1.3f.hdf5", All.CoolTablePath, cooling_redshifts[cooling_redshift_index+1]);
    
    tempfile_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    
    dataset = H5Dopen(tempfile_id, "/Header/Heat_tuning_factor");
    status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &heattune_factor);
    status = H5Dclose(dataset);
    
    /// For normal elements ///
    for (specs = 0; specs<nCoolHeats; specs++)
    {
      sprintf(set_name, "/%s/Heating", CoolHeat_Element_Names[specs]);
      dataset = H5Dopen(tempfile_id, set_name);
      status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, hgrid);
      status = H5Dclose(dataset);
      
      sprintf(set_name, "/%s/Cooling", CoolHeat_Element_Names[specs]);
      dataset = H5Dopen(tempfile_id, set_name);
      status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cgrid);
      status = H5Dclose(dataset);
      
      for (i =0; i <nhdens; i++)
      { 
	for (j = 0; j< ntemps; j++)
	  net_heating[1][specs][i][j] = hgrid[i][j] -cgrid[i][j];
      }
    }

    /// Helium ///
    for (Hes =0; Hes < nHes; Hes++) 
    {
      sprintf(set_name, "/Metal_free/%s/Heating", Henames[Hes]);
      dataset = H5Dopen(tempfile_id, set_name);
      status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, hgrid);
      status = H5Dclose(dataset);
      
      sprintf(set_name, "/Metal_free/%s/Cooling", Henames[Hes]);
      dataset = H5Dopen(tempfile_id, set_name);
      status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cgrid);
      status = H5Dclose(dataset);
      
      sprintf(set_name, "/Metal_free/%s/Temperature", Henames[Hes]);
      dataset = H5Dopen(tempfile_id, set_name);
      status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp);
      status = H5Dclose(dataset);
      
      for (i =0; i <nhdens; i++)
      {
	for (j = 0; j< ntemps; j++)
	{
	  /// multiplied by a tuning factor to account for reionization ///
	  net_heating[1][nCoolHeats+Hes][i][j]  = hgrid[i][j]*heattune_factor -cgrid[i][j]; 
	  u_to_t[1][Hes][i][j] = log10(temp[i][j]);
	}
      }
    }
    status = H5Fclose(tempfile_id);
  }
    
}

void BroadcastCoolingTables()
{
  int position, bufsize, size, i, j, k;
  char *buffer;

  /// get size of the buffer ///
  MPI_Pack_size(2*ntemps*nhdens*(nCoolHeats+nHes), MPI_FLOAT, MPI_COMM_WORLD, &size);
  bufsize = size;
  MPI_Pack_size(2*ntemps*nhdens*nHes, MPI_FLOAT, MPI_COMM_WORLD, &size);
  bufsize += size;
  if (cooling_redshift_index == -1)
  {
    MPI_Pack_size(ntemps*nHes, MPI_FLOAT, MPI_COMM_WORLD, &size);
    bufsize += size;
  }

  /// allocate memory for the buffer ///
  buffer = malloc(bufsize);

  if (ThisTask == 0)
  {
    position = 0;
    for (i = 0; i < 2; i++)
      for (j = 0; j < nCoolHeats+nHes; j++)
	for (k = 0; k < nhdens; k++)
	  MPI_Pack(net_heating[i][j][k], ntemps, MPI_FLOAT, buffer, bufsize, &position, MPI_COMM_WORLD);
    for (i = 0; i < 2; i++)
      for (j = 0; j < nHes; j++)
	for (k = 0; k < nhdens; k++)
	  MPI_Pack(u_to_t[i][j][k], ntemps, MPI_FLOAT, buffer, bufsize, &position, MPI_COMM_WORLD);
    if (cooling_redshift_index == -1)
      for (j = 0; j < nHes; j++)
	MPI_Pack(collisional_electron_abundance[j], ntemps, MPI_FLOAT, buffer, bufsize, &position, MPI_COMM_WORLD);
  }

  MPI_Bcast(buffer, bufsize, MPI_PACKED, 0, MPI_COMM_WORLD);

  if (ThisTask != 0)
  {
    position = 0;
    for (i = 0; i < 2; i++)
      for (j = 0; j < nCoolHeats+nHes; j++)
	for (k = 0; k < nhdens; k++)
	  MPI_Unpack(buffer, bufsize, &position, net_heating[i][j][k], ntemps, MPI_FLOAT, MPI_COMM_WORLD);

    for (i = 0; i < 2; i++)
      for (j = 0; j < nHes; j++)
	for (k = 0; k < nhdens; k++)
	  MPI_Unpack(buffer, bufsize, &position, u_to_t[i][j][k], ntemps, MPI_FLOAT, MPI_COMM_WORLD);

    for (j = 0; j < nHes; j++)
      if (cooling_redshift_index == -1)
	MPI_Unpack(buffer, bufsize, &position, collisional_electron_abundance[j], ntemps, MPI_FLOAT, MPI_COMM_WORLD);
  }

  free(buffer);
  
}

#endif /* closes #ifdef BG_COOLING at the beginning of the file */
