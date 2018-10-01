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


/// Get the redshifts from a text file ///
void GetCoolingRedshifts()
{
  FILE *infile;
  int i = 0;
  char buffer[100], redfilename [100];
  
  sprintf(redfilename, "%sredshifts.dat", All.CoolTablePath);
  infile = fopen(redfilename, "r");
  if (infile == NULL) puts("GetCoolingRedshifts can't open a file");

  fscanf(infile, "%s", buffer);
  nZs = atoi(buffer);
  cooling_redshifts = (float *) malloc(nZs * sizeof(float));
  
  while (fscanf(infile, "%s",  buffer) != EOF )
  {
    cooling_redshifts[i] = atof(buffer); /// fix this to malloc change
    i += 1;
  }
  fclose(infile);

}
 
/// The header file contains all sorts of vital information (the temperature grid, etc.) ///
void ReadCoolingHeader(char *fname)
{
  int i;
  hid_t tempfile_id, dataset, datatype;
  herr_t status;

  /// Fill the constants ///
  tempfile_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  dataset = H5Dopen(tempfile_id, "/Header/Number_of_density_bins");
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nhdens);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Header/Number_of_temperature_bins");
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &ntemps);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Header/Helium/Number_of_helium_fractions");
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nHes);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Header/Number_of_species");
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nCoolHeats);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Header/Abundances/Number_of_abundances");
  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nSolar_Abundances);
  status = H5Dclose(dataset);

  /// Allocate header arrays ///
  tempgrid = (float *) malloc(ntemps * sizeof(float));
  n_hgrid = (float *) malloc(nhdens * sizeof(float));
  ugrid = (float * ) malloc(ntemps * sizeof(float));
  Hefracs = (float *) malloc(nHes * sizeof(float));
  solar_abundances = (float *) malloc(nSolar_Abundances * sizeof(float));

  Henames = (char **) malloc(nHes * sizeof(char*));
  for (i = 0; i < nHes; i++)
    Henames[i] = (char *) malloc(EL_NAME_LENGTH * sizeof(char));

  CoolHeat_Element_Names = (char **) malloc(nCoolHeats * sizeof(char*));
  for (i = 0; i < nCoolHeats; i++)
    CoolHeat_Element_Names[i] = (char *) malloc(EL_NAME_LENGTH * sizeof(char));

  Solar_Abundance_Names = (char **) malloc(nSolar_Abundances * sizeof(char*));
  for (i = 0; i < nSolar_Abundances; i++)
    Solar_Abundance_Names[i] = (char *) malloc(EL_NAME_LENGTH * sizeof(char));

  /// Fill the arrays ///
  dataset = H5Dopen(tempfile_id, "/Header/Temperature_bins");
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, tempgrid);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Header/Hydrogen_density_bins");
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_hgrid);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Header/Energy_density_bins");
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ugrid);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Header/Helium/Helium_fraction_bins");
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, Hefracs);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Header/Abundances/Solar_abundances");
  status = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, solar_abundances);
  status = H5Dclose(dataset);

  datatype = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(datatype, H5T_VARIABLE);
  dataset = H5Dopen(tempfile_id, "/Header/Species_names");
  status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, CoolHeat_Element_Names);  
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Header/Helium/Helium_names");
  status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, Henames);
  status = H5Dclose(dataset);

  dataset = H5Dopen(tempfile_id, "/Header/Abundances/Abund_names");
  status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, Solar_Abundance_Names);
  status = H5Dclose(dataset);

  status = H5Fclose(tempfile_id);
  /// Done ///

}

void BroadCastHeader()
{
  int position, bufsize, size, i;
  char *buffer, *buffer1;

  /// get size of the buffer ///
  MPI_Pack_size(6, MPI_INT, MPI_COMM_WORLD, &size);
  bufsize = size;

  /// allocate memory for the buffer ///
  buffer = malloc(bufsize);

  if (ThisTask == 0)
  {
    position = 0;
    /// Pack the array dimensions ///
    MPI_Pack(&nhdens, 1, MPI_INT, buffer, bufsize, &position, MPI_COMM_WORLD);
    MPI_Pack(&ntemps, 1, MPI_INT, buffer, bufsize, &position, MPI_COMM_WORLD);
    MPI_Pack(&nHes, 1, MPI_INT, buffer, bufsize, &position, MPI_COMM_WORLD);
    MPI_Pack(&nCoolHeats, 1, MPI_INT, buffer, bufsize, &position, MPI_COMM_WORLD);
    MPI_Pack(&nSolar_Abundances, 1, MPI_INT, buffer, bufsize, &position, MPI_COMM_WORLD);
    MPI_Pack(&nZs, 1, MPI_INT, buffer, bufsize, &position, MPI_COMM_WORLD);
  }

  MPI_Bcast(buffer, bufsize, MPI_PACKED, 0, MPI_COMM_WORLD);
  
  if (ThisTask != 0)
  {
    position = 0;

    /// Unpack the array dimensions ///
    MPI_Unpack(buffer, bufsize, &position, &nhdens, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(buffer, bufsize, &position, &ntemps, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(buffer, bufsize, &position, &nHes, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(buffer, bufsize, &position, &nCoolHeats, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(buffer, bufsize, &position, &nSolar_Abundances, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(buffer, bufsize, &position, &nZs, 1, MPI_INT, MPI_COMM_WORLD);

    /// Allocate header arrays ///
    tempgrid = (float *) malloc(ntemps * sizeof(float));
    n_hgrid = (float *) malloc(nhdens * sizeof(float));
    ugrid = (float * ) malloc(ntemps * sizeof(float));
    Hefracs = (float *) malloc(nHes * sizeof(float));
    solar_abundances = (float *) malloc(nSolar_Abundances * sizeof(float));

    cooling_redshifts = (float *) malloc(nZs*sizeof(float));

    Henames = (char **) malloc(nHes * sizeof(char*));
    for (i = 0; i < nHes; i++)
      Henames[i] = (char *) malloc(EL_NAME_LENGTH * sizeof(char));

    CoolHeat_Element_Names = (char **) malloc(nCoolHeats * sizeof(char*));
    for (i = 0; i < nCoolHeats; i++)
      CoolHeat_Element_Names[i] = (char *) malloc(EL_NAME_LENGTH * sizeof(char));

    Solar_Abundance_Names = (char **) malloc(nSolar_Abundances * sizeof(char*));
    for (i = 0; i < nSolar_Abundances; i++)
      Solar_Abundance_Names[i] = (char *) malloc(EL_NAME_LENGTH * sizeof(char));
  }

  /// free allocated memory ///
  free(buffer);

  bufsize = 0;
  /// get size of the buffer ///
  MPI_Pack_size(ntemps, MPI_FLOAT, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(nhdens, MPI_FLOAT, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(ntemps, MPI_FLOAT, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(nHes, MPI_FLOAT, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(nSolar_Abundances, MPI_FLOAT, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(nZs, MPI_FLOAT, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(nHes*EL_NAME_LENGTH, MPI_CHAR, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(nCoolHeats*EL_NAME_LENGTH, MPI_CHAR, MPI_COMM_WORLD, &size);
  bufsize += size;
  MPI_Pack_size(nSolar_Abundances*EL_NAME_LENGTH, MPI_CHAR, MPI_COMM_WORLD, &size);
  bufsize += size;

  /// allocate memory for the buffer ///
  buffer = malloc(bufsize);

  if (ThisTask == 0)
  {
    position = 0;
    MPI_Pack(tempgrid, ntemps, MPI_FLOAT, buffer, bufsize, &position, MPI_COMM_WORLD);
    MPI_Pack(n_hgrid, nhdens, MPI_FLOAT, buffer, bufsize, &position, MPI_COMM_WORLD);
    MPI_Pack(ugrid, ntemps, MPI_FLOAT, buffer, bufsize, &position, MPI_COMM_WORLD);
    MPI_Pack(Hefracs, nHes, MPI_FLOAT, buffer, bufsize, &position, MPI_COMM_WORLD);
    MPI_Pack(solar_abundances, nSolar_Abundances, MPI_FLOAT, buffer, bufsize, &position, MPI_COMM_WORLD);
    MPI_Pack(cooling_redshifts, nZs, MPI_FLOAT, buffer, bufsize, &position, MPI_COMM_WORLD);

    for (i = 0; i < nHes; i++)
      MPI_Pack(Henames[i], EL_NAME_LENGTH, MPI_CHAR, buffer, bufsize, &position, MPI_COMM_WORLD);
    for (i = 0; i < nCoolHeats; i++)
      MPI_Pack(CoolHeat_Element_Names[i], EL_NAME_LENGTH, MPI_CHAR, buffer, bufsize, &position, MPI_COMM_WORLD);
    for (i = 0; i < nSolar_Abundances; i++)
      MPI_Pack(Solar_Abundance_Names[i], EL_NAME_LENGTH, MPI_CHAR, buffer, bufsize, &position, MPI_COMM_WORLD);
  }

  MPI_Bcast(buffer, bufsize, MPI_PACKED, 0, MPI_COMM_WORLD);

  if (ThisTask != 0)
  {
    position = 0;

    MPI_Unpack(buffer, bufsize, &position, tempgrid, ntemps, MPI_FLOAT, MPI_COMM_WORLD);
    MPI_Unpack(buffer, bufsize, &position, n_hgrid, nhdens, MPI_FLOAT, MPI_COMM_WORLD);
    MPI_Unpack(buffer, bufsize, &position, ugrid, ntemps, MPI_FLOAT, MPI_COMM_WORLD);
    MPI_Unpack(buffer, bufsize, &position, Hefracs, nHes, MPI_FLOAT, MPI_COMM_WORLD);
    MPI_Unpack(buffer, bufsize, &position, solar_abundances, nSolar_Abundances, MPI_FLOAT, MPI_COMM_WORLD);
    MPI_Unpack(buffer, bufsize, &position, cooling_redshifts, nZs, MPI_FLOAT, MPI_COMM_WORLD);

    for (i = 0; i < nHes; i++)
      MPI_Unpack(buffer, bufsize, &position, Henames[i], EL_NAME_LENGTH, MPI_CHAR, MPI_COMM_WORLD);
    for (i = 0; i < nCoolHeats; i++)
      MPI_Unpack(buffer, bufsize, &position, CoolHeat_Element_Names[i], EL_NAME_LENGTH, MPI_CHAR, MPI_COMM_WORLD);
    for (i = 0; i < nSolar_Abundances; i++)
      MPI_Unpack(buffer, bufsize, &position, Solar_Abundance_Names[i], EL_NAME_LENGTH, MPI_CHAR, MPI_COMM_WORLD);
  }

  /// free allocated memory ///
  free(buffer);
}

void MakeNamePointers()
{
  int i, j, sili_index;

  SPH_Name_Pointers = (int *) malloc(nCoolHeats * sizeof(int));
  Solar_Name_Pointers = (int *) malloc(nCoolHeats * sizeof(int));

  for (i = 0; i < BG_NELEMENTS; i++)
  {
    if (strcmp(SPH_Element_Name[i], "Silicon") == 0)
      sili_index = i;
  }

  for (i = 0; i < nCoolHeats; i++)
  {
    Solar_Name_Pointers[i] = -999;
    SPH_Name_Pointers[i] = -999;
    for (j = 0; j < nSolar_Abundances; j++)
    {
      if (strcmp(CoolHeat_Element_Names[i], Solar_Abundance_Names[j]) == 0)
	Solar_Name_Pointers[i] = j;
    }

    if(strcmp(CoolHeat_Element_Names[i], "Sulphur") == 0 || strcmp(CoolHeat_Element_Names[i], "Calcium") == 0) /// These elements are tracked!
      SPH_Name_Pointers[i] = -1*sili_index;
    else
    {
      for (j = 0; j < BG_NELEMENTS; j++)
      {
	if (strcmp(CoolHeat_Element_Names[i], SPH_Element_Name[j]) == 0)
	  SPH_Name_Pointers[i] = j;
      }
    }
  }

#ifdef BG_VERBOSE
  if (ThisTask == 0)
    {
      for (i = 0; i < nCoolHeats; i++)
	printf("CoolHeat_Element_Names[%d] is Solar_Abundance_Name %s\n", i, Solar_Abundance_Names[Solar_Name_Pointers[i]]);

      for (i = 0; i < nCoolHeats; i++)
	printf("CoolHeat_Element_Names[%d] is SPH_Element_Name %s\n", i, SPH_Element_Name[SPH_Name_Pointers[i]]);
    }
#endif
}

/// convert abundances, get cooling_redshifts ///
void InitCool()
{
  int i, j, k;
  char tempfname[100];

  /// Initialise element names array ///
  InitChemistry();

  /// Use a general file to get the header ///
  if (ThisTask == 0)
  {
    GetCoolingRedshifts();
    sprintf(tempfname, "%sz_0.000.hdf5", All.CoolTablePath); 
    ReadCoolingHeader(tempfname);
  }

  BroadCastHeader();

  net_heating = (float ****)malloc(2*sizeof(float***));
  u_to_t = (float ****)malloc(2*sizeof(float***));

  for (i = 0; i < 2; i++)
  {
    net_heating[i] = (float ***) malloc((nHes + nCoolHeats) * sizeof(float**));
    u_to_t[i] = (float ***) malloc(nHes * sizeof(float**));
    collisional_electron_abundance = (float **) malloc(nHes * sizeof(float*));
    for (j = 0; j < (nHes + nCoolHeats); j++)
    {
      net_heating[i][j] = (float **) malloc(nhdens * sizeof(float*));
      if (j < nHes)
      {
	u_to_t[i][j] = (float **) malloc(nhdens * sizeof(float*));
	collisional_electron_abundance[j] = (float* ) malloc(ntemps * sizeof(float));
      }
      for (k = 0; k < nhdens; k++)
      {
	net_heating[i][j][k] = (float *) malloc(ntemps * sizeof(float));
	if (j < nHes)
	  u_to_t[i][j][k] = (float *) malloc(ntemps * sizeof(float));
      }
    }
  }

  MakeNamePointers();
} 

void InitChemistry(void)
{
  int i;

  /* copy element names - element names should be read-in from a file */
  strcpy(SPH_Element_Name[0],"Hydrogen");
  strcpy(SPH_Element_Name[1],"Helium");
  strcpy(SPH_Element_Name[2],"Carbon");
  strcpy(SPH_Element_Name[3],"Nitrogen");
  strcpy(SPH_Element_Name[4],"Oxygen ");
  strcpy(SPH_Element_Name[5],"Neon, ");
  strcpy(SPH_Element_Name[6],"Magnesium");
  strcpy(SPH_Element_Name[7],"Silicon");
  strcpy(SPH_Element_Name[8],"Iron");
}

void allocate_header_arrays(void)
{
  int i;

  // bins for interpolation in log10(T [K])
  tempgrid = (float *) malloc(ntemps * sizeof(float));
  
  // bins for interpolation in hydrogen number density, log10(nh  [/cm^3]) 
  n_hgrid = (float *) malloc(nhdens * sizeof(float));
  
  // bins for interpolation in thermal energy per unit mass, log10(u [erg/g])
  ugrid = (float *) malloc(ntemps * sizeof(float));
  
  // bins for interpolation in Helium abundance
  Hefracs = (float *) malloc(nHes * sizeof(float));
  
  // assumed solar abundances used in constructing the tables, and corresponding names
  Solar_Abundance_Names = (char **) malloc(nSolar_Abundances * sizeof(char*));
  for (i = 0; i < nSolar_Abundances; i++)
    Solar_Abundance_Names[i] = (char *) malloc(EL_NAME_LENGTH * sizeof(char));

  solar_abundances = (float *) malloc(nSolar_Abundances * sizeof(float));
  
  // names of groups storing Helium abundances
  Henames = (char **) malloc(nHes * sizeof(char*));
  for (i = 0; i < nHes; i++)
    Henames[i] = (char *) malloc(EL_NAME_LENGTH * sizeof(char));

  // names of chemical elements stored in table
  CoolHeat_Element_Names = (char **) malloc(nCoolHeats * sizeof(char*));
  for(i = 0; i < nCoolHeats; i++)
    CoolHeat_Element_Names[i] = (char *) malloc(EL_NAME_LENGTH * sizeof(char));
}

#endif /* closes #ifdef BG_COOLING at the beginning of the file */
