/* this file is appended to pm_periodic.c in case
 * HPM is selected. 
 */

void hpm_find_eqs(void)
{
  /* call one or the other */
  /*
     hpm_find_eqs_simple(); 
   */
  hpm_find_eqs_selfconsistent();
}

/* This function determines the present equation of state 
 * in case the HPM method is selected. It normally does this
 * by following the ionization history of two fiducial gas
 * elements at the mean density, and at 1.1 times the mean density.
 * From the pressure values, and effective EQS index is then computed.
 */
void hpm_find_eqs_simple(void)
{
  double a3inv;
  double u0, meanWeight, temp;

  /* supppose we want EQS with 20000 K, and alpha=1.2, below z<9 */

  if(All.ComovingIntegrationOn)
    a3inv = 1 / (All.Time * All.Time * All.Time);
  else
    a3inv = 1;


  All.HPM_rho0 = All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

  if(All.Time < 1.0 / (1 + 9))
    {
      temp = 0;
      All.HPM_P0 = 0;
      All.HPM_alpha = 1.2;
    }
  else
    {
      temp = 20000.0;
      meanWeight = 4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC)) * PROTONMASS;	/* note: assuming FULL ionization */
      u0 = temp / (meanWeight / BOLTZMANN * GAMMA_MINUS1) / (All.UnitEnergy_in_cgs / All.UnitMass_in_g);

      All.HPM_entr0 = u0 * GAMMA_MINUS1 / pow(All.HPM_rho0 * a3inv, GAMMA_MINUS1);
      All.HPM_P0 = All.HPM_entr0 * pow(All.HPM_rho0, GAMMA);
      All.HPM_alpha = 1.2;
    }

  if(ThisTask == 0)
    {
      printf("IGM-EQS:  a=%g  T0=%g   P0=%g  alpha=%g\n", All.Time, temp, All.HPM_P0, All.HPM_alpha);
      fflush(stdout);
    }
}


/* This function determines the present equation of state 
 * in case the HPM method is selected. It normally does this
 * by following the ionization history of two fiducial gas
 * elements at the mean density, and at 1.1 times the mean density.
 * From the pressure values, and effective EQS index is then computed.
 */
void hpm_find_eqs_selfconsistent(void)
{
  double dt, dtime, hubble_a = 0, a3inv;
  double time_hubble_a, dmax1, dmax2;
  double u0, u1, meanWeight, temp;

  if(All.ComovingIntegrationOn)
    {
      /* Factors for comoving integration of hydro */
      a3inv = 1 / (All.Time * All.Time * All.Time);
      hubble_a = hubble_function(All.Time);

      time_hubble_a = All.Time * hubble_a;
    }
  else
    a3inv = time_hubble_a = 1;

  dt = (P[0].Ti_endstep - P[0].Ti_begstep) * All.Timebase_interval;	/*  the time-step */

  if(All.ComovingIntegrationOn)
    dtime = All.Time * dt / time_hubble_a;
  else
    dtime = dt;


  All.HPM_rho0 = All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
  All.HPM_rho1 = 1.1 * All.HPM_rho0;


  u0 = All.HPM_entr0 / GAMMA_MINUS1 * pow(All.HPM_rho0 * a3inv, GAMMA_MINUS1);
  u1 = All.HPM_entr1 / GAMMA_MINUS1 * pow(All.HPM_rho1 * a3inv, GAMMA_MINUS1);

  u0 = DoCooling(DMAX(All.MinEgySpec, u0), All.HPM_rho0 * a3inv, dtime, &All.HPM_ne0);
  u1 = DoCooling(DMAX(All.MinEgySpec, u1), All.HPM_rho1 * a3inv, dtime, &All.HPM_ne1);

  All.HPM_entr0 = u0 * GAMMA_MINUS1 / pow(All.HPM_rho0 * a3inv, GAMMA_MINUS1);
  All.HPM_entr1 = u1 * GAMMA_MINUS1 / pow(All.HPM_rho1 * a3inv, GAMMA_MINUS1);

  /* Note: All.HPM_P0 is a "comoving pressure", i.e. computed in gadget's 
   * internal unit convention using the comoving density and the physical entropy.
   */

  All.HPM_P0 = All.HPM_entr0 * pow(All.HPM_rho0, GAMMA);
  All.HPM_P1 = All.HPM_entr1 * pow(All.HPM_rho1, GAMMA);

  All.HPM_alpha = log(All.HPM_P1 / All.HPM_P0) / log(All.HPM_rho1 / All.HPM_rho0);

  /* compute the temperature corresponding to P0, for information purposes only */

  meanWeight = 4.0 / (3 * HYDROGEN_MASSFRAC + 1 + 4 * HYDROGEN_MASSFRAC * All.HPM_ne0) * PROTONMASS;
  temp = meanWeight / BOLTZMANN * GAMMA_MINUS1 * u0 * All.UnitEnergy_in_cgs / All.UnitMass_in_g;

  if(ThisTask == 0)
    {
      printf("IGM-EQS:  a=%g  T0=%g   P0=%g  alpha=%g\n", All.Time, temp, All.HPM_P0, All.HPM_alpha);
      fflush(stdout);
    }
}




/* This function computes the pressure force in case the HPM method is selected
 */
void hpm_pressureforce(void)
{
  double k2, kx, ky, kz, smth;
  double dx, dy, dz;
  double fx, fy, fz, ff;
  double asmth2, facdens, facenthalpy, facdiff, acc_dim;
  int i, j, slab, level, sendTask, recvTask;
  int x, y, z, xl, yl, zl, xr, yr, zr, xll, yll, zll, xrr, yrr, zrr, ip, dim;
  int slab_x, slab_y, slab_z;
  int slab_xx, slab_yy, slab_zz;
  int meshmin[3], meshmax[3], sendmin, sendmax, recvmin, recvmax;
  int rep, ncont, cont_sendmin[2], cont_sendmax[2], cont_recvmin[2], cont_recvmax[2];
  int dimx, dimy, dimz, recv_dimx, recv_dimy, recv_dimz;
  MPI_Status status;


  if(ThisTask == 0)
    {
      printf("Starting HPM pressure force calculation.\n");
      fflush(stdout);
    }

  hpm_find_eqs();

  asmth2 = (2 * M_PI * HPM_SMTH) / PMGRID;
  asmth2 *= asmth2;

  facdens = 1 / (All.BoxSize * All.BoxSize * All.BoxSize);	/* to get density */

  if(All.HPM_alpha > 1.001)
    facenthalpy = All.HPM_alpha / (All.HPM_alpha - 1) * All.HPM_P0 / pow(All.HPM_rho0, All.HPM_alpha);
  else
    facenthalpy = 0;

  facdiff = 1 / (2 * All.BoxSize / PMGRID);	/* for finite differencing */

  /* first, establish the extension of the local patch in the PMGRID  */


  for(j = 0; j < 3; j++)
    {
      meshmin[j] = PMGRID;
      meshmax[j] = 0;
    }

  for(i = 0; i < NumPart; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  slab = to_slab_fac * P[i].Pos[j];
	  if(slab >= PMGRID)
	    slab = PMGRID - 1;

	  if(slab < meshmin[j])
	    meshmin[j] = slab;

	  if(slab > meshmax[j])
	    meshmax[j] = slab;
	}
    }

  MPI_Allgather(meshmin, 3, MPI_INT, meshmin_list, 3, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgather(meshmax, 3, MPI_INT, meshmax_list, 3, MPI_INT, MPI_COMM_WORLD);

  dimx = meshmax[0] - meshmin[0] + 2;
  dimy = meshmax[1] - meshmin[1] + 2;
  dimz = meshmax[2] - meshmin[2] + 2;

  for(i = 0; i < dimx * dimy * dimz; i++)
    workspace[i] = 0;

  for(i = 0; i < N_gas; i++)	/* just bin the gas particles */
    {
      slab_x = to_slab_fac * P[i].Pos[0];
      if(slab_x >= PMGRID)
	slab_x = PMGRID - 1;
      dx = to_slab_fac * P[i].Pos[0] - slab_x;
      slab_x -= meshmin[0];
      slab_xx = slab_x + 1;

      slab_y = to_slab_fac * P[i].Pos[1];
      if(slab_y >= PMGRID)
	slab_y = PMGRID - 1;
      dy = to_slab_fac * P[i].Pos[1] - slab_y;
      slab_y -= meshmin[1];
      slab_yy = slab_y + 1;

      slab_z = to_slab_fac * P[i].Pos[2];
      if(slab_z >= PMGRID)
	slab_z = PMGRID - 1;
      dz = to_slab_fac * P[i].Pos[2] - slab_z;
      slab_z -= meshmin[2];
      slab_zz = slab_z + 1;

      workspace[(slab_x * dimy + slab_y) * dimz + slab_z] += P[i].Mass * (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
      workspace[(slab_x * dimy + slab_yy) * dimz + slab_z] += P[i].Mass * (1.0 - dx) * dy * (1.0 - dz);
      workspace[(slab_x * dimy + slab_y) * dimz + slab_zz] += P[i].Mass * (1.0 - dx) * (1.0 - dy) * dz;
      workspace[(slab_x * dimy + slab_yy) * dimz + slab_zz] += P[i].Mass * (1.0 - dx) * dy * dz;

      workspace[(slab_xx * dimy + slab_y) * dimz + slab_z] += P[i].Mass * (dx) * (1.0 - dy) * (1.0 - dz);
      workspace[(slab_xx * dimy + slab_yy) * dimz + slab_z] += P[i].Mass * (dx) * dy * (1.0 - dz);
      workspace[(slab_xx * dimy + slab_y) * dimz + slab_zz] += P[i].Mass * (dx) * (1.0 - dy) * dz;
      workspace[(slab_xx * dimy + slab_yy) * dimz + slab_zz] += P[i].Mass * (dx) * dy * dz;
    }


  for(i = 0; i < fftsize; i++)	/* clear local density field */
    rhogrid[i] = 0;

  for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;
      if(recvTask < NTask)
	{
	  /* check how much we have to send */
	  sendmin = 2 * PMGRID;
	  sendmax = -1;
	  for(slab_x = meshmin[0]; slab_x < meshmax[0] + 2; slab_x++)
	    if(slab_to_task[slab_x % PMGRID] == recvTask)
	      {
		if(slab_x < sendmin)
		  sendmin = slab_x;
		if(slab_x > sendmax)
		  sendmax = slab_x;
	      }
	  if(sendmax == -1)
	    sendmin = 0;

	  /* check how much we have to receive */
	  recvmin = 2 * PMGRID;
	  recvmax = -1;
	  for(slab_x = meshmin_list[3 * recvTask]; slab_x < meshmax_list[3 * recvTask] + 2; slab_x++)
	    if(slab_to_task[slab_x % PMGRID] == sendTask)
	      {
		if(slab_x < recvmin)
		  recvmin = slab_x;
		if(slab_x > recvmax)
		  recvmax = slab_x;
	      }
	  if(recvmax == -1)
	    recvmin = 0;


	  if((recvmax - recvmin) >= 0 || (sendmax - sendmin) >= 0)	/* ok, we have a contribution to the slab */
	    {
	      recv_dimx = meshmax_list[3 * recvTask + 0] - meshmin_list[3 * recvTask + 0] + 2;
	      recv_dimy = meshmax_list[3 * recvTask + 1] - meshmin_list[3 * recvTask + 1] + 2;
	      recv_dimz = meshmax_list[3 * recvTask + 2] - meshmin_list[3 * recvTask + 2] + 2;

	      if(level > 0)
		{
		  MPI_Sendrecv(workspace + (sendmin - meshmin[0]) * dimy * dimz,
			       (sendmax - sendmin + 1) * dimy * dimz * sizeof(fftw_real), MPI_BYTE, recvTask,
			       TAG_PERIODIC_A, forcegrid,
			       (recvmax - recvmin + 1) * recv_dimy * recv_dimz * sizeof(fftw_real), MPI_BYTE,
			       recvTask, TAG_PERIODIC_A, MPI_COMM_WORLD, &status);
		}
	      else
		{
		  memcpy(forcegrid, workspace + (sendmin - meshmin[0]) * dimy * dimz,
			 (sendmax - sendmin + 1) * dimy * dimz * sizeof(fftw_real));
		}

	      for(slab_x = recvmin; slab_x <= recvmax; slab_x++)
		{
		  slab_xx = (slab_x % PMGRID) - first_slab_of_task[ThisTask];

		  if(slab_xx >= 0 && slab_xx < slabs_per_task[ThisTask])
		    {
		      for(slab_y = meshmin_list[3 * recvTask + 1];
			  slab_y <= meshmax_list[3 * recvTask + 1] + 1; slab_y++)
			{
			  slab_yy = slab_y;
			  if(slab_yy >= PMGRID)
			    slab_yy -= PMGRID;

			  for(slab_z = meshmin_list[3 * recvTask + 2];
			      slab_z <= meshmax_list[3 * recvTask + 2] + 1; slab_z++)
			    {
			      slab_zz = slab_z;
			      if(slab_zz >= PMGRID)
				slab_zz -= PMGRID;

			      rhogrid[PMGRID * PMGRID2 * slab_xx + PMGRID2 * slab_yy + slab_zz] +=
				forcegrid[((slab_x - recvmin) * recv_dimy +
					   (slab_y - meshmin_list[3 * recvTask + 1])) * recv_dimz +
					  (slab_z - meshmin_list[3 * recvTask + 2])];
			    }
			}
		    }
		}
	    }
	}
    }

  /* Do the FFT of the density field */

  rfftwnd_mpi(fft_forward_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);

  /* smooth the density field and deconvolve with CIC  */

  for(y = slabstart_y; y < slabstart_y + nslab_y; y++)
    for(x = 0; x < PMGRID; x++)
      for(z = 0; z < PMGRID / 2 + 1; z++)
	{
	  if(x > PMGRID / 2)
	    kx = x - PMGRID;
	  else
	    kx = x;
	  if(y > PMGRID / 2)
	    ky = y - PMGRID;
	  else
	    ky = y;
	  if(z > PMGRID / 2)
	    kz = z - PMGRID;
	  else
	    kz = z;

	  k2 = kx * kx + ky * ky + kz * kz;

	  if(k2 > 0)
	    {
	      smth = exp(-k2 * asmth2);

	      /* do deconvolution */

	      fx = fy = fz = 1;
	      if(kx != 0)
		{
		  fx = (M_PI * kx) / PMGRID;
		  fx = sin(fx) / fx;
		}
	      if(ky != 0)
		{
		  fy = (M_PI * ky) / PMGRID;
		  fy = sin(fy) / fy;
		}
	      if(kz != 0)
		{
		  fz = (M_PI * kz) / PMGRID;
		  fz = sin(fz) / fz;
		}
	      ff = 1 / (fx * fy * fz);
	      smth *= ff * ff * ff * ff;

	      /* end deconvolution */

	      ip = PMGRID * (PMGRID / 2 + 1) * (y - slabstart_y) + (PMGRID / 2 + 1) * x + z;
	      fft_of_rhogrid[ip].re *= smth;
	      fft_of_rhogrid[ip].im *= smth;
	    }
	}

  /* Do the inverse FFT to get back density field */

  rfftwnd_mpi(fft_inverse_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);

  /* Now rhogrid holds the smoothed density field */

  /* Now compute the pressure field */

  for(i = 0; i < fftsize; i++)
    {
      rhogrid[i] *= facdens;	/* to get comoving baryonic density */

      rhogrid[i] = facenthalpy * pow(rhogrid[i], All.HPM_alpha - 1);
    }

  /* construct the specific enthalpy for the local patch */

  dimx = meshmax[0] - meshmin[0] + 6;
  dimy = meshmax[1] - meshmin[1] + 6;
  dimz = meshmax[2] - meshmin[2] + 6;

  for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
	{

	  /* check how much we have to send */
	  sendmin = 2 * PMGRID;
	  sendmax = -PMGRID;
	  for(slab_x = meshmin_list[3 * recvTask] - 2; slab_x < meshmax_list[3 * recvTask] + 4; slab_x++)
	    if(slab_to_task[(slab_x + PMGRID) % PMGRID] == sendTask)
	      {
		if(slab_x < sendmin)
		  sendmin = slab_x;
		if(slab_x > sendmax)
		  sendmax = slab_x;
	      }
	  if(sendmax == -PMGRID)
	    sendmin = sendmax + 1;


	  /* check how much we have to receive */
	  recvmin = 2 * PMGRID;
	  recvmax = -PMGRID;
	  for(slab_x = meshmin[0] - 2; slab_x < meshmax[0] + 4; slab_x++)
	    if(slab_to_task[(slab_x + PMGRID) % PMGRID] == recvTask)
	      {
		if(slab_x < recvmin)
		  recvmin = slab_x;
		if(slab_x > recvmax)
		  recvmax = slab_x;
	      }
	  if(recvmax == -PMGRID)
	    recvmin = recvmax + 1;

	  if((recvmax - recvmin) >= 0 || (sendmax - sendmin) >= 0)	/* ok, we have a contribution to the slab */
	    {
	      recv_dimx = meshmax_list[3 * recvTask + 0] - meshmin_list[3 * recvTask + 0] + 6;
	      recv_dimy = meshmax_list[3 * recvTask + 1] - meshmin_list[3 * recvTask + 1] + 6;
	      recv_dimz = meshmax_list[3 * recvTask + 2] - meshmin_list[3 * recvTask + 2] + 6;

	      ncont = 1;
	      cont_sendmin[0] = sendmin;
	      cont_sendmax[0] = sendmax;
	      cont_sendmin[1] = sendmax + 1;
	      cont_sendmax[1] = sendmax;

	      cont_recvmin[0] = recvmin;
	      cont_recvmax[0] = recvmax;
	      cont_recvmin[1] = recvmax + 1;
	      cont_recvmax[1] = recvmax;

	      for(slab_x = sendmin; slab_x <= sendmax; slab_x++)
		{
		  if(slab_to_task[(slab_x + PMGRID) % PMGRID] != ThisTask)
		    {
		      /* non-contiguous */
		      cont_sendmax[0] = slab_x - 1;
		      while(slab_to_task[(slab_x + PMGRID) % PMGRID] != ThisTask)
			slab_x++;
		      cont_sendmin[1] = slab_x;
		      ncont++;
		    }
		}

	      for(slab_x = recvmin; slab_x <= recvmax; slab_x++)
		{
		  if(slab_to_task[(slab_x + PMGRID) % PMGRID] != recvTask)
		    {
		      /* non-contiguous */
		      cont_recvmax[0] = slab_x - 1;
		      while(slab_to_task[(slab_x + PMGRID) % PMGRID] != recvTask)
			slab_x++;
		      cont_recvmin[1] = slab_x;
		      if(ncont == 1)
			ncont++;
		    }
		}


	      for(rep = 0; rep < ncont; rep++)
		{
		  sendmin = cont_sendmin[rep];
		  sendmax = cont_sendmax[rep];
		  recvmin = cont_recvmin[rep];
		  recvmax = cont_recvmax[rep];

		  /* prepare what we want to send */
		  if(sendmax - sendmin >= 0)
		    {
		      for(slab_x = sendmin; slab_x <= sendmax; slab_x++)
			{
			  slab_xx = ((slab_x + PMGRID) % PMGRID) - first_slab_of_task[ThisTask];

			  for(slab_y = meshmin_list[3 * recvTask + 1] - 2;
			      slab_y < meshmax_list[3 * recvTask + 1] + 4; slab_y++)
			    {
			      slab_yy = (slab_y + PMGRID) % PMGRID;

			      for(slab_z = meshmin_list[3 * recvTask + 2] - 2;
				  slab_z <= meshmax_list[3 * recvTask + 2] + 4; slab_z++)
				{
				  slab_zz = (slab_z + PMGRID) % PMGRID;

				  forcegrid[((slab_x - sendmin) * recv_dimy +
					     (slab_y - (meshmin_list[3 * recvTask + 1] - 2))) * recv_dimz +
					    slab_z - (meshmin_list[3 * recvTask + 2] - 2)] =
				    rhogrid[PMGRID * PMGRID2 * slab_xx + PMGRID2 * slab_yy + slab_zz];
				}
			    }
			}
		    }

		  if(level > 0)
		    {
		      MPI_Sendrecv(forcegrid,
				   (sendmax - sendmin + 1) * recv_dimy * recv_dimz * sizeof(fftw_real),
				   MPI_BYTE, recvTask, TAG_PERIODIC_B,
				   workspace + (recvmin - (meshmin[0] - 2)) * dimy * dimz,
				   (recvmax - recvmin + 1) * dimy * dimz * sizeof(fftw_real), MPI_BYTE,
				   recvTask, TAG_PERIODIC_B, MPI_COMM_WORLD, &status);
		    }
		  else
		    {
		      memcpy(workspace + (recvmin - (meshmin[0] - 2)) * dimy * dimz,
			     forcegrid, (recvmax - recvmin + 1) * dimy * dimz * sizeof(fftw_real));
		    }
		}
	    }
	}
    }


  dimx = meshmax[0] - meshmin[0] + 2;
  dimy = meshmax[1] - meshmin[1] + 2;
  dimz = meshmax[2] - meshmin[2] + 2;

  recv_dimx = meshmax[0] - meshmin[0] + 6;
  recv_dimy = meshmax[1] - meshmin[1] + 6;
  recv_dimz = meshmax[2] - meshmin[2] + 6;


  for(dim = 0; dim < 3; dim++)	/* Calculate each component of the force. */
    {
      /* get the force component by finite differencing the potential */
      /* note: "workspace" now contains the potential for the local patch, plus a suffiently large buffer region */

      for(x = 0; x < meshmax[0] - meshmin[0] + 2; x++)
	for(y = 0; y < meshmax[1] - meshmin[1] + 2; y++)
	  for(z = 0; z < meshmax[2] - meshmin[2] + 2; z++)
	    {
	      xrr = xll = xr = xl = x;
	      yrr = yll = yr = yl = y;
	      zrr = zll = zr = zl = z;

	      switch (dim)
		{
		case 0:
		  xr = x + 1;
		  xrr = x + 2;
		  xl = x - 1;
		  xll = x - 2;
		  break;
		case 1:
		  yr = y + 1;
		  yl = y - 1;
		  yrr = y + 2;
		  yll = y - 2;
		  break;
		case 2:
		  zr = z + 1;
		  zl = z - 1;
		  zrr = z + 2;
		  zll = z - 2;
		  break;
		}

	      forcegrid[(x * dimy + y) * dimz + z]
		=
		facdiff * ((4.0 / 3) *
			   (workspace[((xl + 2) * recv_dimy + (yl + 2)) * recv_dimz + (zl + 2)]
			    - workspace[((xr + 2) * recv_dimy + (yr + 2)) * recv_dimz + (zr + 2)]) -
			   (1.0 / 6) *
			   (workspace[((xll + 2) * recv_dimy + (yll + 2)) * recv_dimz + (zll + 2)] -
			    workspace[((xrr + 2) * recv_dimy + (yrr + 2)) * recv_dimz + (zrr + 2)]));
	    }

      /* read out the forces, but just for the gas particles */

      for(i = 0; i < N_gas; i++)
	{
	  slab_x = to_slab_fac * P[i].Pos[0];
	  if(slab_x >= PMGRID)
	    slab_x = PMGRID - 1;
	  dx = to_slab_fac * P[i].Pos[0] - slab_x;
	  slab_x -= meshmin[0];
	  slab_xx = slab_x + 1;

	  slab_y = to_slab_fac * P[i].Pos[1];
	  if(slab_y >= PMGRID)
	    slab_y = PMGRID - 1;
	  dy = to_slab_fac * P[i].Pos[1] - slab_y;
	  slab_y -= meshmin[1];
	  slab_yy = slab_y + 1;

	  slab_z = to_slab_fac * P[i].Pos[2];
	  if(slab_z >= PMGRID)
	    slab_z = PMGRID - 1;
	  dz = to_slab_fac * P[i].Pos[2] - slab_z;
	  slab_z -= meshmin[2];
	  slab_zz = slab_z + 1;

	  acc_dim =
	    forcegrid[(slab_x * dimy + slab_y) * dimz + slab_z] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
	  acc_dim += forcegrid[(slab_x * dimy + slab_yy) * dimz + slab_z] * (1.0 - dx) * dy * (1.0 - dz);
	  acc_dim += forcegrid[(slab_x * dimy + slab_y) * dimz + slab_zz] * (1.0 - dx) * (1.0 - dy) * dz;
	  acc_dim += forcegrid[(slab_x * dimy + slab_yy) * dimz + slab_zz] * (1.0 - dx) * dy * dz;

	  acc_dim += forcegrid[(slab_xx * dimy + slab_y) * dimz + slab_z] * (dx) * (1.0 - dy) * (1.0 - dz);
	  acc_dim += forcegrid[(slab_xx * dimy + slab_yy) * dimz + slab_z] * (dx) * dy * (1.0 - dz);
	  acc_dim += forcegrid[(slab_xx * dimy + slab_y) * dimz + slab_zz] * (dx) * (1.0 - dy) * dz;
	  acc_dim += forcegrid[(slab_xx * dimy + slab_yy) * dimz + slab_zz] * (dx) * dy * dz;

	  SphP[i].HydroAccel[dim] = acc_dim;
	}
    }
}
