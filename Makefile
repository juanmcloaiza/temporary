#######################################################################
#  Look at end of file for a brief guide to the compile-time options. #
#######################################################################

####To accrete in the "turbo" way
OPT += -DMYSWITCH #External fixed potential + BH accretion at the centre

#--------------------------------------- Basic operation mode of code
#OPT   +=  -DPERIODIC
#OPT   +=  -DCOOLING
#OPT   +=  -DSFR
#OPT   +=  -DUNEQUALSOFTENINGS


#--------------------------------------- TreePM Options
#OPT   +=  -DPMGRID=64
#OPT   +=  -DASMTH=1.25
#OPT   +=  -DRCUT=4.5
#OPT   +=  -DPLACEHIGHRESREGION=3
#OPT   +=  -DENLARGEREGION=1.2
#OPT   +=  -DONLY_PM
#OPT   +=  -DHPM
#OPT   +=  -DHPM_SMTH=1.5

#--------------------------------------- Single/Double Precision
#OPT   +=  -DDOUBLEPRECISION
OPT   +=  -DDOUBLEPRECISION_FFTW
#OPT   +=  -DFLTROUNDOFFREDUCTION    # enables round off reduction in particle sums
				    # if DOUBLEPRECISION is set, these sums are done in 'long double'
                                    # if single precision is used, they are done in 'double'
                                    # This should in principle allow to make computations
                                    # *exactly* invariant to different numbers of CPUs.

#OPT   +=  -DSOFTDOUBLEDOUBLE       # when this is set, a software implementation of
                                    # 128bit double-double addition is used, implemented as a c++ class.
                                    # Hence, this option requires compilation with a c++ compiler

#--------------------------------------- SFR/feedback model

#OPT   +=  -DSOFTEREQS
#OPT   +=  -DMOREPARAMS
#OPT   +=  -DMETALS
#OPT   +=  -DSTELLARAGE
#OPT   +=  -DWINDS
#OPT   +=  -DQUICK_LYALPHA
#OPT   +=  -DISOTROPICWINDS
#OPT   +=  -DMHM

#--------------------------------------- AGN stuff

#OPT   +=  -DFOF                     # enable FoF output

#OPT   +=  -DBLACK_HOLES             # enable Black-Holes

#OPT   +=  -DBONDI                   # Bondi-Hoyle style accretion model
#OPT   +=  -DENFORCE_EDDINGTON_LIMIT # put a hard limit on the maximum accretion rate
#OPT   +=  -DBH_THERMALFEEDBACK      # couple a fraction of the BH luminosity into surrounding gas
#OPT   +=  -DBH_DRAG                 # Drag on black-holes due to accretion
#OPT   +=  -DSWALLOWGAS              # Enables stochastic accretion of gas particles consistent with growth rate of hole
#OPT   +=  -DEVALPOTENTIAL           # computes gravitational potential
#OPT   +=  -DREPOSITION_ON_POTMIN    # repositions hole on potential minimum

#OPT    +=  -DBUBBLES                # switch on generation of hot bubbles in an
                                     #isolated halo or the the biggest halo in
                                     #the run

#OPT	+=  -DMULTI_BUBBLES 	     #switch on generation of hot bubbles in all
                                     #haloes above certain mass threshold - note
                                     #works only with FOF and BUBBLES have to be
                                     #switched off

#OPT	+=  -DEBUB_PROPTO_BHAR       #Energy content of the bubbles with cosmic
				     #time evolves as an integrated BHAR(z) over a
				     #Salpeter time (Di Matteo 2003 eq. [11])


#--------------------------------------- Things that are always recommended
OPT   +=  -DPEANOHILBERT
OPT   +=  -DWALLCLOCK
#OPT   +=  -DCPUSPEEDADJUSTMENT

#--------------------------------------- Viscous gas treatment 
#OPT   +=  -DNAVIERSTOKES                  #Braginskii-Spitzer parametrization of the shear viscosity: mu = f x T^{5/2}
#OPT   +=  -DNAVIERSTOKES_CONSTANT         #Shear viscosity set constant for all gas particles
#OPT   +=  -DNAVIERSTOKES_BULK             #Bulk viscosity set constant for all gas particles. To run with bulk visocity only one has to set shear viscosity to zero in the parameterfile.
#OPT   +=  -DVISCOSITY_SATURATION          # Both shear and bulk viscosities are saturated, so that unphysical accelerations and entropy increases are avoided. Relevant for the cosmological simulations.
#OPT   +=  -DNS_TIMESTEP                   #Enables timestep criterion based on entropy increase due to internal friction forces
#OPT   +=  -DOUTPUTSTRESS                  #Outputs diagonal and offdiagonal components of viscous shear stress tensor
#OPT   +=  -DOUTPUTBULKSTRESS              #Outputs viscous bulk stress tensor
#OPT   +=  -DOUTPUTSHEARCOEFF              #Outputs variable shear viscosity coefficient in internal code units

#--------------------------------------- Things for special behaviour

#OPT   +=  -DPEDANTIC_MEMORY_HANDLER      # this enables an internal memoru handler (most recently allocated block needs to be freed first)
#OPT   +=  -DPEDANTIC_MEMORY_CEILING=400  # this should be set to define the memory ceiling in MByte
OPT   +=  -DNOGRAVITY
OPT   +=  -DISOTHERM_EQS #JM imported this from the old gadget, use with CAUTION !!!!!
#OPT   +=  -DNOACCEL
#OPT   +=  -DNOISMPRESSURE
#OPT   +=  -DNOFIXEDMASSINKERNEL
#OPT   +=  -DNOGRADHSML
#OPT   +=  -DNOVISCOSITYLIMITER
#OPT   +=  -DNOTREERND
#OPT   +=  -DNOSTOP_WHEN_BELOW_MINTIMESTEP
#OPT   +=  -DNOPMSTEPADJUSTMENT
OPT   +=  -DNOTYPEPREFIX_FFTW
#OPT   +=  -DNO_TREEDATA_IN_RESTART
#OPT   +=  -DNOWINDTIMESTEPPING         # Disable wind reducing timestep
                                        # Not recommended !
#OPT   +=  -DISOTHERM=200  # adds potential of an isothermal sphere
#OPT   +=  -DCOMPUTE_POTENTIAL_ENERGY
#OPT   +=  -DALLOWEXTRAPARAMS
#OPT   +=  -DLONGIDS
#OPT   +=  -DINHOMOG_GASDISTR_HINT   # if the gas distribution is spatially very different
                                     # from collisionless particles, this helps to avoid
                                     # problems in the domain decomposition
#OPT   +=  -DLONG_X=140
#OPT   +=  -DLONG_Y=1
#OPT   +=  -DLONG_Z=1
#OPT   +=  -DTWODIMS
#OPT   +=  -DSPH_BND_PARTICLES
#OPT   +=  -DNEW_RATES               # switches in updated cooling rates from Naoki
#OPT   +=  -DREAD_HSML               # reads hsml from IC file

#OPT    +=  -DVOLUME_CORRECTION
#OPT    +=  -DVOLUME_QUADRIC
#OPT    +=  -DVOLUME_QUINTIC
#OPT    +=  -DVOLUME_BIWEIGHT
#OPT    +=  -DSTOP_AFTER_VOL_CORRECTION

#OPT   +=  -DADAPTIVE_GRAVSOFT_FORGAS  # allows variable softening length for gas particles
                                      # this option require UNEQUALSOFTENINGLENGTH

#OPT   +=  -DADAPTIVE_GRAVSOFT_FORGAS_HSML=0.1  # in case variable gravitational softening for gas
                                            # is activated, this option couples the softening
                                            # length to the SPH smoothing length


#--------------------------------------- Time integration options
OPT   +=  -DSYNCHRONIZATION
#OPT   +=  -DFLEXSTEPS
#OPT   +=  -DPSEUDOSYMMETRIC

OPT   +=  -DAUTO_SWAP_ENDIAN_READIC

#--------------------------------------- Output options
#OPT   +=  -DOUTPUTPOTENTIAL
#OPT   +=  -DOUTPUTACCELERATION
#OPT   +=  -DOUTPUTCHANGEOFENTROPY
#OPT   +=  -DOUTPUTTIMESTEP
#OPT   +=  -DOUTPUTCOOLRATE     # outputs cooling rate, and conduction rate if enabled
#OPT   +=  -DHAVE_HDF5  # needed when HDF5 I/O support is desired
#OPT   +=  -DOUTPUTBSMOOTH
#OPT   +=  -DOUTPUTDENSNORM
#OPT   +=  -DXXLINFO           #  Enables additional output for viscosityand bfield
#OPT   +=  -DOUTPUTLINEOFSIGHT    # enables on-the-fly output of Ly-alpha absorption spectra
#OPT   +=  -DOUTPUTLINEOFSIGHT_SPECTRUM
#OPT   +=  -DOUTPUTLINEOFSIGHT_PARTICLES


#--------------------------------------- Testing and Debugging options
#OPT   +=  -DFORCETEST=0.1
#OPT   +=  -DDEBUG    # enables core-dumps and FPU exceptions
#OPT   +=  -DPARTICLE_DEBUG  # auxiliary communication of IDs


#--------------------------------------- Static NFW Potential
#OPT   +=  -DSTATICNFW
#OPT   +=  -DNFW_C=12
#OPT   +=  -DNFW_M200=100.0
#OPT   +=  -DNFW_Eps=0.01
#OPT   +=  -DNFW_DARKFRACTION=0.87


#--------------------------------------- Static Hernquist Potential
#OPT   +=  -DSTATICHQ
#OPT   +=  -DHQ_M200=1.0
#OPT   +=  -DHQ_C=10
#OPT   +=  -DHQ_DARKFRACTION=0.9


#--------------------------------------- Thermal conduction
#OPT   +=  -DCONDUCTION
#OPT   +=  -DCONDUCTION_CONSTANT
#OPT   +=  -DCONDUCTION_SATURATION


#--------------------------------------- Dark energy
#OPT   +=  -DDARKENERGY # Enables Dark Energy
#OPT   +=  -DTIMEDEPDE  # read w(z) from a file
#OPT   +=  -DRESCALEVINI # rescale v_ini in read_ic / read_ic_cluster
#OPT   +=  -DTIMEDEPGRAV # resacles H and G according to DE model

#--------------------------------------- SPH viscosity options
#OPT   +=  -DCONVENTIONAL_VISCOSITY     # enables the old viscosity
OPT   +=  -DTIME_DEP_ART_VISC          # Enables time dependend viscosity
#OPT   +=  -DNO_SHEAR_VISCOSITY_LIMITER # Turns of the shear viscosity supression
#OPT   +=  -DHIGH_ART_VISC_START=0        # Start with high rather than low viscosity
#OPT   +=  -DALTVISCOSITY               # enables alternative viscosity based on div(v)

#--------------------------------------- Magnetic Field options
#OPT	+=  -DMAGNETIC                   # Turns on B Field (including induction eqn) 
#OPT	+=  -DMAGFORCE                   # Turns on B force
#OPT	+=  -DMAGNETIC_SIGNALVEL         # Extend definition of signal velocity 
                                        # by the magneto sonic waves

#....................................... Induction equation stuff
#OPT	+=  -DDIVBINDUCTION              # Turns on divB term in induction eq.
#OPT	+=  -DCORRECTDB                  # Turns on dW/dh correction in induction eq.

#....................................... Force equation stuff
#OPT	+=  -DMU0_UNITY                  # Sets \mu_0 to unity
#OPT	+=  -DARTBPRES                   # Anti clumping term ala Morris
#OPT	+=  -DDIVBFORCE                  # Subtract div(B) force from M-tensor
#OPT	+=  -DCORRECTBFRC                # Turns on dW/dh correction for B force 

#....................................... Smoothing Options
#OPT	+=  -DBSMOOTH                    # Turns on B field smoothing

#....................................... Artificial magnetic dissipation options
#OPT	+=  -DMAGNETIC_DISSIPATION       # Turns on artificial magnetic dissipation 
#OPT	+=  -DMAGDISSIPATION_PERPEN      # Uses only the perpendicular magnetic 
                                        # field component for dissipation 
#OPT	+=  -DTIME_DEP_MAGN_DISP         # Enables time dependent coefficients
#OPT	+=  -DHIGH_MAGN_DISP_START       # Starts from high coefficient
#OPT	+=  -DROT_IN_MAG_DIS             # Adds the RotB term in dissipation 

#....................................... DivB cleaning options
#OPT	+=  -DDIVBCLEANING_DEDNER        # Turns on hyp/par cleaning (Dedner 2002)
#OPT	+=  -DSMOOTH_PHI                 # Turns on smoothing of Phi

#....................................... magnetic diffusion options
#OPT	+=  -DMAGNETIC_DIFFUSION         # Turns on magnetic diffusion
#OPT	+=  -DSMOOTH_ROTB                # Turns on smoothing of rot(B)

#....................................... Bebugging stuff
#OPT	+=  -DTRACEDIVB                  # Writes div(B) into snapshot
#OPT	+=  -DDBOUTPUT                   # Writes dB/dt into snapshot
#OPT	+=  -DOUTPUT_ROTB                # Writes rotB into snapshot
#OPT	+=  -DOUTPUT_SROTB               # Writes smoothed rotB into snapshot

#....................................... Initil condition stuff
#OPT	+=  -DBINISET                    # Allows to set Bx,By,Bz in parameter file
#OPT	+=  -DBFROMROTA                  # Allows to five vector potential instead of
                                        # B within the IC file
#OPT	+=  -DIGNORE_PERIODIC_IN_ROTA    # Don't use the periodic mapping when 
                                        # calculating rot(A)
                                        # Note A might not be periodic even if B is. 
#OPT	+=  -DBRIOWU                     # Extrapolate A outside simulation in case 
                                        # of Brio & Wu test problem


#--------------------------------------- Glass making
#OPT   +=  -DMAKEGLASS


#--------------------------------------- Cecilia 's model
# All SFR_* options belong to Cecilia
#OPT   +=  -DSFR_METALS
#OPT   +=  -DSFR_FEEDBACK
#OPT   +=  -DSFR_SNI
#OPT   +=  -DSFR_SNII
#OPT   +=  -DSFR_ENRICH
#OPT   +=  -DSFR_DECOUPLING  # Marri-like decoupling
#OPT   +=  -DSFR_PROMOTION
#OPT   +=  -DSFR_DIFFUSION

#-------------------------------------- nonequilibrium proimodal chemisitry
#OPT    += -DNONEQUILIBRIUM
#OPT   +=  -DCHEMISTRY
#OPT   +=  -DCMB
#OPT   +=  -DRADIATION

#--------------------------------------- Cosmic Rays (Martin)
#OPT   +=  -DCOSMIC_RAYS       # Cosmic Rays Master Switch
#OPT   +=  -DCR_IC             # IC files contain CR information
#OPT   +=  -DCR_IC_PHYSICAL
#OPT   +=  -DCR_DISSIPATION    # Catastrophic losses
#OPT   +=  -DCR_THERMALIZATION # Coulomb cooling
#OPT   +=  -DCR_SHOCK=2        # Shock energy is directed into CR
			       # 2 = Mach-Number dependent shocks, Mach-number derived for thermal gas/CR composite
			       # 3 = Mach-Number dependent shocks, Mach-number derived for thermal gas
#OPT   +=  -DCR_DIFFUSION       # Cosmic Ray diffusion
#OPT   +=  -DCR_DIFFUSION_GREEN # alternative diffusion model

#OPT   +=  -DUPDATE_PARANOIA=1 # 1 = Update on every predict
			       # 2 = Update on every energy injection and
			       #     on every predict
#OPT   +=  -DCR_OUTPUT_INJECTION # Output energy injection rate in snapshots
#OPT   +=  -DCR_NO_CHANGE      # Compute changes to CR, but do not execute them
                               # useful for estimating the size of effects
#OPT   +=  -DCOSMIC_RAY_TEST   # starts a test routine instead of the simulation
#OPT   +=  -DCR_NOPRESSURE     # computes CRs as usual, but ignores the pressure in the dynamics

#--------------------------------------- Mach number finder (Christoph)
#OPT  += -DMACHNUM                       # Mach number Master Switch
#OPT  += -DMACHSTATISTIC                 # Dissipated thermal energy at shocks
#OPT  += -DCR_OUTPUT_JUMP_CONDITIONS     # CR: density and thermal energy jump at shocks

#--------------------------------------- Metals from Luca Tornatore

#OPT  += -DLT_METAL_COOLING           # METAL COOLING option
#OPT  += -DLT_MEAN_Z_ESTIMATE
#OPT  += -DLT_SMOOTH_Z

#....................................# STELLAR EVOLUTION options
#OPT  += -DLT_STELLAREVOLUTION        # enable stellar evolution
#OPT  += -DLT_NMet=9                  # number of species
#OPT  += -DLT_NIMFs=1

#OPT  += -DLT_SNII                    # enable snII
#OPT  += -DLT_SNIa                    # enable SnIa
#OPT  += -DLT_Nebulae                 # enable low mass stars

#OPT  += -DLT_ACCOUNT_NONPROC_METALS # don't use so far

#OPT  += -DLT_PM_LIFETIMES            # use Padovani&matteucci 1993 lifetimes
#OPT  += -DLT_MM_LIFETIMES           # use Maeder&Menet 1989 lifetimes

#OPT  +=  LT_STOCHASTIC_IRA

#OPT  += -DLT_LOCAL_IRA
#OPT  += -DLT_AVOID_SELFENRICHMENT
#OPT  += -DLT_AVOID_ENRICH_SFGAS


#OPT  += -DLT_WIND_VELOCITY=500.0   # set the winds' velocity in km/s

#OPT  += -DLT_EJECTA_IN_HOTPHASE
#OPT  += -DLT_SNegy_IN_HOTPHASE
#OPT  += -DLT_HOT_EJECTA
#OPT  += -DLT_HOT_WINDS

#OPT  += -DLT_CharT                  # extra conditions on SF:
                                     #  cooling time < sound cross time
                                     #  ffall time < sound cross time

#.....................................# Other physics options
#OPT  += -DLT_WINDS_EXT               # EXTRA WINDS options
#OPT  += -DLT_WINDS_EXT_NOSF_BRANCH

#OPT  += -DLT_MOD_EFFM                 #

#OPT  += -DLT_USE_TOP_HAT_WEIGHT
#OPT  += -DLT_USE_SOLIDANGLE_WEIGHT

#.....................................# logging/debug options					
#OPT  += -DLT_SEv_INFO                 # infos about enrichment
#OPT  += -DLT_SEv_INFO_DETAILS
#OPT  += -DLT_EXTEGY_INFO              # infos about energy
#OPT  += -DLT_CharT_INFO               # infos about cooling time,
                                      #  sound crossing time,
                                      #  free fall time
#OPT  += -DLT_SEvDbg
#OPT  += -DLT_TRACK_CONTRIBUTES
#.....................................#
#OPT  += -DLT_COST_SE                 # account for cost of Sn in
                                      #  domain dec.

# ====================================
# Blue Gene/L implementation
# ====================================

#OPT += -DBG_SFR
#OPT += -DBG_WINDS
#OPT += -DBG_COOLING
#OPT += -DBG_STELLAR_EVOLUTION
#OPT += -DBG_VERBOSE

#========================================================================


CC       = mpicc        # sets the C-compiler (default)
OPTIMIZE = -Wall -O3 #-g   # optimization and warning flags (default)

#--------------------------------------- Select target computer
#SYSTYPE="ULYSSES"
#SYSTYPE="FERMI"
#SYSTYPE="Dirac"
#SYSTYPE="HG1"
SYSTYPE="OFFICE"

#Added by JM
#--------------------------------------- Adjust settings for target computer
ifeq ($(SYSTYPE),"OFFICE")
CC       =  mpicc   
OPTIMIZE =  -Wall -O3 #-g #-no-multibyte-chars
#GSL_INCL = -I$(GSL_INC)
#GSL_LIBS = -L$(GSL_LIB)
#FFTW_INCL= -I$(FFTW_INC)
#FFTW_LIBS= -L$(FFTW_LIB)
#MPICHLIB =
endif
#

ifeq ($(SYSTYPE),"ULYSSES")
CC       =  mpicc   
OPTIMIZE =  -Wall -O3 #-g #-no-multibyte-chars
GSL_INCL = -I$(GSL_INC)
GSL_LIBS = -L$(GSL_LIB)
FFTW_INCL= -I$(FFTW_INC)
FFTW_LIBS= -L$(FFTW_LIB)
#MPICHLIB =
endif
#End of Added by JM

#MPICHLIB = -lmpich


#GSL_INCL = -I/usr/include
#GSL_LIBS = -L/usr/lib
#FFTW_INCL= -I/opt/fftw/intel_8.1/2.1.5/include
#FFTW_LIBS= -L/opt/fftw/intel_8.1/2.1.5/lib #-ldrfftw_mpi
#MPICHLIB = -L/opt/lam-7.1.2b24-g77/lib -lmpi
#HDF5INCL = -I/opt/hdf5-oscar-1.6.4/include
#HDF5LIB  = -L/opt/hdf5-oscar-1.6.4/lib -lhdf5 -lz


OPTIONS = $(OPTIMIZE) $(OPT)

EXEC   = 2runP-Gadget
#JM added: accrete_particles.o
OBJS   = main.o greenf_diffusion.o  \
	 lineofsight.o kinfb_mhm.o sfr_mhm.o fof.o blackhole.o \
	 run.o predict.o begrun.o endrun.o global.o chemistry_noneq.o \
	 timestep.o init.o restart.o io.o sfr_eff.o \
	 accel.o read_ic.o cooling.o ngb.o  \
	 system.o allocate.o density.o bsmooth.o bubbles.o \
	 gravtree.o hydra.o driftfac.o darkenergy.o \
	 domain.o allvars.o potential.o read_ic_cluster_gas.o \
         forcetree.o read_ic_cluster.o peano.o gravtree_forcetest.o \
	 pm_periodic.o pm_nonperiodic.o longrange.o mymalloc.o \
	 c_metals.o c_enrichment.o c_sfr.o c_hotngbs.o cosmic_rays.o \
	 machfinder.o b_from_rot_a.o smooth_simple.o \
	 lt_getindex.o lt_read_cool_table.o lt_read_metals.o lt_read_sn_data.o \
	 lt_sfr.o lt_sn.o lt_sn_utils.o lt_imf.o accrete_particles.o

INCL   = allvars.h proto.h forcetree.h cooling.h domain.h c_metals.h cosmic_rays.h chemistry.h \
	 machfinder.h dd.h bg_cooling.h Makefile

CFLAGS = $(OPTIONS) $(GSL_INCL) $(FFTW_INCL) $(HDF5INCL)

ifeq (NOTYPEPREFIX_FFTW,$(findstring NOTYPEPREFIX_FFTW,$(OPT)))    # fftw installed with type prefix?
  FFTW_LIBR = $(FFTW_LIBS) -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
else
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(OPT)))
  FFTW_LIBR = $(FFTW_LIBS) -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
else
  FFTW_LIBR = $(FFTW_LIBS) -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
endif
endif


LIBS   = -lm $(HDF5LIB) -g $(MPICHLIB) $(GSL_LIBS) -lgsl -lgslcblas $(FFTW_LIBR)


$(EXEC): $(OBJS)
	$(CC) $(OPTIONS) $(OBJS) $(LIBS) $(RLIBS) -o $(EXEC)

$(OBJS): $(INCL)

pm_periodic.o: pm_hpm.c

clean:
	rm -f $(OBJS) $(EXEC)




###############################################################################
#
# at compile-time. From the list below, please activate/deactivate the
# options that apply to your run. If you modify any of these options,
# make sure that you recompile the whole code by typing "make clean;
# make".
#
# Main code options:
#
#     These affect the physical model that is simulated.
#
#     - PERIODIC:   Set this if you want to have periodic boundary conditions.
#     - COOLING:    This enables radiative cooling and heating. It also enables
#                   an external UV background which is read from a file.
#     - SFR:        This enables star formation using an effective multiphase
#                   models. This option requires cooling.
#     - METALS:     This model activates the tracking of enrichment in gas and
#                   stars. Note that metal-line cooling is not included yet.
#     - STELLARAGE: This stores the formation redshift of each star particle.
#     - WINDS:      This activates galactic winds. Requires star formation.
#     - ISOTROPICWINDS: This makes the wind isotropic. If not set the wind is
#                       spawned in an axial way. Requires winds to be activated.
#     - NOGRAVITY:  This switches off gravity. Makes only sense for pure
#                   SPH simulations in non-expanding space.
#
# Options for SPH:
#
#     - NOFIXEDMASSINKERNEL:  If set, the number of SPH particles in the kernel
#                             is kept constant instead of the mass.
#     - NOGRADHSML:           If actived, an equation of motion without grad(h)
#                             terms is used.
#            Note: To have the default "entropy"-formulation of SPH (Springel &
#                  Hernquist), the switches NOFIXEDMASSINKERNEL and NOGRADHSML
#                  should *not* be set.
#     - NOVISCOSITYLIMITER:   If this is set, there is no explicit upper limit
#                             on the viscosity that tries to prevent particle
#                             'reflection' in case of poor timestepping.
#
# Numerical options:
#
#     - PMGRID:     This enables the TreePM method, i.e. the long-range force
#                   is computed with a PM-algoritthm, and the short range force
#                   with the tree. The parameter has to be set to the size of the
#                   mesh that should be used, (e.g. 64, 96, 128, etc). The mesh
#                   dimensions need not necessarily be a power of two.
#                   Note: If the simulation is not in a periodic box, then a FFT
#                   method for vacuum boundaries is employed, using a mesh with
#                   dimension twice that specified by PMGRID.
#     - PLACEHIGHRESREGION: If this option is set (will only work together
#                   with PMGRID), then the long range force is computed in two
#                   stages: One Fourier-grid is used to cover the whole simulation
#                   volume, allowing the computation of the large-scale force.
#                   A second Fourier mesh is placed on the region occupied by
#                   "high-resolution" particles, allowing the computation of an
#                   intermediate scale force. Finally, the force on very small
#                   scales is supplemented by the tree. This procedure can be useful
#                   for "zoom-simulations", where the majority of particles (the
#                   high-res particles) are occupying only a small fraction of the
#                   volume. To activate this option, the parameter needs to be set
#                   to an integer that encodes the particle types that represent the
#                   high-res particles in the form of a bit mask. For example, if
#                   types 0, 1, and 4 form the high-res particles, set the parameter
#                   to PLACEHIGHRESREGION=1+2+16. The spatial region covered by the
#                   high-res grid is determined automatically from the initial
#                   conditions. Note: If a periodic box is used, the high-res zone
#                   may not intersect the box boundaries.
#     - ENLARGEREGION: The spatial region covered by the high-res zone has a fixed
#                   size during the simulation, which initially is set to the
#                   smallest region that encompasses all high-res particles. Normally, the
#                   simulation will be interrupted, if high-res particles leave this
#                   region in the course of the run. However, by setting this parameter
#                   to a value larger than one, the high-res region can be expanded.
#                   For example, setting it to 1.4 will enlarge its side-length by
#                   40% (it remains centered on the high-res particles). Hence, with
#                   such a setting, the high-res region may expand or move by a
#                   limited amount. If in addition SYNCHRONIZATION is activated, then
#                   the code will be able to continue even if high-res particles
#                   leave the initial high-res grid. In this case, the code will
#                   update the size and position of the grid that is placed onto
#                   the high-resolution region automatically. To prevent that this
#                   potentially happens every single PM step, one should nevertheless
#                   assign a value slightly larger than 1 to ENLARGEREGION.
#     - DOUBLEPRECISION: This makes the code store and compute internal
#                        particle data in double precision. Note that output
#                        files are nevertheless written by converting to single
#                        precision.
#     - NOTREERND:       If this is not set, the tree construction will succeed
#                        even when there are a few particles at identical
#                        locations. This is done by `rerouting' particles once
#                        the node-size has fallen below 1.0e-3 of the softening
#                        length. When this option is activated, this will be
#                        surpressed and the tree construction will always fail
#                        if there are particles at extremely close coordinates.
#     - NOSTOP_WHEN_BELOW_MINTIMESTEP: If this is activated, the code will not
#                        terminate when the timestep falls below the value of
#                        MinSizeTimestep specified in the parameterfile. This
#                        is useful for runs where one wants to enforce a
#                        constant timestep for all particles. This can be done
#                        by activating this option, and by setting Min- and
#                        MaxSizeTimestep to an equal value.
#     - PSEUDOSYMMETRIC: When this option is set, the code will try to "anticipate"
#                        timestep changes by extrapolating the change of the
#                        acceleration into the future. This in general improves the
#                        long-term integration behaviour of periodic orbits.
#     - SYNCHRONIZATION: When this is set, particles may only increase their
#                        timestep if the new timestep will put them into
#                        synchronization with the higher time level. This typically
#                        means that only on half of the timesteps of a particle
#                        an increase of the step may occur.
#     - NOPMSTEPADJUSTMENT: When this is set, the long-range timestep for the
#                        PM force computation is always determined by MaxSizeTimeStep.
#                        Otherwise, it is set to the minimum of MaxSizeTimeStep and
#                        the timestep obtained for the maximum long-range force with
#                        an effective softening scale equal to the PM smoothing-scale.
# - LONG_X/Y/Z:
#     These options can be used together with PERIODIC and NOGRAVITY only.
#     When set, the options define numerical factors that can be used to
#     distorts the periodic simulation cube into a parallelepiped of
#     arbitrary aspect ratio. This can be useful for idealized SPH tests.
#
# - TWODIMS:
#     This effectively switches of one dimension in SPH, i.e. the code
#     follows only 2d hydrodynamics in the xy-, yz-, or xz-plane. This
#     only works with NOGRAVITY, and if all coordinates of the third
#     axis are exactly equal. Can be useful for idealized SPH tests.
#
# - SPH_BND_PARTICLES:
#     If this is set, particles with a particle-ID equal to zero do not
#     receive any SPH acceleration. This can be useful for idealized
#     SPH tests, where these particles represent fixed "walls".
#
#
# Architecture options:
#
#     - T3E:       The code assumes that sizeof(int)=4 holds. A few machines
#                  (like Cray T3E) have sizeof(int)=8. In this case, set the
#                  T3E flag.
#     - NOTYPEPREFIX_FFTW: If this is set, the fftw-header/libraries are accessed
#                  without type prefix (adopting whatever was chosen as default at compile
#                  of fftw). Otherwise, the type prefix 'd' for double is used.
#
# Input options:
#
#     - MOREPARAMS:  Activate this to allow a set of additional parameters in
#                    the parameterfile which control the star formation and
#                    feedback sector. This option must be activated when star
#                    formation is switched on.
#
# Output options:
#
#     - OUTPUTPOTENTIAL: This will force the code to compute gravitational
#                        potentials for all particles each time a snapshot file
#                        is generated. This values are then included in the
#                        snapshot file. Note that the computation of the
#                        values of the potential costs additional time.
#     - OUTPUTACCELERATION: This will include the physical acceleration of
#                        each particle in snapshot files.
#     - OUTPUTCHANGEOFENTROPY: This will include the rate of change of entropy
#                        of gas particles in snapshot files.
#     - OUTPUTTIMESTEP:  This will include an output of the timesteps actually
#                        taken by each particle.
#
# Miscellaneous options:
#
#     - PEANOHILBERT:    This is a tuning option. When set, the code will bring
#                        the particles after each domain decomposition into
#                        Peano-Hilbert order. This improves cache utilization
#                        and performance.
#     - WALLCLOCK:       If set, a wallclock timer is used by the code to
#                        measure internal time consumption (see cpu-log file).
#                        Otherwise a timer that measures consumed processor
#                        ticks is used.
#
# Debugging/testing options:
#
#     - FORCETEST:       This can be set to check the force accuracy of the
#                        code. The option needs to be set to a number between
#                        0 and 1 (e.g. 0.01), which is taken to specify a
#                        random fraction of particles for which at each
#                        timestep forces by direct summation are computed. The
#                        normal tree-forces and the "correct" direct summation
#                        forces are collected in a file. Note that the
#                        simulation itself is unaffected by this option, but it
#                        will of course run much(!) slower
#                        if FORCETEST*NumPart*NumPart >> NumPart. Note: Particle
#                        IDs must be set to numbers >=1 for this to work.
#
###############################################################################




#     - QUICK_LYALPHA:   This only works for cosmological simulations in periodic boxes
#                        with COOLING & SFR. (WINDS, METALS should be deselected).
#                        It will simply convert all gas particles above overdensity
#                        CritPhysOverdensity and with Temperature below 10^5 K to stars.
#                        This should still leave the Ly-Alpha forest largely unaffected,
#                        but should be faster. It is recommended to set GENERATIONS equal
#                        to 1 for maximum speed-up.

# - ISOTHERM_EQS:
#     This special option makes the gas behave like an isothermal gas
#     with equation of state P = cs^2 * rho. The sound-speed cs is set by 
#     the thermal energy per unit mass in the intial conditions, 
#     i.e. cs^2=u. If the value for u is zero, then the initial gas 
#     temperature in the parameter file is used to define the sound speed
#     according to cs^2 = kT/mp, where mp is the proton mass.
#
