Module fiducial

    Implicit none

    save 

    !#############################
    ! PARAMETERS OF FIDUCIAL MODEL
    !#############################

    Real*8,parameter    :: fwhm_arcmin = 5.0d0
!!$    Real*8,parameter    :: prior_bw = -3.23d0
!!$    Real*8,parameter    :: prior_sigma_int = 0.0d0 !0.2d0 !0.d0 !0.1d0!0.20d0 SN Ia hosts
!!$    Real*8,parameter    :: prior_sigma_int_LMC = 0.0d0!0.113d0 ! SAME VALUE AS IN EQUATION (4a) OF EFSTATHIOU'S PAPER
!!$    Real*8,parameter    :: prior_sigma_int_MW = 0.10d0 ! SAME VALUE AS IN SUBSECTION 4.2 OF EFSTATHIOU'S PAPER
!!$    Real*8,parameter    :: prior_alpha_j = 5.d-1
!!$    Real*8,parameter    :: R = 0.410d0                  ! TAKEN FROM PAGE 7 IN R11
!!$    Real*8,parameter    :: a_v = 0.697d0                ! TAKEN FROM PAGE 9 IN R11
!!$    Real*8,parameter    :: a_cal = 0.d0                 ! AUXILIAR PARAMETER TO HAVE DIAGNONAL COVARIANCE MATRIX FOR LMC CEPHEID VARIABLES
!!$    Real*8,parameter    :: NGC4258_distance = 7.60d0    ! TAKEN FROM PAGE 1 IN H13. UNITS : MPC
!!$    Real*8,parameter    :: mu_0_NGC4258 = 5.d0*log10(NGC4258_distance) + 25.d0 ! DEFINITION OF DISTANCE MODULUS
!!$    Real*8,parameter    :: LMC_distance = 49.97d-3       ! TAKEN FROM PAGE 76 IN PIETRZYNSKI. UNITS : MPC
!!$    Real*8,parameter    :: mu_0_LMC = 5.d0*log10(LMC_distance) + 25.d0 
!!$    Real*8,parameter    :: meanOH_LMC = 8.5d0           ! MEAN METALLICITY FOR LMC CEPHEID VARIABLES ASSUMED BY EFSTATHIOU
!!$    Real*8,parameter    :: meanOH_MW = 8.9d0           ! MEAN METALLICITY FOR MW CEPHEID VARIABLES ASSUMED BY EFSTATHIOU
!!$    Real*8,parameter    :: prior_zpH = 28.d0
!!$    Real*8,parameter    :: prior_bH = -2.7d0
!!$    Real*8,parameter    :: prior_mu1 = 30.91d0 ! FROM TABLE 3 IN R11 
!!$    Real*8,parameter    :: prior_mu2 = 31.67d0 ! FROM TABLE 3 IN R11 
!!$    Real*8,parameter    :: prior_mu3 = 32.13d0 ! FROM TABLE 3 IN R11 
!!$    Real*8,parameter    :: prior_mu4 = 31.70d0 ! FROM TABLE 3 IN R11 
!!$    Real*8,parameter    :: prior_mu5 = 32.27d0 ! FROM TABLE 3 IN R11 
!!$    Real*8,parameter    :: prior_mu6 = 32.59d0 ! FROM TABLE 3 IN R11 
!!$    Real*8,parameter    :: prior_mu7 = 31.72d0 ! FROM TABLE 3 IN R11 
!!$    Real*8,parameter    :: prior_mu8 = 31.66d0 ! FROM TABLE 3 IN R11 
!!$    Real*8,parameter    :: prior_mu9 = 24.8d0
!!$    Real*8,parameter    :: prior_mu10 = 18.5d0
!!$    Real*8,parameter    :: prior_zpw = 29.d0
!!$    Real*8,parameter    :: prior_zpw4258 = 30.5d0  ! CENTRAL VALUE OF PRIOR ON zp_{w,4258} 
!!$    Real*8,parameter    :: prior_zpwLMC = 20.98d0
!!$    Real*8,parameter    :: prior_Mw = -5.88d0
!!$    Real*8,parameter    :: prior_Zw = 0.d0         ! CENTRAL VALUE FOR PRIOR ON Zw 
!!$    Real*8,parameter    :: prior_H0 = 70.0d0
!!$    Real*8,parameter    :: prior_bw_from_LMC = -3.35d0 ! CENTRAL VALUE FOR PRIOR ON bw


    !################
    ! SYNFAST PARAMETERS
    !################

    Integer*4,parameter :: simul_type = 1            ! TYPE OF SIMULATION 
    Integer*4,parameter :: nsmax = 2048              ! Nside FOR CMB MAP
    Integer*4,parameter :: nlmax = 2*nsmax           ! HIGHEST MULTIPOLE 
    Integer*4,parameter :: number_of_cmb_simulations = 1000   ! NUMBER OF GAUSSIAN CMB MAPS 
!!$    Integer*4,parameter :: jumping_factor_update = 100           ! NUMBER OF TAKEN STEPS BEFORE UPDATING JUMPING FACTOR (IF NEEDED)
!!$    Integer*4,parameter :: covariance_matrix_update = 10000        ! STEPS TAKEN BEFORE UPDATING COVARIANCE MATRIX (IF NEEDED)
!!$    Integer*4,parameter :: steps_taken_before_definite_run = 1000000 ! STEPS TAKEN BEFORE DEFINITE RUN
!!$    Integer*4,parameter :: number_of_hosts_galaxies = 9 ! TOTAL NUMBER OF HOSTS GALAXIES AS IN R11 (NUMBER INCLUDES NGC4258)
    Integer*4,parameter :: UNIT_EXE_FILE = 90           ! UNIT NUMBER FOR EXECUTION INFORMATION FILE
    Integer*4,parameter :: UNIT_SYNFAST_PAR_FILE = 91           ! UNIT NUMBER FOR FILE
!!$    Integer*4,parameter :: UNIT_PARAMNAMES_FILE = 92           ! UNIT NUMBER FOR PARAMNAMES FILE
!!$    Integer*4,parameter :: UNIT_MCMC_FILE = 93           ! UNIT NUMBER FOR MCMC OUTPUT FILE (CALIBRATING PHASE)
!!$    Integer*4,parameter :: UNIT_MCMC_FINAL_FILE = 94           ! UNIT NUMBER FOR MCMC FINAL OUTPUT FILE
!!$    Integer*4,parameter :: UNIT_HP_FILE = 95           ! UNIT EFFECTIVE HPS FILE
!!$
!!$    Real*8,parameter    :: step_size_changes = 1.d-2             ! CHANGES IN STEP SIZE
!!$    Real*8,parameter    :: cepheid_Period_limit = 205.d0 !205.d0 !60.d0           ! DISREGARD CEPHEID VARIABLES WITH PERIOD GREATER THAN cepheid_Period_limit
!!$    Real*8,parameter    :: cepheid_lower_Period_limit = 0.d0                    ! DISREGARD CEPHEID VARIABLES WITH PERIOD SHORTER THAN cepheid_lower_Period_limit

    Logical,parameter   :: separate_dataA = .false.!.true.               ! INCLUDE DATA SET A AS SINGLE POINTS IF SET IT TRUE
!!$    Logical,parameter   :: separate_dataB = .false.!.true.               ! INCLUDE DATA SET B AS SINGLE POINTS IF SET IT TRUE
!!$    Logical,parameter   :: separate_dataC = .false.!.true.               ! INCLUDE DATA SET C AS SINGLE POINTS IF SET IT TRUE
!!$    Logical,parameter   :: include_dataA = .false.!.true.                ! INCLUDE DATA SET A IF SET IT TRUE
!!$    Logical,parameter   :: include_dataB = .false.!.true.                ! INCLUDE DATA SET B IF SET IT TRUE
!!$    Logical,parameter   :: include_dataC = .false.!.true.                ! INCLUDE DATA SET C IF SET IT TRUE
!!$    Logical,parameter   :: include_table2_R11 = .true.            ! INCLUDE TABLE 2 IN R11 IF SET IT TRUE
!!$    Logical,parameter   :: start_from_fiducial = .true.          ! START MCMC ANALYSIS FROM FIDUCIAL POINT IF SET IT TRUE 
!!$    Logical,parameter   :: testing_Gaussian_likelihood = .false. ! TEST GAUSSIAN LIKELIHOOD IF SET IT TRUE
!!$    Logical,parameter   :: using_hyperparameters = .true.        ! USE HYPER-PARAMETERS IF SET IT TRUE
!!$    Logical,parameter   :: using_jeffreys_prior = .false.        ! USE JEFFREYS PRIOR IF SET IT TRUE, OTHERWISE USE UNIFORM PRIOR [0,1] 
!!$    Logical,parameter   :: hyperparameters_as_mcmc = .false.      ! SET HYPER-PARAMETERS AS MCMC PARAMETERS IF SET IT TRUE
!!$    Logical,parameter   :: use_NGC4258_as_anchor = .false.!.true.       ! USE NFC4258 AS ANCHOR IF SET IT TRUE
!!$    Logical,parameter   :: use_LMC_as_anchor = .true.           ! USE LMC AS ANCHOR IF SET IT TRUE
!!$    Logical,parameter   :: use_MW_as_anchor = .false.!.true.            ! USE MW AS ANCHOR IF SET IT TRUE
!!$    Logical,parameter   :: use_metallicity = .true.             ! USE METALLICITY DEPENDENCE IF SET IT TRUE
!!$    Logical,parameter   :: use_H_band = .false.!.true.                   ! USE H BAND IF SET IT TRUE, OTHERWISE USE W BAND
!!$    Logical,parameter   :: use_HP_in_SNIa = .true.!.false.               ! USE HPs WHEN COMPUTING SNIa CHI2
!!$    Logical,parameter   :: use_HP_in_av = .false.                ! USE HPs WHEN COMPUTINNG av CHI2
!!$    Logical,parameter   :: use_HP_in_anchor = .true.!.false.            ! USE HPs WHEN COMPUTING ANCHOR CHI2
!!$    Logical,parameter   :: use_HP_per_host = .false.              ! USE HPs FOR EACH HOST IN R11 IF SET IT TRUE
!!$    Logical,parameter   :: use_HP_per_cepheid = .true.           ! USE HPs FOR EACH CEPHEID IN R11 IF SET IT TRUE
!!$    Logical,parameter   :: use_HP_per_MW_cepheid = .false.       ! USE HPs FOR EACH CEPHEID IN MW IF SET IT TRUE
!!$    Logical,parameter   :: use_HP_for_MW_dataset = .true.       ! USE HP FOR MW DATASET IF SET IT TRUE (JEFFREY'S PRIOR)
!!$    Logical,parameter   :: doing_R11_analysis = .true.           ! DO R11 ANALYSIS IF SET IT TRUE, OTHERWISE DO EFSTATHIOU'S SECTION 2 (LMC CEPHEIDS ALONE)
!!$    Logical,parameter   :: include_only_cepheids = .false.       ! INCLUDE ONLY CEPHEIDS DATA IF SET IT TRUE
!!$    Logical,parameter   :: all_R11_hosts = .false.             ! INCLUDE ALL CEPHEIDS IN R11 SAMPLE SIMULTANEOUSLY IF SET IT TRUE
!!$    Logical,parameter   :: use_prior_on_zpw4258 = .false. !.true.       ! USE PRIOR ON zp_{w,4258} IS SET IT TRUE
!!$    Logical,parameter   :: use_prior_on_Zw = .false.!.true.              ! USE PRIOR ON Zw IF SET IT TRUE 
!!$    Logical,parameter   :: use_prior_on_bw = .false.!.true.              ! USE PRIOR ON bw IF SET IT TRUE
!!$    Logical,parameter   :: use_HP_in_Zw = .false.                 ! USE HPs WHEN USING PRIOR ON THE METALLICITY IF SET IT TRUE 

    Character(len=*),parameter :: infile = './cmb_angular_power_spectrum/planck2015_lcdm_cl_v2.fits'    ! PATH TO SEED CMB ANGULAR POWER SPECTRUM
    Character(len=*),parameter :: beam_file = " '' "    ! PATH TO BEAM FILE
    Character(len=*),parameter :: almsfile = " '' "     ! PATH ALMS FILE
    Character(len=*),parameter :: plmfile = " '' "      ! PATH TO PLM FILE
!    Character(len=*),parameter :: outfile = './cmb_maps/g_cmb.fits'    ! PATH TO TEMPORARY GAUSSIAN CMB MAP
    Character(len=*),parameter :: outfile_alms = " '' "    ! PATH TO ALMS OUTPUT FILE 
    Character(len=*),parameter :: PATH_TO_SYNFAST_PARAMETER_FILE = './synfast_parameter_files/g_cmb_' ! PATH TO DATA OF TABLE 2 IN R11
    Character(len=*),parameter :: PATH_TO_CMB_MAPS = './cmb_maps/map_g_' ! PATH TO DATA OF TABLE 2 IN R11
    Character(len=*),parameter :: PATH_TO_HEALPIX_DATA = './Healpix_3.30/data' ! PATH TO DATA OF TABLE 2 IN LEEUWEN
!!$    Character(len=5),dimension(number_of_hosts_galaxies), parameter :: host = ['n4536','n4639','n3982','n3370','n3021','n1309',&
!!$    'n4038','n5584','n4258'] ! HOST GALAXIES IN SAME ORDER LISTED IN TABLE 2 OF R11
!!$    Character(len=*),parameter :: EXECUTION_INFORMATION = './output/chains/execution_information.txt' ! PATH TO EXECUTION INFORMATION FILE

End Module fiducial
