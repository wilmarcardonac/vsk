Module fiducial

    Implicit none

    save 

    Real*8,parameter    :: fwhm_arcmin = 5.0d0
    Real*8,parameter    :: miss = -1.63750d-13
    Real*8,parameter :: Pi = 3.141592653589793d0

    Integer*4,parameter :: simul_type = 1            ! TYPE OF SIMULATION 
    Integer*4,parameter :: nsmax = 2048              ! Nside FOR CMB MAP
    Integer*4,parameter :: nlmax = 2*nsmax           ! HIGHEST MULTIPOLE 
    Integer*4,parameter :: number_of_cmb_simulations = 1000   ! NUMBER OF GAUSSIAN CMB MAPS 
    Integer*4,parameter :: nsideC = 2                ! NSIDE FOR V, S, AND K MAPS
    Integer*4 :: npixC                               ! NUMBER OF PIXELS IN V (S OR K) MAPS
    Integer*8 :: n,c,fr                              ! NUMBER OF PIXELS CMB MAPS 
    Integer*4,parameter :: nmasks = 2                ! NUMBER OF MASKS  
    Integer*4,parameter :: ncmbmaps = 2              ! NUMBER OF CMB MAPS IN SMICA/COMMANDER/NILC/SEVEM FITS FILES
    Integer*4,parameter :: UNIT_EXE_FILE = 90           ! UNIT NUMBER FOR EXECUTION INFORMATION FILE
    Integer*4,parameter :: UNIT_SYNFAST_PAR_FILE = 91   ! UNIT NUMBER FOR FILE

    Logical,parameter   :: do_cmb_simulations = .false. ! DO CMB SIMULATIONS IF SET IT TRUE
    Logical,parameter   :: compute_vsk_maps = .true.    ! COMPUTE V,S,K MAPS IF SET IT TRUE
    Logical,parameter   :: compute_vsk_angular_power_spectrum = .true. 

    Character(len=*),parameter :: infile = './cmb_angular_power_spectrum/planck2015_lcdm_cl_v2.fits'    ! PATH TO SEED CMB ANGULAR POWER SPECTRUM
    Character(len=*),parameter :: beam_file = " '' "    ! PATH TO BEAM FILE
    Character(len=*),parameter :: almsfile = " '' "     ! PATH ALMS FILE
    Character(len=*),parameter :: plmfile = " '' "      ! PATH TO PLM FILE
    Character(len=*),parameter :: ORDERING_VSK_MAPS = 'RING'!'NESTED'    ! PATH TO TEMPORARY GAUSSIAN CMB MAP
    Character(len=*),parameter :: outfile_alms = " '' "    ! PATH TO ALMS OUTPUT FILE 
    Character(len=*),parameter :: PATH_TO_SYNFAST_PARAMETER_FILE = './synfast_parameter_files/g_cmb_' ! PATH TO SYNFAST PARAMETERS FILES
    Character(len=*),parameter :: PATH_TO_POLSPICE_PARAMETER_FILE = './polspice_parameter_files/' ! PATH TO POLSPICE PARAMETERS FILES
    Character(len=*),parameter :: PATH_TO_CMB_MAPS = './cmb_maps/map_g_' ! PATH TO SIMULATED CMB MAPS
    Character(len=*),parameter :: PATH_TO_HEALPIX_DATA = './Healpix_3.30/data' ! PATH TO HEALPIX DATA
    Character(len=*),parameter :: SYS_COORD = 'G' ! MAP COORDINATE SYSTEM 
    Character(len=*),parameter :: EXECUTION_INFORMATION = './output/execution_information.txt' ! PATH TO EXECUTION INFORMATION FILE
    character(len=*),parameter :: fmt = '(I4.4)'  ! FORMAT NUMBERS 1 - 1000
    Character(len=*),parameter :: PATH_TO_CMB_MASK = './data/COM_CMB_IQU-common-field-MaskInt_2048_R2.01.fits' ! PLANCK MASK TO BE USED (NESTED ORDERING)
    Character(len=*),parameter :: PATH_TO_PLANCK_CMB_MAP = './data/COM_CMB_IQU-smica-field-Int_2048_R2.01_full.fits' ! PLANCK CMB MAP TO BE USED (NESTED ORDERING)
    Character(len=*),parameter :: PATH_TO_VSK_MASK = './vsk_maps/vsk_mask.fits' ! VSK MASK TO BE USED (RING ORDERING)
    Character(len=*),parameter :: PATH_TO_VSK_SPECTRA = './vsk_angular_power_spectrum/' 

End Module fiducial
