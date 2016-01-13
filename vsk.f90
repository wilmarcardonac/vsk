Program vsk 

  !####################
  ! LOAD NEEDED MODULES 
  !####################

  use healpix_types
  use udgrade_nr, only: udgrade_ring, udgrade_nest
  use pix_tools, only: nside2npix,convert_ring2nest, convert_nest2ring, remove_dipole
  use fitstools, only: getsize_fits, input_map, output_map
  use head_fits
  use fiducial
  use arrays
  use functions 
  use fgsl

  !#################################
  ! DECLARE VARIABLES AND PARAMETERS
  !#################################

  Implicit none
    
  Integer*4 :: m                                   ! INTERGER FOR SHORT LOOPS 

  Real*8,dimension(0:3) :: multipoles              ! SAVES MONOPOLE AND DIPOLE OF CMB MAP 
  Real*8,dimension(1:2) :: zbounds                 ! BOUNDS TO COMPUTE DIPOLE AND MONOPOLE
!!$    Logical,dimension(number_of_parameters) :: plausibility  

  Character(len=4) :: x

  Character(len=80),dimension(1:60) :: header

!  Real*8, allocatable, dimension(:,:) :: maptest

  !############################
  ! INITIALIZATION OF VARIABLES
  !############################

  zbounds(:) = 0.d0 

  npixC = nside2npix(nsideC)

  n = nside2npix(nsmax)

  c = n/npixC     ! NUMBER OF CMB PIXELS PER CELL

  fr = nint(c*0.8) ! FRACTION OF PIXELS TO DETERMINE UNMASKED CELLS

  allocate (cmbmask(0:n-1,1:nmasks),planckmap(0:n-1,1:ncmbmaps), stat = status1)

  call input_map(PATH_TO_CMB_MASK, cmbmask(0:n-1,1:nmasks), n, nmasks) ! READ CMB MASK IN DEFAULT PLANCK ORDERING: NESTED

  call convert_nest2ring(nsmax,cmbmask(0:n-1,1:nmasks))  ! CHANGE ORDERING OF CMB MASK: NESTED->RING

!  call write_minimal_header(header, 'MAP', nside = nsmax, ordering = ORDERING_VSK_MAPS, coordsys = SYS_COORD) ! HEADER OF V, S, K-MAPS

!  call output_map(cmbmask(0:n-1,1:1),header,'ut78.fits')

  call input_map(PATH_TO_PLANCK_CMB_MAP, planckmap(0:n-1,1:ncmbmaps), n, ncmbmaps) ! READ PLANCK CMB MAP IN DEFAULT PLANCK ORDERING: NESTED

  call convert_nest2ring(nsmax,planckmap(0:n-1,1:ncmbmaps))  ! CHANGE ORDERING OF CMB MAP: NESTED->RING

!  call write_minimal_header(header, 'MAP', nside = nsmax, ordering = ORDERING_VSK_MAPS, coordsys = SYS_COORD) ! HEADER OF V, S, K-MAPS

!  call output_map(planckmap(0:n-1,1:1),header,'smica.fits')

  call remove_dipole(nsmax,planckmap(0:n-1,1),1,2,multipoles,zbounds,cmbmask(0:n-1,1))

  open(UNIT_EXE_FILE,file=EXECUTION_INFORMATION)
    
  !#################################################
  ! GENERATE GAUSSIAN CMB SIMULATIONS (IF NECESSARY)
  !#################################################

  If (do_cmb_simulations) then
     
     write(UNIT_EXE_FILE,*) 'COMPUTING GAUSSIAN CMB SIMULATIONS'

     call generate_gaussian_cmb_map()

  Else

     write(UNIT_EXE_FILE,*) 'USING EXISTING GAUSSIAN CMB SIMULATIONS'

  End If

  If (compute_vsk_maps) then

     call compute_variance_skewness_kurtosis_maps(PATH_TO_PLANCK_CMB_MAP,'0000')

     Do m = 1,number_of_cmb_simulations

        write(x,fmt) m

        call compute_variance_skewness_kurtosis_maps(PATH_TO_CMB_MAPS//trim(x)//'.fits',x)

     End Do

  Else

     continue 

  End If

  call compute_vsk_mean_maps()

  ! mkdir vsk_maps_zero_mean

  ! Subroutine : input_map -> subtract mean -> output_map

  ! change target directory in subroutines below

  If (compute_vsk_angular_power_spectrum) then

     call write_parameter_file_polspice(0,'V')

     call write_parameter_file_polspice(0,'S')

     call write_parameter_file_polspice(0,'K')

     call system('./PolSpice_v03-01-06/src/spice -optinfile '//trim(PATH_TO_POLSPICE_PARAMETER_FILE)//&
          ''//trim('vmap_smica')//'.spicerc')

     call system('./PolSpice_v03-01-06/src/spice -optinfile '//trim(PATH_TO_POLSPICE_PARAMETER_FILE)//&
          ''//trim('smap_smica')//'.spicerc')

     call system('./PolSpice_v03-01-06/src/spice -optinfile '//trim(PATH_TO_POLSPICE_PARAMETER_FILE)//&
          ''//trim('kmap_smica')//'.spicerc')

     stop
     call compute_vsk_angular_power_spectra()

  Else

     continue 

  End If

  deallocate(cmbmask,planckmap)

  close(UNIT_EXE_FILE)

End Program vsk




