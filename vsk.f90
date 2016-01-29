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

  Real(kind=DP),dimension(0:DEGREE_REMOVE_DIPOLE*DEGREE_REMOVE_DIPOLE-1) :: multipoles              ! SAVES MONOPOLE AND DIPOLE OF CMB MAP 
  Real(kind=DP),dimension(1:2) :: zbounds                 ! BOUNDS TO COMPUTE DIPOLE AND MONOPOLE

  Logical :: exist

  Character(len=4) :: x

!  Character(len=80),dimension(1:60) :: header

!  Real(kind=DP), dimension(0:2048-1) :: maptest

  !############################
  ! INITIALIZATION OF VARIABLES
  !############################

  zbounds(:) = 0.d0 

  npixC = nside2npix(nsideC)

  n = nside2npix(nsmax)

  c = n/npixC     ! NUMBER OF CMB PIXELS PER CELL

  fr = nint(c*0.9) ! FRACTION OF PIXELS TO DETERMINE UNMASKED CELLS

  allocate (cmbmask(0:n-1,1:nmasks),planckmap(0:n-1,1:ncmbmaps), stat = status1)

  call input_map(PATH_TO_CMB_MASK, cmbmask(0:n-1,1:nmasks), n, nmasks,HPX_DBADVAL) ! READ CMB MASK IN DEFAULT PLANCK ORDERING: NESTED

  call convert_nest2ring(nsmax,cmbmask(0:n-1,1:nmasks))  ! CHANGE ORDERING OF CMB MASK: NESTED->RING

  call input_map(PATH_TO_PLANCK_CMB_MAP, planckmap(0:n-1,1:ncmbmaps), n, ncmbmaps,HPX_DBADVAL) ! READ PLANCK CMB MAP IN DEFAULT PLANCK ORDERING: NESTED. UNITS: K_CMB

  call convert_nest2ring(nsmax,planckmap(0:n-1,1:ncmbmaps))  ! CHANGE ORDERING OF CMB MAP: NESTED->RING

!  call remove_dipole(nsmax,planckmap(0:n-1,1),RING_ORDERING,DEGREE_REMOVE_DIPOLE,multipoles,zbounds,HPX_DBADVAL,cmbmask(0:n-1,1))

  open(UNIT_EXE_FILE,file=EXECUTION_INFORMATION)

  write(UNIT_EXE_FILE,*) 'CMB MASK READ. MONOPOLE AND DIPOLE REMOVED FROM PLANCK CMB MAP'

  write(UNIT_EXE_FILE,*) 'MONOPOLE AND DIPOLE GIVEN BY HEALPIX SUBROUTINE FOR PLANCK MAP ', multipoles
    
  !#################################################
  ! GENERATE GAUSSIAN CMB SIMULATIONS (IF NECESSARY)
  !#################################################

  If (do_cmb_simulations) then
     
     write(UNIT_EXE_FILE,*) 'COMPUTING GAUSSIAN CMB SIMULATIONS. THE CODE DOES NOT OVERWRITE FILES, SO '

     write(UNIT_EXE_FILE,*) 'IF SIMULATIONS ARE FOUND, ONLY THOSE NON-EXISTING WILL BE GENERATED '

     call generate_gaussian_cmb_map()

  Else

     write(UNIT_EXE_FILE,*) 'USING EXISTING GAUSSIAN CMB SIMULATIONS'

  End If

  If (compute_vsk_maps) then

     write(UNIT_EXE_FILE,*) 'COMPUTING VSK MAPS. IF NEW COMPUTATION REQUIRED, THEN CLEAN FOLDER VSK_MAPS'

     inquire(file ='./vsk_maps/vmap_smica.fits',exist=exist)

     If (exist) then
     
        continue

     Else

        call compute_variance_skewness_kurtosis_maps(PATH_TO_PLANCK_CMB_MAP,'0000')

     End If

     Do m = 1,number_of_cmb_simulations

        write(x,fmt) m

        inquire(file ='./vsk_maps/vmap_'//trim(x)//'.fits',exist=exist)

        If (exist) then

           continue
           
        Else

           call compute_variance_skewness_kurtosis_maps(PATH_TO_CMB_MAPS//trim(x)//'.fits',x)

        End If

     End Do

  Else

     continue 

  End If

  If (compute_mean_vsk_maps)  then

     call compute_vsk_mean_maps()

     Do m = 0,npixC-1

        If ( (vsdv(m,1)/vmean(m,1) .gt. 1.d-2) .or. (ssdv(m,1)/smean(m,1) .gt. 1.d-2) .or. &
             ( ksdv(m,1)/kmean(m,1) .gt. 1.d-2) ) then 

           write(UNIT_EXE_FILE,*) 'V, S, K ESTIMATORS HAVE SDV/MEAN RATIO EQUAL TO: ', vsdv(m,1)/vmean(m,1), &
                ssdv(m,1)/smean(m,1), ksdv(m,1)/kmean(m,1), 'FOR PIXEL ', m+1 

           write(UNIT_EXE_FILE,*) 'Nside OF CMB MAP IS ', nsmax

           write(UNIT_EXE_FILE,*) 'Nside OF VSK MAPS IS ', nsideC

        Else

           continue

        End If

     End Do

  Else

     continue

  End If

  If (compute_vsk_angular_power_spectrum) then

     write(UNIT_EXE_FILE,*) 'NOT USING BEAM TO SMOOTH V, S, K MAPS. IF BEAM WANTED, UNCOMMENT CORRESPONDING LINE IN  '

     write(UNIT_EXE_FILE,*) 'SUBROUTINE "write_parameter_file_polspice" '

     call write_parameter_file_polspice(0,'V')

     call write_parameter_file_polspice(0,'S')

     call write_parameter_file_polspice(0,'K')

     call system('./PolSpice_v03-01-06/src/spice -optinfile '//trim(PATH_TO_POLSPICE_PARAMETER_FILE)//&
          ''//trim('vmap_smica')//'.spicerc')

     call system('./PolSpice_v03-01-06/src/spice -optinfile '//trim(PATH_TO_POLSPICE_PARAMETER_FILE)//&
          ''//trim('smap_smica')//'.spicerc')

     call system('./PolSpice_v03-01-06/src/spice -optinfile '//trim(PATH_TO_POLSPICE_PARAMETER_FILE)//&
          ''//trim('kmap_smica')//'.spicerc')

     call compute_vsk_angular_power_spectra()

  Else

     continue 

  End If

  deallocate(cmbmask,planckmap)

  call system('python ./figures/analyze.py')

  close(UNIT_EXE_FILE)

End Program vsk




