Program vsk 

  !####################
  ! LOAD NEEDED MODULES 
  !####################

  use healpix_types
  use udgrade_nr, only: udgrade_ring, udgrade_nest
  use pix_tools, only: nside2npix,convert_ring2nest, convert_nest2ring
  use fitstools, only: getsize_fits, input_map, output_map
  use head_fits
  use fiducial
  use arrays
  use functions 

  !#################################
  ! DECLARE VARIABLES AND PARAMETERS
  !#################################

  Implicit none
    
  Integer*4 :: m                                   ! INTERGER FOR SHORT LOOPS 
  Integer*8 :: i,d,f,j,f2                        ! NUMBER OF PIXELS CMB MAPS 

!!$    Real*8 :: random_uniform                           ! SAVES RANDOM UNIFORM DEVIATE BETWEEN 0 AND 1 
  
!!$    Logical,dimension(number_of_parameters) :: plausibility  

  Character(len=80),dimension(1:60) :: header

  Real*8, allocatable, dimension(:,:) :: maptest

  !############################
  ! INITIALIZATION OF VARIABLES
  !############################

  npixC = nside2npix(nsideC)

  n = nside2npix(nsmax)

  c = n/npixC     ! NUMBER OF CMB PIXELS PER CELL

  fr = nint(c*0.8) ! FRACTION OF PIXELS TO DETERMINE UNMASKED CELLS

  allocate (cmbmask(0:n-1,1:nmasks),planckmap(0:n-1,1:ncmbmaps), stat = status1)

  call input_map(PATH_TO_CMB_MASK, cmbmask(0:n-1,1:nmasks), n, nmasks) ! READ CMB MASK IN DEFAULT PLANCK ORDERING: NESTED

  call convert_nest2ring(nsmax,cmbmask(0:n-1,1:nmasks))  ! CHANGE ORDERING OF CMB MASK: NESTED->RING

!  call input_map(PATH_TO_PLANCK_CMB_MAP, planckmap(0:n-1,1:ncmbmaps), n, ncmbmaps) ! READ PLANCK CMB MAP IN DEFAULT PLANCK ORDERING: NESTED

 ! call convert_nest2ring(nsmax,planckmap(0:n-1,1:ncmbmaps))  ! CHANGE ORDERING OF CMB MAP: NESTED->RING

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

  call compute_variance_skewness_kurtosis_maps(PATH_TO_CMB_MAPS//trim('0001')//'.fits')

!!$
!!$                write(UNIT_RANGES_FILE,*) ''//trim(paramnames(3))//'    20.    40.'
!!$
!!$                write(UNIT_RANGES_FILE,*) ''//trim(paramnames(4))//'    20.    40.'
!!$
!!$                write(UNIT_RANGES_FILE,*) ''//trim(paramnames(5))//'    20.    40.'
!!$
!!$                write(UNIT_RANGES_FILE,*) ''//trim(paramnames(6))//'    20.    40.'
!!$
!!$                write(UNIT_RANGES_FILE,*) ''//trim(paramnames(7))//'    20.    40.'
!!$
!!$                write(UNIT_RANGES_FILE,*) ''//trim(paramnames(8))//'    20.    40.'
!!$
!!$                write(UNIT_RANGES_FILE,*) ''//trim(paramnames(9))//'    20.    30.'
!!$
!!$                write(UNIT_RANGES_FILE,*) ''//trim(paramnames(10))//'    25.    34.'
!!$
!!$                write(UNIT_RANGES_FILE,*) ''//trim(paramnames(11))//'    -3.2    -2.5'
!!$
!!$                write(UNIT_RANGES_FILE,*) ''//trim(paramnames(12))//'    -1.    1.'
!!$
!!$             Else
!!$
!!$                write(UNIT_RANGES_FILE,*) ''//trim(paramnames(1))//'    25.    34.'
!!$
!!$                write(UNIT_RANGES_FILE,*) ''//trim(paramnames(2))//'    -3.2    -2.5'
!!$
!!$                write(UNIT_RANGES_FILE,*) ''//trim(paramnames(3))//'    -1.    1.'
!!$
!!$             End If
!!$
!!$          Else
!!$
!!$             If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor) then
!!$    
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$                    
!!$                   Else
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(1))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(2))//'    20.    40.'
!!$                      
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(3))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(4))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(5))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(6))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(7))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(8))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(9))//'    20.    30.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(10))//'    0.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(11))//'    -7.    -4.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(12))//'    -3.5    -2.5'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(13))//'    55.    95.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(14))//'    -2.    1.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(15))//'    0.    1.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(16))//'    -1.    1.'
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC+MW+NGC4258 AS ANCHORS'
!!$                     
!!$                      stop
!!$ 
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC+MW+NGC4258 AS ANCHORS'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$                    
!!$                   Else
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(1))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(2))//'    20.    40.'
!!$                      
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(3))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(4))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(5))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(6))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(7))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(8))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(9))//'    20.    30.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(10))//'    0.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(11))//'    -7.    -4.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(12))//'    -3.5    -2.5'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(13))//'    55.    95.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(14))//'    -2.    1.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(15))//'    0.    1.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(16))//'    -1.    1.'
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AND NGC4258 AS ANCHORS'
!!$                     
!!$                      stop
!!$ 
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AND NGC4258 AS ANCHORS'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$                    
!!$                   Else
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(1))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(2))//'    20.    40.'
!!$                      
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(3))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(4))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(5))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(6))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(7))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(8))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(9))//'    20.    30.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(10))//'   -7.    -4.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(11))//'    -3.5    -2.5'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(12))//'    55.    95.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(13))//'    -2.    1.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(14))//'    0.    1.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(15))//'    -1.    1.'
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN MW+NGC4258 AS ANCHORS'
!!$                     
!!$                      stop
!!$ 
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN MW+NGC4258 AS ANCHORS'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$                    
!!$                   Else
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(1))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(2))//'    20.    40.'
!!$                      
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(3))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(4))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(5))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(6))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(7))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(8))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(9))//'    20.    30.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(10))//'    0.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(11))//'    -7.    -4.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(12))//'    -3.5    -2.5'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(13))//'    55.    95.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(14))//'    -2.    1.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(15))//'    0.    1.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(16))//'    -1.    1.'
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AND MW AS ANCHORS'
!!$                     
!!$                      stop
!!$ 
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AND MW AS ANCHORS'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( .not.use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$                    
!!$                   Else
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(1))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(2))//'    20.    40.'
!!$                      
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(3))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(4))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(5))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(6))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(7))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(8))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(9))//'    20.    30.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(10))//'   -7.    -4.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(11))//'    -3.5    -2.5'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(12))//'    55.    95.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(13))//'    -2.    1.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(14))//'    0.    1.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(15))//'    -1.    1.'
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN MW AS ANCHOR'
!!$                     
!!$                      stop
!!$ 
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN MW AS ANCHOR'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$                    
!!$                   Else
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(1))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(2))//'    20.    40.'
!!$                      
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(3))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(4))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(5))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(6))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(7))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(8))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(9))//'    20.    30.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(10))//'    0.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(11))//'    15.    25.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(12))//'    -3.5    -2.5'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(13))//'    55.    95.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(14))//'    -2.    1.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(15))//'    0.    1.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(16))//'    -1.    1.'
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AS ANCHOR'
!!$                     
!!$                      stop
!!$ 
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AS ANCHOR'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(1))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(2))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(3))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(4))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(5))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(6))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(7))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(8))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(9))//'    20.    30.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(10))//'    25.    34.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(11))//'    -3.2    -2.5'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(12))//'    55.    95.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(13))//'    -1.    1.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(14))//'    0.    1.'
!!$
!!$                   Else
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(1))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(2))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(3))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(4))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(5))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(6))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(7))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(8))//'    20.    40.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(9))//'    20.    30.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(10))//'    25.    34.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(11))//'    -3.5    -2.5'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(12))//'    55.    95.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(13))//'    -1.    1.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) ''//trim(paramnames(14))//'    0.    1.'
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      write(UNIT_RANGES_FILE,*) 'zpH1    0.    50.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) 'zpH2    0.    50.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) 'zpH3    0.    50.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) 'zpH4    0.    50.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) 'zpH5    0.    50.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) 'zpH6    0.    50.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) 'zpH7    0.    50.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) 'zpH8    0.    50.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) 'zpH4258    0.    50.'
!!$
!!$                      write(UNIT_RANGES_FILE,*) 'bH    -20.    0.'
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             End If ! OF ANCHORS
!!$
!!$          End If ! OF ONLY CEPHEIDS
!!$    
!!$       Else
!!$
!!$          write(UNIT_RANGES_FILE,*) ''//trim(paramnames(1))//'    0.    50.'
!!$
!!$          write(UNIT_RANGES_FILE,*) ''//trim(paramnames(2))//'    -20.    0.'
!!$
!!$          !    write(17,*) 'sigma_int    1.e-10    1 '
!!$
!!$       End If
!!$
!!$       If (hyperparameters_as_mcmc) then
!!$          ! WRITING HARD BOUNDS FOR HYPER-PARAMETERS
!!$          Do m=number_model_parameters+1,number_of_parameters
!!$
!!$             write(string,'(i2.2)') m-number_model_parameters
!!$                    
!!$             If (using_jeffreys_prior) then
!!$
!!$                print *,'MODIFY THIS PART ACCORDING TO PRIOR'
!!$
!!$                stop
!!$
!!$             Else
!!$
!!$                write(UNIT_RANGES_FILE,*) 'alpha_'//trim(string)//'    0.    1.'
!!$
!!$             End If
!!$
!!$          End Do
!!$            
!!$       End If
!!$
!!$       close(UNIT_RANGES_FILE)
!!$
!!$    End If
!!$
!!$    ! OPEN TEMPORARY FILE TO SAVE CHAIN
!!$    open(UNIT_MCMC_FILE,file='./output/mcmc_output.txt')     
!!$
!!$    write(UNIT_EXE_FILE,*) '# NUMBER OF ITERATIONS IN MCMC : ', number_iterations - steps_taken_before_definite_run
!!$
!!$    If (start_from_fiducial .and. .not.testing_Gaussian_likelihood) then
!!$
!!$        write(UNIT_EXE_FILE,*) '# FIDUCIAL MODEL IS (PARAMETERS ARE ORDERED AS IN CHAINS FILES) :', old_point
!!$
!!$        write(UNIT_EXE_FILE,'(a37,es18.10)') '# ln(L/L_max) AT THE FIDUCIAL MODEL :', old_loglikelihood
!!$
!!$    End If
!!$
!!$    If (hyperparameters_as_mcmc) then
!!$
!!$        write(UNIT_EXE_FILE,*) '# WEIGHT   -ln(L/L_{max})    A    bw   ', alpha_string(1:number_hyperparameters)
!!$ 
!!$    Else
!!$
!!$       If (doing_R11_analysis) then
!!$
!!$          If (include_only_cepheids) then
!!$
!!$             write(UNIT_EXE_FILE,*) '# WEIGHT   -ln(L/L_{max})    ', paramnames(1:number_model_parameters) 
!!$
!!$          Else
!!$
!!$             If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor) then
!!$    
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE WHEN LMC+MW+NGC4258 AS ANCHORS'
!!$                     
!!$                      stop
!!$                   
!!$                   Else
!!$
!!$                      write(UNIT_EXE_FILE,*) '# WEIGHT   -ln(L/L_{max})    ', paramnames(1:number_model_parameters) 
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC+NGC4258+MW AS ANCHOR'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC+NGC4258+MW AS ANCHOR'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE WHEN LMC AND NGC4258 AS ANCHORS'
!!$                     
!!$                      stop
!!$                   
!!$                   Else
!!$
!!$                      write(UNIT_EXE_FILE,*) '# WEIGHT   -ln(L/L_{max})    ', paramnames(1:number_model_parameters) 
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC+NGC4258 AS ANCHOR'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC+NGC4258 AS ANCHOR'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE WHEN MW+NGC4258 AS ANCHOR'
!!$                     
!!$                      stop
!!$                   
!!$                   Else
!!$
!!$                      write(UNIT_EXE_FILE,*) '# WEIGHT   -ln(L/L_{max})    ', paramnames(1:number_model_parameters) 
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN MW+NGC4258 AS ANCHOR'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN MW+NGC4258 AS ANCHOR'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE WHEN LMC AND MW AS ANCHORS'
!!$                     
!!$                      stop
!!$                   
!!$                   Else
!!$
!!$                      write(UNIT_EXE_FILE,*) '# WEIGHT   -ln(L/L_{max})    ', paramnames(1:number_model_parameters) 
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC+MW AS ANCHOR'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC+MW AS ANCHOR'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( .not.use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE WHEN MW AS ANCHOR'
!!$                     
!!$                      stop
!!$                   
!!$                   Else
!!$
!!$                      write(UNIT_EXE_FILE,*) '# WEIGHT   -ln(L/L_{max})    ', paramnames(1:number_model_parameters) 
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN MW AS ANCHOR'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN MW AS ANCHOR'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE WHEN LMC AS ANCHOR'
!!$                     
!!$                      stop
!!$                   
!!$                   Else
!!$
!!$                      write(UNIT_EXE_FILE,*) '# WEIGHT   -ln(L/L_{max})    ', paramnames(1:number_model_parameters) 
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AS ANCHOR'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AS ANCHOR'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      write(UNIT_EXE_FILE,*) '# WEIGHT   -ln(L/L_{max})    ', paramnames(1:number_model_parameters) 
!!$                   
!!$                   Else
!!$
!!$                      write(UNIT_EXE_FILE,*) '# WEIGHT   -ln(L/L_{max})    ', paramnames(1:number_model_parameters) 
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      write(UNIT_EXE_FILE,*) '# WEIGHT   -ln(L/L_{max})    zpH1    zpH2    zpH3'//trim(&
!!$                           '    zpH4    zpH5    zpH6    zpH7    zpH8    zpH4258    bH')//'' 
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             End If ! OF ANCHORS
!!$
!!$          End If ! ONLY CEPHEIDS
!!$    
!!$       Else
!!$
!!$          !write(13,*) '# Weight   -ln(L/L_{max})    A    bw    sigma_int ' 
!!$          write(UNIT_EXE_FILE,*) '# WEIGHT   -ln(L/L_{max})   ', paramnames(1:number_model_parameters) 
!!$
!!$       End If
!!$
!!$    End If
!!$
!!$    write(UNIT_EXE_FILE,*)'STARTING SAMPLING OF PARAMETER SPACE'
!!$
    ! LOOP TO EXPLORE PARAMETER SPACE STARTS HERE
!!$    Do m=1,number_iterations
!!$
!!$       ! GENERATE NEW POINT IN PARAMETER SPACE FROM MULTIVARIATE DISTRIBUTION; CODE USES RANLIB LIBRARY (BE CAREFUL WITH X_OLD AND OLD_POINT DEFINITIONS)
!!$       If (testing_Gaussian_likelihood) then
!!$
!!$          call setgmn(x_old,real(Covgauss),number_of_parameters,parm) 
!!$ 
!!$          call genmn(parm,x_new,work)
!!$
!!$       Else
!!$          
!!$          call setgmn(x_old,real(Covguess),number_of_parameters,parm) 
!!$
!!$          call genmn(parm,x_new,work)
!!$
!!$       End If
!!$
!!$       If (doing_R11_analysis) then
!!$
!!$          If (include_only_cepheids) then
!!$
!!$             If (all_R11_hosts) then
!!$
!!$                plausibility(1) = (x_new(1) .le. real(20.d0)) .or. (x_new(1) .ge. real(4.d1))
!!$
!!$                plausibility(2) = (x_new(2) .le. real(20.d0)) .or. (x_new(2) .ge. real(4.d1))
!!$
!!$                plausibility(3) = (x_new(3) .le. real(20.d0)) .or. (x_new(3) .ge. real(4.d1))
!!$
!!$                plausibility(4) = (x_new(4) .le. real(20.d0)) .or. (x_new(4) .ge. real(4.d1))
!!$
!!$                plausibility(5) = (x_new(5) .le. real(20.d0)) .or. (x_new(5) .ge. real(4.d1))
!!$
!!$                plausibility(6) = (x_new(6) .le. real(20.d0)) .or. (x_new(6) .ge. real(4.d1))
!!$
!!$                plausibility(7) = (x_new(7) .le. real(20.d0)) .or. (x_new(7) .ge. real(4.d1))
!!$
!!$                plausibility(8) = (x_new(8) .le. real(20.d0)) .or. (x_new(8) .ge. real(4.d1))
!!$
!!$                plausibility(9) = (x_new(9) .le. real(20.d0)) .or. (x_new(9) .ge. real(30.d0))
!!$
!!$                plausibility(10) =  (x_new(10) .le. real(25.d0)) .or. (x_new(10) .ge. real(34.d0)) 
!!$
!!$                plausibility(11) =  (x_new(11) .le. real(-3.2d0)) .or. (x_new(11) .ge. real(-2.5d0)) 
!!$
!!$                plausibility(12) =  (x_new(12) .le. real(-1.d0)) .or. (x_new(12) .ge. real(1.d0)) 
!!$
!!$             Else
!!$
!!$                plausibility(1) =  (x_new(1) .le. real(25.d0)) .or. (x_new(1) .ge. real(34.d0)) 
!!$
!!$                plausibility(2) =  (x_new(2) .le. real(-3.2d0)) .or. (x_new(2) .ge. real(-2.5d0)) 
!!$
!!$                plausibility(3) =  (x_new(3) .le. real(-1.d0)) .or. (x_new(3) .ge. real(1.d0)) 
!!$
!!$             End If
!!$
!!$          Else
!!$
!!$             If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor) then
!!$    
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      plausibility(1) = (x_new(1) .le. real(20.d0)) .or. (x_new(1) .ge. real(4.d1))
!!$                      
!!$                      plausibility(2) = (x_new(2) .le. real(20.d0)) .or. (x_new(2) .ge. real(4.d1))
!!$
!!$                      plausibility(3) = (x_new(3) .le. real(20.d0)) .or. (x_new(3) .ge. real(4.d1))
!!$
!!$                      plausibility(4) = (x_new(4) .le. real(20.d0)) .or. (x_new(4) .ge. real(4.d1))
!!$
!!$                      plausibility(5) = (x_new(5) .le. real(20.d0)) .or. (x_new(5) .ge. real(4.d1))
!!$
!!$                      plausibility(6) = (x_new(6) .le. real(20.d0)) .or. (x_new(6) .ge. real(4.d1))
!!$
!!$                      plausibility(7) = (x_new(7) .le. real(20.d0)) .or. (x_new(7) .ge. real(4.d1))
!!$
!!$                      plausibility(8) = (x_new(8) .le. real(20.d0)) .or. (x_new(8) .ge. real(4.d1))
!!$
!!$                      plausibility(9) = (x_new(9) .le. real(20.d0)) .or. (x_new(9) .ge. real(30.d0))
!!$
!!$                      plausibility(10) = (x_new(10) .le. real(0.d0)) .or. (x_new(10) .ge. real(40.d0))
!!$
!!$                      plausibility(11) =  (x_new(11) .le. real(-7.d0)) .or. (x_new(11) .ge. real(-4.d0)) 
!!$
!!$                      plausibility(12) =  (x_new(12) .le. real(-3.5d0)) .or. (x_new(12) .ge. real(-2.5d0)) 
!!$
!!$                      plausibility(13) =  (x_new(13) .le. real(55.d0)) .or. (x_new(13) .ge. real(95.d0)) 
!!$
!!$                      plausibility(14) =  (x_new(14) .le. real(-2.d0)) .or. (x_new(14) .ge. real(1.d0)) 
!!$
!!$                      plausibility(15) =  (x_new(15) .le. real(0.d0)) .or. (x_new(15) .ge. real(1.d0)) 
!!$
!!$                      plausibility(16) =  (x_new(16) .le. real(-1.d0)) .or. (x_new(16) .ge. real(1.d0)) 
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC+MW+NGC4258 AS ANCHORS'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC+MW+NGC4258 AS ANCHORS'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      plausibility(1) = (x_new(1) .le. real(20.d0)) .or. (x_new(1) .ge. real(4.d1))
!!$                      
!!$                      plausibility(2) = (x_new(2) .le. real(20.d0)) .or. (x_new(2) .ge. real(4.d1))
!!$
!!$                      plausibility(3) = (x_new(3) .le. real(20.d0)) .or. (x_new(3) .ge. real(4.d1))
!!$
!!$                      plausibility(4) = (x_new(4) .le. real(20.d0)) .or. (x_new(4) .ge. real(4.d1))
!!$
!!$                      plausibility(5) = (x_new(5) .le. real(20.d0)) .or. (x_new(5) .ge. real(4.d1))
!!$
!!$                      plausibility(6) = (x_new(6) .le. real(20.d0)) .or. (x_new(6) .ge. real(4.d1))
!!$
!!$                      plausibility(7) = (x_new(7) .le. real(20.d0)) .or. (x_new(7) .ge. real(4.d1))
!!$
!!$                      plausibility(8) = (x_new(8) .le. real(20.d0)) .or. (x_new(8) .ge. real(4.d1))
!!$
!!$                      plausibility(9) = (x_new(9) .le. real(20.d0)) .or. (x_new(9) .ge. real(30.d0))
!!$
!!$                      plausibility(10) = (x_new(10) .le. real(0.d0)) .or. (x_new(10) .ge. real(40.d0))
!!$
!!$                      plausibility(11) =  (x_new(11) .le. real(-7.d0)) .or. (x_new(11) .ge. real(-4.d0)) 
!!$
!!$                      plausibility(12) =  (x_new(12) .le. real(-3.5d0)) .or. (x_new(12) .ge. real(-2.5d0)) 
!!$
!!$                      plausibility(13) =  (x_new(13) .le. real(55.d0)) .or. (x_new(13) .ge. real(95.d0)) 
!!$
!!$                      plausibility(14) =  (x_new(14) .le. real(-2.d0)) .or. (x_new(14) .ge. real(1.d0)) 
!!$
!!$                      plausibility(15) =  (x_new(15) .le. real(0.d0)) .or. (x_new(15) .ge. real(1.d0)) 
!!$
!!$                      plausibility(16) =  (x_new(16) .le. real(-1.d0)) .or. (x_new(16) .ge. real(1.d0)) 
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AND NGC4258 AS ANCHORS'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AND NGC4258 AS ANCHORS'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      plausibility(1) = (x_new(1) .le. real(20.d0)) .or. (x_new(1) .ge. real(4.d1))
!!$                      
!!$                      plausibility(2) = (x_new(2) .le. real(20.d0)) .or. (x_new(2) .ge. real(4.d1))
!!$
!!$                      plausibility(3) = (x_new(3) .le. real(20.d0)) .or. (x_new(3) .ge. real(4.d1))
!!$
!!$                      plausibility(4) = (x_new(4) .le. real(20.d0)) .or. (x_new(4) .ge. real(4.d1))
!!$
!!$                      plausibility(5) = (x_new(5) .le. real(20.d0)) .or. (x_new(5) .ge. real(4.d1))
!!$
!!$                      plausibility(6) = (x_new(6) .le. real(20.d0)) .or. (x_new(6) .ge. real(4.d1))
!!$
!!$                      plausibility(7) = (x_new(7) .le. real(20.d0)) .or. (x_new(7) .ge. real(4.d1))
!!$
!!$                      plausibility(8) = (x_new(8) .le. real(20.d0)) .or. (x_new(8) .ge. real(4.d1))
!!$
!!$                      plausibility(9) = (x_new(9) .le. real(20.d0)) .or. (x_new(9) .ge. real(30.d0))
!!$
!!$                      plausibility(10) = (x_new(10) .le. real(-7.d0)) .or. (x_new(10) .ge. real(-4.d0))
!!$
!!$                      plausibility(11) =  (x_new(11) .le. real(-3.5d0)) .or. (x_new(11) .ge. real(-2.5d0)) 
!!$
!!$                      plausibility(12) =  (x_new(12) .le. real(55.d0)) .or. (x_new(12) .ge. real(95.d0)) 
!!$
!!$                      plausibility(13) =  (x_new(13) .le. real(-2.d0)) .or. (x_new(13) .ge. real(1.d0)) 
!!$
!!$                      plausibility(14) =  (x_new(14) .le. real(0.d0)) .or. (x_new(14) .ge. real(1.d0)) 
!!$
!!$                      plausibility(15) =  (x_new(15) .le. real(-1.d0)) .or. (x_new(15) .ge. real(1.d0)) 
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN MW+NGC4258 AS ANCHOR'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN MW+NGC4258 AS ANCHOR'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      plausibility(1) = (x_new(1) .le. real(20.d0)) .or. (x_new(1) .ge. real(4.d1))
!!$                      
!!$                      plausibility(2) = (x_new(2) .le. real(20.d0)) .or. (x_new(2) .ge. real(4.d1))
!!$
!!$                      plausibility(3) = (x_new(3) .le. real(20.d0)) .or. (x_new(3) .ge. real(4.d1))
!!$
!!$                      plausibility(4) = (x_new(4) .le. real(20.d0)) .or. (x_new(4) .ge. real(4.d1))
!!$
!!$                      plausibility(5) = (x_new(5) .le. real(20.d0)) .or. (x_new(5) .ge. real(4.d1))
!!$
!!$                      plausibility(6) = (x_new(6) .le. real(20.d0)) .or. (x_new(6) .ge. real(4.d1))
!!$
!!$                      plausibility(7) = (x_new(7) .le. real(20.d0)) .or. (x_new(7) .ge. real(4.d1))
!!$
!!$                      plausibility(8) = (x_new(8) .le. real(20.d0)) .or. (x_new(8) .ge. real(4.d1))
!!$
!!$                      plausibility(9) = (x_new(9) .le. real(20.d0)) .or. (x_new(9) .ge. real(30.d0))
!!$
!!$                      plausibility(10) = (x_new(10) .le. real(0.d0)) .or. (x_new(10) .ge. real(40.d0))
!!$
!!$                      plausibility(11) =  (x_new(11) .le. real(-7.d0)) .or. (x_new(11) .ge. real(-4.d0)) 
!!$
!!$                      plausibility(12) =  (x_new(12) .le. real(-3.5d0)) .or. (x_new(12) .ge. real(-2.5d0)) 
!!$
!!$                      plausibility(13) =  (x_new(13) .le. real(55.d0)) .or. (x_new(13) .ge. real(95.d0)) 
!!$
!!$                      plausibility(14) =  (x_new(14) .le. real(-2.d0)) .or. (x_new(14) .ge. real(1.d0)) 
!!$
!!$                      plausibility(15) =  (x_new(15) .le. real(0.d0)) .or. (x_new(15) .ge. real(1.d0)) 
!!$
!!$                      plausibility(16) =  (x_new(16) .le. real(-1.d0)) .or. (x_new(16) .ge. real(1.d0)) 
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AND MW AS ANCHORS'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AND MW AS ANCHORS'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( .not.use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      plausibility(1) = (x_new(1) .le. real(20.d0)) .or. (x_new(1) .ge. real(4.d1))
!!$                      
!!$                      plausibility(2) = (x_new(2) .le. real(20.d0)) .or. (x_new(2) .ge. real(4.d1))
!!$
!!$                      plausibility(3) = (x_new(3) .le. real(20.d0)) .or. (x_new(3) .ge. real(4.d1))
!!$
!!$                      plausibility(4) = (x_new(4) .le. real(20.d0)) .or. (x_new(4) .ge. real(4.d1))
!!$
!!$                      plausibility(5) = (x_new(5) .le. real(20.d0)) .or. (x_new(5) .ge. real(4.d1))
!!$
!!$                      plausibility(6) = (x_new(6) .le. real(20.d0)) .or. (x_new(6) .ge. real(4.d1))
!!$
!!$                      plausibility(7) = (x_new(7) .le. real(20.d0)) .or. (x_new(7) .ge. real(4.d1))
!!$
!!$                      plausibility(8) = (x_new(8) .le. real(20.d0)) .or. (x_new(8) .ge. real(4.d1))
!!$
!!$                      plausibility(9) = (x_new(9) .le. real(20.d0)) .or. (x_new(9) .ge. real(30.d0))
!!$
!!$                      plausibility(10) = (x_new(10) .le. real(-7.d0)) .or. (x_new(10) .ge. real(-4.d0))
!!$
!!$                      plausibility(11) =  (x_new(11) .le. real(-3.5d0)) .or. (x_new(11) .ge. real(-2.5d0)) 
!!$
!!$                      plausibility(12) =  (x_new(12) .le. real(55.d0)) .or. (x_new(12) .ge. real(95.d0)) 
!!$
!!$                      plausibility(13) =  (x_new(13) .le. real(-2.d0)) .or. (x_new(13) .ge. real(1.d0)) 
!!$
!!$                      plausibility(14) =  (x_new(14) .le. real(0.d0)) .or. (x_new(14) .ge. real(1.d0)) 
!!$
!!$                      plausibility(15) =  (x_new(15) .le. real(-1.d0)) .or. (x_new(15) .ge. real(1.d0)) 
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN MW AS ANCHOR'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN MW AS ANCHOR'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      plausibility(1) = (x_new(1) .le. real(20.d0)) .or. (x_new(1) .ge. real(4.d1))
!!$                      
!!$                      plausibility(2) = (x_new(2) .le. real(20.d0)) .or. (x_new(2) .ge. real(4.d1))
!!$
!!$                      plausibility(3) = (x_new(3) .le. real(20.d0)) .or. (x_new(3) .ge. real(4.d1))
!!$
!!$                      plausibility(4) = (x_new(4) .le. real(20.d0)) .or. (x_new(4) .ge. real(4.d1))
!!$
!!$                      plausibility(5) = (x_new(5) .le. real(20.d0)) .or. (x_new(5) .ge. real(4.d1))
!!$
!!$                      plausibility(6) = (x_new(6) .le. real(20.d0)) .or. (x_new(6) .ge. real(4.d1))
!!$
!!$                      plausibility(7) = (x_new(7) .le. real(20.d0)) .or. (x_new(7) .ge. real(4.d1))
!!$
!!$                      plausibility(8) = (x_new(8) .le. real(20.d0)) .or. (x_new(8) .ge. real(4.d1))
!!$
!!$                      plausibility(9) = (x_new(9) .le. real(20.d0)) .or. (x_new(9) .ge. real(30.d0))
!!$
!!$                      plausibility(10) = (x_new(10) .le. real(0.d0)) .or. (x_new(10) .ge. real(40.d0))
!!$
!!$                      plausibility(11) =  (x_new(11) .le. real(15.d0)) .or. (x_new(11) .ge. real(25.d0)) 
!!$
!!$                      plausibility(12) =  (x_new(12) .le. real(-3.5d0)) .or. (x_new(12) .ge. real(-2.5d0)) 
!!$
!!$                      plausibility(13) =  (x_new(13) .le. real(55.d0)) .or. (x_new(13) .ge. real(95.d0)) 
!!$
!!$                      plausibility(14) =  (x_new(14) .le. real(-2.d0)) .or. (x_new(14) .ge. real(1.d0)) 
!!$
!!$                      plausibility(15) =  (x_new(15) .le. real(0.d0)) .or. (x_new(15) .ge. real(1.d0)) 
!!$
!!$                      plausibility(16) =  (x_new(16) .le. real(-1.d0)) .or. (x_new(16) .ge. real(1.d0)) 
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AS ANCHOR'
!!$                     
!!$                      stop
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AS ANCHOR'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                If (use_metallicity) then 
!!$
!!$                   If (use_H_band) then
!!$
!!$                      plausibility(1) = (x_new(1) .le. real(20.d0)) .or. (x_new(1) .ge. real(4.d1))
!!$
!!$                      plausibility(2) = (x_new(2) .le. real(20.d0)) .or. (x_new(2) .ge. real(4.d1))
!!$
!!$                      plausibility(3) = (x_new(3) .le. real(20.d0)) .or. (x_new(3) .ge. real(4.d1))
!!$                           
!!$                      plausibility(4) = (x_new(4) .le. real(20.d0)) .or. (x_new(4) .ge. real(4.d1))
!!$
!!$                      plausibility(5) = (x_new(5) .le. real(20.d0)) .or. (x_new(5) .ge. real(4.d1))
!!$
!!$                      plausibility(6) = (x_new(6) .le. real(20.d0)) .or. (x_new(6) .ge. real(4.d1))
!!$
!!$                      plausibility(7) = (x_new(7) .le. real(20.d0)) .or. (x_new(7) .ge. real(4.d1))
!!$
!!$                      plausibility(8) = (x_new(8) .le. real(20.d0)) .or. (x_new(8) .ge. real(4.d1))
!!$
!!$                      plausibility(9) = (x_new(9) .le. real(20.d0)) .or. (x_new(9) .ge. real(30.d0))
!!$
!!$                      plausibility(10) =  (x_new(10) .le. real(25.d0)) .or. (x_new(10) .ge. real(34.d0)) 
!!$
!!$                      plausibility(11) =  (x_new(11) .le. real(-3.2d0)) .or. (x_new(11) .ge. real(-2.5d0)) 
!!$
!!$                      plausibility(12) =  (x_new(12) .le. real(55.d0)) .or. (x_new(12) .ge. real(95.d0)) 
!!$
!!$                      plausibility(13) =  (x_new(13) .le. real(-1.d0)) .or. (x_new(13) .ge. real(1.d0)) 
!!$
!!$                      plausibility(14) =  (x_new(14) .le. real(0.d0)) .or. (x_new(14) .ge. real(1.d0)) 
!!$
!!$                   Else
!!$
!!$                      plausibility(1) = (x_new(1) .le. real(20.d0)) .or. (x_new(1) .ge. real(4.d1))
!!$
!!$                      plausibility(2) = (x_new(2) .le. real(20.d0)) .or. (x_new(2) .ge. real(4.d1))
!!$
!!$                      plausibility(3) = (x_new(3) .le. real(20.d0)) .or. (x_new(3) .ge. real(4.d1))
!!$                           
!!$                      plausibility(4) = (x_new(4) .le. real(20.d0)) .or. (x_new(4) .ge. real(4.d1))
!!$
!!$                      plausibility(5) = (x_new(5) .le. real(20.d0)) .or. (x_new(5) .ge. real(4.d1))
!!$
!!$                      plausibility(6) = (x_new(6) .le. real(20.d0)) .or. (x_new(6) .ge. real(4.d1))
!!$
!!$                      plausibility(7) = (x_new(7) .le. real(20.d0)) .or. (x_new(7) .ge. real(4.d1))
!!$
!!$                      plausibility(8) = (x_new(8) .le. real(20.d0)) .or. (x_new(8) .ge. real(4.d1))
!!$
!!$                      plausibility(9) = (x_new(9) .le. real(20.d0)) .or. (x_new(9) .ge. real(30.d0))
!!$
!!$                      plausibility(10) =  (x_new(10) .le. real(25.d0)) .or. (x_new(10) .ge. real(34.d0)) 
!!$
!!$                      plausibility(11) =  (x_new(11) .le. real(-3.5d0)) .or. (x_new(11) .ge. real(-2.5d0)) 
!!$
!!$                      plausibility(12) =  (x_new(12) .le. real(55.d0)) .or. (x_new(12) .ge. real(95.d0)) 
!!$
!!$                      plausibility(13) =  (x_new(13) .le. real(-1.d0)) .or. (x_new(13) .ge. real(1.d0)) 
!!$
!!$                      plausibility(14) =  (x_new(14) .le. real(0.d0)) .or. (x_new(14) .ge. real(1.d0)) 
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   If (use_H_band) then
!!$
!!$                      plausibility(1) = (x_new(1) .le. real(0.d0)) .or. (x_new(1) .ge. real(5.d1))
!!$
!!$                      plausibility(2) = (x_new(2) .le. real(0.d0)) .or. (x_new(2) .ge. real(5.d1))
!!$
!!$                      plausibility(3) = (x_new(3) .le. real(0.d0)) .or. (x_new(3) .ge. real(5.d1))
!!$
!!$                      plausibility(4) = (x_new(4) .le. real(0.d0)) .or. (x_new(4) .ge. real(5.d1))
!!$
!!$                      plausibility(5) = (x_new(5) .le. real(0.d0)) .or. (x_new(5) .ge. real(5.d1))
!!$
!!$                      plausibility(6) = (x_new(6) .le. real(0.d0)) .or. (x_new(6) .ge. real(5.d1))
!!$
!!$                      plausibility(7) = (x_new(7) .le. real(0.d0)) .or. (x_new(7) .ge. real(5.d1))
!!$
!!$                      plausibility(8) = (x_new(8) .le. real(0.d0)) .or. (x_new(8) .ge. real(5.d1))
!!$
!!$                      plausibility(9) = (x_new(9) .le. real(0.d0)) .or. (x_new(9) .ge. real(50.d0))
!!$
!!$                      plausibility(10) =  (x_new(10) .le. real(-20.d0)) .or. (x_new(10) .ge. real(0.d0)) 
!!$
!!$                   Else
!!$
!!$                      print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
!!$                     
!!$                      stop
!!$                        
!!$                   End If
!!$
!!$                End If
!!$
!!$             Else
!!$
!!$                print *, 'USER MUST SET TRUE AT LEAST ONE ANCHOR DISTANCE IN FIDUCIAL MODULE'
!!$
!!$                stop
!!$
!!$             End If ! OF ANCHOR
!!$
!!$          End If ! OF CEPHEIDS
!!$    
!!$       Else
!!$
!!$          plausibility(1) = (x_new(1) .le. real(0.d0)) .or. (x_new(1) .ge. real(5.d1))
!!$          plausibility(2) = (x_new(2) .le. real(-2.d1)) .or. (x_new(2) .ge. real(0.d0))
!!$          !plausibility(3) =  (x_new(3) .gt. real(0.d0)) .or. (x_new(3) .lt. real(-10.d0))    ! limit log10(sigma_int)
!!$
!!$       End If
!!$
!!$       If (hyperparameters_as_mcmc) then
!!$          ! CHECKING PLAUSIBILITY FOR HYPER-PARAMETERS
!!$          Do n=number_model_parameters+1,number_of_parameters
!!$                   
!!$             plausibility(n) =  (x_new(n) .gt. real(1.d0)) .or. (x_new(n) .lt. real(0.d0))    ! limit alpha_j
!!$
!!$          End Do
!!$            
!!$       End If
!!$
!!$       Do n=1,number_of_parameters
!!$
!!$          If (plausibility(n)) then
!!$
!!$             non_plausible_parameters = .true.
!!$
!!$             exit
!!$
!!$          Else 
!!$
!!$             non_plausible_parameters = .false.
!!$
!!$          End If
!!$
!!$       End Do
!!$
!!$       ! NEW POINT GENERATED 
!!$
!!$       Do n=1,number_of_parameters
!!$
!!$          If (n .gt. number_model_parameters) then
!!$
!!$             If (using_jeffreys_prior) then
!!$
!!$                current_point(n) = 10**(dble(x_new(n))) ! CONVERTING LOG10(alpha_j) to alpha_j 
!!$
!!$             Else
!!$
!!$                current_point(n) = dble(x_new(n))
!!$              
!!$             End If
!!$
!!$          Else
!!$
!!$             current_point(n) = dble(x_new(n))
!!$
!!$          End If
!!$
!!$       End Do
!!$
!!$       ! EVALUATE LOG_LIKELIHOOD FOR CURRENT POINT IN PARAMETER SPACE
!!$       If (testing_Gaussian_likelihood) then
!!$
!!$          current_loglikelihood = log_Gaussian_likelihood(current_point)
!!$
!!$       Else
!!$
!!$          If (non_plausible_parameters) then
!!$
!!$             current_loglikelihood = -1.d10
!!$
!!$          Else
!!$
!!$             If (using_hyperparameters) then    
!!$
!!$                If (doing_R11_analysis) then
!!$
!!$                   If (include_only_cepheids) then
!!$
!!$                      If (all_R11_hosts) then
!!$                        
!!$                         current_loglikelihood = log_R11_likelihood_W_cepheids(current_point(1:number_model_parameters-3),&
!!$                              current_point(number_model_parameters-2),current_point(number_model_parameters-1),&
!!$                              current_point(number_model_parameters),prior_sigma_int)
!!$
!!$                      Else
!!$
!!$                         current_loglikelihood = log_likelihood_only_cepheids(galaxy,current_point(1),&
!!$                              current_point(2),current_point(3),prior_sigma_int)
!!$
!!$                      End If
!!$
!!$                   Else
!!$
!!$                      If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor) then !!!HERE
!!$    
!!$                         If (use_metallicity) then 
!!$
!!$                            If (use_H_band) then
!!$                     
!!$                               print *,'H BAND NOT IMPLEMENTED USING LMC+MW+NGC4258 AS ANCHOR AND INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                               stop
!!$                     
!!$                            Else
!!$                          
!!$                               current_loglikelihood = log_R11_likelihood_W_LMC_MW_NGC4258(current_point(1:&
!!$                                    number_model_parameters-6),&
!!$                                    current_point(number_model_parameters-5),current_point(number_model_parameters-4),&
!!$                                    current_point(number_model_parameters-3),current_point(number_model_parameters-2),&
!!$                                    current_point(number_model_parameters-1),current_point(number_model_parameters),&
!!$                                    prior_sigma_int,prior_sigma_int_LMC)
!!$
!!$                            End If
!!$
!!$                         Else
!!$
!!$                            If (use_H_band) then
!!$                               
!!$                               print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC+MW+NGC4258 AS ANCHORS'
!!$                     
!!$                               stop
!!$
!!$                            Else
!!$
!!$                               print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC+MW+NGC4258 AS ANCHORS'
!!$                     
!!$                               stop
!!$
!!$                            End If
!!$
!!$                         End If
!!$
!!$                      Else If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                         If (use_metallicity) then 
!!$
!!$                            If (use_H_band) then
!!$                     
!!$                               print *,'H BAND NOT IMPLEMENTED USING LMC AND NGC4258 AS ANCHOR AND INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                               stop
!!$                     
!!$                            Else
!!$                          
!!$                               current_loglikelihood = log_R11_likelihood_W_LMC_NGC4258(current_point(1:number_model_parameters-6),&
!!$                                    current_point(number_model_parameters-5),current_point(number_model_parameters-4),&
!!$                                    current_point(number_model_parameters-3),current_point(number_model_parameters-2),&
!!$                                    current_point(number_model_parameters-1),current_point(number_model_parameters),&
!!$                                    prior_sigma_int,prior_sigma_int_LMC)
!!$
!!$                            End If
!!$
!!$                         Else
!!$
!!$                            If (use_H_band) then
!!$                               
!!$                               print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AND NGC4258 AS ANCHORS'
!!$                     
!!$                               stop
!!$
!!$                            Else
!!$
!!$                               print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AND NGC4258 AS ANCHORS'
!!$                     
!!$                               stop
!!$
!!$                            End If
!!$
!!$                         End If
!!$
!!$                      Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                         If (use_metallicity) then 
!!$
!!$                            If (use_H_band) then
!!$                     
!!$                               print *,'H BAND NOT IMPLEMENTED USING MW AS ANCHOR AND INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                               stop
!!$                     
!!$                            Else
!!$                          
!!$                               current_loglikelihood = log_R11_likelihood_W_MW_NGC4258(current_point(1:number_model_parameters-6),&
!!$                                    current_point(number_model_parameters-5),current_point(number_model_parameters-4),&
!!$                                    current_point(number_model_parameters-3),current_point(number_model_parameters-2),&
!!$                                    current_point(number_model_parameters-1),current_point(number_model_parameters),&
!!$                                    prior_sigma_int,prior_sigma_int_MW)
!!$
!!$                            End If
!!$
!!$                         Else
!!$
!!$                            If (use_H_band) then
!!$                               
!!$                               print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR MW+NGC4258 AS ANCHOR'
!!$                     
!!$                               stop
!!$
!!$                            Else
!!$
!!$                               print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR MW+NGC4258 AS ANCHOR'
!!$                     
!!$                               stop
!!$
!!$                            End If
!!$
!!$                         End If
!!$
!!$                      Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                         If (use_metallicity) then 
!!$
!!$                            If (use_H_band) then
!!$                     
!!$                               print *,'H BAND NOT IMPLEMENTED USING LMC AND MW AS ANCHOR AND INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                               stop
!!$                     
!!$                            Else
!!$                          
!!$                               current_loglikelihood = log_R11_likelihood_W_LMC_MW(current_point(1:number_model_parameters-6),&
!!$                                    current_point(number_model_parameters-5),current_point(number_model_parameters-4),&
!!$                                    current_point(number_model_parameters-3),current_point(number_model_parameters-2),&
!!$                                    current_point(number_model_parameters-1),current_point(number_model_parameters),&
!!$                                    prior_sigma_int,prior_sigma_int_LMC)
!!$
!!$                            End If
!!$
!!$                         Else
!!$
!!$                            If (use_H_band) then
!!$                               
!!$                               print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AND MW AS ANCHORS'
!!$                     
!!$                               stop
!!$
!!$                            Else
!!$
!!$                               print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AND MW AS ANCHORS'
!!$                     
!!$                               stop
!!$
!!$                            End If
!!$
!!$                         End If
!!$
!!$                      Else If ( ( .not.use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$                         If (use_metallicity) then 
!!$
!!$                            If (use_H_band) then
!!$                     
!!$                               print *,'H BAND NOT IMPLEMENTED USING MW AS ANCHOR AND INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                               stop
!!$                     
!!$                            Else
!!$                          
!!$                               current_loglikelihood = log_R11_likelihood_W_MW(current_point(1:number_model_parameters-6),&
!!$                                    current_point(number_model_parameters-5),current_point(number_model_parameters-4),&
!!$                                    current_point(number_model_parameters-3),current_point(number_model_parameters-2),&
!!$                                    current_point(number_model_parameters-1),current_point(number_model_parameters),&
!!$                                    prior_sigma_int,prior_sigma_int_MW)
!!$
!!$                            End If
!!$
!!$                         Else
!!$
!!$                            If (use_H_band) then
!!$                               
!!$                               print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR MW AS ANCHOR'
!!$                     
!!$                               stop
!!$
!!$                            Else
!!$
!!$                               print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR MW AS ANCHOR'
!!$                     
!!$                               stop
!!$
!!$                            End If
!!$
!!$                         End If
!!$
!!$                      Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                         If (use_metallicity) then 
!!$
!!$                            If (use_H_band) then
!!$                     
!!$                               print *,'H BAND NOT IMPLEMENTED USING LMC AS ANCHOR AND INCLUDING METALLICITY DEPENDENCE'
!!$                     
!!$                               stop
!!$                     
!!$                            Else
!!$                          
!!$                               current_loglikelihood = log_R11_likelihood_W_LMC(current_point(1:number_model_parameters-6),&
!!$                                    current_point(number_model_parameters-5),current_point(number_model_parameters-4),&
!!$                                    current_point(number_model_parameters-3),current_point(number_model_parameters-2),&
!!$                                    current_point(number_model_parameters-1),current_point(number_model_parameters),&
!!$                                    prior_sigma_int,prior_sigma_int_LMC)
!!$
!!$                            End If
!!$
!!$                         Else
!!$
!!$                            If (use_H_band) then
!!$                               
!!$                               print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AS ANCHOR'
!!$                     
!!$                               stop
!!$
!!$                            Else
!!$
!!$                               print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE FOR LMC AS ANCHOR'
!!$                     
!!$                               stop
!!$
!!$                            End If
!!$
!!$                         End If
!!$
!!$                      Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$                         If (use_metallicity) then 
!!$
!!$                            If (use_H_band) then
!!$                     
!!$                               current_loglikelihood = log_R11_likelihood_H_NGC4258(current_point(1:number_model_parameters-5),&
!!$                                    current_point(number_model_parameters-4),current_point(number_model_parameters-3),&
!!$                                    current_point(number_model_parameters-2),current_point(number_model_parameters-1),&
!!$                                    current_point(number_model_parameters),prior_sigma_int)
!!$
!!$                            Else
!!$
!!$                               current_loglikelihood = log_R11_likelihood_W(current_point(1:number_model_parameters-5),&
!!$                                    current_point(number_model_parameters-4),current_point(number_model_parameters-3),&
!!$                                    current_point(number_model_parameters-2),current_point(number_model_parameters-1),&
!!$                                    current_point(number_model_parameters),prior_sigma_int)
!!$
!!$                            End If
!!$
!!$                         Else
!!$
!!$                            If (use_H_band) then
!!$
!!$                               current_loglikelihood = log_R11_likelihood_H(current_point(1:number_model_parameters-1),&
!!$                                    current_point(number_model_parameters),prior_sigma_int)
!!$
!!$                            Else
!!$
!!$                               print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
!!$                     
!!$                               stop
!!$
!!$                            End If
!!$
!!$                         End If
!!$
!!$                      End If ! OF ANCHORS 
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$!                   current_loglikelihood = log_Efstathiou_likelihood_hyperparameters(current_point(1),&
!!$ !                       current_point(2),prior_sigma_int)
!!$
!!$                   current_loglikelihood = log_Efstathiou_likelihood(current_point(1),current_point(2),prior_sigma_int_LMC)
!!$              
!!$                End If
!!$
!!$             Else
!!$
!!$                current_loglikelihood = log_Efstathiou_likelihood(current_point(1),current_point(2),prior_sigma_int_LMC)
!!$
!!$             End If
!!$
!!$          End If
!!$
!!$       End If
!!$       ! LOG_LIKELIHOOD FOR CURRENT POINT COMPUTED
!!$
!!$       !MAKE DECISION ABOUT CURRENT POINT : ACCEPT OR REJECT IT
!!$       If (current_loglikelihood .ge. old_loglikelihood) then ! ACCEPT CURRENT POINT
!!$
!!$          If (m .gt. steps_taken_before_definite_run) then
!!$
!!$             number_accepted_points = number_accepted_points + 1         
!!$
!!$          End If
!!$
!!$          ! COMPUTING ACCEPTANCE PROBABILITY FOR CURRENT POINT
!!$          acceptance_probability(m) = min(1.d0,exp(current_loglikelihood - old_loglikelihood))    
!!$        
!!$          If (m .le. steps_taken_before_definite_run) then ! WRITE OUT INFORMATION IN TEMPORARY FILE
!!$               
!!$             write(UNIT_MCMC_FILE,*) weight,-old_loglikelihood,old_point(1:number_of_parameters)
!!$
!!$          Else ! WRITE OUT INFORMATION IN DEFINITE FILE
!!$
!!$             write(UNIT_MCMC_FINAL_FILE,*) weight,-old_loglikelihood,old_point(1:number_of_parameters)
!!$
!!$          End If
!!$       
!!$          weight = 1    
!!$
!!$          old_loglikelihood = current_loglikelihood
!!$        
!!$          Do i=1,number_of_parameters 
!!$
!!$             old_point(i) = current_point(i)
!!$
!!$             If (i .gt. number_model_parameters) then
!!$
!!$                If (using_jeffreys_prior) then
!!$
!!$                   x_old(i) = log10(real(old_point(i))) ! CONVERTING alpha_j TO  log10(alpha_j)
!!$
!!$                Else
!!$
!!$                   x_old(i) = real(old_point(i))
!!$
!!$                End If
!!$
!!$             Else
!!$
!!$                x_old(i) = real(old_point(i))
!!$
!!$             End If
!!$
!!$          End Do
!!$   
!!$       Else ! ACCEPT OR REJECT THE CURRENT POINT ACCORDING TO :
!!$
!!$          random_uniform = dble(genunf(real(0.),real(1.)))
!!$
!!$          If ( random_uniform .le. exp(current_loglikelihood-old_loglikelihood)) then ! ACCEPT CURRENT POINT
!!$
!!$             If (m .gt. steps_taken_before_definite_run) then
!!$
!!$                number_accepted_points = number_accepted_points + 1         
!!$
!!$             End If
!!$                
!!$             acceptance_probability(m) = min(1.d0,dexp(current_loglikelihood - old_loglikelihood))    
!!$
!!$             If (m .le. steps_taken_before_definite_run) then ! WRITE OUT INFORMATION TO TEMPORARY FILE
!!$
!!$                write(UNIT_MCMC_FILE,*) weight,-old_loglikelihood,old_point(1:number_of_parameters)
!!$
!!$             Else ! WRITE OUT INFORMATION TO DEFINITE FILE
!!$                   
!!$                write(UNIT_MCMC_FINAL_FILE,*) weight,-old_loglikelihood,old_point(1:number_of_parameters)
!!$
!!$             End If
!!$
!!$             weight = 1
!!$
!!$             old_loglikelihood = current_loglikelihood
!!$
!!$             Do i=1,number_of_parameters 
!!$
!!$                old_point(i) = current_point(i)
!!$
!!$                If (i .gt. number_model_parameters) then
!!$                        
!!$                   If (using_jeffreys_prior) then
!!$
!!$                      x_old(i) = real(log10(old_point(i))) ! CONVERTING alpha_j TO log10(alpha_j)
!!$
!!$                   Else
!!$
!!$                      x_old(i) = real(old_point(i))
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   x_old(i) = real(old_point(i))
!!$
!!$                End If
!!$
!!$             End Do
!!$
!!$          Else   ! REJECT CURRENT POINT 
!!$
!!$             If (m .gt. steps_taken_before_definite_run) then
!!$
!!$                number_rejected_points = number_rejected_points + 1            
!!$
!!$             End If
!!$
!!$             acceptance_probability(m) = min(1.d0,exp(current_loglikelihood - old_loglikelihood))    
!!$
!!$             weight = weight + 1
!!$
!!$             Do i=1,number_of_parameters 
!!$
!!$                If (i .gt. number_model_parameters) then
!!$
!!$                   If (using_jeffreys_prior) then
!!$
!!$                      x_old(i) = real(log10(old_point(i))) ! CONVERT alpha_j TO log10(alpha_j)
!!$
!!$                   Else
!!$
!!$                      x_old(i) = real(old_point(i))
!!$
!!$                   End If
!!$
!!$                Else
!!$
!!$                   x_old(i) = real(old_point(i))
!!$
!!$                End If
!!$
!!$             End Do
!!$
!!$          End If
!!$
!!$       End If
!!$       ! DECISION ABOUT CURRENT POINT MADE
!!$
!!$       ! COMPUTING AVERAGE ACCEPTANCE PROBABILITY AND UPDATING BOTH COVARIANCE MATRIX AND JUMPING FACTOR (IF NEEDED)
!!$       If ((mod(m,jumping_factor_update) .eq. 0) .and. (m .le. steps_taken_before_definite_run) ) then
!!$
!!$          average_acceptance_probability = sum(acceptance_probability(m-jumping_factor_update+1:m))&
!!$               /real(jumping_factor_update)
!!$
!!$          !            write(15,*) 'CURRENT AVERAGE ACCEPTANCE PROBABILITY = ',average_acceptance_probability
!!$            
!!$          ! UPDATE JUMPING FACTOR IF NEEDED        
!!$          If (average_acceptance_probability .lt. 0.1) then 
!!$               
!!$             jumping_factor = (1.d0 - step_size_changes)    !    Decreasing step size
!!$
!!$             If (testing_Gaussian_likelihood) then
!!$
!!$                Covgauss = jumping_factor*Covgauss
!!$
!!$             Else
!!$
!!$                Covguess = jumping_factor*Covguess
!!$
!!$             End If
!!$
!!$          Else if (average_acceptance_probability .gt. 0.4) then
!!$
!!$             jumping_factor = (1.d0 + step_size_changes)    !    Increasing step size 
!!$
!!$             If (testing_Gaussian_likelihood) then
!!$
!!$                Covgauss = jumping_factor*Covgauss
!!$
!!$             Else
!!$
!!$                Covguess = jumping_factor*Covguess
!!$
!!$             End If
!!$
!!$          End If
!!$          ! JUMPING FACTOR UPDATED (IF IT WAS NEEDED)
!!$             
!!$          not_good_aap = (average_acceptance_probability .lt. 0.1) .or. (average_acceptance_probability .gt. 0.4)
!!$
!!$          If ( (mod(m,covariance_matrix_update) .eq. 0) .and. not_good_aap) then
!!$                
!!$             call stat('./output/mcmc_output.txt',buff,status1)
!!$
!!$             If ((status1 .eq. 0) .and. (buff(8) .gt. 0)) then
!!$              
!!$                If (testing_Gaussian_likelihood) then
!!$
!!$                   call system('cd output; python compute_covariance_matrix_Gaussian.py')
!!$
!!$                   call read_covariance_matrix_mcmc(Covgauss)
!!$
!!$                   close(UNIT_MCMC_FILE)
!!$
!!$                   call system('rm ./output/mcmc_output.txt')
!!$
!!$                   open(UNIT_MCMC_FILE,file='./output/mcmc_output.txt')
!!$
!!$                Else
!!$
!!$                   If (using_hyperparameters) then
!!$
!!$                      If (hyperparameters_as_mcmc) then
!!$                                
!!$                         call system('cd output; python compute_covariance_matrix_HP.py')
!!$                                
!!$                      Else
!!$
!!$                         call system('cd output; python compute_covariance_matrix.py')
!!$                                
!!$                      End If
!!$
!!$                   Else
!!$                           
!!$                      call system('cd output; python compute_covariance_matrix.py')
!!$
!!$                   End If
!!$
!!$                   call read_covariance_matrix_mcmc(Covguess)
!!$
!!$                   close(UNIT_MCMC_FILE)
!!$
!!$                   call system('rm ./output/mcmc_output.txt')
!!$
!!$                   open(UNIT_MCMC_FILE,file='./output/mcmc_output.txt')
!!$
!!$                End If
!!$
!!$             End If
!!$
!!$          End If
!!$
!!$       End If
!!$
!!$    End Do
!!$    ! LOOP TO EXPLORE PARAMETER SPACE ENDED
!!$
!!$    write(UNIT_EXE_FILE,*) 'NUMBER OF REJECTED POINTS = ', number_rejected_points
!!$
!!$    write(UNIT_EXE_FILE,*) 'ACCEPTANCE RATIO = ', dble(number_iterations - steps_taken_before_definite_run&
!!$    - number_rejected_points)/dble(number_iterations - steps_taken_before_definite_run)
!!$ 
!!$    ! CLOSE FILE STORING CHAIN
!!$    close(UNIT_MCMC_FINAL_FILE)
!!$    ! CLOSE TEMPORARY FILE FOR CHAINS
!!$    close(UNIT_MCMC_FILE)
!!$
!!$    !ANALYZE SAMPLES, MAKE FIGURES, COMPUTE BESTFIT AND HYPER-PARAMETERS (IF NEEDED)
!!$    If (testing_Gaussian_likelihood) then
!!$
!!$        call system('cd analyzer; python analyze.py')
!!$
!!$    Else
!!$
!!$        If (using_hyperparameters) then
!!$
!!$            If (hyperparameters_as_mcmc) then
!!$
!!$                call system('cd analyzer; python analyze_HP_as_MCMC.py')
!!$
!!$            Else
!!$               
!!$               If (doing_R11_analysis) then  !MUST IMPLEMENT OTHER OPTIONS LATER!!!!!!!!!!!!!!!!!
!!$
!!$                  If (include_only_cepheids) then
!!$
!!$                     call system('cd analyzer; python analyze_HP.py')
!!$
!!$                  Else
!!$
!!$                     If (use_H_band) then
!!$
!!$                        call system('cd analyzer; python analyze_HP_R11_H.py')
!!$
!!$                     Else
!!$
!!$                        call system('cd analyzer; python analyze_HP_R11_W.py')
!!$
!!$                     End If
!!$
!!$                  End If
!!$
!!$               Else
!!$
!!$                  call system('cd analyzer; python analyze_HP.py')
!!$
!!$               End If
!!$
!!$            End If
!!$
!!$        Else
!!$
!!$            call system('cd analyzer; python analyze.py')
!!$
!!$        End If    
!!$
!!$    End If
!!$    
!!$    call read_bestfit_mcmc(bestfit)
!!$
!!$    call read_means_mcmc(means)
!!$
!!$    If (doing_R11_analysis) then
!!$
!!$       If (include_only_cepheids) then 
!!$
!!$          write(UNIT_EXE_FILE,*) 'BESTFIT IS : '
!!$
!!$          Do m=1,number_model_parameters
!!$
!!$             write(UNIT_EXE_FILE,*) ''//trim(paramnames(m))//' = ', bestfit(m)
!!$
!!$          End Do
!!$
!!$       Else
!!$
!!$          If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor) then
!!$    
!!$             write(UNIT_EXE_FILE,*) 'BESTFIT IS : '
!!$
!!$             Do m=1,number_model_parameters
!!$
!!$                write(UNIT_EXE_FILE,*) ''//trim(paramnames(m))//' = ', bestfit(m)
!!$
!!$             End Do
!!$
!!$             chi2SNIabestfit = 0.d0
!!$
!!$             Do n=1,number_of_hosts_galaxies-1 ! SN Ia
!!$
!!$                If (use_HP_in_SNIa) then
!!$                   
!!$                   write(UNIT_EXE_FILE,*) 'NEED TO IMPLEMENT CHI2 VALUE WHEN HPs IN SNIa'
!!$
!!$                   stop
!!$!                   chi2SNIabestfit = log(new_chi2(chi2R11_SNIa(bestfit(n),bestfit(13),bestfit(15),n))) + &
!!$ !                       log(N_tilde_R11_SNIa(n)) + chi2SNIabestfit
!!$
!!$                Else
!!$
!!$                   chi2SNIabestfit = chi2R11_SNIa(bestfit(n),bestfit(13),bestfit(15),n) + chi2SNIabestfit
!!$
!!$                End If
!!$              
!!$             End Do
!!$
!!$             write(UNIT_EXE_FILE,*) 'CHI2 CONTRIBUTION OF SNIa IS', chi2SNIabestfit
!!$
!!$             print *,'USE OF THREE ANCHORS SIMULTANEOUSLY NOT IMPLEMENTED YET'
!!$
!!$             stop
!!$
!!$          Else If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$             write(UNIT_EXE_FILE,*) 'BESTFIT IS : '
!!$
!!$             Do m=1,number_model_parameters
!!$
!!$                write(UNIT_EXE_FILE,*) ''//trim(paramnames(m))//' = ', bestfit(m)
!!$
!!$             End Do
!!$
!!$             chi2SNIabestfit = 0.d0
!!$
!!$             Do n=1,number_of_hosts_galaxies-1 ! SN Ia
!!$
!!$                If (use_HP_in_SNIa) then
!!$                   
!!$                   write(UNIT_EXE_FILE,*) 'NEED TO IMPLEMENT CHI2 VALUE WHEN HPs IN SNIa'
!!$
!!$                   stop
!!$!                   chi2SNIabestfit = log(new_chi2(chi2R11_SNIa(bestfit(n),bestfit(13),bestfit(15),n))) + &
!!$ !                       log(N_tilde_R11_SNIa(n)) + chi2SNIabestfit
!!$
!!$                Else
!!$
!!$                   chi2SNIabestfit = chi2R11_SNIa(bestfit(n),bestfit(13),bestfit(15),n) + chi2SNIabestfit
!!$
!!$                End If
!!$              
!!$             End Do
!!$
!!$             write(UNIT_EXE_FILE,*) 'CHI2 CONTRIBUTION OF SNIa IS', chi2SNIabestfit
!!$
!!$             print *,'NGC4258+LMC NOT IMPLEMENTED YET'
!!$
!!$             stop
!!$
!!$          Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$             print *,'NGC4258+MW NOT IMPLEMENTED YET'
!!$
!!$             stop
!!$          
!!$          Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$             print *,'MW+LMC NOT IMPLEMENTED YET'
!!$
!!$             stop
!!$
!!$          Else If ( ( .not.use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$             print *,'MW NOT IMPLEMENTED YET'
!!$
!!$             stop
!!$
!!$          Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$
!!$             write(UNIT_EXE_FILE,*) 'BESTFIT IS : '
!!$
!!$             Do m=1,number_model_parameters
!!$
!!$                write(UNIT_EXE_FILE,*) ''//trim(paramnames(m))//' = ', bestfit(m)
!!$
!!$             End Do
!!$
!!$          Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$             If (use_metallicity) then 
!!$
!!$                If (use_H_band) then
!!$
!!$                   write(UNIT_EXE_FILE,*) 'BESTFIT IS : '
!!$
!!$                   Do m=1,number_model_parameters
!!$
!!$                      write(UNIT_EXE_FILE,*) ''//trim(paramnames(m))//' = ', bestfit(m)
!!$
!!$                   End Do
!!$
!!$                Else
!!$
!!$                   write(UNIT_EXE_FILE,*) 'BESTFIT IS : '
!!$
!!$                   Do m=1,number_model_parameters
!!$
!!$                      write(UNIT_EXE_FILE,*) ''//trim(paramnames(m))//' = ', bestfit(m)
!!$
!!$                   End Do
!!$
!!$                End If
!!$
!!$                chi2SNIabestfit = 0.d0
!!$
!!$                Do n=1,number_of_hosts_galaxies-1 ! SN Ia
!!$
!!$                   If (use_HP_in_SNIa) then
!!$                   
!!$                      write(UNIT_EXE_FILE,*) 'NEED TO IMPLEMENT CHI2 VALUE WHEN HPs IN SNIa'
!!$
!!$                      continue
!!$                      !                   chi2SNIabestfit = log(new_chi2(chi2R11_SNIa(bestfit(n),bestfit(13),bestfit(15),n))) + &
!!$                      !                       log(N_tilde_R11_SNIa(n)) + chi2SNIabestfit
!!$
!!$                   Else
!!$
!!$                      chi2SNIabestfit = chi2R11_SNIa(bestfit(n),bestfit(12),bestfit(14),n) + chi2SNIabestfit
!!$
!!$                      print *, chi2R11_SNIa(bestfit(n),bestfit(12),bestfit(14),n)
!!$
!!$                   End If
!!$              
!!$                End Do
!!$
!!$                write(UNIT_EXE_FILE,*) 'CHI2 CONTRIBUTION OF SNIa IS', chi2SNIabestfit
!!$
!!$             Else
!!$
!!$                If (use_H_band) then
!!$
!!$                   write(UNIT_EXE_FILE,*) 'BESTFIT IS : '
!!$
!!$                   write(UNIT_EXE_FILE,*) 'bH = ', bestfit(number_of_parameters)
!!$
!!$                Else
!!$
!!$                   print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE'
!!$                     
!!$                   stop
!!$                        
!!$                End If
!!$
!!$             End If
!!$
!!$          Else
!!$
!!$             print *, 'USER MUST SET TRUE AT LEAST ONE ANCHOR DISTANCE IN FIDUCIAL MODULE'
!!$
!!$             stop
!!$
!!$          End If
!!$    
!!$       End If
!!$    
!!$    Else
!!$
!!$       write(UNIT_EXE_FILE,*) 'BESTFIT IS : '
!!$
!!$       Do m=1,number_model_parameters
!!$
!!$          write(UNIT_EXE_FILE,*) ''//trim(paramnames(m))//' = ', bestfit(m)
!!$
!!$       End Do
!!$
!!$       !write(15,*) 'sigma_int = ', bestfit(3)
!!$
!!$    End If
!!$
!!$    If (hyperparameters_as_mcmc .and. using_hyperparameters) then
!!$    ! WRITING BESTFIT FOR HYPER-PARAMETERS
!!$        Do m=number_model_parameters+1,number_of_parameters
!!$
!!$            write(string,'(i2)') m-number_model_parameters
!!$
!!$            write(UNIT_EXE_FILE,*) 'alpha_'//trim(string)//' = ', bestfit(m)
!!$
!!$        End Do
!!$            
!!$    End If
!!$
!!$    If (doing_R11_analysis) then
!!$
!!$       If (include_only_cepheids) then
!!$
!!$          write(UNIT_EXE_FILE,*) 'MEANS FOR THE SAMPLES ARE : '
!!$
!!$          Do m=1,number_model_parameters
!!$
!!$             write(UNIT_EXE_FILE,*) ''//trim(paramnames(m))//' = ', means(m)
!!$
!!$          End Do
!!$
!!$       Else
!!$
!!$          If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor) then
!!$    
!!$             print *,'USE OF THREE ANCHORS SIMULTANEOUSLY NOT IMPLEMENTED YET'
!!$
!!$             stop
!!$
!!$          Else If ( ( use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$             print *,'NGC4258+LMC NOT IMPLEMENTED YET'
!!$
!!$             stop
!!$
!!$          Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$             print *,'NGC4258+MW NOT IMPLEMENTED YET'
!!$
!!$             stop
!!$          
!!$          Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$             print *,'MW+LMC NOT IMPLEMENTED YET'
!!$
!!$             stop
!!$
!!$          Else If ( ( .not.use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. use_MW_as_anchor ) then
!!$
!!$             print *,'MW NOT IMPLEMENTED YET'
!!$
!!$             stop
!!$
!!$          Else If ( ( .not.use_NGC4258_as_anchor .and. use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$             If (use_metallicity) then 
!!$
!!$                If (use_H_band) then
!!$
!!$                   print *,'H BAND NOT IMPLEMENTED INCLUDING METALLICITY DEPENDENCE WHEN LMC AS ANCHOR'
!!$                     
!!$                   stop
!!$                
!!$                Else
!!$
!!$                   write(UNIT_EXE_FILE,*) 'MEANS FOR THE SAMPLES ARE : '
!!$
!!$                   Do m=1,number_model_parameters
!!$
!!$                      write(UNIT_EXE_FILE,*) ''//trim(paramnames(m))//' = ', means(m)
!!$
!!$                   End Do
!!$
!!$                End If
!!$
!!$             Else
!!$
!!$                If (use_H_band) then
!!$
!!$                   print *,'H BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AS ANCHOR'
!!$                     
!!$                   stop
!!$
!!$                Else
!!$
!!$                   print *,'W BAND NOT IMPLEMENTED WITHOUT METALLICITY DEPENDENCE WHEN LMC AS ANCHOR'
!!$                     
!!$                   stop
!!$                        
!!$                End If
!!$
!!$             End If
!!$
!!$          Else If ( ( use_NGC4258_as_anchor .and. .not.use_LMC_as_anchor ) .and. (.not.use_MW_as_anchor) ) then
!!$
!!$             If (use_metallicity) then 
!!$
!!$                If (use_H_band) then
!!$
!!$                   write(UNIT_EXE_FILE,*) 'MEANS FOR THE SAMPLES ARE : '
!!$
!!$                   Do m=1,number_model_parameters
!!$
!!$                      write(UNIT_EXE_FILE,*) ''//trim(paramnames(m))//' = ', means(m)
!!$
!!$                   End Do
!!$
!!$                Else
!!$
!!$                   write(UNIT_EXE_FILE,*) 'MEANS FOR THE SAMPLES ARE : '
!!$
!!$                   Do m=1,number_model_parameters
!!$
!!$                      write(UNIT_EXE_FILE,*) ''//trim(paramnames(m))//' = ', means(m)
!!$
!!$                   End Do
!!$
!!$                End If
!!$

  deallocate(cmbmask,planckmap)

  close(UNIT_EXE_FILE)

End Program vsk




