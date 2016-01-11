Module functions

  Implicit none
     
contains

  Subroutine generate_gaussian_cmb_map()

    use healpix_types
    use pix_tools, only: nside2npix
    use fitstools, only: getsize_fits, input_map, output_map
    use arrays
    use fiducial

    Implicit none

    Integer*4 :: p

    Character(len=4) :: x

    Logical :: exist

    Do p=1,number_of_cmb_simulations

       write(x,fmt) p

       call write_parameter_file_synfast(p)

       inquire(file = PATH_TO_CMB_MAPS//trim(x)//'.fits',exist=exist)

       If (exist) then
          
          continue

       Else

          call system('synfast -d '//trim(PATH_TO_SYNFAST_PARAMETER_FILE)//''//trim(x)//'.par')

       End If

    End Do

  End subroutine generate_gaussian_cmb_map

  Subroutine write_parameter_file_synfast(iseed)

    use fiducial
    Implicit none

    Integer*4 :: iseed
    Character(len=4) :: x

    write(x,fmt) iseed

    open(UNIT_SYNFAST_PAR_FILE,file=PATH_TO_SYNFAST_PARAMETER_FILE//trim(x)//'.par')

    write(UNIT_SYNFAST_PAR_FILE,*) 'simul_type = ', simul_type

    write(UNIT_SYNFAST_PAR_FILE,*) 'nsmax = ', nsmax

    write(UNIT_SYNFAST_PAR_FILE,*) 'nlmax = ', nlmax

    write(UNIT_SYNFAST_PAR_FILE,*) 'infile = ', infile

    write(UNIT_SYNFAST_PAR_FILE,*) 'iseed = ', iseed

    write(UNIT_SYNFAST_PAR_FILE,*) 'fwhm_arcmin = ', fwhm_arcmin

    write(UNIT_SYNFAST_PAR_FILE,*) 'beam_file = ', beam_file

    write(UNIT_SYNFAST_PAR_FILE,*) 'almsfile = ', almsfile

    write(UNIT_SYNFAST_PAR_FILE,*) 'plmfile = ', plmfile

    write(UNIT_SYNFAST_PAR_FILE,*) 'outfile = ', PATH_TO_CMB_MAPS//trim(x)//'.fits'

    write(UNIT_SYNFAST_PAR_FILE,*) 'outfile_alms = ', outfile_alms

    write(UNIT_SYNFAST_PAR_FILE,*) 'winfiledir = ', PATH_TO_HEALPIX_DATA

    close(UNIT_SYNFAST_PAR_FILE)

  End subroutine write_parameter_file_synfast

  Subroutine compute_variance_skewness_kurtosis_maps(PATH_TO_CMB_MAP)

    use healpix_types
    use udgrade_nr, only: udgrade_ring, udgrade_nest
    use pix_tools, only: nside2npix,convert_ring2nest, convert_nest2ring
    use fitstools, only: getsize_fits, input_map, output_map
    use head_fits
    use fiducial
    use arrays

    Implicit none

    Integer*4 :: m
    Integer*8 :: d,i,f,j,f2 ! NUMBERS OF PIXELS CMB MAPS
    Integer*8, allocatable, dimension(:) :: indexc 

    Real*8, allocatable, dimension(:,:) :: vmap,kmap,smap,cmbmapa,cmbmapb,kmap2,kmap3
    Real*8, allocatable, dimension(:) :: pcmb,pcmb2

    Character(len=*) :: PATH_TO_CMB_MAP
    Character(len=80),dimension(1:60) :: header,headercmb

    allocate (kmap(0:npixC-1,1:1), smap(0:npixC-1,1:1), vmap(0:npixC-1,1:1),& 
         cmbmapa(0:n-1,1:1), cmbmapb(0:n-1,1:1),&
         kmap2(0:n-1,1:1),kmap3(0:npixC-1,1:1), stat = status1)

    call write_minimal_header(header, 'MAP', nside = nsideC, ordering = ORDERING_VSK_MAPS, coordsys = SYS_COORD) ! HEADER OF V, S, K-MAPS

    call write_minimal_header(headercmb, 'MAP', nside = nsmax, ordering = ORDERING_VSK_MAPS, coordsys = SYS_COORD) ! HEADER OF CMB MAP

    If (PATH_TO_CMB_MAP .eq. PATH_TO_PLANCK_CMB_MAP) then

       cmbmapa(0:,1:1) = planckmap(0:,1:1)*cmbmask(0:,1:1)  ! USES MASK UT78

    Else

       call input_map(PATH_TO_CMB_MAP, cmbmapa(0:n-1,1:1), n, 1) 

       cmbmapa(0:,1:1) = cmbmapa(0:,1:1)*cmbmask(0:,1:1)  ! USES MASK UT78

    End If

    kmap(0:npixC-1,1:1) = miss ! INITIALIZATION OF K MAP
    smap(0:npixC-1,1:1) = miss ! INITIALIZATION OF S MAP
    vmap(0:npixC-1,1:1) = miss ! INITIALIZATION OF V MAP
    kmap3(0:npixC-1,1:1) = miss ! INITIALIZATION OF AUXILIAR MAP

    Do m=0,npixC-1

       print *, 'CURRENT CELL :', m+1

       kmap3(m,1:1) = 1.d0  ! One assigned to current cell 

       call udgrade_nest(kmap3,2,kmap2,2048)  ! CHANGES RESOLUTION OF K AUXILIARY MAP TO THAT OF CMB MAP 

       !       call output_map(kmap3,header,'./output/kmap3.fits')

       !      call output_map(kmap2,headercmb,'./output/kmap2.fits')

       allocate (indexc(0:c-1), stat = status4) ! USED TO STORE INDICES OF CURRENT CELL

       d = 0   ! COUNTS PIXELS WHICH EQUAL ONE IN CURRENT CELL

       Do i = 0,n-1                                                
                                          
          If (kmap2(i,1:1) .eq. miss) then 

             continue

          Else  ! ASSIGN VALUES OF INDICES IN CURRENT CELL (CORRESPONDING TO CMB NSIDE)

             indexc(d) = i

             d = d +1

          End if

       End Do

       allocate (pcmb(0:c-1), stat = status5)

       f = 0   ! COUNTS NON-ZERO PIXELS IN CURRENT CELL

       Do i = 0,c-1
 
          pcmb(i) = cmbmapa(indexc(i),1:1)

          If (pcmb(i) .eq. 0.d0) then 

             continue

          Else

             f = f + 1

          End if

       End Do
 
       If (f .lt. fr) then ! REJECT CURRENT CELL 

!!$          kmap(m,1:1) = 0.d0

!!$          smap(m,1:1) = 0.d0

          vmap(m,1:1) = miss 

!!$          !    mask(m,1:1) = 0.d0
!!$          !  continue

       Else

          allocate (pcmb2(0:f-1),stat = status6)

          j = 0

          Do i=0,c-1

             If (pcmb(i) .eq. 0.d0) then

                continue

             Else

                pcmb2(j) = pcmb(i)

                j = j + 1

             End If

          End Do

!!$          print *, 'Computing sample kurtosis, sample skewness and unbiased sample variance in cell ', m+1,'...'
!!$          kmap(m,1:1) = kurtosis(pcmb2)
!!$          smap(m,1:1) = skewness(pcmb2)
          call variance(pcmb2,vmap(m,1:1))

          deallocate(pcmb2,stat = status6)

       End if
!!$
!!$       !##########################################################################################
!!$       !##########################################################################################
!!$       ! We deallocate memory of cmb map pixels used to compute current cell and  array of indices
!!$       !##########################################################################################
!!$       !##########################################################################################
!!$
!!$       deallocate (indexc, stat = status4)
!!$
!!$       deallocate (pcmb, stat = status5)
!!$
!!$       !deallocate (kmap2)
!!$
!!$       print *,'deallocate status ', status4, status5
!!$       !####################################################################
!!$       !####################################################################
!!$       ! We recover the initial value to the map used to extract the indices
!!$       !####################################################################
!!$       !####################################################################
!!$
!!$       kmap3(m,1:1) = miss
!!$
!!$       kmap2(:,1:1) = 0.d0
!!$
!!$       !#######################################################################################
!!$       !#######################################################################################
!!$       ! The loop to compute sample kurtosis, sample skewness and unbiased sample variance ends
!!$       !#######################################################################################
!!$       !#######################################################################################
!!$
    End DO
!!$
!!$    If (p == 1) then
!!$       f2 = 0                                                                                 
!!$       Do m=0,npixC-1
!!$          If (kmap(m,1) == 0.d0) then                                                   
!!$             f2 = f2 + 1                                                              
!!$          Else                                                                          
!!$             continue                                                                  
!!$          End if
!!$       End Do
!!$    Else
!!$       Go to 1
!!$    End if
!!$
!!$1   Do m=0,npixC-1                                                                           
!!$       If (kmap(m,1) == 0.d0) then                                                          
!!$          vmap(m,1:1) = mean_2(vmap(:,1),f2)                                               
!!$          smap(m,1:1) = mean_2(smap(:,1),f2)                                               
!!$          kmap(m,1:1) = mean_2(kmap(:,1),f2)                                               
!!$       Else                                                                                 
!!$          continue                                                                         
!!$       End if
!!$    End Do
!!$
!!$    !###################################################
!!$    !###################################################
!!$    ! We convert K-S-V maps from NESTED to RING ordering
!!$    !###################################################
!!$    !###################################################
!!$
!!$    !Go to 2
!!$
!!$    !#####################################################
!!$    !#####################################################
!!$    ! We write the K,S,V maps and the mask into fits files
!!$    !#####################################################
!!$    !#####################################################
!!$
!!$    call output_map(kmap,header,'/space/wilmar.cardona/projects/B/ksv-maps/planck-nside/gaussian/48c/nm/kmap_'//trim(x)//'.fits')
!!$
!!$    call output_map(smap,header,'/space/wilmar.cardona/projects/B/ksv-maps/planck-nside/gaussian/48c/nm/smap_'//trim(x)//'.fits')
!!$
!!$    call output_map(vmap,header,'/space/wilmar.cardona/projects/B/ksv-maps/planck-nside/gaussian/48c/nm/vmap_'//trim(x)//'.fits')


    deallocate(vmap,kmap,smap,cmbmapa,cmbmapb,kmap2,kmap3)

  End Subroutine compute_variance_skewness_kurtosis_maps

  Subroutine variance(array,var)

    use fgsl
    use mod_unit

    Implicit none

    Integer(fgsl_size_t),parameter :: nsize = size(array)

    Real(fgsl_double) :: array(nsize)
    Real(fgsl_double) :: var
    
    call unit_init(200)

    var = fgsl_stats_variance(array,1_fgsl_size_t,nsize)

    call unit_finalize()
    
  End Subroutine variance

End module functions
