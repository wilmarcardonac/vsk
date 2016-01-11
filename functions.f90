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

  Subroutine compute_variance_skewness_kurtosis_maps(PATH_TO_CMB_MAP,x)

    use healpix_types
    use udgrade_nr, only: udgrade_ring, udgrade_nest
    use pix_tools, only: nside2npix,convert_ring2nest, convert_nest2ring
    use fitstools, only: getsize_fits, input_map, output_map
    use head_fits
    use fiducial
    use arrays

    Implicit none

    Integer*4 :: m
    Integer*8 :: d,i,f,j ! NUMBERS OF PIXELS CMB MAPS
    Integer*8, allocatable, dimension(:) :: indexc 

    Real*8, allocatable, dimension(:,:) :: vmap,kmap,smap,cmbmapa,cmbmapb,kmap2,kmap3,vskmask
    Real*8, allocatable, dimension(:) :: pcmb,pcmb2

    Character(len=*) :: PATH_TO_CMB_MAP
    Character(len=80),dimension(1:60) :: header
    Character(len=4) :: x

    Logical :: exist,computing_data

    computing_data = .false.

    allocate (kmap(0:npixC-1,1:1), smap(0:npixC-1,1:1), vmap(0:npixC-1,1:1),& 
         cmbmapa(0:n-1,1:1), cmbmapb(0:n-1,1:1),vskmask(0:npixC-1,1:1),stat = status1)

    call write_minimal_header(header, 'MAP', nside = nsideC, ordering = ORDERING_VSK_MAPS, coordsys = SYS_COORD) ! HEADER OF V, S, K-MAPS

    If (PATH_TO_CMB_MAP .eq. PATH_TO_PLANCK_CMB_MAP) then

       cmbmapa(0:,1:1) = planckmap(0:,1:1)*cmbmask(0:,1:1)  ! USES MASK UT78 

       computing_data = .true.

    Else

       call input_map(PATH_TO_CMB_MAP, cmbmapa(0:n-1,1:1), n, 1) 

       cmbmapa(0:,1:1) = cmbmapa(0:,1:1)*cmbmask(0:,1:1)  ! USES MASK UT78

    End If

    kmap(0:npixC-1,1) = miss ! INITIALIZATION OF K MAP
    smap(0:npixC-1,1) = miss ! INITIALIZATION OF S MAP
    vmap(0:npixC-1,1) = miss ! INITIALIZATION OF V MAP

    Do m=0,npixC-1

       allocate(kmap2(0:n-1,1:1),kmap3(0:npixC-1,1:1))

       kmap2 = cmbmapa

       call udgrade_ring(kmap2,nsmax,kmap3,nsideC)

       kmap3(0:npixC-1,1) = miss ! INITIALIZATION OF AUXILIAR MAP

       kmap3(m,1) = 1.d0  ! One assigned to current cell 

       call udgrade_ring(kmap3,nsideC,kmap2,nsmax)  ! CHANGES RESOLUTION OF K AUXILIARY MAP TO THAT OF CMB MAP 

       allocate (indexc(0:c-1), stat = status4) ! USED TO STORE INDICES OF CURRENT CELL

       d = 0   ! COUNTS PIXELS WHICH EQUAL ONE IN CURRENT CELL

       Do i = 0,n-1                                                
                  
          If ( kmap2(i,1) .eq. miss) then 

             continue

          Else  ! ASSIGN VALUES OF INDICES IN CURRENT CELL (CORRESPONDING TO CMB NSIDE)

             indexc(d) = i

             d = d +1

          End if

       End Do

       deallocate(kmap2,kmap3)

       allocate (pcmb(0:c-1), stat = status5)

       f = 0   ! COUNTS NON-ZERO PIXELS IN CURRENT CELL OF THE CMB MASK 

       Do i = 0,c-1
 
          pcmb(i) = cmbmapa(indexc(i),1)

          If (cmbmask(indexc(i),1) .eq. 0.d0) then 

             continue

          Else

             f = f + 1

          End if

       End Do
 
       If (f .lt. fr) then ! MASK CURRENT CELL 

          vmap(m,1) = miss 

          smap(m,1) = miss

          kmap(m,1) = miss

          vskmask(m,1) = 0.d0

       Else

          allocate (pcmb2(0:f-1),stat = status6)

          j = 0

          Do i=0,c-1

             If (cmbmask(indexc(i),1) .eq. 0.d0) then

                continue

             Else

                pcmb2(j) = pcmb(i)

                j = j + 1

             End If

          End Do

          call variance(pcmb2,vmap(m,1))

          call skewness(pcmb2,smap(m,1))

          call kurtosis(pcmb2,kmap(m,1))

          vskmask(m,1) = 1.d0

          deallocate(pcmb2,stat = status6)

       End if

       deallocate (indexc,pcmb, stat = status5)

    End Do

    If (computing_data) then 

       call output_map(vmap,header,'./vsk_maps/vmap_'//trim('smica')//'.fits')

       call output_map(smap,header,'./vsk_maps/smap_'//trim('smica')//'.fits')

       call output_map(kmap,header,'./vsk_maps/kmap_'//trim('smica')//'.fits')

    Else

       call output_map(vmap,header,'./vsk_maps/vmap_'//trim(x)//'.fits')

       call output_map(smap,header,'./vsk_maps/smap_'//trim(x)//'.fits')

       call output_map(kmap,header,'./vsk_maps/kmap_'//trim(x)//'.fits')

    End If

    inquire(file = './vsk_maps/vsk_mask.fits',exist=exist)

    If (exist) then

       continue

    Else

       call output_map(vskmask,header,'./vsk_maps/vsk_mask.fits')

    End If

    deallocate(vmap,kmap,smap,cmbmapa,cmbmapb,vskmask)

  End Subroutine compute_variance_skewness_kurtosis_maps

  Subroutine variance(array,var)

    use fgsl

    Implicit none

    Integer(fgsl_size_t) :: nsize 

    Real(fgsl_double) :: array(:),var

    nsize = size(array)
    
    var = fgsl_stats_variance(array,1_fgsl_size_t,nsize)

  End Subroutine variance

  Subroutine skewness(array,var)

    use fgsl

    Implicit none

    Integer(fgsl_size_t) :: nsize 

    Real(fgsl_double) :: array(:),var

    nsize = size(array)
    
    var = fgsl_stats_skew(array,1_fgsl_size_t,nsize)

  End Subroutine skewness

  Subroutine kurtosis(array,var)

    use fgsl

    Implicit none

    Integer(fgsl_size_t) :: nsize 

    Real(fgsl_double) :: array(:),var

    nsize = size(array)
    
    var = fgsl_stats_kurtosis(array,1_fgsl_size_t,nsize)

  End Subroutine kurtosis

End module functions
