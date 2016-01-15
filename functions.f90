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

  Subroutine compute_vsk_angular_power_spectra()

    use healpix_types
    use pix_tools, only: nside2npix
    use fitstools, only: getsize_fits, input_map, output_map
    use arrays
    use fiducial

    Implicit none

    Integer*4 :: p

    Character(len=4) :: x

    Do p=1,number_of_cmb_simulations

       write(x,fmt) p

       call write_parameter_file_polspice(p,'V')

       call write_parameter_file_polspice(p,'S')

       call write_parameter_file_polspice(p,'K')

       call system('./PolSpice_v03-01-06/src/spice -optinfile '//trim(PATH_TO_POLSPICE_PARAMETER_FILE)//&
            ''//trim('vmap')//'_'//trim(x)//'.spicerc')

       call system('./PolSpice_v03-01-06/src/spice -optinfile '//trim(PATH_TO_POLSPICE_PARAMETER_FILE)//&
            ''//trim('smap')//'_'//trim(x)//'.spicerc')


       call system('./PolSpice_v03-01-06/src/spice -optinfile '//trim(PATH_TO_POLSPICE_PARAMETER_FILE)//&
            ''//trim('kmap')//'_'//trim(x)//'.spicerc')

    End Do

  End subroutine compute_vsk_angular_power_spectra

  Subroutine write_parameter_file_polspice(iseed,map_type)

    use fiducial
    Implicit none

    Integer*4 :: iseed

    Character(len=4) :: x,y
    Character(len=1) :: map_type   ! IT CAN ONLY TAKE VALUES 'V','S',AND 'K'

    write(x,fmt) iseed

    If (iseed .eq. 0) then

       If (map_type .eq. 'V') then

          open(UNIT_SYNFAST_PAR_FILE,file=PATH_TO_POLSPICE_PARAMETER_FILE//trim('vmap_smica')//'.spicerc')

       Else If (map_type .eq. 'S') then

          open(UNIT_SYNFAST_PAR_FILE,file=PATH_TO_POLSPICE_PARAMETER_FILE//trim('smap_smica')//'.spicerc')

       Else If (map_type .eq. 'K') then

          open(UNIT_SYNFAST_PAR_FILE,file=PATH_TO_POLSPICE_PARAMETER_FILE//trim('kmap_smica')//'.spicerc')

       Else

          print *,'"map_type" VARIABLE CAN ONLY TAKE VALUES "V","S","K"'

          stop

       End If

    Else

       If (map_type .eq. 'V') then

          open(UNIT_SYNFAST_PAR_FILE,file=PATH_TO_POLSPICE_PARAMETER_FILE//trim('vmap')//'_'//trim(x)//'.spicerc')

       Else If (map_type .eq. 'S') then

          open(UNIT_SYNFAST_PAR_FILE,file=PATH_TO_POLSPICE_PARAMETER_FILE//trim('smap')//'_'//trim(x)//'.spicerc')

       Else If (map_type .eq. 'K') then

          open(UNIT_SYNFAST_PAR_FILE,file=PATH_TO_POLSPICE_PARAMETER_FILE//trim('kmap')//'_'//trim(x)//'.spicerc')

       Else

          print *,'"map_type" VARIABLE CAN ONLY TAKE VALUES "V","S","K"'

          stop

       End If

    End If

!    write(UNIT_SYNFAST_PAR_FILE,*) 'beam = ', sqrt(3.d0/Pi)*3600.d0/dble(nsideC)*3.d0  

    If (iseed .eq. 0) then

       If (map_type .eq. 'V') then

          write(UNIT_SYNFAST_PAR_FILE,*) 'clfile = ', PATH_TO_VSK_SPECTRA//trim('vl_smica.cl')//''

          write(UNIT_SYNFAST_PAR_FILE,*) 'corfile = ', PATH_TO_VSK_SPECTRA//trim('cor_vsmica.cor')//''

          write(UNIT_SYNFAST_PAR_FILE,*) 'mapfile = ', './zero_mean_vsk_maps/vmap_smica.fits'

       Else if (map_type .eq. 'S') then

          write(UNIT_SYNFAST_PAR_FILE,*) 'clfile = ', PATH_TO_VSK_SPECTRA//trim('sl_smica.cl')//''

          write(UNIT_SYNFAST_PAR_FILE,*) 'corfile = ', PATH_TO_VSK_SPECTRA//trim('cor_ssmica.cor')//''

          write(UNIT_SYNFAST_PAR_FILE,*) 'mapfile = ', './zero_mean_vsk_maps/smap_smica.fits'

       Else

          write(UNIT_SYNFAST_PAR_FILE,*) 'clfile = ', PATH_TO_VSK_SPECTRA//trim('kl_smica.cl')//''

          write(UNIT_SYNFAST_PAR_FILE,*) 'corfile = ', PATH_TO_VSK_SPECTRA//trim('cor_ksmica.cor')//''

          write(UNIT_SYNFAST_PAR_FILE,*) 'mapfile = ', './zero_mean_vsk_maps/kmap_smica.fits'

       End If

    Else

       If (map_type .eq. 'V') then

          write(UNIT_SYNFAST_PAR_FILE,*) 'clfile = ', PATH_TO_VSK_SPECTRA//trim('vl')//'_'//trim(x)//'.cl'

          write(UNIT_SYNFAST_PAR_FILE,*) 'corfile = ', PATH_TO_VSK_SPECTRA//trim('corv')//'_'//trim(x)//'.cor'

          write(UNIT_SYNFAST_PAR_FILE,*) 'mapfile = ', './zero_mean_vsk_maps/vmap_'//trim(x)//'.fits'

       Else if (map_type .eq. 'S') then

          write(UNIT_SYNFAST_PAR_FILE,*) 'clfile = ', PATH_TO_VSK_SPECTRA//trim('sl')//'_'//trim(x)//'.cl'

          write(UNIT_SYNFAST_PAR_FILE,*) 'corfile = ', PATH_TO_VSK_SPECTRA//trim('cors')//'_'//trim(x)//'.cor'

          write(UNIT_SYNFAST_PAR_FILE,*) 'mapfile = ', './zero_mean_vsk_maps/smap_'//trim(x)//'.fits'

       Else

          write(UNIT_SYNFAST_PAR_FILE,*) 'clfile = ', PATH_TO_VSK_SPECTRA//trim('kl')//'_'//trim(x)//'.cl'

          write(UNIT_SYNFAST_PAR_FILE,*) 'corfile = ', PATH_TO_VSK_SPECTRA//trim('cork')//'_'//trim(x)//'.cor'

          write(UNIT_SYNFAST_PAR_FILE,*) 'mapfile = ', './zero_mean_vsk_maps/kmap_'//trim(x)//'.fits'

       End If

    End If

    write(UNIT_SYNFAST_PAR_FILE,*) 'maskfile = ', PATH_TO_VSK_MASK

    write(UNIT_SYNFAST_PAR_FILE,*) 'subav = YES'

    write(y,fmt) nsideC 

    write(UNIT_SYNFAST_PAR_FILE,*) 'pixelfile = ', PATH_TO_HEALPIX_DATA//trim('/pixel_window_')//'n'//trim(y)//'.fits'

    close(UNIT_SYNFAST_PAR_FILE)

  End subroutine write_parameter_file_polspice

  Subroutine compute_variance_skewness_kurtosis_maps(PATH_TO_CMB_MAP,x)

    use healpix_types
    use udgrade_nr, only: udgrade_ring, udgrade_nest
    use pix_tools, only: nside2npix,convert_ring2nest, convert_nest2ring,remove_dipole
    use fitstools, only: getsize_fits, input_map, output_map
    use head_fits
    use fiducial
    use arrays

    Implicit none

    Integer*4 :: m
    Integer(kind=I8B) :: d,i,f,j ! NUMBERS OF PIXELS CMB MAPS
    Integer(kind=I8B), allocatable, dimension(:) :: indexc 

    Real(kind=DP), allocatable, dimension(:,:) :: vmap,kmap,smap,cmbmapa,kmap2,kmap3,vskmask
    Real(kind=DP), allocatable, dimension(:) :: pcmb,pcmb2
    Real(kind=DP),dimension(0:DEGREE_REMOVE_DIPOLE*DEGREE_REMOVE_DIPOLE-1) :: multipoles              ! SAVES MONOPOLE AND DIPOLE OF CMB MAP 
    Real(kind=DP),dimension(1:2) :: zbounds                 ! BOUNDS TO COMPUTE DIPOLE AND MONOPOLE

    Character(len=*) :: PATH_TO_CMB_MAP
    Character(len=80),dimension(1:60) :: header
    Character(len=4) :: x

    Logical :: exist,computing_data

    zbounds(:) = 0.d0 

    computing_data = .false.

    allocate (kmap(0:npixC-1,1:1), smap(0:npixC-1,1:1), vmap(0:npixC-1,1:1),& 
         cmbmapa(0:n-1,1:1), vskmask(0:npixC-1,1:1),stat = status1)

    call write_minimal_header(header, 'MAP', nside = nsideC, ordering = ORDERING_VSK_MAPS, coordsys = SYS_COORD) ! HEADER OF V, S, K-MAPS

    If (PATH_TO_CMB_MAP .eq. PATH_TO_PLANCK_CMB_MAP) then

       cmbmapa(0:,1:1) = planckmap(0:,1:1)*cmbmask(0:,1:1)  ! USES MASK UT78. MONOPOLE AND DIPOLE REMOVED ALREADY

       computing_data = .true.

    Else

       call input_map(PATH_TO_CMB_MAP, cmbmapa(0:n-1,1:1), n, 1) 

       cmbmapa = cmbmapa*1.d-6 ! CONVERSION OF UNITS IN CMB MAP: \mu K_CMB -> K_CMB AS IN PLANCK MAP

       call remove_dipole(nsmax,cmbmapa(0:n-1,1),RING_ORDERING,DEGREE_REMOVE_DIPOLE,multipoles,zbounds,HPX_DBADVAL,cmbmask(0:n-1,1))

       cmbmapa(0:,1:1) = cmbmapa(0:,1:1)*cmbmask(0:,1:1)  ! USES MASK UT78

    End If

    kmap(0:npixC-1,1) = HPX_DBADVAL ! INITIALIZATION OF K MAP
    smap(0:npixC-1,1) = HPX_DBADVAL ! INITIALIZATION OF S MAP
    vmap(0:npixC-1,1) = HPX_DBADVAL ! INITIALIZATION OF V MAP

    Do m=0,npixC-1

       allocate(kmap2(0:n-1,1:1),kmap3(0:npixC-1,1:1))

       kmap2 = cmbmapa  ! RING ORDERING 

       call udgrade_ring(kmap2,nsmax,kmap3,nsideC)

       kmap3(0:npixC-1,1) = HPX_DBADVAL ! INITIALIZATION OF AUXILIAR MAP IN RING ORDERING

       kmap3(m,1) = 1.d0  ! One assigned to current cell 

       call udgrade_ring(kmap3,nsideC,kmap2,nsmax)  ! CHANGES RESOLUTION OF K AUXILIARY MAP TO THAT OF CMB MAP 

       allocate (indexc(0:c-1), stat = status4) ! USED TO STORE INDICES OF CURRENT CELL

       d = 0   ! COUNTS PIXELS WHICH EQUAL ONE IN CURRENT CELL

       Do i = 0,n-1                                                
                  
          If ( kmap2(i,1) .eq. HPX_DBADVAL) then 

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

          vmap(m,1) = HPX_DBADVAL 

          smap(m,1) = HPX_DBADVAL

          kmap(m,1) = HPX_DBADVAL

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

    If (computing_data) then  ! VARIANCE, KURTOSIS, SKEWNESS MAPS OF UNMASKED PIXELS ARE WRITTEN TO FILES

       call output_map(vmap,header,'./vsk_maps/vmap_smica.fits')

       call output_map(smap,header,'./vsk_maps/smap_smica.fits')

       call output_map(kmap,header,'./vsk_maps/kmap_smica.fits')

    Else                      ! VARIANCE, KURTOSIS, SKEWNESS MAPS OF UNMASKED PIXELS ARE WRITTEN TO FILES

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

    deallocate(vmap,kmap,smap,cmbmapa,vskmask)

  End Subroutine compute_variance_skewness_kurtosis_maps

  Subroutine compute_vsk_mean_maps()

    use healpix_types
    use fitstools, only: input_map, output_map
    use head_fits
    use fiducial
    use arrays

    Implicit none

    Integer*4 :: m,t

    Real(kind=DP), allocatable, dimension(:,:) :: vmap,kmap,smap
!    Real(kind=DP) :: fmissval 

    Character(len=4) :: x    
    Character(len=80),dimension(1:60) :: header

    call write_minimal_header(header, 'MAP', nside = nsideC, ordering = ORDERING_VSK_MAPS, coordsys = SYS_COORD) ! HEADER OF V, S, K-MAPS

    allocate(vsim(0:npixC-1,1:number_of_cmb_simulations),ssim(0:npixC-1,1:number_of_cmb_simulations),&
       ksim(0:npixC-1,1:number_of_cmb_simulations),vmean(0:npixC-1,1:1),smean(0:npixC-1,1:1),kmean(0:npixC-1,1:1),&
    vsdv(0:npixC-1,1:1),ssdv(0:npixC-1,1:1),ksdv(0:npixC-1,1:1),vmap(0:npixC-1,1:1),smap(0:npixC-1,1:1),&
    kmap(0:npixC-1,1:1),stat = status1)

    If (status1 .eq. 0) then

       continue

    Else

       print *, 'PROBLEM WITH MEMORY IN SUBROUTINE compute_vsk_mean_maps'

       stop

    End If

    Do m=1,number_of_cmb_simulations

       write(x,fmt) m

       call input_map('./vsk_maps/vmap_'//trim(x)//'.fits', vsim(0:npixC-1,m:m), npixC, 1, fmissval = HPX_DBADVAL)

       call input_map('./vsk_maps/smap_'//trim(x)//'.fits', ssim(0:npixC-1,m:m), npixC, 1, fmissval = HPX_DBADVAL)

       call input_map('./vsk_maps/kmap_'//trim(x)//'.fits', ksim(0:npixC-1,m:m), npixC, 1, fmissval = HPX_DBADVAL)

    End Do

    Do m=0,npixC-1

       call mean(vsim(m,1:number_of_cmb_simulations),vmean(m,1))

       call sdv(vsim(m,1:number_of_cmb_simulations),vsdv(m,1))

       call mean(ssim(m,1:number_of_cmb_simulations),smean(m,1))

       call sdv(ssim(m,1:number_of_cmb_simulations),ssdv(m,1))

       call mean(ksim(m,1:number_of_cmb_simulations),kmean(m,1))

       call sdv(ksim(m,1:number_of_cmb_simulations),ksdv(m,1))

    End Do
    
    call output_map(vmean,header,'./vsk_maps/vmap_mean.fits')

    call output_map(smean,header,'./vsk_maps/smap_mean.fits')

    call output_map(kmean,header,'./vsk_maps/kmap_mean.fits')

    call output_map(vsdv,header,'./vsk_maps/vmap_sdv.fits')

    call output_map(ssdv,header,'./vsk_maps/smap_sdv.fits')

    call output_map(ksdv,header,'./vsk_maps/kmap_sdv.fits')

    call input_map('./vsk_maps/vmap_smica.fits',vmap(0:npixC-1,1:1),npixC,1, fmissval = HPX_DBADVAL)

    call input_map('./vsk_maps/smap_smica.fits',smap(0:npixC-1,1:1),npixC,1, fmissval = HPX_DBADVAL)

    call input_map('./vsk_maps/kmap_smica.fits',kmap(0:npixC-1,1:1),npixC,1, fmissval = HPX_DBADVAL)

    Do m=0,npixC-1

       If (vmap(m,1) .eq. HPX_DBADVAL) then 

          continue

       Else

          vmap(m,1) = vmap(m,1) - vmean(m,1) 

          smap(m,1) = smap(m,1) - smean(m,1) 

          kmap(m,1) = kmap(m,1) - kmean(m,1) 

       End If

    End Do

    call output_map(vmap,header,'./zero_mean_vsk_maps/vmap_smica.fits')

    call output_map(smap,header,'./zero_mean_vsk_maps/smap_smica.fits')

    call output_map(kmap,header,'./zero_mean_vsk_maps/kmap_smica.fits')

    Do m=1,number_of_cmb_simulations

       write(x,fmt) m

       Do t=0,npixC-1

          If (vmap(t,1) .eq. HPX_DBADVAL) then

             continue 

          Else

             vsim(t,m:m) = vsim(t,m:m) - vmean(t,1:1) 

             ssim(t,m:m) = ssim(t,m:m) - smean(t,1:1) 
 
             ksim(t,m:m) = ksim(t,m:m) - kmean(t,1:1) 

          End If

       End Do

       call output_map(vsim(0:npixC-1,m:m),header,'./zero_mean_vsk_maps/vmap_'//trim(x)//'.fits')

       call output_map(ssim(0:npixC-1,m:m),header,'./zero_mean_vsk_maps/smap_'//trim(x)//'.fits')

       call output_map(ksim(0:npixC-1,m:m),header,'./zero_mean_vsk_maps/kmap_'//trim(x)//'.fits')

    End Do

    deallocate(vsim,ssim,ksim,vmap,smap,kmap)

  End subroutine compute_vsk_mean_maps

  Subroutine mean(array,var)

    use fgsl

    Implicit none

    Integer(fgsl_size_t) :: nsize 

    Real(fgsl_double) :: array(:),var

    nsize = size(array)
    
    var = fgsl_stats_mean(array,1_fgsl_size_t,nsize)

  End Subroutine mean

  Subroutine variance(array,var)

    use fgsl

    Implicit none

    Integer(fgsl_size_t) :: nsize 

    Real(fgsl_double) :: array(:),var

    nsize = size(array)
    
    var = fgsl_stats_variance(array,1_fgsl_size_t,nsize)

  End Subroutine variance

  Subroutine sdv(array,var)

    use fgsl

    Implicit none

    Integer(fgsl_size_t) :: nsize 

    Real(fgsl_double) :: array(:),var

    nsize = size(array)
    
    var = fgsl_stats_sd(array,1_fgsl_size_t,nsize)

  End Subroutine sdv

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
