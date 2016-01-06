! edited to run with Healpix 3.*
module deal_with_files
! subroutine check_headers
! subroutine assert_read_piofile
! subroutine set_bad_pixels_to_0
! subroutine input_files
! subroutine output_results
! subroutine check_header
! subroutine set_nlmax
! subroutine check_header_asc
! subroutine check_open_output_files
! subroutine write_corcl_file
! subroutine read_corcl_file
! subroutine check_mylread
! subroutine generate_gaussian_beam
! subroutine my_generate_beam
! function my_isfits
! subroutine get_pixwin_filename
! subroutine read_twod_fits_DP
! subroutine write_twod_fits_DP
! subroutine write_nd_fits_DP
! subroutine check_twod_fits_header
! subroutine check_nd_fits_header
! subroutine deletefile
!!!!!!!! function my_getnumext_fits
! subroutine write_spice_header

contains
!=======================================================================
subroutine check_headers
!=======================================================================
  use spice_common
  use misc_utils, only: fatal_error, string
#ifdef PIO
  use piolib
#endif
  implicit none
  integer(I4B) :: my_npixtot
  logical(LGT) :: head_defined=.false.!must be .false. to get nside, npixtot, nlmax from file
  ! then subsequent calls with .true. will compare previsou nside and npixtot to file value
  integer(I8B) :: myerr,my_group,my_nx8,firstline,lastline
  character(len=FILENAMELEN) :: grpname='' !SP
  integer(i4b) :: imap, jstokes, k

  do imap=1,2
     do jstokes=1,3
        do k=1,nlmaps(jstokes,imap)
           call check_header(listmapfile(k,jstokes,imap), verbose, megaverbose, .true., &
                &            nside, npixtot, nlmax, ordering_tqu(k,jstokes,imap), dry,&
                &            polarization, head_defined) 
        enddo
     enddo
  enddo


  if (clmapinput) call check_header_asc(cl_inmap_file, nside, nlmax, ncor,&
 &                                      head_defined, verbose, megaverbose, npixtot)
  if (masks_present) then
     call check_header(maskfile, verbose, megaverbose, .false.,&
          &            nside, npixtot, nlmax, ordering_m1, dry,&
          &            .false., head_defined)
  endif
  if (masks2_present) then
     call check_header(maskfile2, verbose, megaverbose, .false.,&
 &                     nside, npixtot, nlmax, ordering_m2, dry, &
 &                     .false., head_defined)
  endif
  if (clmaskinput) then
     ncmask = -1 ! read from file actual number of columns
     call check_header_asc(cl_inmask_file, nside, nlmax, ncmask,&
 &                         head_defined, verbose, megaverbose, npixtot)
  endif
  if (weights_present) then
     call check_header(weightfile, verbose, megaverbose, .false.,&
 &                     nside, npixtot, nlmax, ordering_w1, dry, &
 &                     .false., head_defined)
  endif
  if (weights2_present) then
     call check_header(weightfile2, verbose, megaverbose, .false.,&
 &                     nside, npixtot, nlmax, ordering_w2, dry, &
 &                     .false., head_defined)
  endif

  if (windowinput) then
     if (megaverbose) write(*,*) 'Check header for file '//trim(window_in_file)
#ifdef PIO
     myerr=PIOgetgrpname(grpname,window_in_file)
     my_group=PIOopentab2dgrp(grpname,"r")
     myerr=PIOgettab2dlines(firstline,lastline,window_in_file,my_group)
     my_nx8=PIOgettab2dcolumngrp(my_group)
     if (myerr /= 0) then
        write(*,*) 'ERROR in check_headers'
        write(*,*) 'Could not check header for file'
        write(*,*) trim(window_in_file)
        write(*,*) PIOerrmess(myerr)
     else
        windowfile_ny=lastline-firstline ! +1 !SP
        windowfile_nx=my_nx8
     endif
#else 
     call check_twod_fits_header(window_in_file,my_npixtot,windowfile_nx,windowfile_ny)
#endif
     if (windowfile_nx < nlmax+1 .or. windowfile_ny < 2*nlmax+1) then
        write(*,*) 'ERROR in check_headers'
        write(*,*) 'The input window function dimensions are too small for current study'
        write(*,*) 'nx found =',windowfile_nx,',   nx needed =',nlmax+1
        write(*,*) 'ny found =',windowfile_ny,',   ny needed =',2*nlmax+1
        CALL FATAL_ERROR
     endif
     if ((windowfile_nx > nlmax+1 .or. windowfile_ny > 2*nlmax+1) .and. verbose) then
        write(*,*) 'WARNING in check_headers'
        write(*,*) 'The input window function dimensions are too large for current study'
        write(*,*) 'only a subset will be used'
        write(*,*) 'nx found =',windowfile_nx,',   nx needed =',nlmax+1
        write(*,*) 'ny found =',windowfile_ny,',   ny needed =',2*nlmax+1
     endif
  endif

  ! This can be put only here because the default name of pixwin_file depends
  ! on the value of nside
  if (correct_pix) call get_pixwin_filename(pixwin_file_name_default, &
 &                                          pixwin_file_name_input, &
 &                                          pixwin_file_name, &
 &                                          default_pix_file,nside,healpix_data,stringlength)

end subroutine check_headers

!--------------------------------------------------------------------
subroutine assert_read_piofile(nread, nexpected, file, code)
  use healpix_types
  use misc_utils, only: fatal_error
  implicit none
  integer(I8B) :: nread
  integer(I4B) :: nexpected
  character(len=*) :: file, code

  if (nread /= nexpected) then
     write(*,*) 'ERROR in '//code
     write(*,*) 'Inconsistent number of data for file'
     write(*,*) trim(file)
     write(*,*) 'Number of pixels detected :',nread
     write(*,*) 'Npixtot should be :',nexpected
     CALL FATAL_ERROR
  endif
  return
end subroutine assert_read_piofile
!--------------------------------------------------------------------
subroutine set_bad_pixels_to_0(map, npix, nmaps, bad_value, verbose, name)
  use spice_parameters, only: KMAP
  use healpix_types
  implicit none

  !integer, parameter :: KMAP = SP
  real(KMAP), dimension(0:npix-1,1:nmaps), intent(inout) :: map
  integer(i4b),               intent(in)    :: npix, nmaps
  real(KMAP),                   intent(in)    :: bad_value
  logical(lgt),               intent(in)    :: verbose
  character(len=*),           intent(in)    :: name

  integer(i4b) :: n_nan, i, j

  n_nan = 0
  if (verbose) write(*,*) 'Setting bad pixels to 0 ('//trim(name)//')'
  do j=1,nmaps
     do i=0,npix-1
        if (abs(map(i,j)/bad_value-1.0) < 1.e-6) map(i,j) = 0.
        if (map(i,j) /= map(i,j)) then ! to detect NaN
           n_nan = n_nan + 1
           map(i,j) = 0.
        endif
     enddo
  enddo
  if (n_nan > 0) then
     write(*,*)          '=============================================================================='
     write(*,'(a,i0,a)') ' **WARNING**: Detected ',n_nan,' NaN-valued pixels in input '//trim(name)//', replaced with 0.'
     write(*,*)          '             Note that the mask/weight were NOT updated.'
     write(*,*)          '=============================================================================='
  endif

  return
end subroutine set_bad_pixels_to_0
!--------------------------------------------------------------------
subroutine myfitsreadmap(mapfile, map, npix, nmaps, fmissval, code)
  use spice_parameters, only: KMAP
  use healpix_types
  use fitstools, only: input_map, getsize_fits
  implicit none

  !integer, parameter :: KMAP = SP
  character(len=*),                        intent(in)  :: mapfile
  real(KMAP), dimension(0:npix-1,1:nmaps), intent(out) :: map
  integer(I4B),                            intent(in)  :: npix, nmaps
  real(KMAP),       optional,              intent(in)  :: fmissval
  character(len=*), optional,              intent(in)  :: code
  integer(I4B) :: njunk, map1type, istokes

  njunk= getsize_fits(mapfile, type=map1type)
  if (map1type == 3 .and. nmaps == 3) then
     ! cut sky polarized map
     do istokes = 1, 3
        call input_map(mapfile, map(0:,istokes:istokes), npix, 1, fmissval=fmissval, extno= istokes-1)
     enddo
  else
     ! full sky map (polarized or not) or unpolarized cut sky map
     call input_map(mapfile, map, npix, nmaps, fmissval=0._KMAP)
  endif

return
end subroutine myfitsreadmap

!=======================================================================
subroutine input_files
!=======================================================================
  use spice_common
  use misc_utils, only: assert_alloc, fatal_error, string
  use pix_tools,  only: convert_nest2ring
#ifdef PIO
  use piolib
  use my_pio_routines
  use fitstools, only: read_dbintab
#else
  use fitstools, only: read_dbintab,input_map, getsize_fits
#endif
  implicit none

  integer(I4B) :: ipixbad,i,j,my_npixtot, status,n_nan
  integer(I4B) :: istokes, map1type, map2type, njunk, kem
  real(DP), allocatable, dimension(:) :: my_kcross
  real(DP) :: dnullval
  logical(LGT) :: anynull

  integer(I8B) :: myerr,mygroup,nbdata, p
  character(len=FILENAMELEN) :: grpname
  real(SP),pointer,dimension(:) :: tmpmap
  real(DP),pointer,dimension(:) :: tmpmapDP
  character(len=*), parameter :: code = 'input_files'
  real(KMAP), allocatable, dimension(:,:) :: tmp_buffer
  integer(I4B) :: nmap_buffer
  integer(I4B) :: jmap, kstokes, k

  nmap_buffer = nmap

! read all the maps (if applicable), re-order them and combine them

  if (maxval(nlmaps(1:3,1)) > 0) then
     allocate(map_in(0:npixtot-1,1:nmap),stat=status)
     call assert_alloc(status, code, 'map_in')
     map_in = 0. ! 2015-03-02
  endif
  if (maxval(nlmaps(1:3,2)) > 0) then
     allocate(map2_in(0:npixtot-1,1:nmap),stat=status)
     call assert_alloc(status, code, 'map2_in')
     map2_in = 0. ! 2015-03-02
  endif

  if (allocated(map_in)) then
     allocate(tmp_buffer(0:npixtot-1,1:nmap_buffer),stat=status) 
     call assert_alloc(status, code, 'tmp_buffer')

     do jmap=1,2 ! 1st (and optional second) map(s)
        do i=1,nlmaps(1, jmap)
           ! read (polarized) map, and set bad pixels to 0
#ifdef PIO
           do kstokes=1,nmap ! T,Q,U
              if (megaverbose) write(*,'(a,i1,a,a)') &
                   & ' Importing (',jmap,')' ,trim(listmapfile(i,kstokes,jmap))
              call mypioreadmap (listmapfile(i,kstokes,jmap), tmp_buffer(0:,kstokes), npixtot, code)
           enddo
           call set_bad_pixels_to_0(tmp_buffer, npixtot, nmap, HPX_SBADVAL, .false., 'map')
#else
           kstokes = 1
           if (megaverbose) write(*,'(a,i1,a,a)') &
                & ' Importing (',jmap,') ',trim(listmapfile(i,kstokes,jmap))
           call myfitsreadmap(listmapfile(i,kstokes,jmap), tmp_buffer, npixtot, nmap, fmissval=0._KMAP)
#endif
           ! reorder to RING if necessary
           if (ordering_tqu(i,1,jmap) == 2) then
              if (megaverbose) print*,'... and reordering ...'
              call convert_nest2ring(nside, tmp_buffer)  
           endif
           ! combine
           if (megaverbose) print*,'... and combining '//trim(string(listmapw8(i,1,jmap)))
           if (jmap == 1) map_in  = map_in  + listmapw8(i,1,jmap) * tmp_buffer
           if (jmap == 2) map2_in = map2_in + listmapw8(i,1,jmap) * tmp_buffer
        enddo
     enddo
     deallocate(tmp_buffer)
  endif

  if (masks_present) then
     if (megaverbose) write(*,*) 'Read input mask file'
     allocate(mask_map(0:npixtot-1,1:nmask1))
#ifdef PIO
     ! do explicit loop to avoid overloading memory stack
     ! mask_map is in stack because it has the 'target' attribute
!      myerr=PIOreadmapobject(tmpmap,maskfile,' ')
!      call assert_read_piofile(myerr, npixtot, maskfile, code)
!      do i=1, npixtot
!         mask_map(i-1,1) = tmpmap(i)
!      enddo
!      myerr=PIOdeletemaptable(tmpmap)
     call mypioreadmap(maskfile, mask_map(0:,1), npixtot, code)
     if (maskfilep_toread) then
!         myerr=PIOreadmapobject(tmpmap,maskfilep,' ')
!         call assert_read_piofile(myerr, npixtot, maskfilep, code)
!         do i=1, npixtot
!            mask_map(i-1,2) = tmpmap(i)
!         enddo
!         myerr=PIOdeletemaptable(tmpmap)
        call mypioreadmap(maskfilep, mask_map(0:,2), npixtot, code)
     endif
#else
     call input_map(maskfile, mask_map(:,1:1), npixtot, 1, fmissval=0._KMAP)
     if (maskfilep_toread) then
        call input_map(maskfilep, mask_map(:,2:2), npixtot, 1, fmissval=0._KMAP)
     endif
#endif	
     ! guaranty that the masks are zero or one
     ipixbad=0
     do i=0,npixtot-1
        if (mask_map(i,1) /= 0. .and. mask_map(i,1) /= 1.) then
           mask_map(i,1)=1.
           ipixbad=ipixbad+1
        endif
     enddo
     if (ipixbad > 0 .and. verbose) then
        write(*,*)
        write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*) 'WARNING in input_files :'
        write(*,*) 'The mask file has pixels different from 0 or 1 '
        write(*,*) 'Pixels different from 0 are set equal to 1'
        write(*,*) 'Fraction of pixels changed :',real(ipixbad)/real(npixtot)
        write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*)
     endif
  endif

  if (masks2_present) then
     if (megaverbose) write(*,*) 'Read input mask file (2)'
     allocate(mask2_map(0:npixtot-1,1:nmask2))
#ifdef PIO
     call mypioreadmap(maskfile2, mask2_map(0:,1), npixtot, code)
     if (maskfilep2_toread) then
        call mypioreadmap(maskfilep2, mask2_map(0:,2), npixtot, code)
     endif
#else
     call input_map(maskfile2, mask2_map, npixtot, 1, fmissval=0._KMAP)
     if (maskfilep2_toread) then
        call input_map(maskfilep2, mask2_map(:,2:2), npixtot, 1, fmissval=0._KMAP)
     endif
#endif	
     ! guaranty that the masks are zero or one
     ipixbad=0
     do i=0,npixtot-1
        if (mask2_map(i,1) /= 0. .and. mask2_map(i,1) /= 1.) then
           mask2_map(i,1)=1.
           ipixbad=ipixbad+1
        endif
     enddo
     if (ipixbad > 0 .and. verbose) then
        write(*,*)
        write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*) 'WARNING in input_files :'
        write(*,*) 'The mask2 file has pixels different from 0 or 1 '
        write(*,*) 'Pixels different from 0 are set equal to 1'
        write(*,*) 'Fraction of pixels changed :',real(ipixbad)/real(npixtot)
        write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*)
     endif
  endif

  if (weights_present) then
     if (megaverbose) write(*,*) 'Read input weight file'
     allocate(weight_map(0:npixtot-1,1:nweight1))
#ifdef PIO
     call mypioreadmap(weightfile, weight_map(0:,1), npixtot, code)
     if (weightfilep_toread) then
        call mypioreadmap(weightfilep, weight_map(0:,2), npixtot, code)
     endif     
#else
     call input_map(weightfile, weight_map(:,1:1), npixtot, 1, fmissval=0._KMAP)
     if (weightfilep_toread) then
        call input_map(weightfilep, weight_map(:,2:2), npixtot, 1, fmissval=0._KMAP)
     endif
#endif
     if (minval(weight_map) < 0.) then
        write(*,*) 'ERROR in input_files :'
        write(*,*) 'The weight map has negative values.'
        CALL FATAL_ERROR
     endif
  endif

  if (weights2_present) then
     if (megaverbose) write(*,*) 'Read input weight file (2)'
     allocate(weight2_map(0:npixtot-1,1:nweight2))
#ifdef PIO
     call mypioreadmap(weightfile2, weight2_map(0:,1), npixtot, code)
     if (weightfilep2_toread) then
        call mypioreadmap(weightfilep2, weight2_map(0:,2), npixtot, code)
     endif
#else
     call input_map(weightfile2, weight2_map, npixtot, 1, fmissval=0._KMAP)
     if (weightfilep2_toread) then
        call input_map(weightfilep2, weight2_map(:,2:2), npixtot, 1, fmissval=0._KMAP)
     endif
#endif
     if (minval(weight2_map) < 0.) then
        write(*,*) 'ERROR in input_files :'
        write(*,*) 'The weight map (2) has negative values.'
        CALL FATAL_ERROR
     endif
  endif

  if (clmapinput) then 
     if (megaverbose) write(*,*) 'Read input precomputed Cls for the map'
     allocate(cl_map_precomp (0:nlmax,1:ncor))
     call read_corcl_file(cl_map_precomp,nlmax,cl_inmap_file,.false.,ncor)
  endif
  
  if (clmaskinput) then
     if (megaverbose) write(*,*) 'Read input precomputed Cls for the masks/weights'
     allocate(cl_mask_precomp(0:nlmax,1:ncmask))
     call read_corcl_file(cl_mask_precomp,nlmax,cl_inmask_file,.false.,ncmask)
  endif

  if (windowinput) then
     if (megaverbose) write(*,*) 'Read input window function'
#ifdef PIO
     nbdata=PIOreadtab2dobject(tmpmapDP,window_in_file,' ')
     if (nbdata < 0_I8B) then
        write(*,*) 'ERROR in input_files :'
        write(*,*) 'Something went wrong while reading file'
        write(*,*) trim(window_in_file)
        write(*,*) PIOerrmess(nbdata)
     endif
     allocate(kcross(0:nlmax,0:2*nlmax))
     do j=0,2*nlmax
        do i=0,nlmax
           kcross(i,j)=tmpmapDP(1 + i + j*windowfile_nx)
        enddo
     enddo
     myerr=PIOdeletetab2dtable(tmpmapDP)
#else
     my_npixtot= windowfile_nx * windowfile_ny
     
     allocate(my_kcross(0:my_npixtot-1))
     call read_twod_fits_DP(window_in_file,my_npixtot,windowfile_nx,windowfile_ny,my_kcross)
     allocate(kcross(0:nlmax,0:2*nlmax))
     do j=0,2*nlmax
        do i=0,nlmax
           kcross(i,j)=my_kcross(i + j*windowfile_nx)
        enddo
     enddo
     deallocate(my_kcross)
#endif
  endif
    
  if (correct_pix) then
     if (megaverbose) write(*,*) 'Read pixel window function'
     allocate(wl(0:4*nside,1)) 
     call read_dbintab(pixwin_file_name,wl,nside*4+1,1,dnullval,anynull)
  endif

  if (subtract_noise_cor.or.subtract_noise_cl) then
     allocate(xi_noise(0:nlmax,1:ncor))
     allocate(cl_noise(0:nlmax,1:ncor))
  endif

  if (subtract_noise_cor) then
     if (megaverbose) write(*,*) 'Read input noise correlation function'
     call read_corcl_file(xi_noise,nlmax,noisecor_file_name,.true.,ncor)
  endif

  if (subtract_noise_cl) then
     if (megaverbose) write(*,*) 'Read input noise Cls'
     call read_corcl_file(cl_noise,nlmax,noisecl_file_name,.false.,ncor)
  endif

  if (correct_transfer_function) then 
     if (megaverbose) write(*,*) 'Read transfer function'
     allocate(transfer_function(0:nlmax,1:ncor))
     call read_corcl_file(transfer_function, nlmax, tf_file, .false., ncor)
  endif

  if (correct_beam .or. correct_beam2) then
     allocate(gb(0:nlmax,1:1)) !<<<<<<<<<<<<<<<<
     gb = 1.0_dp
     allocate(gb2(0:nlmax,1:1)) !<<<<<<<<<<<<<<<<
     gb2 = 1.0_dp
     if (map2_present .and. (fwhm /= fwhm2 .or. trim(beam_file) /= trim(beam_file2))) then
        write(*,*) 'Warning: the 2 beams are different'
        if (fwhm /= fwhm2) print*,fwhm,fwhm2
        if (trim(beam_file) /= trim(beam_file2)) print*,trim(beam_file),trim(beam_file2)
     endif
  endif

  if (correct_beam) then
     if (beam_present) then
        if (megaverbose) write(*,*) 'Read input beam'
        call my_generate_beam(nlmax, gb, megaverbose, beam_file)
     else
        call generate_gaussian_beam(fwhm, nlmax, gb)
     endif
  endif

  if (correct_beam2) then
     if (beam2_present) then
        if (megaverbose) write(*,*) 'Read input beam (2)'
        call my_generate_beam(nlmax, gb2, megaverbose, beam_file2)
     else
        call generate_gaussian_beam(fwhm2, nlmax, gb2)
     endif
  endif

end subroutine input_files

!=======================================================================
subroutine output_results
!=======================================================================
  use spice_common
  use misc_utils, only: fatal_error
#ifdef PIO
  use piolib
#else
  use spice_parameters, only : decouple, version, thetamax, &
       & apodize, apodizesigma, apodizetype, map_present, map2_present, &
       masks_present, masks2_present, weights_present, weights2_present, &
       weightpower, weightpower2, &
       weightpowerp, weightpowerp2, &
       fwhm, fwhm2, correct_beam, correct_beam2, beam_present, beam2_present,&
       beam_file, beam_file2, &
       correct_transfer_function, subtract_average
#endif
  use head_fits, only: add_card
  use misc_utils, only: assert_alloc
  implicit none
  integer(I4B) :: i,j,k,my_npixtot,my_nx,my_ny,status
  real(DP), allocatable, dimension(:) :: my_kcross

  integer(I8B) :: my_npixtot8,my_nx8,my_ny8, my_nz8, myerr,mygroup, inc !SP
  ! integer(I8B), allocatable, dimension(:) :: index1,index2 !SP
  character(len=FILENAMELEN) :: grpname='' !SP
  character(len=FILENAMELEN) :: command='' !SP
  character(len=10) :: str !SP
  character(len=80), dimension(1:120) :: header
  integer(i4b), dimension(1:3) :: dims

  if (coroutput) then
     if (megaverbose) write(*,*) 'Output correlation function file'
     call write_corcl_file(xi_final,mu,nlmax,cor_file_name,.true.,ncor,nside)
  endif

  if (cloutput) then
     if (megaverbose) write(*,*) 'Output Cl file'
!!!!     print*,size(cl_final,1),size(cl_final,2),size(mu),nlmax,ncor,nside
     call write_corcl_file(cl_final,mu,nlmax,cl_file_name,.false.,ncor,nside)
  endif

  if (clmapoutput) then
     if (megaverbose) write(*,*) 'Output raw Cls for the map'
     call write_corcl_file(cl,mu,nlmax,cl_outmap_file,.false.,ncor,nside)
  endif

  if (clmaskoutput) then
     if (megaverbose) write(*,*) 'Output raw Cls for the weights/masks'
     call write_corcl_file(cl_mask,mu,nlmax,cl_outmask_file,.false.,ncmask,nside)
  endif

  !---------------------------------------------------------------------
  if (windowoutput) then
     if (.not. allocated(kcross)) then
        print*,'Window function not computed. Skipping output'
     else
        if (megaverbose) write(*,*) 'Output the window function'
#ifdef PIO
     my_npixtot=(nlmax+1)*(2*nlmax+1)

     allocate(my_kcross(0:my_npixtot-1),stat=status)
     call assert_alloc(status,'output_results','my_kcross')
     my_npixtot8=my_npixtot
     inc=0
     do j=0,2*nlmax
        do i=0,nlmax
           my_kcross(inc)=kcross(i,j)
           inc=inc+1
        enddo
     enddo     
     myerr=PIOgetgrpname(grpname,window_out_file)
     my_nx8 = INT(nlmax+1,kind=I8B) !SP
     myerr=PIOcreatetab2dgrp(grpname,my_nx8) !SP
     myerr=PIOcreatetab2dobject(window_out_file,'PIODOUBLE') !SP
     mygroup=PIOopentab2dgrp(grpname,'w')
!!$     write(str,*) 2*nlmax !SP
!!$     command = 'tab=*,0:'//trim(adjustl(str)) !SP
     write(command,'(''tab=*,'',I6,'':'',I6)') 0,2*nlmax !SP
     myerr=PIOwritetab2dobject(my_kcross(0:my_npixtot-1), & !SP
          & window_out_file,trim(command),mygroup) !SP
     myerr=PIOclosetab2dgrp(mygroup)
     deallocate(my_kcross) !SP
#else
     my_npixtot=(nlmax+1)*(2*nlmax+1)
     my_nx=nlmax+1
     my_ny=2*nlmax+1
     allocate(my_kcross(0:my_npixtot-1),stat=status)
     call assert_alloc(status,'output_results','my_kcross')
     inc=0
     do j=0,2*nlmax
        do i=0,nlmax
           my_kcross(inc)=kcross(i,j)
           inc=inc+1
        enddo
     enddo     
     status=0
     call deletefile(window_out_file,status)
     call write_twod_fits_DP(window_out_file,my_npixtot,my_nx,my_ny,my_kcross)
     deallocate(my_kcross)
#endif
  endif
  endif

  !---------------------------------------------------------------------
  if (kernelsoutput) then
     if (megaverbose) write(*,*) 'Output the kernels '//trim(kernels_out_file)
     if (nkernels /= 1 .and. nkernels /= 4) then
        print*,nkernels
        print*,'ERROR: Expected either 1 or 4 kernels'
        call fatal_error
     endif
     dims(1:3) = (/ nlmax+1, 2*nlmax+1, nkernels /)
     my_npixtot8 = product(dims(1:3))
     allocate(my_kcross(0:my_npixtot8-1))
     inc=0
     do k=1, nkernels
        do j=0,2*nlmax
           do i=0,nlmax
              my_kcross(inc)=kernels(i,j, k)
              inc=inc+1
           enddo
        enddo
     enddo
#ifdef PIO
     myerr=PIOGetGrpname(grpname, kernels_out_file)
     my_nx8 = INT(  dims(1), kind=I8B)
     my_ny8 = INT(  dims(2), kind=I8B)
     my_nz8 = INT(  dims(3), kind=I8B)
     write(command,'(''tab=*,*,'',I6,'':'',I6)') 0, my_nz8
     myerr=PIOCreateTAB3DGrp(grpname, my_nx8, my_ny8)
     mygroup=PIOOpenTAB3DGrp(grpname,'w')
     myerr=PIOCreateTAB3DObject(kernels_out_file,'PIODOUBLE')
     myerr=PIOWriteTAB3DObject(my_kcross, kernels_out_file, trim(command), mygroup)
     myerr=PIOWriteKeywordObject(my_nx8, '1st_dim', 'nx', kernels_out_file, mygroup)
     myerr=PIOWriteKeywordObject(my_ny8, '2nd_dim', 'ny', kernels_out_file, mygroup)
     myerr=PIOWriteKeywordObject(my_nz8, '3rd_dim', 'nz', kernels_out_file, mygroup)
     myerr=PIOCloseTAB3DGrp(mygroup)
#else
     header = ''
     call add_card(header,'COMMENT','------------------------------------------------------------')
     call add_card(header,'COMMENT','          Coupling Kernels for Spice                        ')
     call add_card(header,'COMMENT','   These 4 kernels relates the average Spice C(l) estimator ')
     call add_card(header,'COMMENT','       to the underlying ''true'' CMB power spectra.        ')
     call add_card(header,'COMMENT','                                                            ')
     call add_card(header,'COMMENT',' These kernels are not required internally by Spice         ')
     call add_card(header,'COMMENT','           (except for Kern(*,*,4))                         ')
     call add_card(header,'COMMENT','  and are only provided for convenience                     ')
     call add_card(header,'COMMENT','  (eg, for cosmological interpretation of the Spice C(l))   ')
     call add_card(header,'COMMENT','                                                            ')
     call add_card(header,'COMMENT',' These kernels depend on the apodization scheme, thetamax   ')
     call add_card(header,'COMMENT','  and Lmax.                             ')
     if (nkernels == 4) then
        call add_card(header,'COMMENT',' The way to use them depends on the estimators              ')
        call add_card(header,'COMMENT','    (decoupled or not) being considered                     ')
        call add_card(header,'COMMENT','   (option ''decouple'' of Spice).                          ')
        call add_card(header,'COMMENT','                                                            ')
        call add_card(header,'COMMENT','     E-B Decoupled Estimators                               ')
        call add_card(header,'COMMENT','  <C_TT(l1)>           = Sum_l2 Kern(l1, l2, 1) C_TT(l2)_true    ')
        call add_card(header,'COMMENT','  <C_EE(l1)>           = Sum_l2 Kern(l1, l2, 3) C_EE(l2)_true    ')
        call add_card(header,'COMMENT','  <C_BB(l1)>           = Sum_l2 Kern(l1, l2, 3) C_BB(l2)_true    ')
        call add_card(header,'COMMENT','  <C_TE(l1)>           = Sum_l2 Kern(l1, l2, 4) C_TE(l2)_true    ')
        call add_card(header,'COMMENT','  <C_TB(l1)>           = Sum_l2 Kern(l1, l2, 4) C_TB(l2)_true    ')
        call add_card(header,'COMMENT','  <C_EB(l1)>           = Sum_l2 Kern(l1, l2, 3) C_EB(l2)_true    ')
        call add_card(header,'COMMENT','          with l1 in {0,lmax}, l2 in {0,2*lmax}                  ')
        call add_card(header,'COMMENT','     NB: In the decoupled case, Kern(*,*,2) is *NOT* used   ')
        call add_card(header,'COMMENT','       (see Eq. (91) of Chon et al. 2004)                   ')
        call add_card(header,'COMMENT','                                                            ')
        call add_card(header,'COMMENT','                                                            ')
        call add_card(header,'COMMENT','                                                            ')
        call add_card(header,'COMMENT','     Coupled estimators                                              ')
        call add_card(header,'COMMENT','  <C_TT(l1)>           = Sum_l2 Kern(l1, l2, 1)  C_TT(l2)_true       ')
        call add_card(header,'COMMENT','  <C_EE(l1)+C_BB(l1)>  = Sum_l2 Kern(l1, l2, 2) (C_EE(l2)+C_BB(l2))  ')
        call add_card(header,'COMMENT','  <C_EE(l1)-C_BB(l1)>  = Sum_l2 Kern(l1, l2, 3) (C_EE(l2)-C_BB(l2))  ')
        call add_card(header,'COMMENT','  <C_TE(l1)>           = Sum_l2 Kern(l1, l2, 4)  C_TE(l2)_true       ')
        call add_card(header,'COMMENT','  <C_TB(l1)>           = Sum_l2 Kern(l1, l2, 4)  C_TB(l2)_true       ')
        call add_card(header,'COMMENT','  <C_EB(l1)>           = Sum_l2 Kern(l1, l2, 3)  C_EB(l2)_true       ')
        call add_card(header,'COMMENT','          with l1 in {0,lmax}, l2 in {0,2*lmax}                      ')
        call add_card(header,'COMMENT','       (see Eqs. (56)-(60) of Chon et al. 2004)                      ')
        call add_card(header,'COMMENT','                                                                     ')
        call add_card(header,'COMMENT','                                                             ')
        call add_card(header,'COMMENT','                                                             ')
        call add_card(header,'COMMENT','   Sanity checks:                                            ')
        call add_card(header,'COMMENT','   1) The Kernels have dimension (Lmax+1)*(2*Lmax+1)*4       ')
        call add_card(header,'COMMENT','   2) The Kernels each have the form                         ')
        call add_card(header,'COMMENT','   Kern(l1,l2,i) = S(l1,l2,i) * (2*l2+1)                     ')
        call add_card(header,'COMMENT','   where S(l1,l2,i) = S(l2,l1,i) is symmetric for i=1,2,3,4  ')
        call add_card(header,'COMMENT','        (see Eqs. (58) and (62) of Chon et al. 2004)         ')
        call add_card(header,'COMMENT','                                                             ')
     else
        call add_card(header,'COMMENT','  <C_TT(l1)>           = Sum_l2 Kern(l1, l2, 1) C_TT(l2)_true    ')
        call add_card(header,'COMMENT','          with l1 in {0,lmax}, l2 in {0,2*lmax}                  ')
        call add_card(header,'COMMENT','                        ')
        call add_card(header,'COMMENT','   Sanity checks:                                            ')
        call add_card(header,'COMMENT','   1) The Kernel has dimension (Lmax+1)*(2*Lmax+1)       ')
        call add_card(header,'COMMENT','   2) The Kernel has the form                            ')
        call add_card(header,'COMMENT','   Kern(l1,l2) = S(l1,l2) * (2*l2+1)                     ')
        call add_card(header,'COMMENT','   where S(l1,l2) = S(l2,l1) is symmetric                    ')
        call add_card(header,'COMMENT','        (see Eqs. (58) and (62) of Chon et al. 2004)         ')
        call add_card(header,'COMMENT','                                                             ')
     endif
     call add_card(header,'COMMENT','                        ')
     call add_card(header,'COMMENT','                        ')
     call add_card(header,'COMMENT','                        ')
     call add_card(header,'COMMENT','                        ')
     call add_card(header,'COMMENT','                        ')
     call add_card(header,'COMMENT','  ref: Chon et al, 2004, MNRAS 350, 914                      ')
     
     call write_spice_header(header, size(header), nlmax, ncor, nside, .false.)
     status=0
     call deletefile(kernels_out_file,status)
     call write_nd_fits_DP(kernels_out_file, my_npixtot8, 3, dims, my_kcross, &
          &  header, size(header))
#endif
     deallocate(my_kcross)
  endif

  !---------------------------------------------------------------------
  if (do_cov) then
     if (megaverbose) write(*,*) 'Output the covariance matrix '//trim(cov_out_file)
     if (ncov /= 1) then
        print*,ncov
        print*,'ERROR: Expected 1'
        call fatal_error
     endif
     dims(1:3) = (/ nlmax+1, nlmax+1, ncov /)
     my_npixtot8 = product(dims(1:3))
     allocate(my_kcross(0:my_npixtot8-1))
     inc=0
     do k=1, dims(3)
        do j=0, dims(2)-1
           do i=0, dims(1)-1
              my_kcross(inc)=cov_matrix(i,j, k)
              inc=inc+1
           enddo
        enddo
     enddo
     print*,'matrix dimensions: ',dims
     print*,'matrix min and max: ', minval(my_kcross), maxval(my_kcross)
#ifdef PIO
     myerr=PIOGetGrpname(grpname, cov_out_file)
     my_nx8 = INT(  dims(1), kind=I8B)
     my_ny8 = INT(  dims(2), kind=I8B)
     my_nz8 = INT(  dims(3), kind=I8B)
     write(command,'(''tab=*,*,'',I6,'':'',I6)') 0, my_nz8-1
     myerr=PIOCreateTAB3DGrp(grpname, my_nx8, my_ny8)
     if (myerr /= 0) write(*,*), myerr, 'creategrp ',PIOerrmess(myerr)
     mygroup=PIOOpenTAB3DGrp(grpname,'w')
     myerr=PIOCreateTAB3DObject(cov_out_file,'PIODOUBLE')
     if (myerr /= 0) write(*,*), myerr, 'createobj ',PIOerrmess(myerr)
     myerr=PIOWriteTAB3DObject(my_kcross, cov_out_file, trim(command), mygroup)
     if (myerr < 0) write(*,*), myerr, 'writeobj ',PIOerrmess(myerr)
     myerr=PIOWriteKeywordObject(my_nx8, '1st_dim', 'nx', cov_out_file, mygroup)
     myerr=PIOWriteKeywordObject(my_ny8, '2nd_dim', 'ny', cov_out_file, mygroup)
     myerr=PIOWriteKeywordObject(my_nz8, '3rd_dim', 'nz', cov_out_file, mygroup)
     if (myerr /= 0) write(*,*), myerr, 'writekw ',PIOerrmess(myerr)
     myerr=PIOCloseTAB3DGrp(mygroup)
     if (myerr /= 0) write(*,*), myerr, 'closegrp ',PIOerrmess(myerr)
#else
     header = ''
     call add_card(header,'COMMENT','------------------------------------------------------------')
     call add_card(header,'COMMENT','          Covariance Matrix for Spice                        ')
     call add_card(header,'COMMENT','                        ')
     call add_card(header,'COMMENT','                        ')
     call add_card(header,'COMMENT','                        ')
     call add_card(header,'COMMENT','                        ')
     call add_card(header,'COMMENT','                        ')
     
     call write_spice_header(header, size(header), nlmax, ncor, nside, .false.)
     status=0
     call deletefile(cov_out_file,status)
     call write_nd_fits_DP(cov_out_file, my_npixtot8, 3, dims, my_kcross, &
          &  header, size(header))
#endif
     deallocate(my_kcross)
  endif

end subroutine output_results

!=======================================================================
subroutine check_header(file_name,verbose,megaverbose,mapfile, &
 &                      nside,npixtot,nlmax,ordering,dry,polarization, &
 &                      head_defined)
!=======================================================================
  use healpix_types
  use misc_utils, only: fatal_error
#ifdef PIO
  use piolib  
#else
  use fitstools, only: getsize_fits, getnumext_fits
#endif
  implicit none

  character(len=*), intent(IN)    :: file_name
  logical,          intent(IN)    :: verbose,megaverbose,mapfile
  integer(I4B),     intent(INOUT) :: nside,npixtot
  integer(I4B),     intent(OUT)   :: nlmax,ordering
  logical,          intent(IN)    :: dry,polarization
  logical,          intent(INOUT) :: head_defined
  ! if head_defined=.false. get nside, nlmax and npixtot from file, and switch head_defined to .true.
  ! if head_defined=.true. compare provided nside and npixtot to those read in file

  integer(I4B) :: nside_in,nmap_in,npixtot_in,type
  logical :: first_time=.true.
  save first_time 
  
  integer(I8B) :: mygroup,myerr
  character(len=FILENAMELEN) :: grpname='' !SP
  character(len=FILENAMELEN), dimension(:), pointer :: flgnamePIO,maptypePIO,mapnamePIO
  character(len=FILENAMELEN) :: coordsys_in='',ordering_in='' !SP
  character(len=FILENAMELEN), save :: coordsys='' !SP


#ifdef PIO
  myerr   = PIOgetgrpname(grpname,file_name)
  mygroup = PIOopenmapgrp(grpname,"r")
  nside_in= PIOnsidegrp(mygroup)
  myerr   = PIOcoordsysgrp(coordsys_in,mygroup)
  myerr   = PIOorderinggrp(ordering_in,mygroup)
  myerr   = PIOclosemapgrp(mygroup) 
  npixtot_in=12*nside_in**2
  nmap_in=1

  if (ordering_in=='RING') then
     ordering=1
  else
     ordering=2
  endif
#else  
  npixtot_in=getsize_fits(file_name,nmaps=nmap_in,ordering=ordering, &
 &                        nside=nside_in,type=type) !EH 2003-06
  if (type == 3) then
     npixtot_in = 12*nside_in**2 ! cut sky file  !EH 2003-06
     nmap_in = getnumext_fits(file_name) ! polarized cut sky map, EH 2007-04-03
  endif
#endif


  if (dry) then
     nside_in=nint(sqrt(dble(npixtot_in)/12.d0))
     if (12*nside_in**2 /= npixtot_in) then
        write(*,*) 'ERROR in check header for file'
        write(*,*) trim(file_name)
        write(*,*) 'npixtot is not a multiple of 12*nside**2'
        write(*,*) 'Even with the dry option, I REFUSE to go further.'
        CALL FATAL_ERROR
     endif
     ordering=1
     nmap_in=1
     if (mapfile .and. polarization) nmap_in=3

     if (first_time) then
        if (verbose) then
           write(*,*)
           write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           write(*,*) 'WARNING in check_header:'
           write(*,*) 'You wrote the fits file like a pig so you are using'
           write(*,*) 'the dry option at your own risk.'
           write(*,*) 'Only LOUSY consistency check is done.'
           write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           write(*,*)
        endif
        first_time=.false.
     endif
  endif

  if (megaverbose) then
     if (dry) then
        write(*,*) 'Check approximately header for file '//trim(file_name)
     else
        write(*,*) 'Check header for file '//trim(file_name)
     endif
  endif

#ifdef PIO
     coordsys=coordsys_in
#endif

  if (.not. head_defined) then
     nside=nside_in
     npixtot=12*nside**2
     call set_nlmax(nlmax,nside)
!!$#ifdef PIO !SP
!!$     coordsys=coordsys_in !SP
!!$#endif !SP
     if (megaverbose) write(*,*) 'nside for input map file =',nside
     head_defined=.true.
  endif

  if (nmap_in > 1 .and. type==2 .and. verbose .and. ((mapfile .and. .not. polarization) &
 &    .or.(.not. mapfile))) then
     write(*,*)
     write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(*,*) 'WARNING in check_header:'
     write(*,*) 'I detect multiple maps in file'
     write(*,*) trim(file_name)
     write(*,*) 'Only the first map will be read.'
     write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(*,*)
  endif

#ifndef PIO
!  if (nmap_in /= 3 .and. mapfile .and. polarization) then
  if (mapfile .and. polarization) then
     if (nmap_in < 3) then
          write(*,*) 'ERROR in check_header for file'
	  write(*,*) trim(file_name)
	  write(*,*) 'You had requested polarization analysis'
	  write(*,*) 'while the input file does not seem to contain'
	  write(*,*) 'Stokes parameters maps'
	  CALL FATAL_ERROR
     endif
     if (nmap_in > 3) then
          write(*,*) 'WARNING in check_header for file'
	  write(*,*) trim(file_name)
	  write(*,*) 'has more than 3 columns.'
	  write(*,*) 'Will assume that first 3 are Stokes parameters maps.'
     endif
  endif
#endif

  if (nside_in /= nside) then
     write(*,*) 'ERROR in check_header for file'
     write(*,*) trim(file_name)
     write(*,*) 'The input file has wrong value of nside :'
     write(*,*) 'nside_in =',nside_in
     write(*,*) 'Other files have Nside = ',nside
     CALL FATAL_ERROR
  endif

  if (npixtot_in /= npixtot) then
     write(*,*) 'ERROR in check_header for file'
     write(*,*) trim(file_name)
     write(*,*) 'mismatch between npixtot and nside'
     write(*,*) 'npixtot_in =',npixtot_in
     write(*,*) 'expected npixtot =',npixtot,' for Nside =',nside_in
     CALL FATAL_ERROR
  endif

  if (ordering < 1 .or. ordering > 2) then
     write(*,*) 'ERROR in check_header for file'
     write(*,*) trim(file_name)
     write(*,*) 'pixel ordering is not RING nor NESTED :'
     write(*,*) 'ordering =',ordering
     CALL FATAL_ERROR
  endif

  if (ordering == 2) then
     write(*,*) 'WARNING in check_header for file'
     write(*,*) trim(file_name)
     write(*,*) 'pixel ordering is NESTED :'
     write(*,*) 'data will be reordered before usage'
  endif

#ifdef PIO
  if (coordsys_in /= coordsys) then
     write(*,*) 'ERROR in check_header for file'
     write(*,*) trim(file_name)
     write(*,*) 'mismatch between coordinate systems'
     CALL FATAL_ERROR
  endif
#endif
end subroutine check_header

!=======================================================================
subroutine set_nlmax(nlmax,nside)
!=======================================================================
  use healpix_types
  implicit none
  integer(I4B), intent(inout) :: nlmax
  integer(I4B), intent(in)    :: nside

  if (nlmax > 0) then
     nlmax = min(nlmax,3*nside-1)
  else
     nlmax = 3*nside-1
  endif
end subroutine set_nlmax

!=======================================================================
subroutine check_header_asc(file_name, nside, nlmax, ncor, head_defined, &
 &                          verbose, megaverbose, npixtot)
!=======================================================================
  ! check header of ASCII file  or  FITS ASCII table containing C(l)
  use healpix_types
  use misc_utils, only: fatal_error
#ifdef PIO
  use piolib
  use my_pio_routines, only: pio_objecttype
#else
  use fitstools, only: getsize_fits
#endif
  implicit none

  character(len=*), intent(in)    :: file_name
  logical,          intent(in)    :: verbose,megaverbose
  logical,          intent(inout) :: head_defined
  integer(I4B),     intent(inout) :: nside,nlmax,ncor
  integer(I4B),     intent(out)   :: npixtot

  integer(I4B) :: lunit=10
  character(len=FILENAMELEN) :: header
  character(len=15) :: stringjunk
  integer(I4B) :: nlmax_in,ncor_in,nside_in

  integer(I4B) :: junk, ncol, polarisation, type
  character(len=8) :: starter
!  logical :: my_isfits
#ifdef PIO
  character(len=20) :: object_type
  integer(I8B) :: myerr,mygroup,nlmaxD,ncorD,nsideD
  character(len=FILENAMELEN) :: comment='',grpname='' !SP
#endif

#ifdef PIO
  myerr   = PIOgetgrpname(grpname,file_name)
  object_type = trim(pio_objecttype(file_name))
  if (trim(object_type) == 'VECT') mygroup = PIOopenVECTGrp(grpname,'r')
  if (trim(object_type) == 'CL'  ) mygroup = PIOopenClGrp  (grpname,'r')
  myerr   = PIOreadkeywordgrp(nlmaxD,comment,'nlmax',mygroup)
  myerr   = PIOreadkeywordgrp(ncorD, comment,'ncor', mygroup)
  myerr   = PIOreadkeywordgrp(nsideD,comment,'nside',mygroup)
  if (trim(object_type) == 'VECT') myerr   = PIOcloseVECTGrp(mygroup)
  if (trim(object_type) == 'CL'  ) myerr   = PIOcloseClGrp  (mygroup)
  nlmax_in = nlmaxD
  ncor_in  = ncorD
  nside_in = nsideD
#else
  if (my_isfits(file_name)) then
     junk = getsize_fits(file_name, nside=nside_in, mlpol=nlmax_in, &
          &             nmaps=ncol, polarisation=polarisation, type=type)
     if (type /= 1) then
        write(*,*) trim(file_name)//' is not an ASCII FITS file.'
        call fatal_error
     endif
     if (polarisation == 1) then
        ncor_in = ncol
        ! correlation file : one extra column
!        if (ncol == 2) ncor_in = 1 ! unpolarized
        if (ncol == 5) ncor_in = 4  ! polarized
!        if (ncol == 7) ncor_in = 6 ! polarized + B coupling
     else
        ncor_in = 1
     endif
  else
     open(file=file_name,unit=lunit,form='formatted',status='old',action='read')

     read(lunit,'(A)',err=1,end=1) header
     read(header,'(A22,3(I12))',err=1,end=1) stringjunk,nlmax_in,ncor_in,nside_in
     close(lunit)
  endif
#endif

  if (head_defined) then
     if (nside_in /= nside) then
        write(*,*) 'ERROR in check_header_asc'
        write(*,*) 'file : '//trim(file_name)
        write(*,*) 'The value of nside corresponding to this file is inconsistent'
        write(*,*) 'with the one of the input map.'
        write(*,*) 'nside of the file :',nside_in
        write(*,*) 'nside of the input map :',nside
        CALL FATAL_ERROR
     endif
     if (nlmax_in < nlmax) then
        write(*,*) 'ERROR in check_head_asc'
        write(*,*) 'file : '//trim(file_name)
        write(*,*) 'This file has the following value of nlmax :',nlmax_in
        write(*,*) 'This value is inconsistent with the current value of nlmax.'
        write(*,*) 'We indeed have nlmax =',nlmax
        CALL FATAL_ERROR
     elseif (nlmax_in > nlmax .and. verbose) then
        write(*,*) 'WARNING in check_head_asc'
        write(*,*) 'file : '//trim(file_name)
        write(*,*) 'This file has the following value of nlmax :',nlmax_in
        write(*,*) 'This value is larger than the chosen value of nlmax.'
        write(*,*) 'We indeed have nlmax =',nlmax
     endif
     ! if ncor <= 0, adapt it to the one in file
     if (ncor <= 0) ncor = ncor_in
     ! otherwise, be without mercy
     if (ncor_in < ncor) then
        if (ncor_in /= 1) then
           write(*,*) 'ERROR in check_head_asc'           
           write(*,*) 'Inconsitent header in '
           write(*,*) 'file : '//trim(file_name)
           write(*,*) 'ncor read :',ncor
           CALL FATAL_ERROR
        endif
        if (ncor == 4) then
           write(*,*) 'ERROR in check_head_asc'           
           write(*,*) 'This file corresponds to unpolarized data, while you have'
           write(*,*) 'asked for a calculation including polarization.'
           CALL FATAL_ERROR
        endif
     elseif (ncor_in > ncor) then
        if (ncor_in /= 4) then
           write(*,*) 'ERROR in check_head_asc'           
           write(*,*) 'Inconsitent header in '
           write(*,*) 'file : '//trim(file_name)
           write(*,*) 'ncor read :',ncor
           CALL FATAL_ERROR
        endif
        if (verbose) then
           write(*,*) 'WARNING in check_head_asc'
           write(*,*) 'file : '//trim(file_name)
           write(*,*) 'This file corresponds to polarized data, while you have'
           write(*,*) 'asked for a calculation without polarization.'
        endif
     endif
  else
     nside=nside_in
     npixtot=12*nside**2
     call set_nlmax(nlmax,nside)
     if (megaverbose) write(*,*) 'nside for input map file =',nside
     head_defined=.true.
  endif
  return
1 write(*,*) 'ERROR while reading/parsing '//trim(file_name)
  write(*,*) 'Is it really an ASCII file generated by Spice?'
  call fatal_error
end subroutine check_header_asc

!=======================================================================
subroutine check_open_output_files
!=======================================================================
  use spice_common
  use misc_utils, only: fatal_error
  implicit none

  integer(I4B) :: lunit=10
  character(len=filenamelen) :: xx_file_name
  character(len=10) :: status

  if (overwrite) then
     status = 'unknown'
  else
     status='new'
  endif
  if (coroutput) then
     xx_file_name = trim(adjustl(cor_file_name))
     if (xx_file_name(1:1) == '!') then
        xx_file_name=xx_file_name(2:filenamelen)
        status = 'unknown'
     endif
     open(file=xx_file_name,unit=lunit,form='formatted',status=status,err=1)
     close(lunit)
  endif

  if (cloutput) then
     xx_file_name = trim(adjustl(cl_file_name))
     if (xx_file_name(1:1) == '!') then
        xx_file_name=xx_file_name(2:filenamelen)
        status = 'unknown'
     endif
     open(file=xx_file_name,unit=lunit,form='formatted',status=status,err=2)
     close(lunit)
  endif

  if (clmapoutput) then
     xx_file_name = trim(adjustl(cl_outmap_file))
     if (xx_file_name(1:1) == '!') then
        xx_file_name=xx_file_name(2:filenamelen)
        status = 'unknown'
     endif
     open(file=xx_file_name,unit=lunit,form='formatted',status=status,err=3)
     close(lunit)        
  endif

  if (clmaskoutput) then
     xx_file_name = trim(adjustl(cl_outmask_file))
     if (xx_file_name(1:1) == '!') then
        xx_file_name=xx_file_name(2:filenamelen)
        status = 'unknown'
     endif
     open(file=xx_file_name,unit=lunit,form='formatted',status=status,err=4)
     close(lunit)        
  endif

  return

1 write(*,*) 'ERROR in check_open_output_files'
  write(*,*) 'Impossible to create file'
  write(*,*) trim(cor_file_name)
  CALL FATAL_ERROR

2 write(*,*) 'ERROR in check_open_output_files'
  write(*,*) 'Impossible to create file'
  write(*,*) trim(cl_file_name)
  CALL FATAL_ERROR

3 write(*,*) 'ERROR in check_open_output_files'
  write(*,*) 'Impossible to create file'
  write(*,*) trim(cl_outmap_file)
  CALL FATAL_ERROR

4 write(*,*) 'ERROR in check_open_output_files'
  write(*,*) 'Impossible to create file'
  write(*,*) trim(cl_outmask_file)
  CALL FATAL_ERROR
end subroutine check_open_output_files

!=======================================================================
subroutine write_corcl_file(data, mu, nlmax, file_name, cor_file, ncor, nside)
!=======================================================================
  use healpix_types
#ifdef PIO
  use piolib
  use piolib_tuple
  use spice_parameters, only : decouple, mytype_cl_tuple, mytype_vect_tuple
  use spice_common, only : out_cl_grp ! EH
  use misc_utils, only: string, assert, fatal_error
#else
  use spice_parameters, only : decouple
  use spice_common, only : fits_out !, &
       !mapfile, maskfile, weightfile, maskfilep, weightfilep, &
       !mapfile2, maskfile2, weightfile2, maskfilep2, weightfilep2, &
       !tf_file
  use fitstools, only: write_asctab
  use head_fits, only: add_card
  use misc_utils, only: fatal_error
#endif
  implicit none

  integer(I4B), intent(in) :: nlmax,ncor,nside
  real(DP), intent(in), dimension(0:nlmax,1:ncor), target :: data
!  data(0:nlmax,1:ncor),mu(nlmax+1)
  real(DP), intent(in), dimension(1:nlmax+1) :: mu
  character(len=*), intent(in) :: file_name
  logical, intent(in) :: cor_file

  integer(I4B) :: l,i,nh,iw
  integer(I4B) :: lunit=10
  character(len=FILENAMELEN) :: fn_bang, fn_nobang
!  character(len=10) :: status
  character(len=filenamelen) :: header

  character(len=80), dimension(1:120)            :: headfits
  real(DP), allocatable, dimension(:,:), target :: tmpdata

#ifdef PIO
  character(len=FILENAMELEN) :: command='',grpname='',object='' !SP
  integer(I8B) :: nlmaxD,ncorD,nsideD  !,nlmaxDp1
  integer(I8B) :: myerr,mygroup
  real(DP),allocatable,dimension(:) :: cos_tab
  character*2 :: theta_exten='TH'
  character*2 :: npairs_exten='NP'
  character*2, dimension(1:9) :: all_exten
  character(len=DMCPIOSTRINGMAXLEN) :: decoupleS

  type(PIOPtrOnArrayDble), dimension(:), pointer        :: InData 
  character(len=512) :: line
  integer(PIOSHORT) :: n2_short 
  integer(I4B) :: mytype
  character(len=8) :: stype
#endif

  !-----------------------------------------------------------------------
  if (ncor > 9) then
     print*,'ncor (number of C(l) or C(theta) to write) is larger than 9: ', ncor
     print*,'Aborting'
     call fatal_error
  endif

  fn_bang = trim(adjustl(file_name))
  if (fn_bang(1:1) == '!') then 
     fn_nobang = fn_bang(2:filenamelen)  !  no leading '!' (for PIOlib/DMC)
  else
     fn_nobang = fn_bang            !  no leading '!' (for PIOlib/DMC)
     fn_bang = '!'//trim(fn_bang)   !  with leading '!' (for FITS)
  endif

#ifdef PIO  
  nlmaxD=nlmax
  ncorD=ncor
  nsideD=nside
!  nlmaxDp1=nlmaxD+1_I8B


  if (PIOCheckTemporary(fn_nobang) == 0) then
     ! if temporary file, create default group
     print*,trim(fn_nobang)
     print*,'temporary output file. Use default group'
     grpname = out_cl_grp
  else
     ! otherwise, get group from full path name
     myerr=PIOgetgrpname(grpname,fn_nobang)
  endif
  print*, '<'//trim(grpname)//'>'


  mytype = mytype_cl_tuple
  write(command,'(''begin='',I6,'';end='',I6)') 0, nlmax

  if (mytype == mytype_cl_tuple .or. mytype == mytype_vect_tuple) then
     ! ----------------  write in Tuple C(l) ----------------
     if (mytype == mytype_cl_tuple) then
        myerr   = PIOCreateClGrp(grpname)
        mygroup = PIOOpenClGrp(grpname,'w')
        stype = 'CL'
     else
        myerr   = PIOCreateVectGrp(grpname)
        mygroup = PIOOpenVectGrp(grpname,'w')
        stype = 'VECT'
     endif
     ! group keywords
     myerr=PIOWriteKeywordGrp(nlmaxD,'Max_Multipole',                'nlmax',mygroup)
     myerr=PIOWriteKeywordGrp(ncorD, 'Number_of_Fields',             'ncor', mygroup)
     myerr=PIOWriteKeywordGrp(nsideD,'Nside_Parameter_of_Input_Map', 'nside',mygroup)
     ! object keywords
     decoupleS = 'FALSE'
     if (decouple) decoupleS = 'TRUE'
     myerr=PIOWriteKeywordObject(decoupleS  ,'EB_Decoupled','decouple', trim(fn_nobang),mygroup)

     ! prepare data
     if (cor_file) then
        print*,'write cor_file'
        ! write theta and power spectra
        allocate(InData(1:ncor+1))
        allocate(tmpdata(0:nlmax,1:1))
        tmpdata(0:nlmax,1)       = dacos(mu(1:nlmax+1)) ! theta
        InData(1)%IdxSple        => tmpdata(0:nlmax,1)
        do iw = 1, ncor
           InData(iw+1)%IdxSple  => data(0:nlmax,iw)
        enddo
        n2_short = ncor + 1
     else
        print*,'write cl_file'
        ! write power spectra
        allocate(InData(ncor))
        do iw = 1, ncor
           InData(iw)%IdxSple => data(0:nlmax, iw)
        enddo
        n2_short = ncor
     endif

     ! create tuple object
     myerr = PIOCreateTupleObject(trim(fn_nobang), 'PIODOUBLE', n2_short, mygroup, stype, stype) 
     line = 'creating tuple: '//PIOErrMess(myerr)
     call assert(myerr >= 0, line)

     ! write tuple object
     if (mytype == mytype_cl_tuple) then
        myerr = PIOWriteClTupleObject(InData, trim(fn_nobang), command, mygroup)
     else
        myerr = PIOWriteVectTupleObject(InData, trim(fn_nobang), command, mygroup)
     endif
     line = 'writing tuple: '//PIOErrMess(myerr)
     call assert(myerr >= 0, line)

     if (allocated(tmpdata)) deallocate(tmpdata)
     deallocate(InData)
     ! close group
     if (mytype == mytype_cl_tuple) then
        myerr = PIOcloseClGrp(mygroup)
     else
        myerr = PIOcloseVectGrp(mygroup)
     endif

  else
     ! ---------------- write in flat Vect -------------------
     if (decouple) then
        all_exten = (/ '  ','E ','B ','TE','TB','EB', 'ET', 'BT', 'BE' /)
     else
        all_exten = (/ '  ','Q ','U ','TQ','TU','QU', 'QT', 'UT', 'UQ' /)
     endif
     write(command,'(''begin='',I6,'';end='',I6)') 0,nlmax

     myerr=PIOcreateVECTGrp(grpname) 
     mygroup=PIOopenVECTGrp(grpname,'w')

     myerr=PIOwritekeywordgrp(nlmaxD,' ','nlmax',mygroup)
     myerr=PIOwritekeywordgrp(ncorD, ' ','ncor', mygroup)
     myerr=PIOwritekeywordgrp(nsideD,' ','nside',mygroup)

     if (cor_file) then
        print*,'write cor_file'
        ! write correlation functions
        do iw = 1, ncor
           object=trim(fn_nobang)//trim(all_exten(iw))
           myerr=PIOcreateVECTObject(object,'PIODOUBLE',mygroup)
           myerr=PIOwriteVECTObject(data(0:nlmax,iw),object,command,mygroup)
        enddo
        ! write cos(theta)
        object=trim(fn_nobang)//theta_exten
        myerr=PIOcreateVECTObject(object,'PIODOUBLE',mygroup)
        myerr=PIOwriteVECTObject(mu(1:nlmax+1),object,command,mygroup)
        ! write theta
        object=trim(fn_nobang)//npairs_exten
        myerr=PIOcreateVECTObject(object,'PIODOUBLE',mygroup)
        allocate(cos_tab(0:nlmax)) !SP
        do l=0,nlmax
           cos_tab(l)=dacos(mu(l+1))
        enddo
        myerr=PIOwriteVECTObject(cos_tab(0:nlmax),object,command,mygroup)
        deallocate(cos_tab)
     else
        print*,'write cl_file'
        ! write power spectra
        do iw = 1, ncor
           object=trim(fn_nobang)//trim(all_exten(iw))
           myerr=PIOcreateVECTObject(object,'PIODOUBLE',mygroup)
           myerr=PIOwriteVECTObject(data(0:nlmax,iw),object,command,mygroup)
        enddo
     endif
     ! close group
     myerr = PIOcloseVECTGrp(mygroup)
  endif
#else
  !----------------- Non Piolib ---------------------
  if (fits_out) then 
     ! FITS file
     ! ---- create FITS header --------
     headfits = ''
     if (cor_file) then
        call add_card(headfits,'EXTNAME',"'ANGULAR CORRELATION'")
        call add_card(headfits,'TTYPE1',"ANGLE","Angular Separation")
        call add_card(headfits,'TUNIT1',"RAD")
        call add_card(headfits,'TTYPE2',"TT","Temperature Angular Correlation")
        if (ncor > 1) then
           call add_card(headfits,'TTYPE3',"QQ")
           call add_card(headfits,'TTYPE4',"UU")
           call add_card(headfits,'TTYPE5',"TQ")
        endif
        if (ncor > 4) then
           call add_card(headfits,'TTYPE6',"TU")
           call add_card(headfits,'TTYPE7',"QU")
        endif
        if (ncor > 6) then
           call add_card(headfits,'TTYPE8',"QT")
           call add_card(headfits,'TTYPE9',"UT")
           call add_card(headfits,'TTYPE10',"UQ")
        endif
     else
        call add_card(headfits,'EXTNAME',"'ANGULAR POWER SPECTRUM'")
        call add_card(headfits,'TTYPE1',"TT")
        if (ncor > 1) then
           call add_card(headfits,'TTYPE2',"EE")
           call add_card(headfits,'TTYPE3',"BB")
           call add_card(headfits,'TTYPE4',"TE")
        endif
        if (ncor > 4) then
           call add_card(headfits,'TTYPE5',"TB")
           call add_card(headfits,'TTYPE6',"EB")
        endif
        if (ncor > 6) then
           call add_card(headfits,'TTYPE7',"ET")
           call add_card(headfits,'TTYPE8',"BT")
           call add_card(headfits,'TTYPE9',"BE")
        endif
     endif
     call write_spice_header(headfits, size(headfits), nlmax, ncor, nside, .true.)

     ! ---- write FITS file --------
     nh = size(headfits)
     if (cor_file) then
        allocate(tmpdata(0:nlmax,1:ncor+1))
        tmpdata(0:nlmax,1)        = dacos(mu(1:nlmax+1)) ! theta
        tmpdata(0:nlmax,2:ncor+1) = data(0:nlmax,1:ncor)
        call write_asctab(tmpdata, nlmax, ncor+1, headfits, nh, fn_bang)
        deallocate(tmpdata)
     else
        call write_asctab(data, nlmax, ncor, headfits, nh, fn_bang)
     endif

  else
     ! ordinary ASCII file
     open(file=fn_nobang,unit=lunit,form='formatted',status='unknown',action='write')

     write(header,'(A22,3(I12))') '# nlmax, ncor, nside =',nlmax,ncor,nside
     write(lunit,'(A)') trim(header)
     do l=0,nlmax
        if (cor_file) then
           if (ncor == 1) then
              write(lunit, '(3(E24.16))') dacos(mu(l+1)), mu(l+1), data(l,1)
           else
              write(lunit, '(11(E24.16))') dacos(mu(l+1)), mu(l+1), (data(l,i),i=1,ncor)
           endif
        else
           if (ncor == 1) then
              write(lunit,'(I5,1X,E24.16)') l, data(l,1)
           else
              write(lunit,'(I5,1X,9(E24.16))') l, (data(l,i),i=1,ncor)
           endif
        endif
     enddo
  
     close(lunit)
  endif
#endif
end subroutine write_corcl_file

!=======================================================================
subroutine read_corcl_file(data, nlmax, file_name, cor_file, ncor)
!=======================================================================
!   routine to read C(l) files, C(theta) (correlation) files
!   or transfer function files
!
!=======================================================================

  use healpix_types
  use misc_utils, only: fatal_error
#ifdef PIO
  use piolib
  use piolib_tuple
  use my_pio_routines, only: pio_objecttype
  use spice_parameters, only : decouple
#else
  use fitstools, only: fits2cl
#endif
  implicit none

  integer(I4B),     intent(in)  :: nlmax, ncor
  real(DP),         intent(out) :: data(0:nlmax,1:ncor)
  character(len=*), intent(in)  :: file_name
  logical,          intent(in)  :: cor_file

  real(DP), dimension(1:ncor) ::  data_in
  real(DP)     :: xjunk1,xjunk2
  integer(I4B) :: lin,i,nlmax_in,ncor_in,nside_in,iw
  integer(I4B) :: lunit=10
  character(len=FILENAMELEN) :: header
  character(len=15) :: stringjunk

  integer(I8B) :: lread

  character(len=80), dimension(1:60) :: headfits
  real(DP), allocatable, dimension(:,:) :: tmpdata
!  logical :: my_isfits
#ifdef PIO
  character(len=FILENAMELEN) :: object
  character(len=20) :: object_type
  integer(PIOLONG) :: myerr, mygroup
  integer(PIOSHORT) :: ncor_short
  character(len=DMCPIOSTRINGMAXLEN) :: grpname
  character(len=20) :: command
  type(PIOPtrOnArrayDble), dimension(:), pointer  :: OutData 

  real(DP),pointer,dimension(:) :: data_in_arr
  character(len=2) :: theta_exten='TH'
  character(len=2) :: npairs_exten='NP'
  character(len=2), dimension(1:9) :: all_exten
#endif

  !-----------------------------------------------------------------------
  if (ncor > 9) then
     print*,'ncor (number of C(l) or C(theta) to read) is larger than 9: ', ncor
     print*,'Aborting'
     call fatal_error
  endif

  data(0:nlmax,1:ncor)=0._DP

#ifdef PIO
  if (decouple) then
     all_exten = (/ '  ','E ','B ','TE','TB','EB', 'ET', 'BT', 'BE' /)
  else
     all_exten = (/ '  ','Q ','U ','TQ','TU','QU', 'QT', 'UT', 'UQ' /)
  endif
  command = ' '

  ! identify object type (VECT vs CL, TUPLE vs FLAT) and open it
  myerr       = PIOGetGRPName(grpname,file_name)
  object_type = trim(pio_objecttype(file_name))
  if (trim(object_type) == 'VECT') mygroup     = PIOOpenVECTGrp(grpname,'r')
  if (trim(object_type) == 'CL'  ) mygroup     = PIOOpenCLGrp(  grpname,'r')
  myerr = PIOGetNbVecTuple(ncor_short, file_name, mygroup)

  if (ncor_short > 1) then
     ! ---- read tuple CL or VECT object ----
     if (ncor /= ncor_short) then
        print*,'Was expecting ',ncor,' fields in '//trim(file_name)
        print*,'Found ',ncor_short,' in tuple.'
        call fatal_error
     endif
     if (trim(object_type) == 'CL'  ) myerr = PIOReadClTupleObject  (OutData, trim(file_name), command, mygroup)
     if (trim(object_type) == 'VECT') myerr = PIOReadVectTupleObject(OutData, trim(file_name), command, mygroup)
     lread = myerr - 1
     call check_mylread(lread,trim(file_name),nlmax)
     do iw = 1, ncor
        data(0:lread, iw) = OutData(iw)%IdxSple
     enddo
     myerr = PIODeleteTupleTable(OutData, mygroup)
  else
     ! ---- read flat CL or VECT object ----
     do iw = 1, ncor
        object=trim(file_name)//trim(all_exten(iw))
        if (trim(object_type) == 'VECT') myerr=PIOreadVECTObject(data_in_arr, object, command)
        if (trim(object_type) == 'CL')   myerr=PIOreadCLObject  (data_in_arr, object, command)
        lread=myerr-1
        call check_mylread(lread,object,nlmax)
        data(0:lread,iw)=data_in_arr(1:lread+1)
        myerr=PIOdeleteVECTTable(data_in_arr)
        !      if (object_type == 'VECT') myerr=PIOdeleteVECTTable(data_in_arr)
        !      if (object_type == 'CL'  ) myerr=PIOdeleteCLTable(  data_in_arr)
     enddo
  endif

#else
  ! Non PIOlib
  if (my_isfits(file_name)) then 
     ! FITS file
     if (cor_file) then
        allocate(tmpdata(0:nlmax, 1:ncor+1))
        call fits2cl(file_name, tmpdata, nlmax, ncor+1, headfits)
        data(0:nlmax, 1:ncor) = tmpdata(0:nlmax, 2:ncor+1)
        deallocate(tmpdata)
     else
        call fits2cl(file_name, data, nlmax, ncor, headfits)
     endif
  else
     ! ordinary ASCII file
     open(file=file_name,unit=lunit,form='formatted',status='old')
     read(lunit,'(A)',err=1,end=2) header
     read(header,'(A22,3(I12))',err=1,end=2) stringjunk,nlmax_in,ncor_in,nside_in
     
     do lread=0,nlmax
        if (cor_file) then
           read(lunit,*,err=1,end=2) xjunk1,xjunk2,(data_in(i),i=1,ncor)
        else
           read(lunit,*,err=1,end=2) lin,(data_in(i),i=1,ncor)
        endif
        data(lread,1:ncor)=data_in(1:ncor)
     enddo

     close(lunit)
  endif
#endif

  return

1 write(*,*) 'ERROR in read_corcl_file'
  write(*,*) 'file '//trim(file_name)
  write(*,*) 'seems to be corrupted.'
  close(lunit)
  CALL FATAL_ERROR

2 write(*,*) 'ERROR in read_corcl_file'
  write(*,*) 'file '//trim(file_name)
  write(*,*) 'seems to be truncated.'
  close(lunit)
  CALL FATAL_ERROR

end subroutine read_corcl_file

!=======================================================================
subroutine check_mylread(lread,file_name,nlmax, warn_only)
!=======================================================================
  ! checks number of multipoles read from PIOLIB/DMC C(l) file
  !---------------------------------------------------------------------
  use healpix_types
  use misc_utils, only: fatal_error
  implicit none
  integer(I8B) :: lread
  integer(I4B) :: nlmax
  character(len=*) :: file_name
  logical(LGT), optional :: warn_only
  logical(LGT) :: crash

  crash = .true.
  if (present(warn_only)) crash = .not. warn_only

  if (lread < 0) then
     write(*,*) 'ERROR while reading C(l)/correlation/beam file: '
     write(*,*) 'file '//trim(file_name)
     write(*,*) 'This file is not valid'
     CALL FATAL_ERROR
  elseif (lread > nlmax) then
     if (crash) then
        write(*,*) 'ERROR while reading C(l)/correlation/beam file: '
     else
        write(*,*) 'WARNING while reading C(l)/correlation/beam file: '
     endif
     write(*,*) 'file '//trim(file_name)
     write(*,'(i12,a3,i12)') lread, ' > ',nlmax
     write(*,*) 'This file corresponds to a larger nlmax'
     if (crash) CALL FATAL_ERROR
  endif
end subroutine check_mylread

!=======================================================================
subroutine generate_gaussian_beam(fwhm,nlmax,gb)
!=======================================================================
  ! this routine was moved from deal_with_xi_and_cl.f90 to here, for clarity
  ! 2013-03-28: also returns polarized B(l)
  !======================================================================
  use healpix_types
  implicit none
  real(DP),     intent(in) :: fwhm
  integer(I4B), intent(in) :: nlmax
  real(DP), dimension(0:,1:), intent(out) :: gb

  real(DP) :: arcmin2rad,sigma2fwhm,sigma, fact_pol
  integer(I4B) :: l, nd
  !======================================================================

  nd   = size(gb,2)
  arcmin2rad = PI / (180.0_dp * 60.0_dp)
  sigma2fwhm = sqrt(8.0_dp * log(2.0_dp))

  sigma    = fwhm * arcmin2rad / sigma2fwhm ! in radians
  fact_pol = exp(2.0_dp*sigma**2) ! correction for polarised fields

  ! temperature
  do l=0,nlmax
     gb(l,1) = exp(-0.5_dp * l*(l+1.0_dp) * sigma**2)
  enddo    
  ! electric or gradient
  if (nd > 1) gb(0:nlmax,2) = gb(0:nlmax,1) * fact_pol
  ! magnetic or curl
  if (nd > 2) gb(0:nlmax,3) = gb(0:nlmax,1) * fact_pol

end subroutine generate_gaussian_beam

!=======================================================================
subroutine my_generate_beam(lmax, gb, megaverbose, beam_file)
!=======================================================================
! Hacked from Healpix 1.2, modified to take into account verbose mode
! This routine should be cleaned up.
! adapted to Healpix 2.0, can read DMC file
! 2013-03-28: adapted to Healpix 3.0
!    bref is not the same as in generate_beam of Healpix 3.10
!==========================================================================
  use healpix_types
  use misc_utils, only: fatal_error
#ifdef PIO
  use piolib
  use my_pio_routines, only: pio_objecttype
#else
  use fitstools, only : fits2cl, getsize_fits
#endif
  implicit none
  real(kind=DP), dimension(0:,1:), intent(out) :: gb
  integer(kind=I4B), intent(in) :: lmax
  character(len=*),  intent(in) :: beam_file
  logical(LGT),      intent(in) :: megaverbose
  
  integer(kind=i4b) :: nsize, type, nlheader, nl, nd, lunit, il, i, junk
  character(len=80), dimension(1:180) :: header
  character(len=1600) :: str
  character(len=80) :: card
!  logical :: my_isfits

  real(DP),pointer,dimension(:) :: tmpvectDP
  integer(I8B) :: myerr
  character(len=20) :: object_type
  integer(i4b) :: l100
  real(DP):: bref
  !==========================================================================
  ! test if name of external is given and valid

  nl = size(gb, 1)
  nd = size(gb, 2)
  gb = 0.0_dp
  
  if (nl <= lmax.and.megaverbose) then
     write(*,*) 'WARNING in Generate_beam:'
     write(*,*) 'beam array only available up to ',nl
  endif
  nl = min(nl, lmax+1)
  
#ifdef PIO
 ! ------- piolib case ------------
  myerr = -1
  object_type = trim(pio_objecttype(beam_file))
  if (object_type == 'VECT') myerr=PIOreadVECTObject(tmpvectDP,beam_file,' ')
  if (object_type == 'CL')   myerr=PIOreadCLObject  (tmpvectDP,beam_file,' ')
  call check_mylread(myerr, beam_file, nl, warn_only=.true.)
  gb(0:nl-1,1) = tmpvectDP(1:nl)
  myerr=PIODeleteVECTTable(tmpvectDP)
!   if (object_type == 'VECT') myerr=PIOdeleteVECTTable(tmpvectDP)
!   if (object_type == 'CL'  ) myerr=PIOdeleteCLTable(  tmpvectDP)
#else
 !------------ Non piolib case ----------
  lunit = 15
  
  ! read file according to its type
  if (.not.my_isfits(beam_file)) then 
     ! ordinary ascii file ?
     lunit = 32
     open(unit=lunit,file=beam_file,status='old', &
          &          form='formatted',action='read')
     do
        read(lunit,'(a)', end=100, err=100) str
        if (str(1:1) /= '#') then
           read(str,*) il, gb(il,1)
           if (il == nl-1) exit
        endif
     enddo
100  continue
     close(lunit)
     if (il < (nl-1).and.megaverbose) then
        write(*,*) 'WARNING in generate beam :'
        write(*,*) 'Beam transfer function only available up to l= ',il
        write(*,*) 'The larger multipoles will be set to 0'
     endif
       
  else 
     junk = getsize_fits(beam_file,type=type)
     if (type == 1 .or. type == 2) then ! ASCII or BINARY fits table
        ! FITS file with ascii table
        call fits2cl(beam_file, gb, nl-1, nd, header, fmissval=0.0_dp)
        ! if Grad absent, replicate Temperature; if Curl absent, replicate Grad
        l100 = min(100, nl-1)
        bref = 1.e-3*sum(abs(gb(0:l100,1)))
        do i=2,nd
           print*,sum(abs(gb(0:l100,i))) ,  bref
           if ( sum(abs(gb(0:l100,i))) <  bref ) then
              print 9002,' column #',i,' empty, fill in with column #',i-1
              gb(:,i) = gb(:,i-1)
           endif
        enddo
!         do i=2,nd
!            ! if Grad and/or Curl absent, replicate Temperature
!            if ( sum(abs(gb(:,i))) < 1.e-7 ) gb(:,i) = gb(:,1)
!         enddo
     else
        write(*,'(a)') 'ERROR in generate_beam'
        write(*,'(a)') 'the file '//trim(beam_file) &
             &                //' is of unknown type,'
        write(*,'(a)') ' or does not exist.'
        call fatal_error
     endif
  endif
#endif
 !--------- end ----------
9002 format(a,i3,a,i3)

  return
end subroutine my_generate_beam


!=======================================================================
function my_isfits(filename) result(fits)
  !---------------------------------------------------------------------
  use misc_utils, only: file_present, fatal_error

  character(len=*), intent(in) :: filename
  logical :: fits, found_unix, found_extended
  !
  integer :: lunit
  character(len=8) :: card

  found_extended = file_present(filename) ! test existence of physical or virtual file
  inquire(file=filename, exist=found_unix) ! test existence of physical file

  if (.not.found_extended) then
     call fatal_error(trim(filename)//' not found.')
  endif

  if (found_unix) then
     ! real file: try to read its first line
     lunit = 100
     open(unit=lunit, file=filename, status='old', form='formatted', action='read')
     read(lunit,'(a)') card
     close(lunit)

     card = adjustl(card)
     fits = (card(1:8) == 'SIMPLE  ' .OR. card(1:8) == 'XTENSION')
  else
     ! virtual file: assume it to be FITS
     fits = .true.
  endif

  return
end function my_isfits
!=======================================================================
subroutine get_pixwin_filename(pixwin_file_name_default, &
 &                             pixwin_file_name_input, &
 &                             pixwin_file_name, &
 &                             default,nside,healpix_data,stringlength)
!=======================================================================
! pixwin_file_name_default : default name for the pixel window file
! pixwin_file_name_input  : pixel window file if default is not taken
! default     (T/F) : whether or not default name is taken for pixel
!                     window file.
!=======================================================================
  use healpix_types
  use misc_utils, only: fatal_error
  implicit none

  integer(I4B), intent(in) :: stringlength
  character(len=*), intent(in) :: pixwin_file_name_default
  character(len=*), intent(in) :: pixwin_file_name_input
  character(len=*), intent(in) :: healpix_data
  character(len=*), intent(out) :: pixwin_file_name
  logical :: default
  integer(I4B) :: nside

  character(len=8) :: string='00000000'
  integer(I4B) :: leng
  logical :: exist

! Default name for input window file
  if (default) then
     call convert_to_ascii(string,nside)
     leng=len(trim(healpix_data)//trim(pixwin_file_name_default))
     if (leng > stringlength-9) then
        write(*,*) 'ERROR in get_pixwin_filename :'
        write(*,*) 'string pixwin_file_name_default is too long.'
        CALL FATAL_ERROR
     endif
     pixwin_file_name=trim(healpix_data)//trim(pixwin_file_name_default) &
 &                    //string(5:8)//'.fits'
! Customized file name
  else
     leng=len(trim(pixwin_file_name_input))
     if (leng > stringlength) then
        write(*,*) 'ERROR in get_pixwin_filename :'
        write(*,*) 'string pixwin_file_name_input is too long.'
        CALL FATAL_ERROR
     endif
     pixwin_file_name=trim(pixwin_file_name_input)
  endif

  inquire(file=pixwin_file_name,exist=exist)
  if (.not.exist) then
     write(*,*) 'ERROR in get_pixwin_filename'
     write(*,*) 'INPUT file : '//trim(pixwin_file_name)
     write(*,*) 'does not exist.'
     CALL FATAL_ERROR
  endif
end subroutine get_pixwin_filename

!===============================================================================
subroutine read_twod_fits_DP(file_name,npixtot,nx,ny,map)
!===============================================================================
  use healpix_types
  use misc_utils, only: fatal_error
  implicit none
  character*(*), intent(in) :: file_name
  integer(i4b),  intent(in) :: npixtot
  integer(i4b),  intent(in) :: nx,ny ! use-less
  real(DP),      intent(out) :: map(0:npixtot-1)

  integer(I4B) :: naxes(2)
  integer(I4B) :: status,unit,readwrite,blocksize,nfound,nbuffer
  integer(I4B) :: group,firstpix,naxis1,naxis2
  logical :: anynull
  real(DP) :: rnullval
  integer(I4B) :: i,j,inc
  character*30 :: errtext

  status=0
  call ftgiou(unit,status)
  readwrite=0

  call ftopen(unit,file_name,readwrite,blocksize,status)
  call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)
  
  naxis1 = naxes(1)
  naxis2 = naxes(2)

  nbuffer=naxis1*naxis2
  if (nbuffer > npixtot) then
     print*,'array too small in read_twod_fits'
     print*,nbuffer,npixtot
  endif
  group=1
  firstpix=1
  rnullval=max_sp

  call ftgpvd(unit,group,firstpix,nbuffer,rnullval,map,anynull,status)

  call ftgerr(status,errtext)
  call ftclos(unit, status)
  call ftfiou(unit, status)
  if (status > 0) then
     write(*,*) 'ERROR in read_twod_fits :',errtext
     CALL FATAL_ERROR
  endif

end subroutine read_twod_fits_DP


!===============================================================================
subroutine write_twod_fits_DP(file_name,npixtot,nx,ny,map)
!===============================================================================
  use healpix_types
  use misc_utils, only: fatal_error
  implicit none
  character*(*), intent(in) :: file_name
  integer(I4B),  intent(in) :: npixtot,nx,ny
  real(DP),      intent(in) :: map(0:npixtot-1)

  integer(I4B) :: status,unit,blocksize,bitpix,naxis,naxes(2)
  integer(I4B) :: i,j,group,fpixel,nelements
  logical :: simple,extend
  character*30 :: errtext

  status=0
  call ftgiou(unit,status)
  blocksize=1
  call ftinit(unit,file_name,blocksize,status)
  simple=.true.
  bitpix=-64
  naxis=2
  naxes(1)=nx
  naxes(2)=ny
  extend=.true.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  group=1
  fpixel=1
  nelements=npixtot
  call ftpprd(unit,group,fpixel,nelements,map,status)
  call ftclos(unit, status)
  call ftfiou(unit, status)
  if (status > 0) then
     write(*,*) 'ERROR in write_twod_fits :',errtext
     CALL FATAL_ERROR
  endif

end subroutine write_twod_fits_DP

!===============================================================================
subroutine write_nd_fits_DP(file_name, npixtot, naxes, naxis, map, header, nlheader)
!===============================================================================
  use healpix_types
  use misc_utils, only: fatal_error
  use fitstools, only: printerror, putrec
  implicit none
  character(len=*), intent(in) :: file_name
  integer(I8B),     intent(in) :: npixtot
  integer(I4B),     intent(in) :: naxes
  integer(I4B),     intent(in) :: naxis(1:naxes)
  real(DP),         intent(in) :: map(0:npixtot-1)
  character(len=80),intent(in) :: header(1:nlheader)
  integer(i4b),     intent(in) :: nlheader

  integer(I4B) :: status, unit, blocksize, bitpix
  integer(I4B) :: i, j, group, fpixel, nelements
!  logical      :: simple, extend
!  character*30 :: errtext

  ! open file
  status=0
  call ftgiou(unit,status)
  call printerror(status)
  blocksize=1
  call ftinit(unit, file_name, blocksize, status)
  call printerror(status)

  if (npixtot /= product(naxis(1:naxes)) ) then
     print*,'Error ',naxes, npixtot, naxis(1:naxes)
     call fatal_error
  endif
  ! write minimal header for image
!  simple=.true.
  bitpix=-64
!  extend=.true.
  call ftphps(unit,         bitpix, naxes, naxis,                        status)
!  call ftphpr(unit, simple, bitpix, naxes, naxis(1:naxes), 0, 1, extend, status)
  call printerror(status)

  ! write user provided header
  do i=1, nlheader
     if (trim(header(i)) /= '') call putrec(unit, header(i), status)
  enddo
  call printerror(status)

  ! write data for image
  group=1
  fpixel=1
  nelements=npixtot
  call ftpprd(unit, group, fpixel, nelements, map, status)
  call printerror(status)

  ! close and exit
  call ftclos(unit, status)
  call ftfiou(unit, status)
  call printerror(status)

end subroutine write_nd_fits_DP

!===============================================================================
subroutine check_twod_fits_header(file_name,npixtot,nx,ny)
!===============================================================================
  use healpix_types
  use misc_utils, only: fatal_error
  implicit none
  character*(*), intent(IN)  :: file_name
  integer(I4B),  intent(OUT) :: npixtot,nx,ny
  
  logical(LGT) :: ultraverbose=.false.
  integer(I4B) :: naxes(2)
  integer(I4B) :: status,unit,readwrite,blocksize,nfound
   
  status=0
  call ftgiou(unit,status)
  readwrite=0
  call ftopen(unit,file_name,readwrite,blocksize,status)
  call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)
  
  if (nfound .ne. 2) then
     write(*,*) 'ERROR in check_twod_fits_header :'
     write(*,*) 'Failure in reading the NAXIS keywords'
     CALL FATAL_ERROR
  endif

  nx=naxes(1)
  ny=naxes(2)
  npixtot=nx*ny

  if (ultraverbose) write(*,*) 'Dimensions of the image detected ',nx,ny
  
  call ftclos(unit, status)
  call ftfiou(unit, status)

end subroutine check_twod_fits_header

!===============================================================================
subroutine check_nd_fits_header(file_name, npixtot, naxes, naxis)
!===============================================================================
  ! determine image dimension, EH, 2009-01-27
  !
  use healpix_types
  use misc_utils, only: fatal_error
  implicit none
  character(len=*), intent(IN)  :: file_name
  integer(I4B),  intent(OUT) :: npixtot ! number of pixels
  integer(I4B),  intent(OUT) :: naxes   ! number of dimensions
  integer(I4B),  intent(OUT), dimension(1: ) :: naxis ! size in each dimension
  
  logical(LGT) :: ultraverbose=.false.
  integer(I4B) :: status, unit, readwrite, blocksize, nfound
  character(len=80) :: comment
   
  status=0
  call ftgiou(unit,status)
  readwrite=0
  call ftopen(unit, file_name, readwrite, blocksize, status)

  ! number of dimension
  call ftgkyj(unit,'NAXIS', naxes, comment, status)
  if (size(naxis) < naxes) then
     write(*,*) 'ERROR in check_nd_fits_header :'
     write(*,*) 'Found ',naxes,' dimension'
     write(*,*) 'Receiving array has size ',size(naxis)
     call fatal_error
  endif

  ! size in each dimension
  naxis(1:) = 1
  call ftgknj(unit,'NAXIS', 1, naxes, naxis, nfound, status)
  if (nfound .ne. naxes) then
     write(*,*) 'ERROR in check_nd_fits_header :'
     write(*,*) 'Failure in reading the NAXIS keywords'
     CALL FATAL_ERROR
  endif

  ! total number of elements
  npixtot = product(naxis(1:naxes))

  if (ultraverbose) write(*,*) 'Dimensions of the image detected ',naxis(1:naxes)
  
  call ftclos(unit, status)
  call ftfiou(unit, status)

end subroutine check_nd_fits_header

!===============================================================================
subroutine deletefile(filename,status)
!===============================================================================
! Extracted from cookbood.f of the cfitsio library
! A simple little routine to delete a FITS file
!===============================================================================
  use healpix_types
  integer(I4B) :: status,unit,blocksize
  character*(*) :: filename

! Simply return if status is greater than zero
  if (status .gt. 0)return

! Get an unused Logical Unit Number to use to open the FITS file
  call ftgiou(unit,status)

! Try to open the file, to see if it exists
  call ftopen(unit,filename,1,blocksize,status)

  if (status .eq. 0)then
!    file was opened;  so now delete it 
     call ftdelt(unit,status)
  else if (status .eq. 103)then
!    file doesn't exist, so just reset status to zero and clear errors
     status=0
     call ftcmsg
  else
!    there was some other error opening the file; delete the file anyway
     status=0
     call ftcmsg
     call ftdelt(unit,status)
  end if

!  Free the unit number for later reuse
  call ftfiou(unit, status)
end subroutine deletefile

! !=======================================================================
! function my_getnumext_fits(filename)
!   !=======================================================================
!   ! patch to bugged Healpix 2.01 getnumext_fits
!   !
!   !  result = my_getnumext_fits(filename)
!   !    returns the number of extensions present in FITS file 'filename'
!   !
!   ! EH, Nov 2004
!   ! April 2007: close file on exit
!     !=======================================================================
!   use healpix_types
!   implicit none
!   character(LEN=*), intent(IN)             :: filename
!   integer(i4b)                             :: my_getnumext_fits
!   !
!   integer(i4b) :: status, unit, readwrite, blocksize, nhdu
!   !-----------------------------------------------------------------------
!   status         = 0
!   unit           = 149
!   my_getnumext_fits = 0
!   readwrite      = 0 ! Read only
!   call ftopen(unit, filename, readwrite, blocksize, status)
!   if (status > 0) then
!      !       call printerror(status)
!      print*,'FITS Error: ',status
!      call ftclos(unit, status)
!      return
!   endif
    
!   call ftthdu(unit, nhdu, status)
!   my_getnumext_fits = nhdu - 1

!   call ftclos(unit, status)
!   return
! end function my_getnumext_fits

!====================================================================
subroutine write_spice_header(headfits, nlh, nlmax, ncor, nside, aboutmap)
  !====================================================================
  use healpix_types
#ifndef PIO
  use spice_parameters, only : decouple, version, thetamax, &
       & apodize, apodizesigma, apodizetype, map_present, map2_present, &
       masks_present, masks2_present, weights_present, weights2_present, &
       weightpower, weightpower2, &
       weightpowerp, weightpowerp2, &
       fwhm, fwhm2, correct_beam, correct_beam2, beam_present, beam2_present,&
       beam_file, beam_file2, &
       correct_transfer_function, subtract_average, &
       subtract_dipole, &
       pairsthresholding, npairsthreshold, tolerance
  use spice_common, only : fits_out, &
       mapfile, maskfile, weightfile, maskfilep, weightfilep, &
       mapfile2, maskfile2, weightfile2, maskfilep2, weightfilep2, &
       tf_file, listmapfile, listmapw8
  use head_fits, only: add_card
  implicit none
  character(len=80), intent(inout), dimension(1:nlh) :: headfits
  integer(I4B),     intent(in)                  :: nlh, nlmax, ncor, nside
  logical(LGT),     intent(in)                  :: aboutmap

  call add_card(headfits)
  call add_card(headfits,'COMMENT','---------------------------------')
  call add_card(headfits)
  call add_card(headfits,"CREATOR","Spice", "Software creating the FITS file")
  call add_card(headfits,"VERSION",version, "Version of the simulation software")
  call add_card(headfits)
  call add_card(headfits,"POLAR",  (ncor>1),"Polarisation included (True/False)")
  call add_card(headfits,"BCROSS", (ncor>4),"Magnetic cross terms included (True/False)")
  call add_card(headfits,"ASYMCL", (ncor>6),"Asymmetric pol cross terms (XY vs YX) included")
  call add_card(headfits,'APODIZE',apodize, "Apodization of Xi  (True/False)")
  call add_card(headfits,'DECOUPLE',decouple,"Decouple E and B polarization (True/False)")
  call add_card(headfits,'SUBAV',subtract_average .or. subtract_dipole &
       & ,"Remove average from T map (True/False)")
  call add_card(headfits,'SUBDIPOL',subtract_dipole,"Remove dipole from T map (True/False)")
  call add_card(headfits,'TF_CORRE',correct_transfer_function,'Tranfer funct. correction (True/False)')
  call add_card(headfits)
  call add_card(headfits,'MAX-LPOL',nlmax,  "Maximum L multipole order")
  call add_card(headfits,'NCOR',   ncor,    "Number of fields (1 or 4)")
  call add_card(headfits,'NSIDE',  nside,   "Resolution Parameter of Input Map")
  call add_card(headfits,'THETAMAX',thetamax, "[Deg] Largest lag in angul. correl. fct Xi")
  call add_card(headfits,'TOLERANC',tolerance, "Relative tolerance on E/B decoupling integrals")
  if (apodize) then
     call add_card(headfits,'APOTYPE',apodizetype, "Apodization type: 0=Gaussian, 1=Cosine")
     call add_card(headfits,'APOSIGMA',apodizesigma, "[Deg] Apodization parameter")
  endif
  if (pairsthresholding) then
     call add_card(headfits,'NPTHRESH',npairsthreshold,"Minimal # of pixel pairs used in correlation")
  endif
  if (aboutmap) then
     if (map_present)  call add_card(headfits,'MAPFILE1',trim(listmapfile(1,1,1)),'Input map file 1')
     if (map2_present) call add_card(headfits,'MAPFILE2',trim(listmapfile(1,1,2)),'Input map file 2')
  endif
  if (masks_present) then
     call add_card(headfits,'MASKFIL1',trim(maskfile),'Input mask file 1')
     if (maskfilep /= maskfile) &
          & call add_card(headfits,'MASKF_P1',trim(maskfilep),'Input pol mask file 1')
  endif
  if (masks2_present) then
     call add_card(headfits,'MASKFIL2',trim(maskfile2),'Input mask file 2')
     if (maskfilep2 /= maskfile2) &
          & call add_card(headfits,'MASKF_P2',trim(maskfilep2),'Input pol mask file 2')
  endif
  if (weights_present) then
     call add_card(headfits,'W8FILE1',trim(weightfile),'Input weight file 1')
     call add_card(headfits,'W8POWER1',weightpower,'power index applied to weight file data 1')
     if (weightpower /= weightpowerp .or. weightfilep /= weightfile) then
        call add_card(headfits,'W8FIL_P1',trim(weightfilep),'Input pol weight file 1')
        call add_card(headfits,'W8POW_P1',weightpowerp,'power index applied to weight pol file data 1')
     endif
  endif
  if (weights2_present) then
     call add_card(headfits,'W8FILE2',trim(weightfile2),'Input weight file 2')
     call add_card(headfits,'W8POWER2',weightpower2,'power index applied to weight file data 2')
     if (weightpower2 /= weightpowerp2 .or. weightfilep2 /= weightfile2) then
        call add_card(headfits,'W8FIL_P2',trim(weightfilep),'Input pol weight file 2')
        call add_card(headfits,'W8POW_P2',weightpowerp2,'power index applied to weight pol file data 2')
     endif
  endif
  if (correct_beam) then
     if (beam_present) then 
        call add_card(headfits,'BEAMFIL1',trim(beam_file),'Beam file 1')
     else
        call add_card(headfits,'FWHM1',fwhm/60.0_dp,'[Deg] beam FWHM 1')
     endif
  endif
  if (correct_beam2) then
     if (beam2_present) then 
        call add_card(headfits,'BEAMFIL2',trim(beam_file2),'Beam file 2')
     else
        call add_card(headfits,'FWHM2',fwhm2/60.0_dp,'[Deg] beam FWHM 2')
     endif
  endif
  if (correct_transfer_function) then
     call add_card(headfits,'TF_FILE',trim(tf_file),'Transfer function file')
  endif
  ! missing:  maskfilep, weightfilep, ...
  call add_card(headfits)
  call add_card(headfits,'COMMENT','---------------------------------')
  call add_card(headfits)
#endif

  return
end subroutine write_spice_header

end module deal_with_files
