module deal_with_xi_and_cl
  public :: do_cl_with_alm
  public :: do_xi_from_cl, do_cl_from_xi
  public :: correct_final_cl, correct_final_xi, correct_xi_from_mask
  public :: alm2cl_spice
  private

contains

  !========================================================
  subroutine sub_alm2cl(alm1, i1, alm2, i2, cl, i3)
    !========================================================
    use healpix_types
    use spice_parameters, only: KMAPC
    use misc_utils, only: fatal_error
    integer(I4B),                        intent(in) :: i1, i2, i3
    complex(KMAPC), dimension(1:,0:,0:), intent(in) :: alm1, alm2
    real(DP),       dimension(0:,1:),    intent(out):: cl

    integer(I4B) :: nlmax, nmmax, l, mm
    complex(DPC) :: dc
    real(DP), parameter :: two = 2.000000000000000000_dp
    real(DP), parameter :: one = 1.000000000000000000_dp

    nlmax = min( size(alm1, 2), size(alm2, 2), size(cl, 1)) - 1
    nmmax = min( size(alm1, 3), size(alm2, 3)) - 1
    j1 = size(alm1, 1)
    j2 = size(alm2, 1)
    j3 = size(cl,   2)

    if (i1 > j1 .or. i2 > j2 .or. i3 > j3) then
       call fatal_error('invalid index in alm -> C(l)')
    endif
    do l = 0, nlmax
       mm = min(l, nmmax)
       dc =          sum(      alm1(i1,l,1:mm)*conjg(alm2(i2,l,1:mm)))
       dc = (dc + conjg(dc)) + alm1(i1,l,0)   *      alm2(i2,l,0)
       cl(l,i3) = real(dc, kind=DP) / (two*l + one)
    enddo

    return
  end subroutine sub_alm2cl

  !========================================================
  subroutine alm2cl_spice(nlmax, nmmax, alm1, alm2, cl, symmetric)
  !========================================================
    ! adapted from Healpix alm2cl with SPC alm and DP cl
    !
    ! computes C(l) from a_lm, in the order
    ! TT, [EE, TE, BB, [TB, EB, [ET, BT, BE]]]
    !
    ! TE= alm1_T * alm2_E 
    ! unless symmetric is set : TE = (alm1_T * alm2_E + alm1_E * alm2_T)/2
    !=======================================================
    use healpix_types
    use spice_parameters, only: KMAPC
    implicit none
    integer(I4B),                        intent(in) :: nlmax, nmmax
    complex(KMAPC), dimension(1:,0:,0:), intent(in) :: alm1, alm2
    real(DP)    ,   dimension(0:, 1: ),  intent(out):: cl
    logical(LGT),   optional,            intent(in):: symmetric
    ! 
    integer(I4B) :: l, ncl, na1, na2, mm, k1, k2
    !complex(DPC)     :: dc, dcs
    !real(DP), parameter :: two = 2.000000000000000000_dp
    !real(DP), parameter :: one = 1.000000000000000000_dp
    real(DP), parameter :: half = 0.500000000000000000_dp
    logical(LGT) :: polarisation, bcoupling, do_sym, asympol

    real(DP), allocatable, dimension(:,:) :: cl_work
    
    !========================================================

    ncl = size(cl, 2)
    na1 = size(alm1, 1)
    na2 = size(alm2, 1)
    polarisation = (na1 >= 3 .and. na2 >= 3 .and. ncl >=4)
    bcoupling    = (ncl >=6) .and. polarisation
    asympol      = (ncl >=9) .and. polarisation
    do_sym = .false.
    cl = 0.0_DP
    if (present(symmetric)) do_sym = symmetric
    if (polarisation .and. do_sym) then
       print*,'Symmetric TE C(l)'
    endif
    if (polarisation) then
       allocate(cl_work(0:nlmax, 1:2))
    endif

     ! TT power spectrum
    call sub_alm2cl(alm1, 1, alm2, 1, cl, 1)
    if (polarisation) then
       ! GG or EE power spectrum
       call sub_alm2cl(alm1, 2, alm2, 2, cl, 2)
       ! CC or BB power spectrum
       call sub_alm2cl(alm1, 3, alm2, 3, cl, 3)

       ! TG or TE power spectrum
       call sub_alm2cl(alm1, 1, alm2, 2, cl_work, 1)
       call sub_alm2cl(alm1, 2, alm2, 1, cl_work, 2)
       k1 = 4  ;   k2 = k1 +3
                    cl(0:,k1) =  cl_work(0:,1)
       if (asympol) cl(0:,k2) =                cl_work(0:,2)
       if (do_sym ) cl(0:,k1) = (cl_work(0:,1)+cl_work(0:,2)) * half

       if (bcoupling) then
          ! TC or TB power spectrum
          call sub_alm2cl(alm1, 1, alm2, 3, cl_work, 1)
          call sub_alm2cl(alm1, 3, alm2, 1, cl_work, 2)
          k1 = 5  ;   k2 = k1 +3
                       cl(0:,k1) =  cl_work(0:,1)
          if (asympol) cl(0:,k2) =                cl_work(0:,2)
          if (do_sym ) cl(0:,k1) = (cl_work(0:,1)+cl_work(0:,2)) * half


          ! GC or EB power spectrum
          call sub_alm2cl(alm1, 2, alm2, 3, cl_work, 1)
          call sub_alm2cl(alm1, 3, alm2, 2, cl_work, 2)
          k1 = 6  ;   k2 = k1 +3
                       cl(0:,k1) =  cl_work(0:,1)
          if (asympol) cl(0:,k2) =                cl_work(0:,2)
          if (do_sym ) cl(0:,k1) = (cl_work(0:,1)+cl_work(0:,2)) * half

       endif
    endif

    if (allocated(cl_work)) deallocate(cl_work)

!     ! TT power spectrum
!     do l = 0, nlmax
!        mm = min(l, nmmax)
!        dc =          sum(      alm1(1,l,1:mm)*conjg(alm2(1,l,1:mm)))
!        dc = (dc + conjg(dc)) + alm1(1,l,0)   *      alm2(1,l,0)
!        cl(l,1) = real(dc, kind=DP) / (two*l + one)
!     enddo

!     if (polarisation) then
!        ! GG or EE power spectrum
!        do l = 0, nlmax
!           mm = min(l, nmmax)
!           dc =          sum(      alm1(2,l,1:mm)*conjg(alm2(2,l,1:mm)))
!           dc = (dc + conjg(dc)) + alm1(2,l,0)   *      alm2(2,l,0)
!           cl(l,2) = real(dc, kind=DP) / (two*l + one)
!        enddo
       
!        ! CC or BB power spectrum
!        do l = 0, nlmax
!           mm = min(l, nmmax)
!           dc =          sum(      alm1(3,l,1:mm)*conjg(alm2(3,l,1:mm)))
!           dc = (dc + conjg(dc)) + alm1(3,l,0)   *      alm2(3,l,0)
!           cl(l,3) = real(dc, kind=DP) / (two*l + one)
!        enddo
       
!        ! TG or TE power spectrum
!        do l = 0, nlmax
!           mm = min(l, nmmax)
!           dc =          sum(      alm1(1,l,1:mm)*conjg(alm2(2,l,1:mm)))
!           dc = (dc + conjg(dc)) + alm1(1,l,0)   *      alm2(2,l,0)
!           if (do_sym) then
!              dcs =          sum(        alm2(1,l,1:mm)*conjg(alm1(2,l,1:mm)))
!              dcs = (dcs + conjg(dcs)) + alm2(1,l,0)   *      alm1(2,l,0)
!              dc = (dc + dcs) * half
!           endif
!           cl(l,4) = real(dc, kind=DP) / (two*l + one)             
!        enddo
!     endif 

!     if (bcoupling) then
!        ! TC or TB power spectrum
!        do l = 0, nlmax
!           mm = min(l, nmmax)
!           dc =          sum(      alm1(1,l,1:mm)*conjg(alm2(3,l,1:mm)))
!           dc = (dc + conjg(dc)) + alm1(1,l,0)   *      alm2(3,l,0)
!           if (do_sym) then
!              dcs =          sum(        alm2(1,l,1:mm)*conjg(alm1(3,l,1:mm)))
!              dcs = (dcs + conjg(dcs)) + alm2(1,l,0)   *      alm1(3,l,0)
!              dc  = (dc + dcs) * half
!           endif
!           cl(l,5) = real(dc, kind=DP) / (two*l + one)
!        enddo
       
!        ! GC or EB power spectrum
!        do l = 0, nlmax
!           mm = min(l, nmmax)
!           dc =          sum(      alm1(2,l,1:mm)*conjg(alm2(3,l,1:mm)))
!           dc = (dc + conjg(dc)) + alm1(2,l,0)   *      alm2(3,l,0)
!           if (do_sym) then
!              dcs =          sum(        alm2(2,l,1:mm)*conjg(alm1(3,l,1:mm)))
!              dcs = (dcs + conjg(dcs)) + alm2(2,l,0)   *      alm1(3,l,0)
!              dc  = (dc + dcs) * half
!           endif
!           cl(l,6) = real(dc, kind=DP) / (two*l + one)
!        enddo
!     endif

!     if (asympol) then
!        ! GT or ET power spectrum
!        do l = 0, nlmax
!           mm = min(l, nmmax)
!           dc =          sum(      alm1(2,l,1:mm)*conjg(alm2(1,l,1:mm)))
!           dc = (dc + conjg(dc)) + alm1(2,l,0)   *      alm2(1,l,0)
!           cl(l,7) = real(dc, kind=DP) / (two*l + one)             
!        enddo
!        ! CT or BT power spectrum
!        do l = 0, nlmax
!           mm = min(l, nmmax)
!           dc =          sum(      alm1(3,l,1:mm)*conjg(alm2(1,l,1:mm)))
!           dc = (dc + conjg(dc)) + alm1(3,l,0)   *      alm2(1,l,0)
!           cl(l,8) = real(dc, kind=DP) / (two*l + one)             
!        enddo
!        ! CG or BE power spectrum
!        do l = 0, nlmax
!           mm = min(l, nmmax)
!           dc =          sum(      alm1(3,l,1:mm)*conjg(alm2(2,l,1:mm)))
!           dc = (dc + conjg(dc)) + alm1(3,l,0)   *      alm2(2,l,0)
!           cl(l,9) = real(dc, kind=DP) / (two*l + one)
!        enddo

!     endif

    return
  end subroutine alm2cl_spice

!=======================================================================
subroutine do_cl_with_alm(nside, npixtot, nlmax, map, cl, &
     &                        verbose, message, nmap, ncor, polarization, &
     &                        map2, cl_auto1, cl_auto2, symmetric)
!=======================================================================
!   map1 -> alm1,  [map2 -> alm2]
!   output cl = alm1*alm1  or alm1*alm2
  use healpix_types
  use spice_parameters, only: KMAP, KMAPC
  use misc_utils, only: fatal_error
  use alm_tools, only: map2alm !, alm2cl
  use fitstools, only: read_dbintab
  implicit none
    
  integer(I4B), intent(in)  :: nside,nlmax,npixtot,nmap,ncor
  real(KMAP),   intent(in)  :: map(0:npixtot-1,1:nmap)
  real(DP),     intent(out) :: cl(0:nlmax,1:ncor)
  logical,      intent(in)  :: verbose,polarization
  character(len=*), intent(in) :: message
  real(KMAP), optional,intent(in)  :: map2(0:npixtot-1,1:nmap)
  real(DP),   optional,intent(out)  :: cl_auto1(0:nlmax,1:ncor)
  real(DP),   optional,intent(out)  :: cl_auto2(0:nlmax,1:ncor)
  logical(LGT), optional,intent(in)  :: symmetric
  
  integer(I4B) :: l,m
  real(DP), allocatable, dimension(:,:) :: w8
  !complex(SPC), allocatable, dimension(:,:,:) :: alm, alm2, alm3
  complex(KMAPC), allocatable, dimension(:,:,:) :: alm, alm2, alm3
  logical(LGT) :: do_cross, simple, doublemask, do_cross_and_auto
  real(DP), dimension(0:1) :: zrange = (/ 0.0_dp, 0.0_dp /)
  logical(LGT) :: my_symmetric, asym_cl
  real(DP), allocatable, dimension(:,:) :: cl_temp

  my_symmetric = .false.
  if (present(symmetric)) my_symmetric = symmetric
  do_cross = present(map2)
  do_cross_and_auto = do_cross .and. present(cl_auto1) .and. present(cl_auto2)

  !print*,map(0,:), map(npixtot-1,:)
  !if (do_cross) print*,map2(0,:),map2(npixtot-1,:)

  if (verbose) write(*,*) 'Compute Cl '//trim(message)
  if (.not. polarization) then
     ! unpolarized (either T map,   or T or P mask)
     simple     = (nmap == 1 .and. ncor == 1)
     doublemask = (nmap == 2 .and. ncor == 3) ! dealing with masks
     asym_cl    = (nmap == 2 .and. ncor == 4) ! dealing with masks
     if (.not.(simple .or. doublemask  .or. asym_cl)) then
        write(*,*) 'do_cl_with_alm: expects nmap and ncor either 1,1 or 2,3 or 2,4.'
        write(*,*) '                got :', nmap,ncor
        call fatal_error
     endif
     allocate(w8(1:2*nside,1:1))
     w8 = 1.0_DP
     allocate(alm(1:1,0:nlmax,0:nlmax))
     alm = 0
     if (doublemask .or. asym_cl) then
        ! case of different mask for T and P
        allocate(alm2(1:1,0:nlmax,0:nlmax))
        if (do_cross) then
           if (my_symmetric .or. asym_cl) then ! 2013-10-07: symmetric mask C(l)
              if (my_symmetric) then
                 print*,'Symmetric mask'
              else
                 print*,'Asymmetric masks'
              endif
              allocate(alm3(1:1,0:nlmax,0:nlmax))
              allocate(cl_temp(0:nlmax,1:2))
              call map2alm(nside,nlmax,nlmax,map (0:,1),alm  , zrange,w8) ! T_1
              call map2alm(nside,nlmax,nlmax,map2(0:,1),alm2 , zrange,w8) ! T_2
              call map2alm(nside,nlmax,nlmax,map2(0:,2),alm3 , zrange,w8) ! P_2
              call alm2cl_spice(nlmax, nlmax, alm , alm2,  cl(:,1:1))      ! T_1 x T_2
              call alm2cl_spice(nlmax, nlmax, alm , alm3,  cl_temp(:,1:1)) ! T_1 x P_2
              call map2alm(nside,nlmax,nlmax,map (0:,2),alm   , zrange,w8) ! P_1
              call alm2cl_spice(nlmax, nlmax, alm , alm3,  cl(:,2:2))      ! P_1 x P_2
              call alm2cl_spice(nlmax, nlmax, alm , alm2,  cl_temp(:,2:2))  ! P_1 x T_2
              if (my_symmetric) then
                 cl(:,3:3) = 0.5*(cl_temp(:,1:1) + cl_temp(:,2:2))         ! ( P_1 x T_2 + T_1 x P_2 )/2
              else
                 cl(:,3:3) = cl_temp(:,1:1)
                 cl(:,4:4) = cl_temp(:,2:2)
              endif
              deallocate(cl_temp, alm3)
           else
              call map2alm(nside,nlmax,nlmax,map (:,1),alm  , zrange,w8) ! T_1
              call map2alm(nside,nlmax,nlmax,map2(:,1),alm2 , zrange,w8) ! T_2
              call alm2cl_spice(nlmax, nlmax, alm , alm2,  cl(:,1:1)) ! T_1 x T_2
              call map2alm(nside,nlmax,nlmax,map2(:,2),alm2 , zrange,w8) ! P_2
              call alm2cl_spice(nlmax, nlmax, alm , alm2,  cl(:,3:3)) ! T_1 x P_2
              call map2alm(nside,nlmax,nlmax,map (:,2),alm  , zrange,w8) ! P_1
              call alm2cl_spice(nlmax, nlmax, alm , alm2,  cl(:,2:2)) ! P_1 x P_2
              !           print*,minval(map(:,1)),maxval(map(:,1)),minval(map2(:,1)),maxval(map2(:,1))
              !           print*,'cl mask 1:',cl(0:10,1)
              !           print*,minval(map(:,2)),maxval(map(:,2))
              !           print*,'cl mask 2:',cl(0:10,2)
              !           print*,minval(map2(:,2)),maxval(map2(:,2))
              !           print*,'cl mask 3:',cl(0:10,3)
           endif
        else
           call map2alm(nside,nlmax,nlmax,map (:,1),alm  , zrange,w8) ! T
           call map2alm(nside,nlmax,nlmax,map (:,2),alm2 , zrange,w8) ! P
           call alm2cl_spice(nlmax, nlmax, alm , alm,  cl(:,1:1)) ! TT
           call alm2cl_spice(nlmax, nlmax, alm2, alm2, cl(:,2:2)) ! PP
           call alm2cl_spice(nlmax, nlmax, alm , alm2, cl(:,3:3)) ! TP
        endif
        deallocate(alm2)
     else
        ! same mask for T and P: either T only signal map, or T mask(s)
        if (do_cross) then
           allocate(alm2(1:1,0:nlmax,0:nlmax))
           alm2 = 0
           call map2alm(nside,nlmax,nlmax,map (:,1),alm,  zrange,w8)
           call map2alm(nside,nlmax,nlmax,map2(:,1),alm2, zrange,w8)
           call alm2cl_spice(nlmax, nlmax, alm, alm2, cl)
           if (do_cross_and_auto) then
              call alm2cl_spice(nlmax, nlmax, alm,  alm,  cl_auto1)
              call alm2cl_spice(nlmax, nlmax, alm2, alm2, cl_auto2)
           endif
           deallocate(alm2)
        else
           call map2alm(nside,nlmax,nlmax,map(:,1),alm, zrange,w8)
           call alm2cl_spice(nlmax, nlmax, alm, alm, cl)
        endif
     endif
     deallocate(w8)
     deallocate(alm)
          
  else
     ! polarized (signal map)
     allocate(w8(1:2*nside,1:3))
     allocate(alm(1:3,0:nlmax,0:nlmax))
     w8 = 1.0_DP
     if (do_cross) then
        allocate(alm2(1:3,0:nlmax,0:nlmax))
        call map2alm(nside,nlmax,nlmax,map, alm, zrange,w8)
        call map2alm(nside,nlmax,nlmax,map2,alm2,zrange,w8)
        call alm2cl_spice(nlmax, nlmax, alm, alm2, cl, symmetric=symmetric)
        if (do_cross_and_auto) then
           call alm2cl_spice(nlmax, nlmax, alm,  alm,  cl_auto1)
           call alm2cl_spice(nlmax, nlmax, alm2, alm2, cl_auto2)
        endif
        deallocate(alm2)
     else
        call map2alm(nside,nlmax,nlmax,map, alm, zrange,w8)
        call alm2cl_spice(nlmax, nlmax, alm, alm, cl)
     endif
     deallocate(alm)
     deallocate(w8)
     
  endif
  
end subroutine do_cl_with_alm


!=============================================================================
subroutine do_xi_from_cl(xi,cl,Pl,nlmax,verbose,message,ncor,polarization)
  ! ----------------------------------------
  ! C(l) -> xi(theta)
  ! ----------------------------------------
  ! June 2008: produces TU and QU as well
  ! July 2008: correction on QU ( new_value = (-2) * old_value )
  !
  ! d^l_00 = Pl(*,l,0)
  ! d^l_{2 2} = 2(Pl(*,l,1)+Pl(*,l,2))
  ! d^l_{2-2} = 2(Pl(*,l,1)-Pl(*,l,2))
  ! d^l{{20} = Pl(*,l,3) =  = Pl(*,l,4)
  !
  ! Note: Equations implemented here differ from those in Chon et al (2004) in
  ! in 2 respects
  ! - there is a typo in the paper, and Eq. (27) reads
  !     C^{TE} - i C^{TB} = 2 \pi \int \xi_X d_{20} d\cos\beta
  !  and is therefore consistent with Eq. (44);
  ! - to be consistent with Healpix
  !     Q_spice = Q_healpix = + Q_chon 
  !     U_spice = U_healpix = - U_chon 
  !
!=============================================================================
  use healpix_types
  implicit none
    
  integer(I4B), intent(in) :: nlmax,ncor
  real(DP), intent(out) :: xi(0:nlmax,1:ncor) ! TT, QQ, UU, TQ, TU, QU, [QT, UT, UQ]
  real(DP), intent(in)  :: cl(0:nlmax,1:ncor) ! TT, EE, BB, TE, TB, EB, [ET, BT, BE]
  real(DP), intent(in)  :: Pl(0:nlmax,0:nlmax,0:ncor-1)
  logical, intent(in)  :: verbose,polarization
  character(len=*), intent(in)  :: message
  
  integer(I4B) :: itheta,l
  logical(LGT) :: bcoupling, asym_cl
  real(DP) :: twolp1
  
  bcoupling    = (ncor >=6) .and. polarization
  asym_cl      = (ncor > 6) .and. polarization
  if (verbose) write(*,*) 'Compute xi from Cl '//trim(message)
    
  xi=0.0_DP
  if (.not. polarization) then     
     do itheta = 0,nlmax
        do l = 0,nlmax
           xi(itheta,1) = xi(itheta,1) + cl(l,1)*Pl(itheta,l,0)*dble(2*l+1)
        enddo
        xi(itheta,1)=xi(itheta,1)/FOURPI
     enddo     
  else
     do itheta = 0,nlmax
        do l = 0,nlmax
           twolp1 = dble(2*l+1)
           xi(itheta,1) = xi(itheta,1) + cl(l,1)*Pl(itheta,l,0)*twolp1
           xi(itheta,2) = xi(itheta,2) &
                &           + 2.0_DP*(cl(l,2)*Pl(itheta,l,1)+cl(l,3)*Pl(itheta,l,2))*twolp1
           xi(itheta,3) = xi(itheta,3) &
                &           + 2.0_DP*(cl(l,3)*Pl(itheta,l,1)+cl(l,2)*Pl(itheta,l,2))*twolp1
           xi(itheta,4) = xi(itheta,4) + cl(l,4)*Pl(itheta,l,3)*twolp1
           if (bcoupling) then
              xi(itheta,5) = xi(itheta,5) + cl(l,5)*Pl(itheta,l,4)*twolp1
              !               xi(itheta,6) = xi(itheta,6) + cl(l,6)*(Pl(itheta,l,2)-Pl(itheta,l,1))*twolp1! July 1st 2008
              xi(itheta,6) = xi(itheta,6) &
                   &           + 2.0_DP*cl(l,6)*(Pl(itheta,l,1)-Pl(itheta,l,2))*twolp1 ! July 1st, 2008
           endif
           if (asym_cl) then
              xi(itheta,7) = xi(itheta,7) + cl(l,7)*Pl(itheta,l,3)*twolp1
              xi(itheta,8) = xi(itheta,8) + cl(l,8)*Pl(itheta,l,4)*twolp1
              xi(itheta,9) = xi(itheta,9) &
                   &           + 2.0_DP*cl(l,9)*(Pl(itheta,l,1)-Pl(itheta,l,2))*twolp1 ! July 1st, 2008
           endif
        enddo
        xi(itheta,:) = xi(itheta,:)/FOURPI
     enddo
  endif

end subroutine do_xi_from_cl

  
!====================================================================================================
subroutine do_cl_from_xi(xi,cl,w,Pl,nlmax,verbose,message,ncor, &
     &                       polarization,decouple,mu,apodizesigma,thetamax,apodizetype)
  ! ----------------------------------------
  ! xi(theta) -> C(l)
  ! ----------------------------------------
  ! June 2008: produces TB and EB spectra as well
  ! July 2008: correction on input QU ( new_value = (-2) * old_value )
  ! Jan  2009: renormalize TE even in non 'decouple' configuration
!====================================================================================================
  use healpix_types
  use apodize_mod, only : apodizefunction
  use windows, only : get_TE_factor
  implicit none

  integer(I4B), intent(in)  :: nlmax,ncor
  real(DP),     intent(in)  :: w(nlmax+1)
  real(DP),     intent(in)  :: Pl(0:,0:,0:)
  real(DP),     intent(in)  :: xi(0:nlmax,1:ncor) ! TT, QQ, UU, TQ, TU, QU or TT, EE, BB, TQ, TU, QU
  real(DP),     intent(out) :: cl(0:nlmax,1:ncor) ! TT, EE, BB, TE, TB, EB
  real(DP),     intent(in), OPTIONAL :: apodizesigma
  real(DP),     intent(in), OPTIONAL :: thetamax
  real(DP),     intent(in), OPTIONAL :: mu(1:nlmax+1) !Tabulated roots of LP
  integer(I4B), intent(in), OPTIONAL :: apodizetype
  logical,          intent(in) :: verbose,polarization,decouple
  character(len=*), intent(in) :: message
  
  integer(I4B) :: i,l
  real(DP) :: wxi(ncor)
  real(DP) :: Fl(0:nlmax)     ! Normalization factor for EE and BB spectra, decoupling case
  real(DP) :: Crossl(0:nlmax) ! Normalization factor for TE spectra
  real(DP) :: nfact = 2.0_DP !SZ(cmbfast), 1 for KKS
  real(DP) :: tempo
  logical(LGT) :: bcoupling, asym_cl

  if (verbose) write(*,*) 'Compute Cl from xi '//trim(message)
  
  if (decouple) Fl=0.0_DP
  bcoupling    = (ncor >=6) .and. polarization
  asym_cl      = (ncor > 6) .and. polarization

  cl=0.0_DP
  if (.not. polarization) then
     !------------------------
     ! ***** unpolarized *****
     !------------------------
     do i=0,nlmax
        wxi(1)=w(i+1)*xi(i,1)*twopi
        do l=0,nlmax
           cl(l,1) = cl(l,1) + Pl(i,l,0)*wxi(1)
        enddo
     enddo
  else
     if (.not. bcoupling) then
        !---------------------------------
        ! ***** polarized (4 fields) *****
        !---------------------------------
        do i=0,nlmax
           wxi(1) = w(i+1) * xi(i,1) * twopi
           wxi(2) = nfact**2 * w(i+1)* xi(i,2) * pi
           wxi(3) = nfact**2 * w(i+1)* xi(i,3) * pi
           wxi(4) = nfact * w(i+1) * xi(i,4) * pi
           if (decouple) then ! compute renormalization factor
              call apodizefunction(mu(i+1), apodizesigma, thetamax, tempo, apodizetype)
              tempo = tempo / sin(dacos(mu(i+1))/2.0_dp)**2 !csc^2(beta/2) factor, eq 65 of Chon et al.
              tempo = tempo * w(i+1) * 2.0_dp  ! d^l_{2-2} = 2(Pl(*,l,1)-Pl(*,l,2))
           endif
           do l=0,nlmax
              cl(l,1) = cl(l,1) + wxi(1)*Pl(i,l,0)
              if (.not. decouple) then
                 cl(l,2) = cl(l,2) + (wxi(2)*Pl(i,l,1) + wxi(3)*Pl(i,l,2))
                 cl(l,3) = cl(l,3) + (wxi(3)*Pl(i,l,1) + wxi(2)*Pl(i,l,2))
              else
                 cl(l,2) = cl(l,2) + wxi(2)*(Pl(i,l,1)-Pl(i,l,2)) !\int 1/2(C-QQ+UU)d_{2-2}
                 cl(l,3) = cl(l,3) + wxi(3)*(Pl(i,l,1)-Pl(i,l,2)) !\int 1/2(C+QQ-UU)d_{2-2}
                 Fl(l) = Fl(l) + tempo*(Pl(i,l,1)-Pl(i,l,2))
              endif
              cl(l,4) = cl(l,4) + wxi(4)*Pl(i,l,3)
           enddo
        enddo

        ! **** bug correction 2009-01-23: renormalization of TE must be done in any case
        !     not only in 'decouple' case
        ! renormalize the TE spectra by the sum of the window {}_{x}K_{ll'}
        ! Computed from the 3js
        call get_TE_factor(Crossl)
        cl(2:nlmax,4) = cl(2:nlmax,4) / Crossl(2:nlmax)

        if (decouple) then
           ! renormalize the polarized Cls by the sum of the window {}_{-2}K_{ll'}
           cl(2:nlmax,2) = cl(2:nlmax,2) / Fl(2:nlmax)
           cl(2:nlmax,3) = cl(2:nlmax,3) / Fl(2:nlmax)
        endif

     else

        !---------------------------------
        ! ***** polarized (6 or 9 fields) *****
        !---------------------------------
        do i=0,nlmax
           wxi(1) =            w(i+1) * xi(i,1) * twopi
           wxi(2) = nfact**2 * w(i+1) * xi(i,2) * pi
           wxi(3) = nfact**2 * w(i+1) * xi(i,3) * pi
           wxi(4) = nfact    * w(i+1) * xi(i,4) * pi
           wxi(5) = nfact    * w(i+1) * xi(i,5) * pi
!            wxi(6) = nfact**2 * w(i+1) * xi(i,6) * twopi ! July 1st, 2008
           wxi(6) = nfact**2 * w(i+1) * xi(i,6) * pi      ! July 1st, 2008
           if (asym_cl) then
              wxi(7) = nfact    * w(i+1) * xi(i,7) * pi
              wxi(8) = nfact    * w(i+1) * xi(i,8) * pi
              wxi(9) = nfact**2 * w(i+1) * xi(i,9) * pi
           endif
           if (decouple) then ! compute renormalization factor
              call apodizefunction(mu(i+1),apodizesigma,thetamax,tempo,apodizetype)
              tempo = tempo / sin(dacos(mu(i+1))/2.0_dp)**2 !csc^2(beta/2) factor, eq 65 of Chon et al.
              tempo = tempo * w(i+1) * 2.0_dp  ! d^l_{2-2} = 2(Pl(*,l,1)-Pl(*,l,2))
           endif
           do l=0,nlmax
              cl(l,1) = cl(l,1) + wxi(1)*Pl(i,l,0)
              if (.not. decouple) then
                 cl(l,2) = cl(l,2) + (wxi(2)*Pl(i,l,1) + wxi(3)*Pl(i,l,2))
                 cl(l,3) = cl(l,3) + (wxi(3)*Pl(i,l,1) + wxi(2)*Pl(i,l,2))
              else
                 cl(l,2) = cl(l,2) + wxi(2)*(Pl(i,l,1)-Pl(i,l,2)) !\int 1/2(C-QQ+UU)d_{2-2}
                 cl(l,3) = cl(l,3) + wxi(3)*(Pl(i,l,1)-Pl(i,l,2)) !\int 1/2(C+QQ-UU)d_{2-2}
                 Fl(l)   = Fl(l)   + tempo *(Pl(i,l,1)-Pl(i,l,2))
              endif
              cl(l,4) = cl(l,4) + wxi(4)*Pl(i,l,3)
              cl(l,5) = cl(l,5) + wxi(5)*Pl(i,l,4)
!               cl(l,6) = cl(l,6) + wxi(6)*(Pl(i,l,2)-Pl(i,l,1)) ! July 1st, 2008
              cl(l,6) = cl(l,6) + wxi(6)*(Pl(i,l,1)-Pl(i,l,2))   ! July 1st, 2008
              if (asym_cl) then 
                 cl(l,7) = cl(l,7) + wxi(7)*Pl(i,l,3)
                 cl(l,8) = cl(l,8) + wxi(8)*Pl(i,l,4)
                 cl(l,9) = cl(l,9) + wxi(9)*(Pl(i,l,1)-Pl(i,l,2))
              endif
           enddo
        enddo

        ! **** bug correction 2009-01-23: renormalization of TE must be done in any case
        !     not only in 'decouple' case
        ! renormalize the TE spectra by the sum of the window {}_{x}K_{ll'}
        ! Computed from the 3js
        call get_TE_factor(Crossl)
        cl(2:nlmax,4) = cl(2:nlmax,4) / Crossl(2:nlmax)
        cl(2:nlmax,5) = cl(2:nlmax,5) / Crossl(2:nlmax)
        if (asym_cl) then
           ! same Crossl ????
           cl(2:nlmax,7) = cl(2:nlmax,7) / Crossl(2:nlmax)
           cl(2:nlmax,8) = cl(2:nlmax,8) / Crossl(2:nlmax)
        endif

        if (decouple) then
           ! renormalize the polarized Cls by the sum of the window {}_{-2}K_{ll'}
           cl(2:nlmax,2) = cl(2:nlmax,2) / Fl(2:nlmax)
           cl(2:nlmax,3) = cl(2:nlmax,3) / Fl(2:nlmax)
        endif
     endif
  endif

end subroutine do_cl_from_xi

!=======================================================================
!subroutine correct_final_cl(only_non_trivial)
subroutine correct_final_cl(cl_local, only_non_trivial, verbose_local)
!=======================================================================
! corrects for the beam, pixel window and transfer functions
! also apply scaling.
! NB:
! Beam and pixel windows are assumed to be the same for temperature and polarization
!
! 2012-07-13: replaced allocatable 'tmparray' with automatic array
!-----------------------------------------------------------------------
  use healpix_types
  use misc_utils, only: fatal_error
  use spice_common, only: nlmax, ncor, &
       & correct_beam, correct_beam2, correct_pix, correct_transfer_function, normalize, &
       & gb,           gb2,           wl,          transfer_function,         quad_uK
  implicit none
  
  real(dp), dimension(0:,1:), intent(inout) :: cl_local
  logical,                    intent(in)    :: only_non_trivial, verbose_local
    
  real(DP) :: sigb2,quad_uK2,rl
  integer(I4B) :: l,i,bad, status
!  real(DP), dimension(:), allocatable :: tmparray
  real(DP), dimension(0:nlmax) :: tmparray
  integer(i4b) :: ncor_local, nlmax_local

  nlmax_local = ubound(cl_local,1)
  ncor_local = ubound(cl_local,2)
  if (verbose_local) then
     if (nlmax /= nlmax_local) print*, 'NLmax = ', nlmax, nlmax_local
     if (ncor /= ncor_local)   print*, 'NCor  = ', ncor,  ncor_local
  endif

  if (correct_beam.or.correct_beam2) then
     if (verbose_local) write(*,*) 'Correcting Cl for the beam window function'
!     allocate(tmparray(0:nlmax_local))

     if (correct_beam2) then
        tmparray(0:nlmax_local) = gb(0:nlmax_local,1) * gb2(0:nlmax_local,1) ! cross correlate maps with different beams
     else
        tmparray(0:nlmax_local) = gb(0:nlmax_local,1) * gb(0:nlmax_local,1) ! auto-correlation, or same beam for both maps
     endif

     bad = count( tmparray(0:nlmax_local) == 0 )
     if (bad > 0) then
        write(*,*) 'ERROR in correct_final_cl :'
        write(*,*) 'Provided beam window function vanishes at ',bad,' multipoles'
        write(*,*) 'Can not correct for beam'
        CALL FATAL_ERROR
     endif
     do i=1,ncor_local ! TT [, EE, BB, TE, TB, EB]
        cl_local(0:nlmax_local,i) = cl_local(0:nlmax_local,i) / tmparray(0:nlmax_local)
     enddo
!     deallocate(tmparray)
  Endif
       
  if (correct_pix) then
     if (verbose_local) write(*,*) 'Correcting Cl for the pixel window function'
!!!     allocate(tmparray(0:nlmax_local))
     tmparray(0:nlmax_local) = wl(0:nlmax_local,1) * wl(0:nlmax_local,1)

     do i=1,ncor_local
        cl_local(0:nlmax_local,i) = cl_local(0:nlmax_local,i) / tmparray(0:nlmax_local)
     enddo
!!!!     deallocate(tmparray)
!!!!!     print*,"done correct. of pix wind funct"
  endif

  if (correct_transfer_function) then
     bad = count( transfer_function(0:nlmax,1:ncor) == 0 )
     if (bad > 0) then
        write(*,*) 'ERROR in correct_final_cl :'
        write(*,*) 'Provided transfer function vanishes at ',bad,' multipoles'
        write(*,*) 'Can not correct for transfer function'
        CALL FATAL_ERROR
     endif
     do i=1,ncor_local ! TT [, EE, BB, TE, TB, EB]
        cl_local(0:,i) = cl_local(0:,i) / transfer_function(0:nlmax_local,i)
     enddo
  endif

  if (only_non_trivial) return
  
  if (normalize) then
     if (verbose_local) write(*,*) 'Normalizing Cl as requested'
     quad_uK2=quad_uK**2
     cl_local=cl_local*quad_uK2
  endif

end subroutine correct_final_cl

!=======================================================================
subroutine correct_final_xi(xi_local)
!=======================================================================
  use spice_common, only: DP, normalize, megaverbose, quad_uK ! ,xi_final
  implicit none
  real(DP), dimension(0:,1:), intent(inout) :: xi_local
  
  real(DP) :: quad_uK2
    
  if (normalize) then
     if (megaverbose) write(*,*) 'Normalizing xi as requested'
     quad_uK2=quad_uK**2
!     xi_final=xi_final*quad_uK2
     xi_local = xi_local * quad_uK2
  endif
end subroutine correct_final_xi

!=======================================================================
subroutine correct_xi_from_mask(nlmax, ncor, ncmask, xi, xi_mask, xi_final, verbose)
!=======================================================================
  use healpix_types
  implicit none
  integer(i4b), intent(in) :: nlmax, ncor, ncmask
  logical(LGT), intent(in) :: verbose
  real(DP), dimension(0:nlmax,1:ncor),   intent(in)  :: xi
  real(DP), dimension(0:nlmax,1:ncmask), intent(in)  :: xi_mask
  real(DP), dimension(0:nlmax,1:ncor),   intent(out) :: xi_final

  logical(LGT):: warning
  integer(i4b), dimension(1:9) :: ncountbinrm, kindex
  integer(i4b) :: j, l, kmask

  warning=.false.
  ncountbinrm(:)=0

  !!maskthreshold = -1.d30

  kindex = (/ 1, 1, 1, 1, 1, 1, 1, 1, 1 /)
  if (ncmask == 3) kindex = (/  1, 2, 2, 3, 3, 2, 3, 3, 2 /)
  if (ncmask == 4) kindex = (/  1, 2, 2, 3, 3, 2, 4, 4, 2 /)
     do j= 1, ncor
        kmask = kindex(j)
        do l=0,nlmax
           if (xi_mask(l,kmask) /= 0.d0) then
               xi_final(l,j)=xi(l,j)/xi_mask(l,kmask)
           else
              warning=.true.
              xi_final(l,j)=0.d0
              ncountbinrm(j)=ncountbinrm(j)+1
           endif
        enddo
     enddo
!!!     if (allocated(sintheta)) deallocate(sintheta)

     if (warning.and.verbose) then
        write(*,*)
        write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*) 'WARNING in compute'
!         if (pairsthresholding) then
!            write(*,*) 'Some of the values of xi_mask are below npairthreshold = ',npairsthreshold
!         else
           write(*,*) 'Some of the values of xi_mask are equal to zero'
!        endif
        write(*,*) 'Corresponding values of xi_final are set to zero'
        write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*)
        write(*,'(a,6(i5))') 'number of bins removed =',ncountbinrm(1:ncor)
        write(*,'(a,i6)') 'out of ',nlmax+1
     endif

   end subroutine correct_xi_from_mask
end module deal_with_xi_and_cl
