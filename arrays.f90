Module arrays

  use healpix_types

  Implicit none
  
  Integer :: status1,status2,status3,status4,status5,status6
  !    Integer*4,dimension(13) :: buff
  !   Integer*4,allocatable,dimension(:) :: ID

  Real(kind=DP), allocatable, dimension(:,:) :: map, gmap,cmbmask,planckmap,vsim,ssim,ksim,vmean,smean,kmean
  Real(kind=DP), allocatable, dimension(:,:) :: vsdv,ssdv,ksdv
  Real(kind=DP), allocatable, dimension(:,:) :: kl, sl, vl

!!$    Character(len=6),allocatable,dimension(:) :: FieldHipp

End module arrays
