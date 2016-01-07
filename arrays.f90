Module arrays

    Integer :: status1,status2,status3,status4,status5,status6
!    Integer*4,dimension(13) :: buff
 !   Integer*4,allocatable,dimension(:) :: ID

    Real*8, allocatable, dimension(:,:) :: map, gmap
!!$    Real*8, allocatable, dimension(:) :: PeriodB,HB,Sigma_mB,VB,IIB,current_point
!!$    Real*8, allocatable, dimension(:) :: PeriodC,HC,Sigma_mC,VC,IIC,mvi5av,Sigma_mvi5av,logP,Mw,sigmaMw
!!$    Real*8, allocatable, dimension(:) :: Period,H,Sigma_m,V,II,PeriodR11,VIR11,F160WR11,eF160WR11,OHR11
!!$    Real*8, allocatable, dimension(:,:,:,:,:) :: cov,inv_cov
!!$    Real*4, allocatable, dimension(:) :: acceptance_probability
!!$    Real*8, allocatable, dimension(:,:) :: Covgauss,Covguess
!!$    Real*8 :: jumping_factor                           ! SAVES JUMPING FACTOR FOR MCMC (INCREASE IF WANT BIGGER STEP SIZE, DECREASE OTHERWISE) 

!!$    Character(len=10),allocatable,dimension(:) :: Name,NameA,NameB,NameC
!!$    Character(len=5),allocatable,dimension(:) :: Field,Fieldmvi
!!$    Character(len=6),allocatable,dimension(:) :: FieldHipp

End module arrays
