module utilities

  implicit none 

  public :: utility_recip_lattice 
  public :: utility_cart_to_frac
  public :: utility_frac_to_cart
  public :: utility_diag
  public :: heavyside
  public :: fermi_dirac
  public :: trace
  public :: gauss

  external  :: cheev

  contains

  function heavyside(x) result (hside)
    
    real*8, intent(in)  :: x
    real*8              :: hside

    if (x < 0) then
       hside = 0.0
    else if (x >= 1) then
       hside = 1.0
    end if
 
  end function heavyside

  function fermi_dirac(x,temp) result (fd)

    use constants,    only : kboltz

    real*8, intent(in)  :: x, temp
    real*8              :: fd

    fd = 1.0/(exp(x/(temp*kboltz)) + 1)

  end function fermi_dirac

  function gauss(x,x0,sigma) result (gaussfunc)

    use constants,    only : twopi

    real*8, intent(in)  :: x, x0, sigma
    real*8              :: gaussfunc

    if ((x-x0) > -6*sigma .and. (x-x0) < 6*sigma) then
       gaussfunc = exp(-((x-x0)**2)/(2*sigma**2))/sqrt(twopi*sigma**2)
    else
       gaussfunc = 0.0
    end if

  end function gauss

  function trace(A,N) result (tr)

    integer    , intent(in) :: N
    complex*16 , intent(in) :: A(N,N)
    complex*16              :: tr
    integer                 :: i  

    tr = 0.0
    do i = 1, N
       tr = tr + A(i,i)
    end do
  
  end function trace

    

  subroutine utility_recip_lattice (real_lat,recip_lat,volume)  !
    !==================================================================!
    !                                                                  !
    !!  Calculates the reciprical lattice vectors and the cell volume
    !                                                                  !
    !===================================================================

    implicit none
    real*8, intent(in)  :: real_lat (3, 3)
    real*8, intent(out) :: recip_lat (3, 3)
    real*8, intent(out) :: volume

    recip_lat(1,1)=real_lat(2,2)*real_lat(3,3)-real_lat(3,2)*real_lat(2,3)
    recip_lat(2,1)=real_lat(3,2)*real_lat(1,3)-real_lat(3,3)*real_lat(1,2)
    recip_lat(3,1)=real_lat(1,2)*real_lat(2,3)-real_lat(2,2)*real_lat(1,3)
    recip_lat(1,2)=real_lat(2,3)*real_lat(3,1)-real_lat(3,3)*real_lat(2,1)
    recip_lat(2,2)=real_lat(3,3)*real_lat(1,1)-real_lat(1,3)*real_lat(3,1)
    recip_lat(3,2)=real_lat(2,1)*real_lat(1,3)-real_lat(2,3)*real_lat(1,1)
    recip_lat(1,3)=real_lat(2,1)*real_lat(3,2)-real_lat(3,1)*real_lat(2,2)
    recip_lat(2,3)=real_lat(3,1)*real_lat(1,2)-real_lat(1,1)*real_lat(3,2)
    recip_lat(3,3)=real_lat(2,2)*real_lat(1,1)-real_lat(2,1)*real_lat(1,2)

    volume=real_lat(1,1)*recip_lat(1,1) + &
         real_lat(2,1)*recip_lat(2,1) + &
         real_lat(3,1)*recip_lat(3,1)

    recip_lat=2.*acos(-1.0)*recip_lat/volume
    volume=abs(volume)

  end subroutine utility_recip_lattice

  subroutine utility_cart_to_frac(cart,frac,recip_lat)
    !==================================================================!
    !                                                                  !
    !!  Convert from Cartesian to fractional coordinates
    !                                                                  !
    !===================================================================
    implicit none

    real*8, intent(in)  :: recip_lat(3,3)
    real*8, intent(out) :: frac(3)
    real*8, intent(in)  :: cart(3)

    integer :: i

    do i=1,3
       frac(i)=recip_lat(i,1)*cart(1) + recip_lat(i,2)*cart(2) + recip_lat(i,3)*cart(3)
    end do

    frac=frac/(2.0*acos(-1.0))


  end subroutine utility_cart_to_frac

  subroutine utility_frac_to_cart(frac,cart,real_lat)
    !==================================================================!
    !                                                                  !
    !!  Convert from fractional to Cartesian coordinates
    !                                                                  !
    !===================================================================
    implicit none

    real*8, intent(in)  :: real_lat(3,3)
    real*8, intent(in)  :: frac(3)
    real*8, intent(out) :: cart(3)

    integer :: i

    do i=1,3
       cart(i)=real_lat(i,1)*frac(1) + real_lat(i,2)*frac(2) + real_lat(i,3)*frac(3)
    end do

    return

  end subroutine utility_frac_to_cart

  subroutine utility_diag(MAT,EIG,N)

    implicit none

    integer, intent(in) :: N
    integer :: INF, LWORK
    integer, parameter :: LWMAX = 1000
    real*8, dimension(N), intent(out) :: EIG
    complex*16, dimension(N,N), intent(inout) :: MAT
    complex*16 :: W(LWMAX)
    real*8, dimension(3*N-2) :: RW

    call zheev('V','U',N,MAT,N,EIG,W,-1,RW,INF)
    LWORK = min(LWMAX, int(W(1)))
    call zheev('V','U',N,MAT,N,EIG,W,LWORK,RW,INF)

  end subroutine utility_diag

end module utilities
