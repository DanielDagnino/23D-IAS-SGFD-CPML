  !*********************************************************************/
  !** This code has been done in the Barcelona Center for Subsurface 
  !** Imaging (BCSI)
  !** Authors: Daniel Dagnino.
  !*********************************************************************/
  !** Short resume of this code/file:
  !** Compute of the initial and last point that include a specific 
  !** percent of the total weight of a function (area bellow the function).
  !** The support of a function.
  !** support_funct1: Area of |f(x)|.
  !** support_funct2: Area of f(x)^2.
  !** Meaning for support_funct2:
  !** t0_{sou} = { max t0 | p0_cut = \int_{t=0}^{t0}    f(t)^2 }
  !** t2_{sou} = { min t1 | p1_cut = \int_{t=t1}^{\inf} f(t)^2 }
  !** Example:
  !** f(t):  0  0  0  0  0  1  2  3  2  1  0  0  0  0  0
  !** t   :  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
  !** In case: p0_cut = p1_cut = 0., then:
  !** t0_cut= 6
  !** t1_cut= 10
  !*********************************************************************/
  module m_support_funct
    implicit none
    
    !**********************************************************!
    contains
    
    !**********************************************************!
    !**********************************************************!
    !**********************************************************!
    subroutine support_funct1( f, nt, dt, p0_cut, p1_cut, t0_cut, t1_cut )
      use m_data_kind
      implicit none
      ! The variables which are passed to the function.
      real(sp), intent(out) :: t0_cut, t1_cut
      integer, intent(in)  :: nt
      real(sp), intent(in) :: f(:)
      real(sp), intent(in) :: dt, p0_cut, p1_cut
      ! The variables which are generated inside the function.
      integer :: it
      real(sp) :: integ
      real(sp) :: p0, p1
!dir$ assume_aligned f(1):64
      
      !----------------------------------------------------------!
      ! Total weight.
      integ = 0._sp
      do it=1,nt
        integ = integ + abs(f(it))
      end do
      
      ! Check.
      if ( integ == 0._sp ) stop '***** ERROR - support_funct1: integ == 0._sp *****'
      
      ! Weight to be integrated.
      p0 = (1._sp-p0_cut)*integ
      p1 = (1._sp-p1_cut)*integ
      
      !----------------------------------------------------------!
      ! t1
      integ = 0._sp
      do it=1,nt
        integ = integ + abs(f(it))
        if ( integ > p1 ) then
          t1_cut = dt*real(it-1,sp)
          exit
        end if
      end do
      
      ! t0
      integ = 0._sp
      do it=nt,1,-1
        integ = integ + abs(f(it))
        if ( integ > p0 ) then
          t0_cut = dt*real(it-1,sp)
          exit
        end if
      end do
      
      !----------------------------------------------------------!
    end subroutine support_funct1
    
    
    
    !**********************************************************!
    !**********************************************************!
    !**********************************************************!
    subroutine support_funct2( f, nt, dt, p0_cut, p1_cut, t0_cut, t1_cut )
      use m_data_kind
      implicit none
      ! The variables which are passed to the function.
      real(sp), intent(out) :: t0_cut, t1_cut
      integer, intent(in)  :: nt
      real(sp), intent(in) :: f(:)
      real(sp), intent(in) :: dt, p0_cut, p1_cut
      ! The variables which are generated inside the function.
      integer :: it
      real(sp) :: integ
      real(sp) :: p0, p1
!dir$ assume_aligned f(1):64
      
      !----------------------------------------------------------!
      ! Total weight.
      integ = 0._sp
      do it=1,nt
        integ = integ + f(it)**2
      end do
      
      ! Check.
      if ( integ == 0._sp ) stop '***** ERROR - support_funct2: integ == 0._sp *****'
      
      ! Weight to be integrated.
      p0 = (1._sp-p0_cut)*integ
      p1 = (1._sp-p1_cut)*integ
      
      !----------------------------------------------------------!
      ! t1
      integ = 0._sp
      do it=1,nt
        integ = integ + f(it)**2
        if ( integ > p1 ) then
          t1_cut = dt*real(it-1,sp)
          exit
        end if
      end do
      
      ! t0
      integ = 0._sp
      do it=nt,1,-1
        integ = integ + f(it)**2
        if ( integ > p0 ) then
          t0_cut = dt*real(it-1,sp)
          exit
        end if
      end do
      
      !----------------------------------------------------------!
    end subroutine support_funct2
    
  !**********************************************************!
  end module m_support_funct
      
      