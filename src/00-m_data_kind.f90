  !*********************************************************************/
  !** This code has been done in the Barcelona Center for Subsurface 
  !** Imaging (BCSI)
  !** Authors: Daniel Dagnino.
  !*********************************************************************/
  !** Short resume of this code/file:
  !** 
  !** 
  !** 
  !*********************************************************************/
  module m_data_kind
    implicit none
!     ! At least 6 significant digits of precision and an exponent range of at least 37.
!     integer, parameter :: sp = selected_real_kind(6, 37)
!     ! At least 15 significant digits of precision and an exponent range of at least 307.
!     integer, parameter :: rl = selected_real_kind(15, 307)
    
    ! Double
    integer, parameter :: rl = kind(1.0d0)
    
    ! Single
    integer, parameter :: sp = kind(1.0d0)
    
    !----------------------------------------------------------!
  end module m_data_kind




