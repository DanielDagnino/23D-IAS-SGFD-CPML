  !*********************************************************************/
  !** This code has been done in the Barcelona Center for Subsurface 
  !** Imaging (BCSI).
  !** Authors: Daniel Dagnino.
  !*********************************************************************/
  !** Short resume of this code/file:
  !** Storage variables definition for the solver code.
  !*********************************************************************/
  module m_solv3Dahvc_time
    use m_data_kind
    implicit none
! #ifdef valgrind
! !     integer, parameter :: save_every = 10
!     integer, parameter :: save_every = 100
! #else
! !     integer, parameter :: save_every = 3500
! !     integer :: save_every = 3200
!     integer :: save_every = 2500
    integer :: save_every = 8000
! #endif
    integer :: nt_solver = -1
    real(rl) :: dt_solver
  end module m_solv3Dahvc_time



