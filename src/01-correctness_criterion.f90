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
  module m_solv3Dahvc_correctness_criterion
    use m_data_kind
    implicit none
    real(rl), private, parameter :: courant_stable = 1.25_rl
    
    !*********************************************************************/
    contains
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    pure real(rl) function dt_stability_CFL( c_max, dx )
      implicit none
      ! The variables which are passed to the function.
      real(rl), intent(in) :: c_max, dx
      !----------------------------------------------------------!
      ! Stability.
      ! courant_number = c_max*dt/dx \le courant_stable
      ! => dt_stability_CFL \leq dx*courant_stable/c_max
      dt_stability_CFL = (1._rl/courant_stable)*dx/c_max
    end function dt_stability_CFL
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    pure real(rl) function dx_correct( c_min, freq_max, npw )
      implicit none
      ! The variables which are passed to the function.
      integer, intent(in) :: npw
      real(rl), intent(in) :: c_min, freq_max
      !----------------------------------------------------------!
      ! Minimum of npw points per wave length.
      ! wave length = c_min/freq_max
      dx_correct = c_min/(freq_max*dble(npw))
    end function dx_correct
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    pure subroutine time_adapt( t_fin, store, dt, nt )
      implicit none
      ! The variables which are passed to the function.
      real(rl), intent(in):: t_fin
      integer, intent(inout) :: nt
      integer, intent(in) :: store
      real(rl), intent(inout) :: dt
      ! The variables which are generated inside the function.
      integer :: k
      
      !----------------------------------------------------------!
      ! New time sampling calculation.
      nt = 1 + ceiling(t_fin/dt)   ! (nt-1)*dt >= t_fin
      
      ! Increase nt up to become nt-1 divisible by store.
      do k=0,store-1
        if ( mod(nt-1+k,store) == 0 ) then
          nt = nt+k
          exit
        end if
      end do
      
      dt = t_fin/dble(nt-1)   ! (nt-1)*dt = t_fin
      
    end subroutine time_adapt
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
#ifdef testing
    integer function nt_store( nt, store )
#else
    pure integer function nt_store( nt, store )
#endif
      implicit none
      ! The variables which are passed to the function.
      integer, intent(in) :: nt
      integer, intent(in) :: store
      
      !----------------------------------------------------------!
      ! 
#ifdef testing
      if ( mod(nt-1,store) /= 0 ) stop 'ERROR function nt_store'
#endif
      
      nt_store = 1 + (nt-1)/store
      
    end function nt_store
    
    !***********************************************************************/
  end module m_solv3Dahvc_correctness_criterion
  
  
  
  
  