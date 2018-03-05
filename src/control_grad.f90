  !*********************************************************************/
  !** This code has been done in the Barcelona Center for Subsurface 
  !** Imaging (BCSI).
  !** Authors: Daniel Dagnino.
  !*********************************************************************/
  !** Short resume of this code/file:
  !** 
  !*********************************************************************/
  module m_control_grad
    implicit none
    
    !*********************************************************************/
    contains
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine control_grad( model, grad, mean, sigma_o_mean, z_mean )
      use m_data_kind
      use m_geo
      implicit none
      ! The variables which are passed to the function.
      real(sp), intent(out) :: mean, z_mean, sigma_o_mean
      real(sp), allocatable, intent(in) :: model(:,:,:), grad(:,:,:)
      ! The variables which are generated inside the function.
      integer :: ix, iy, iz, cont
      real(sp) :: raux
      real(sp) :: mean2
!dir$ assume_aligned model(1,1,1):64,grad(1,1,1):64
      
      !----------------------------------------------------------!
      sigma_o_mean = 0.
      mean   = 0.
      mean2  = 0.
      z_mean = 0.
      cont = 0
      
      do iz=1,nz
      do ix=1,nx
      do iy=grid%ind_bath(ix,iz),ny
        cont   = cont + 1
        raux   = abs(grad(ix,iy,iz)/model(ix,iy,iz))
        mean   = mean + raux
        mean2  = mean2 + raux*raux
        z_mean = z_mean + dx*dble(iy-1)*raux
      end do
      end do
      end do
      
      ! 
      z_mean       = z_mean/mean
      mean         = mean /dble(cont)
      mean2        = mean2/dble(cont)
      sigma_o_mean = ( mean2 - mean**2 )/(mean**2)
      
    end subroutine control_grad
    
  !**********************************************************!
  end module m_control_grad




