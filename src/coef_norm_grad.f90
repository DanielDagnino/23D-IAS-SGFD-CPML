  !*********************************************************************/
  !** This code has been done in the Barcelona Center for Subsurface 
  !** Imaging (BCSI).
  !** Authors: Daniel Dagnino.
  !*********************************************************************/
  !** Short resume of this code/file:
  !** 
  !*********************************************************************/
  module m_coef_norm_grad
    implicit none
    
    !*********************************************************************/
    contains
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine coef_norm_grad( var, n, norm_var, norm_grad, coef_norm, grad )
      use m_data_kind
      implicit none
      ! The variables which are passed to the function.
      real(sp), intent(out)  :: coef_norm
      real(sp), intent(out)  :: norm_var, norm_grad
      real(sp), intent(in)   :: n
      real(sp), allocatable, intent(in)    :: var(:,:,:)
      real(sp), allocatable, intent(inout) :: grad(:,:,:,:)
      ! The variables which are generated inside the function.
      real(sp) :: sig
!dir$ assume_aligned var(1,1,1):64,grad(1,1,1,1):64
      
      !----------------------------------------------------------!
      ! We rescale the gradient to avoid grad=0.
      norm_var  = 1./maxval(abs(var))
      norm_grad = 1./maxval(abs(grad(1,:,:,:)))
      
      if ( n == 1. ) then
        sig = -1.
      else
        sig = -(1.-n)/abs(1.-n)
      end if
      
      grad(1,:,:,:) = sig * (norm_grad*grad(1,:,:,:)) * ( (norm_var*var)**(n/(1.-n)) )
!       grad(1,:,:,:) = sig * (norm_grad*grad(1,:,:,:)) * ( (norm_var*var)**((2.*n-1.)/(1.-n)) )
      
      ! 
      coef_norm   = 1./maxval(abs(grad(1,:,:,:)))
      grad(1,:,:,:) = coef_norm*grad(1,:,:,:)
      
    end subroutine coef_norm_grad
    
  !**********************************************************!
  end module m_coef_norm_grad





