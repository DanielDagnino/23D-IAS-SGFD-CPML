  !*********************************************************************/
  !** This code has been done in the Barcelona Center for Subsurface 
  !** Imaging (BCSI).
  !** Authors: Daniel Dagnino.
  !** Goal: This is a 2D/3D Acoustic Solver with the following 
  !**       characteristics:
  !**       - Inhomogenous density and compressibility.
  !**       - Time Domain.
  !**       - Spacial discretization in a staggered grid.
  !**       - Runge-Kutta-4th time integration method.
  !**       - 6th order spacial derivatives.
  !**       - Convolutional Perfectly Matched Layers in the boundaries.
  !**       - Free-Surface in the top boundary.
  !*********************************************************************/
  module m_solv3Dahvc_RK_loop
    implicit none
    
    contains
    !*********************************************************************/
    !*********************************************************************/
    elemental subroutine RK_loop( it, shgat, vf, ktmp, k1, k2, k3, k4 )
      use m_solv3Dahvc_data_kind
      use m_sg_type
      use m_solv3Dahvc_solver, only: dt_o_2, dt_solver, apply_FS
      use m_solv3Dahvc_RK_dif_op_split
      use m_solv3Dahvc_RK_dif_op_FS
      use m_solv3Dahvc_RK_step_split
      use m_solv3Dahvc_RK_step_FS
      use m_solv3Dahvc_RK_new_step_split
      use m_solv3Dahvc_RK_new_step_FS
      use m_solv3Dahvc_RK_sou
      implicit none
      ! The variables which are passed to the function.
      integer, intent(in) :: it
      type(shot_gather), intent(in) :: shgat
      type(field), intent(inout) :: vf, ktmp, k1, k2, k3, k4   ! Auxiliar fields for RK4.
      
      !*********************************************************************/
      ! k1
      !*********************************************************************/
      ! \dot \vec F = \hat L \vec F
      call RK_dif_op_split( vf, k1, shgat%mg )
      
      ! Operator for the FS condition.
      if ( apply_FS ) call RK_dif_op_FS( vf, k1, shgat%mg )
      
      ! \vec source: s(t) \delta(x-x_R) \delta(y-y_R) \delta(z-z_R)
      call RK_sou( k1%u, .true., it, shgat )
      
      ! FS Mirror.
      if ( apply_FS ) call RK_mirror( k1, shgat%mg )
      
      !*********************************************************************/
      ! k2
      !*********************************************************************/
      ! RK intermediated step
      call RK_step_split( dt_o_2, vf, k1, ktmp, shgat%mg )
      
      ! RK intermediated step for the FS condition.
      if ( apply_FS ) call RK_step_FS( dt_o_2, vf, k1, ktmp, shgat%mg )
      
      ! \dot \vec F = \hat L \vec F
      call RK_dif_op_split( ktmp, k2, shgat%mg )
      
      ! Operator for the FS condition.
      if ( apply_FS ) call RK_dif_op_FS( k2, ktmp, shgat%mg )
      
      ! \vec source: s(t) \delta(x-x_R) \delta(y-y_R) \delta(z-z_R)
      call RK_sou( k2%u, .false., it, shgat )
      
      ! FS Mirror.
      if ( apply_FS ) call RK_mirror( k2, shgat%mg )
      
      !*********************************************************************/
      ! k3
      !*********************************************************************/
      ! RK intermediated step
      call RK_step_split( dt_o_2, vf, k2, ktmp, shgat%mg )
      
      ! RK intermediated step for the FS condition.
      if ( apply_FS ) call RK_step_FS( dt_o_2, vf, k2, ktmp, shgat%mg )
      
      ! Operator 
      call RK_dif_op_split( ktmp, k3, shgat%mg )

      ! FS condition.
      if ( apply_FS ) call RK_dif_op_FS( ktmp, k3, shgat%mg )
      
      ! \vec source: s(t) \delta(x-x_R) \delta(y-y_R) \delta(z-z_R)
      call RK_sou( k3%u, .false., it, shgat )
      
      ! FS Mirror.
      if ( apply_FS ) call RK_mirror( k3, shgat%mg )
      
      !*********************************************************************/
      ! k4
      !*********************************************************************/
      ! RK intermediated step
      call RK_step_split( dt_solver, vf, k3, ktmp, shgat%mg )
      
      ! RK intermediated step for the FS condition.
      if ( apply_FS ) call RK_step_FS( dt_solver, vf, k3, ktmp, shgat%mg )
      
      ! Operator 
      call RK_dif_op_split( ktmp, k4, shgat%mg )
      
      ! FS condition.
      if ( apply_FS ) call RK_dif_op_FS( ktmp, k4, shgat%mg )
      
      ! \vec source: s(t) \delta(x-x_R) \delta(y-y_R) \delta(z-z_R)
      call RK_sou( k4%u, .true., it+1, shgat )
      
      ! FS Mirror.
      if ( apply_FS ) call RK_mirror( k4, shgat%mg )
      
      !*********************************************************************/
      ! n+1: wave eq
      !*********************************************************************/      
      ! RK final step.
      call RK_new_step_split( vf, k1, k2, k3, k4, shgat%mg )
      
      ! FS condition.
      if ( apply_FS ) call RK_new_step_FS( vf, k1, k2, k3, k4, shgat%mg )
      
    end subroutine RK_loop
  
  !***********************************************************************/
  end module m_solv3Dahvc_RK_loop


