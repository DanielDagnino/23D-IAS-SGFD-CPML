  !*********************************************************************/
  !** This code has been done in the Barcelona Center for Subsurface 
  !** Imaging (BCSI).
  !** Authors: Daniel Dagnino.
  !*********************************************************************/
  !** Short resume of this code/file:
  !** Calibration term to adjust tr_o to tr_s (except the ones with inv)
  !** base in minimize the L2 norm adding a constant term. We apply:
  !** tr_s \approx calibration*tr_o
  !** 
  !** Name sufix meaning of the functions:
  !** Global calibartion: coef_calib 
  !** Receiver calibartion: array_calib 
  !** Calibration by maximum: max
  !** Calibration by product: L2
  !**                         if auto: <tr_s*tr_s>/<tr_o*tr_o>
  !**                         if cros: <tr_o*tr_s>/<tr_o*tr_o>
  !** With data preconditioning: prec
  !** For the inverse calibration (adjust tr_s to tr_o): inv
  !*********************************************************************/
  module m_calib
    implicit none
    
    !**********************************************************!
    contains
    
    
    
    !**********************************************************!
    !**********************************************************!
    !**********************************************************!
    subroutine coef_calib( shgat )
      use m_data_kind
      use m_sg_type
      use m_sg_func
      use m_support_funct, only: support_funct1
      use m_type_inversion, only: offset_cal_min, offset_cal_max, p0_cut, p1_cut
      use m_phys, only: c_water
      implicit none
      ! The variables which are passed to the function.
      type(shot_gather), intent(inout) :: shgat
      ! The variables which are generated inside the function.
      integer :: ir, is, it
      real(sp) :: denom
      real(sp) :: dt, dt_sou, t0_sou_cut, t1_sou_cut
      real(sp) :: offset_SR
      real(sp) :: coef
      
      !----------------------------------------------------------!
      ! Time grid.
      dt = shgat%t1_obs/real(shgat%nt_tr_obs-1,sp)
      
      ! Support of the observed source wavelet.
      do is=1,shgat%n_sou
        dt_sou = shgat%t1_sou/real(shgat%nt_sou-1,sp)
        call support_funct1( shgat%sou(is,:), shgat%nt_sou, dt_sou, p0_cut, p1_cut, t0_sou_cut, t1_sou_cut )
      end do
      
      !----------------------------------------------------------!
      ! Calibration R per R.
      coef = 0._sp
      denom = 0._sp
      do ir=1,shgat%n_rec
        offset_SR = offset_sou_rec( shgat, ir )
        if ( offset_SR > offset_cal_min .and. offset_SR < offset_cal_max ) then
          
          ! Coef. for the calibration.
          do it=1,shgat%nt_tr_obs
            coef = coef + shgat%pr_dat(ir,it)*shgat%tr_obs(ir,it)*shgat%tr_syn(ir,it)
            denom = denom + shgat%pr_dat(ir,it)*shgat%tr_obs(ir,it)*shgat%tr_obs(ir,it)
          end do
          
        end if
      end do
      
      !----------------------------------------------------------!
      ! Chek error.
      if ( denom == 0._sp ) stop '***** ERROR - coef_calib: denom == 0. *****'
      if ( coef == 0._sp ) stop '***** ERROR - coef_calib: coef == 0. *****'
      
      !----------------------------------------------------------!
      shgat%calib = coef/denom
      
      !----------------------------------------------------------!
    end subroutine coef_calib
    
    
    
    !**********************************************************!
    !**********************************************************!
    !**********************************************************!
    subroutine coef_calib2( shgat )
      use m_data_kind
      use m_sg_type
      use m_sg_func
      use m_support_funct, only: support_funct1
      use m_type_inversion, only: offset_cal_min, offset_cal_max, p0_cut, p1_cut
      use m_phys, only: c_water
      implicit none
      ! The variables which are passed to the function.
      type(shot_gather), intent(inout) :: shgat
      ! The variables which are generated inside the function.
      integer :: ir, is, it
      real(sp) :: denom
      real(sp) :: dt, dt_sou, t0_sou_cut, t1_sou_cut
      real(sp) :: offset_SR
      real(sp) :: coef
      
      !----------------------------------------------------------!
      ! Time grid.
      dt = shgat%t1_obs/real(shgat%nt_tr_obs-1,sp)
      
      ! Support of the observed source wavelet.
      do is=1,shgat%n_sou
        dt_sou = shgat%t1_sou/real(shgat%nt_sou-1,sp)
        call support_funct1( shgat%sou(is,:), shgat%nt_sou, dt_sou, p0_cut, p1_cut, t0_sou_cut, t1_sou_cut )
      end do
      
      !----------------------------------------------------------!
      ! Calibration R per R.
      coef = 0._sp
      denom = 0._sp
      do ir=1,shgat%n_rec
        offset_SR = offset_sou_rec( shgat, ir )
        if ( offset_SR > offset_cal_min .and. offset_SR < offset_cal_max ) then
          
          ! Coef. for the calibration.
          do it=1,shgat%nt_tr_obs
            coef = coef + shgat%pr_dat(ir,it)*shgat%tr_syn(ir,it)*shgat%tr_syn(ir,it)
            denom = denom + shgat%pr_dat(ir,it)*shgat%tr_obs(ir,it)*shgat%tr_obs(ir,it)
          end do
          
        end if
      end do
      
      !----------------------------------------------------------!
      ! Chek error.
      if ( denom == 0._sp ) stop '***** ERROR - coef_calib2: denom == 0. *****'
      if ( coef == 0._sp ) stop '***** ERROR - coef_calib2: coef == 0. *****'
      
      !----------------------------------------------------------!
      shgat%calib = sqrt(coef/denom)
      
      !----------------------------------------------------------!
    end subroutine coef_calib2
    
    
    
    !**********************************************************!
    !**********************************************************!
    !**********************************************************!
    subroutine array_calib_L2_cros( shgat )
      use m_data_kind
      use m_sg_type
      use m_sg_func
      use m_support_funct, only: support_funct1
      use m_type_inversion, only: p0_cut, p1_cut
      use m_phys, only: c_water
      implicit none
      ! The variables which are passed to the function.
      type(shot_gather), intent(inout) :: shgat
      ! The variables which are generated inside the function.
      integer :: ir, is, it
      real(sp) :: denom
      real(sp) :: dt, dt_sou, t_cut, t0_sou_cut, t1_sou_cut
      real(sp) :: dist_SR
      
      !----------------------------------------------------------!
      ! Time grid.
      dt = shgat%t1_obs/real(shgat%nt_tr_obs-1,sp)
      
      ! Support of the observed source wavelet.
      do is=1,shgat%n_sou
        dt_sou = shgat%t1_sou/real(shgat%nt_sou-1,sp)
        call support_funct1( shgat%sou(is,:), shgat%nt_sou, dt_sou, p0_cut, p1_cut, t0_sou_cut, t1_sou_cut )
      end do
      
      !----------------------------------------------------------!
      ! Calibration R per R.
      do ir=1,shgat%n_rec
      
        ! Some info of SR distance in time and space.
        dist_SR = dist_sou_rec( shgat, ir )
        t_cut = t1_sou_cut + dist_SR/c_water
        
        ! Coef. for the calibration.
        shgat%calib(ir) = 0._sp
        denom = 0._sp
        do it=1,shgat%nt_tr_obs
          shgat%calib(ir) = shgat%calib(ir) + shgat%tr_obs(ir,it)*shgat%tr_syn(ir,it)
          denom = denom + shgat%tr_obs(ir,it)*shgat%tr_obs(ir,it)
        end do
        
        shgat%calib(ir) = shgat%calib(ir)/denom
        
        !----------------------------------------------------------!
        ! Chek error.
        if ( denom == 0._sp ) stop '***** ERROR - array_calib_L2_cros: denom == 0. *****'
        if ( shgat%calib(ir) == 0._sp ) stop '***** ERROR - array_calib_L2_cros: shgat%calib(ir) == 0. *****'
        
      end do
      
      !----------------------------------------------------------!
    end subroutine array_calib_L2_cros   
    
    
    
    !**********************************************************!
    !**********************************************************!
    !**********************************************************!
    subroutine array_calib_L2_prec_auto( shgat )
      use m_data_kind
      use m_sg_type
      implicit none
      ! The variables which are passed to the function.
      type(shot_gather), intent(inout) :: shgat
      ! The variables which are generated inside the function.
      integer :: ir, it
      real(sp) :: denom
      
      !----------------------------------------------------------!
      ! Calibration R per R.
      do ir=1,shgat%n_rec
        
        ! Coef. for the calibration.
        shgat%calib(ir) = 0._sp
        denom = 0._sp
        do it=1,shgat%nt_tr_obs
          shgat%calib(ir) = shgat%calib(ir) + shgat%pr_dat(ir,it)*shgat%tr_syn(ir,it)*shgat%tr_syn(ir,it)
          denom = denom + shgat%pr_dat(ir,it)*shgat%tr_obs(ir,it)*shgat%tr_obs(ir,it)
        end do
        
        shgat%calib(ir) = sqrt(shgat%calib(ir)/denom)
        
        !----------------------------------------------------------!
        ! Fix error.
        ! Remove trace if pr_dat*tr_o=0
        if ( denom == 0._sp ) then
          if ( shgat%alg /= 2 ) shgat%kill_tr(ir) = .true.
          shgat%calib(ir) = 0._sp
        end if
        
      end do
      
      !----------------------------------------------------------!
    end subroutine array_calib_L2_prec_auto
    
    
    
    !**********************************************************!
    !**********************************************************!
    !**********************************************************!
    subroutine array_inv_calib_L2_prec_auto( shgat )
      use m_data_kind
      use m_sg_type
      implicit none
      ! The variables which are passed to the function.
      type(shot_gather), intent(inout) :: shgat
      ! The variables which are generated inside the function.
      integer :: ir, it
      real(sp) :: denom
      
      !----------------------------------------------------------!
      ! Calibration R per R.
      do ir=1,shgat%n_rec
        
        ! Coef. for the calibration.
        shgat%calib(ir) = 0._sp
        denom = 0._sp
        do it=1,shgat%nt_tr_obs
          shgat%calib(ir) = shgat%calib(ir) + shgat%pr_dat(ir,it)*shgat%tr_obs(ir,it)*shgat%tr_obs(ir,it)
          denom = denom + shgat%pr_dat(ir,it)*shgat%tr_syn(ir,it)*shgat%tr_syn(ir,it)
        end do
        
        shgat%calib(ir) = sqrt(shgat%calib(ir)/denom)
        
        !----------------------------------------------------------!
        ! Fix error.
        ! Remove trace if pr_dat*tr_s=0
        if ( denom == 0._sp ) then
          if ( shgat%alg /= 2 ) shgat%kill_tr(ir) = .true.
          shgat%calib(ir) = 0._sp
        end if
        
      end do
      
      !----------------------------------------------------------!
    end subroutine array_inv_calib_L2_prec_auto
    
    
    
    !**********************************************************!
    !**********************************************************!
    !**********************************************************!
    subroutine array_calib_L2_prec_cros( shgat )
      use m_data_kind
      use m_sg_type
      implicit none
      ! The variables which are passed to the function.
      type(shot_gather), intent(inout) :: shgat
      ! The variables which are generated inside the function.
      integer :: ir, it
      real(sp) :: denom
      
      !----------------------------------------------------------!
      ! Calibration R per R.
      do ir=1,shgat%n_rec
        
        ! Coef. for the calibration.
        shgat%calib(ir) = 0._sp
        denom = 0._sp
        do it=1,shgat%nt_tr_obs
          shgat%calib(ir) = shgat%calib(ir) + shgat%pr_dat(ir,it)*shgat%tr_obs(ir,it)*shgat%tr_syn(ir,it)
          denom = denom + shgat%pr_dat(ir,it)*shgat%tr_obs(ir,it)*shgat%tr_obs(ir,it)
        end do
        
        shgat%calib(ir) = shgat%calib(ir)/denom
        
        !----------------------------------------------------------!
        ! Fix error.
        ! Remove trace if pr_dat*tr_o=0
        if ( denom == 0._sp ) then
          if ( shgat%alg /= 2 ) shgat%kill_tr(ir) = .true.
          shgat%calib(ir) = 0._sp
        end if
        
      end do
      
      !----------------------------------------------------------!
    end subroutine array_calib_L2_prec_cros
    
    
    
    !**********************************************************!
    !**********************************************************!
    !**********************************************************!
    subroutine array_calib_max_prec( shgat )
      use m_data_kind
      use m_sg_type
      implicit none
      ! The variables which are passed to the function.
      type(shot_gather), intent(inout) :: shgat
      ! The variables which are generated inside the function.
      integer :: ir
      real(sp) :: denom
      
      !----------------------------------------------------------!
      ! Calibration R per R.
      do ir=1,shgat%n_rec
        
        ! Coef. for the calibration.
        shgat%calib(ir) = maxval(abs( shgat%pr_dat(ir,:)*shgat%tr_syn(ir,:) ))
        denom = maxval(abs( shgat%pr_dat(ir,:)*shgat%tr_obs(ir,:) ))
        
        shgat%calib(ir) = shgat%calib(ir)/denom
        
        !----------------------------------------------------------!
        ! Fix error.
        ! Remove trace if pr_dat*tr_o=0
        if ( denom == 0._sp ) then
          if ( shgat%alg /= 2 ) shgat%kill_tr(ir) = .true.
          shgat%calib(ir) = 0._sp
        end if
        
      end do
      
      !----------------------------------------------------------!
    end subroutine array_calib_max_prec
    
    
    
    !**********************************************************!
    !**********************************************************!
    !**********************************************************!
    subroutine array_calib_max( shgat )
      use m_data_kind
      use m_sg_type
      implicit none
      ! The variables which are passed to the function.
      type(shot_gather), intent(inout) :: shgat
      ! The variables which are generated inside the function.
      integer :: ir
      real(sp) :: denom
      
      !----------------------------------------------------------!
      ! Calibration R per R.
      do ir=1,shgat%n_rec
        
        ! Coef. for the calibration.
        shgat%calib(ir) = maxval(abs( shgat%pr_dat(ir,:)*shgat%tr_syn(ir,:) ))
        denom = maxval(abs( shgat%pr_dat(ir,:)*shgat%tr_obs(ir,:) ))
        
        shgat%calib(ir) = shgat%calib(ir)/denom
        
        !----------------------------------------------------------!
        ! Fix error.
        ! Remove trace if pr_dat*tr_o=0
        if ( denom == 0._sp ) then
          if ( shgat%alg /= 2 ) shgat%kill_tr(ir) = .true.
          shgat%calib(ir) = 0._sp
        end if
        
      end do
      
      !----------------------------------------------------------!
    end subroutine array_calib_max
    
    
    
    !**********************************************************!
    !**********************************************************!
    !**********************************************************!
    elemental real(sp) function croscor( shgat, ir )
      use m_data_kind
      use m_sg_type
      implicit none
      ! The variables which are passed to the function.
      type(shot_gather), intent(in) :: shgat
      integer, intent(in) :: ir
      ! The variables which are generated inside the function.
      integer :: it
      real(sp) :: denom1, denom2
      real(sp) :: dt
      
      !----------------------------------------------------------!
      ! Time grid.
      dt = shgat%t1_obs/real(shgat%nt_tr_obs-1,sp)
      
      !----------------------------------------------------------!
!       do ir=1,shgat%n_rec
        
        ! Coef. for the calibration.
        croscor = 0._sp
        denom1 = 0._sp
        denom2 = 0._sp
        do it=1,shgat%nt_tr_obs
          croscor = croscor + shgat%tr_obs(ir,it)*shgat%tr_syn(ir,it)
          denom1 = denom1 + shgat%tr_syn(ir,it)*shgat%tr_syn(ir,it)
          denom2 = denom2 + shgat%tr_obs(ir,it)*shgat%tr_obs(ir,it)
        end do
        
        croscor = croscor/(sqrt(denom1)*sqrt(denom2))
        
        !----------------------------------------------------------!
        ! Chek error.
        if ( denom1 == 0._sp .or. denom2 == 0._sp ) croscor = 0._sp
        if ( croscor < 0._sp ) then
          croscor = 0._sp
        else
          croscor = max( croscor, 0.25_sp )
        end if
        
!       end do
      
      !----------------------------------------------------------!
    end function croscor
    
    
    
  !**********************************************************!
  end module m_calib



