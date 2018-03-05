  !*********************************************************************/
  !** This code has been done in the Barcelona Center for Subsurface 
  !** Imaging (BCSI).
  !** Authors: Daniel Dagnino.
  !*********************************************************************/
  !** Short resume of this code/file:
  !** This subroutine does the following:
  !** 
  !** Loop over shots:
  !**   Forward propagation > Direct field + trace + Illumination.
  !**   Adjoint funtion > Adjoint source + Misfit.
  !**   Backward propagation > Adjoint field + Gradient.
  !** Repeat loop.
  !** 
  !** Sum the contribution of all the shots.
  !** 
  !*********************************************************************/
  module m_solv3Dahv_dir_adj_all
    implicit none
    
    !*********************************************************************/
    contains
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine solv3Dahv_dir_adj_all( n_sg, shgat, buoy, cmpr, freq_0_inv, freq_inv, error_sum, mf_sg_tmp, grad_sum, illu_sum, &
                                      start_grad_sg, grad_sg, mf_sg_a, mf_sg_b, wei_grad )
      use m_mpi
      use m_mat
      use m_data_kind
      use m_work, only: n_shot
      use m_sg_type
      use m_sg_func
      use m_sg_pr
      use m_sg_filt
      use m_sg_pos
      use m_sg_get_picks
      use m_sg_data, only: folder_tmp, folder_pos_wv
      use m_geo
      use m_mygeo
      use m_type_inversion, only: source_inver, prec_d, meth_OF, dyw_Wd, offset_pr_min, offset_pr_max, c_mean
      use m_dyw, only: dyw_read, dyw_read_pos, dyw_reread_pr, ck_mem, dyw_read_J, dyw_reread_tr
      use m_adj_sou_L1
      use m_adj_sou_L2
      use m_adj_sou_L2_array
      use m_adj_sou_p
      use m_adj_sou_e
      use m_adj_sou_L2_water
      use m_adj_sou_L2_water_Wd
      use m_adj_sou_L2_water_Wd_p
      use m_adj_sou_L2_EA
      use m_adj_sou_L2_dip_EA
      use m_adj_sou_L2_dip_EA_AGC
      use m_adj_sou_L2_dip_EA_AGC_env
      use m_adj_sou_L2_dip_EA_AGC_env_2
      use m_adj_sou_L2_dip_EA_AGC_env_3
      use m_adj_sou_L2_EA_cc
      use m_adj_sou_CC_W
      use m_adj_sou_WCC_EA
      use m_solv3Dahvc_data_kind
      use m_solv3Dahvc_time
      use m_solv3Dahv_dir
      use m_solv3Dahv_adj
!       use m_clear_div_sou
      use m_apply_weight
      use m_dag_time
      use m_work, only: n_shot, mean_grad0_sg, sig_grad0_sg, mean_grad0_it0_sg, first_normalization, coef_sea, mf_i0_f0
      use m_type_inversion, only: use_H0
      use m_Hessian
      use m_Hessian_red
      use m_compute_J
      use m_compute_J_nmo
      use m_smooth_op
#ifdef usempi
      use mpi
#endif
      implicit none
      ! The variables which are passed to the function.
      integer, allocatable, intent(in) :: n_sg(:)
      type(shot_gather), allocatable, intent(inout) :: shgat(:)
      real(rl), intent(in) :: freq_0_inv, freq_inv
      real(rl), intent(out) :: error_sum
      real(rl), allocatable, intent(inout) :: grad_sum(:,:,:), illu_sum(:,:,:)
      real(rl), allocatable, intent(in)  :: cmpr(:,:,:), buoy(:,:,:)
      real(rl), allocatable, intent(inout)  :: grad_sg(:,:,:,:)
      logical, intent(inout) :: start_grad_sg
      real(rl), allocatable, intent(inout)  :: mf_sg_tmp(:)
      real(rl), allocatable, intent(in)  :: mf_sg_a(:), mf_sg_b(:)
      logical, intent(in) :: wei_grad
      ! The variables which are generated inside the function.    
      integer :: ix, iy, iz, k, i
      integer :: isg, ir, cont
      real(rl) :: error, mean, mean2, raux, raux1, raux2, norm, g_max, coef
      real(rl) :: freq_cnt
      real(rl), allocatable :: mat_aux_mg(:,:,:)
      real(rl) :: r_max_mpi
      real(rl), allocatable :: v_aux(:)
      real(rl), allocatable :: v_aux2(:)
      real(4), allocatable :: sr_aux(:), sr_aux2(:)
      type(field) :: last_vf
      integer     :: i_save
      integer, allocatable :: ind_save(:)
      real(8) :: time_solver, time_all
      integer :: stal
      character(200) :: memory
      character(100) :: str_tmp
      integer :: ix_shift, iy_shift, iz_shift
      logical :: initialize
      
      logical :: dyw_read_Wd, clean
      
      real(4) :: y
      real(4), allocatable, target :: voffset(:)
      real(4), allocatable, target :: vp_0(:)
      real(4), allocatable, target :: tt_0(:), nvp_0(:)
      
      real(rl) :: vel_z_mean
      real(rl), allocatable :: vel_x(:)
      
      real :: A_random1, A_random2, k_norm
      
!dir$ assume_aligned mf_sg_tmp(1):64
!dir$ assume_aligned grad_sg(1,1,1,1):64
!dir$ assume_aligned illu_sum(1,1,1):64,grad_sum(1,1,1):64
!dir$ assume_aligned cmpr(1,1,1):64,buoy(1,1,1):64
      
      !**********************************************************!
      ! 
      time_all = rtc()
      
      grad_sum  = 0._rl
      illu_sum  = 0._rl
      error_sum = 0._rl
      mf_sg_tmp = 0._rl
      
      allocate( v_aux(nx), sr_aux(nx), sr_aux2(nx), stat=stal ); if ( stal/=0 ) stop 'AE solv3Dahv_dir_adj_all 1'
      v_aux = 0._rl
      sr_aux  = 0.
      sr_aux2 = 0.
      
      initialize = .true.
      dyw_read_Wd = .true.
      
      !**********************************************************!
      if ( .not. start_grad_sg ) then
        grad_sg = 0._rl
      end if
      
      !**********************************************************!
      ! Clean ...
      do k=1,n_shot
        ! Is the shot in this group.
        clean = .true.
        do isg=1,n_sg(rank+1)
          if ( shgat(isg)%ord == k ) then
            clean = .false.
            exit
          end if
        end do
        ! If not, remove.
        if ( clean ) then
          mean_grad0_it0_sg(k) = 0.
          mean_grad0_sg(k) = 0.
          sig_grad0_sg(k) = 0.
        end if
      end do
      
      !**********************************************************!
      ! 
      coef_sea = 0.
      
      if ( rank==0 ) then
        
        ! 
        do i=1,10
          k_norm = 2.*pi*(0.5*dble(i))
          call random_number(A_random1)
          call random_number(A_random2)
          A_random1 = 2.*A_random1-1.
          A_random2 = 2.*A_random2-1.
          do k=1,n_shot
            coef_sea(k) = coef_sea(k) + &
              A_random1*sin( k_norm*dble(k-1)/dble(n_shot-1) ) + &
              A_random2*cos( k_norm*dble(k-1)/dble(n_shot-1) )      
          end do
        end do
        
        coef_sea = coef_sea - minval(coef_sea)
        coef_sea = coef_sea/maxval(abs(coef_sea))
        coef_sea = 0.1 + coef_sea
        
      end if
      
#ifdef usempi
      allocate( v_aux2(n_shot), stat=stal ); if ( stal/=0 ) stop 'AE solv3Dahv_dir_adj_all 2'
      v_aux2 = coef_sea
      call MPI_ALLREDUCE( v_aux2, coef_sea, n_shot, mpi_fwi, MPI_SUM, MPI_COMM_WORLD, ierr )
      deallocate( v_aux2, stat=stal ); if ( stal/=0 ) stop 'dAE solv3Dahv_dir_adj_all 2'
#else
      coef_sea = 1.
#endif
      
      !**********************************************************!
      !**********************************************************!
      !**********************************************************!
      ! Shot.
      do isg=1,n_sg(rank+1)
        !----------------------------------------------------------!
        ! Label.
#ifdef testing
        write(*,'(i4,1x,i4,1x,i4,1x,i4)') rank, shgat(isg)%id, isg, n_sg(rank+1)
#endif
        
        !----------------------------------------------------------!
        ! Direct.
!         write(*,*) 'xxx -2'
        if ( dyw_read_pos ) then
          
          call get_one_sg_pos( folder_pos_wv, shgat(isg), 0._sp, offset_pr_max )
          call get_my_geometry( shgat(isg) )
          call get_one_sg_bath( shgat(isg) )
          
          call allo_one_sg_kill( shgat(isg) )
          call allo_one_sg_calib( shgat(isg) )
          
          if ( meth_OF == 8 .or. meth_OF == 80 .or. meth_OF == 81 .or. meth_OF == 82 .or. meth_OF == 83 .or. meth_OF == 84 ) then
            call allo_one_sg_picks_s( shgat(isg) )
            call allo_one_sg_dly( shgat(isg) )
          else if ( meth_OF == 10 ) then
            call allo_one_sg_picks_s( shgat(isg) )
            call allo_one_sg_dly( shgat(isg) )
            call get_one_picks( folder_pos_wv, shgat(isg) )
          else if ( meth_OF == 11 ) then
            call allo_one_sg_picks_both( shgat(isg) )
            call allo_one_sg_dly( shgat(isg) )
          end if
          
        end if
        
        ! All the models have the same geometry.
        if ( initialize ) then
          allocate( mat_aux_mg(shgat(isg)%mg%iLx,shgat(isg)%mg%iLy,shgat(isg)%mg%iLz), stat=stal ); if ( stal/=0 ) stop 'AE solv3Dahv_dir_adj_all 1'
          mat_aux_mg = 0._rl
          initialize = .false.
        end if
        
        !----------------------------------------------------------!
        ! 
        if ( start_grad_sg ) then
          allocate( grad_sg(n_sg(rank+1),shgat(isg)%mg%iLx,shgat(isg)%mg%iLy,shgat(isg)%mg%iLz), stat=stal ); if ( stal/=0 ) stop 'AE solv3Dahv_dir_adj_all 1'
          start_grad_sg = .false.
          grad_sg = 0._rl
        end if
        
        !----------------------------------------------------------!
        ! 
!         write(*,*) 'xxx -1'
        if ( dyw_read .or. dyw_reread_tr ) then
          call get_one_filt_tr( shgat(isg), freq_0_inv, freq_inv, shgat(isg)%t1_inv )
          if ( .not. source_inver ) call get_one_filt_sou( shgat(isg), freq_inv, shgat(isg)%t1_inv )
        end if
        
        !----------------------------------------------------------!
        ! 
!         freq_cnt = 0.5_rl*( freq_0_inv + freq_inv )
!         write(*,*) 'xxx 0'
        freq_cnt = freq_inv
        call solv3Dahv_dir( time_solver, shgat(isg), buoy, cmpr, freq_cnt, &
                            i_save, ind_save, last_vf, mat_aux_mg )
        
!         ! Weight by distance soure-point*receiver-point.
!         if ( use_H0 ) then
!           call apply_weight( shgat(isg), mat_aux_mg, shgat(isg)%mg )
!         end if
        
        ! 
        coef = 1./maxval(mat_aux_mg)
        do iz=1,shgat(isg)%mg%iLz
        iz_shift = iz+(shgat(isg)%mg%iz-1)
        if ( iz_shift<=nz ) then
        do iy=1,shgat(isg)%mg%iLy
        iy_shift = iy+(shgat(isg)%mg%iy-1)
        if ( iy_shift<=ny ) then
        do ix=1,shgat(isg)%mg%iLx
        ix_shift = ix+(shgat(isg)%mg%ix-1)
        if ( ix_shift<=nx ) then
          illu_sum(ix_shift,iy_shift,iz_shift) = illu_sum(ix_shift,iy_shift,iz_shift) + coef*mat_aux_mg(ix,iy,iz)
        end if
        end do
        end if
        end do
        end if
        end do
        
        !----------------------------------------------------------!
        ! 
        if ( ck_mem .and. rank == 0 ) then
          write(str_tmp,*) shgat(isg)%id
          
          memory = 'du -ch ' // trim(adjustl( folder_tmp )) // '/shot_'
          memory = trim( memory ) // adjustl(str_tmp)
          memory = trim( memory ) // '_k*_step_*.txt'
          
!           memory = 'du -c ' // trim(adjustl( folder_tmp )) // '/shot_'
!           memory = trim( memory ) // adjustl(str_tmp)
!           memory = trim( memory ) // '_k*_step_*.txt'
          
          write(*,*) '----- Memory requiered for one shot (no RAM) ----'
          call system( memory )
          write(*,*) '-------------------------------------------------'
          
          ck_mem = .false.
        end if
        
        !----------------------------------------------------------!
        ! 
        if ( shgat(isg)%alg /= 2 ) shgat(isg)%kill_tr = .false.
        
        !----------------------------------------------------------!
        ! 
!         write(*,*) 'xxx 1'
        if ( dyw_read .or. dyw_reread_pr ) then
          call get_one_sg_pr( shgat(isg), prec_d, freq_inv, offset_pr_min, offset_pr_max )
        end if
        
        !----------------------------------------------------------!
        !----------------------------------------------------------!
        !----------------------------------------------------------!
        if ( dyw_Wd .and. dyw_read_Wd ) then
          ! 
          ! Allocate.
          allocate( voffset(shgat(isg)%n_rec), stat=stal ); if ( stal/=0 ) stop "AE0 solv3Dahv_dir_adj_all"
          
          ! Allocate first.
          if ( isg == 1 ) then
            allocate( vp_0(ny), stat=stal ); if ( stal/=0 ) stop "AE5 solv3Dahv_dir_adj_all"
            allocate( nvp_0(ny), stat=stal ); if ( stal/=0 ) stop "AE7 solv3Dahv_dir_adj_all"
            allocate( tt_0(ny), stat=stal ); if ( stal/=0 ) stop "AE9 solv3Dahv_dir_adj_all"
          end if
          
          ! From m to km.
          do ir=1,shgat(isg)%n_rec
            voffset(ir) = 1000.*offset_sou_rec( shgat(isg), ir )
          end do
          
          ! Mean (in the horizontal axis) velocity.
          vp_0  = 0.
          do iz=1,shgat(isg)%mg%iLz
          iz_shift = iz+(shgat(isg)%mg%iz-1)
          if ( iz_shift<=nz ) then
          do iy=1,shgat(isg)%mg%iLy
          iy_shift = iy+(shgat(isg)%mg%iy-1)
          if ( iy_shift<=ny ) then
          do ix=1,shgat(isg)%mg%iLx
          ix_shift = ix+(shgat(isg)%mg%ix-1)
          if ( ix_shift<=nx ) then
            vp_0(iy) = vp_0(iy) + 1000.*sqrt(buoy(ix_shift,iy_shift,iz_shift)*cmpr(ix_shift,iy_shift,iz_shift))
          end if
          end do
          end if
          end do
          end if
          end do
          vp_0 = vp_0/real(shgat(isg)%mg%iLx*shgat(isg)%mg%iLz)
          
          ! Time of the NMO.
          tt_0 = 0.
          nvp_0(1) = vp_0(1)
          do iy=2,ny
            y = 1000.*dx*dble(iy-1)
            tt_0(iy) = tt_0(iy-1) + 2.*1000.*dx/vp_0(iy-1) 
            nvp_0(iy) = 2.*y/tt_0(iy)
          end do
          
          !----------------------------------------------------------!
          ! Wd.
          call get_cal_Wd( voffset, tt_0, nvp_0, ny, shgat(isg)%t1_obs, shgat(isg)%n_rec, isg )
        end if
        
        !----------------------------------------------------------!
        ! 
!         write(*,*) 'xxx 2'
        if ( meth_OF == 1 ) then
          call adj_sou_L1( shgat(isg), error )
        else if ( meth_OF == 2 ) then
          call adj_sou_L2( shgat(isg), error, isg )
        else if ( meth_OF == 21 ) then
          call adj_sou_L2_array( shgat(isg), error, isg )
        else if ( meth_OF == 3 ) then
          call adj_sou_e( shgat(isg), error, freq_inv )
        else if ( meth_OF == 4 ) then
          call adj_sou_p( shgat(isg), error, freq_inv )
        else if ( meth_OF == 5 ) then
          call adj_sou_L2_water( shgat(isg), error, isg )
        else if ( meth_OF == 6 ) then
          call adj_sou_L2_water_Wd( shgat(isg), error, isg )
        else if ( meth_OF == 7 ) then
          call adj_sou_L2_water_Wd_p( shgat(isg), error, isg )
        else if ( meth_OF == 8 ) then
          call adj_sou_L2_EA( shgat(isg), error, isg )
        else if ( meth_OF == 80 ) then
          call adj_sou_L2_dip_EA( shgat(isg), error, isg )
        else if ( meth_OF == 81 ) then
          call adj_sou_L2_dip_EA_AGC( shgat(isg), error, isg )
        else if ( meth_OF == 82 ) then
          call adj_sou_L2_dip_EA_AGC_env( shgat(isg), error, isg )
        else if ( meth_OF == 83 ) then
          call adj_sou_L2_dip_EA_AGC_env_2( shgat(isg), error, isg )
        else if ( meth_OF == 84 ) then
          call adj_sou_L2_dip_EA_AGC_env_3( shgat(isg), error, isg )
        else if ( meth_OF == 9 ) then
          call adj_sou_L2_EA_cc( shgat(isg), error, isg )
        else if ( meth_OF == 10 ) then
          call adj_sou_CC_W( shgat(isg), error, freq_inv )
        else if ( meth_OF == 11 ) then
          call adj_sou_WCC_EA( shgat(isg), error, freq_inv )
        else
          stop '***** ERROR solv3Dahv_dir_adj_all: Wrong meth_OF value *****'
        end if
        
        if ( error == 0._sp ) then
          write(*,*) " ***** WARNING solv3Dahv_dir_adj_all: error == 0._sp ***** ", shgat(isg)%id
        end if
        
          if ( isnan(real(error)) ) then
            write(*,*) " ***** WARNING solv3Dahv_dir_adj_all: isnan(real(error)) ***** "
            stop " ***** ERROR solv3Dahv_dir_adj_all: isnan(real(error)) ***** "
          end if
          
        !----------------------------------------------------------!
        ! Adjoint
!         write(*,*) 'xxx 3'
        call solv3Dahv_adj( time_solver, shgat(isg), i_save, ind_save, last_vf, mat_aux_mg )
        
        ! 
        if ( error <= 0._sp .or. isnan(real(error)) ) then
          error = 0._sp
          mat_aux_mg = 0._sp
        end if
        
        !----------------------------------------------------------!
        if ( maxval(abs(mat_aux_mg)) == 0._sp ) then
          write(*,*) " ***** WARNING solv3Dahv_dir_adj_all: maxval(abs(mat_aux_mg)) == 0._sp ***** ", shgat(isg)%id
        end if
        
        if ( any(isnan(real(mat_aux_mg))) ) then
          write(*,*) " ***** WARNING solv3Dahv_dir_adj_all: any(isnan(real(mat_aux_mg))) ***** "
          stop " ***** ERROR solv3Dahv_dir_adj_all: any(isnan(real(mat_aux_mg))) ***** "
        end if
        
        !----------------------------------------------------------!
        ! Velocity gradient.
!         write(*,*) 'xxx 4'
        do iz=1,shgat(isg)%mg%iLz
        iz_shift = iz+(shgat(isg)%mg%iz-1)
        if ( iz_shift<=nz ) then
        do iy=1,shgat(isg)%mg%iLy
        iy_shift = iy+(shgat(isg)%mg%iy-1)
        if ( iy_shift<=ny ) then
        do ix=1,shgat(isg)%mg%iLx
        ix_shift = ix+(shgat(isg)%mg%ix-1)
        if ( ix_shift<=nx ) then
          mat_aux_mg(ix,iy,iz) = -mat_aux_mg(ix,iy,iz)*sqrt(cmpr(ix_shift,iy_shift,iz_shift))
        end if
        end do
        end if
        end do
        end if
        end do
        
        !----------------------------------------------------------!
        ! Grad manipulation.
!         ! 
!         call clear_div_sou( shgat(isg), mat_aux_mg )
        
        ! Hess.
        if ( use_iH ) then
          if ( .not. full_H ) call mult_inv_H_v_diag( shgat(isg)%mg%iLx, shgat(isg)%mg%iLy, shgat(isg)%mg%iLz, mat_aux_mg )
        end if
        
#ifdef valgrind
      mat_aux_mg = 1._sp
#endif
        
        !write(*,*) 'xxx 5'
        ! Weight
        if ( use_H0 ) then
          call apply_weight( shgat(isg), mat_aux_mg, shgat(isg)%mg )
        end if
        
        !----------------------------------------------------------!
        ! 
        allocate( vel_x(shgat(isg)%mg%iLy), stat=stal ); if ( stal/=0 ) stop 'AE solv3Dahv_dir_adj_all 5'
        
        !----------------------------------------------------------!
        ! Smooth operator.
!         call apply_smooth_op_small( shgat(isg), mat_aux_mg, freq_inv )
        
        vel_x = 0._rl
        do iz=1,shgat(isg)%mg%iLz
        iz_shift = iz+(shgat(isg)%mg%iz-1)
        if ( iz_shift<=nz ) then
        do iy=1,shgat(isg)%mg%iLy
        iy_shift = iy+(shgat(isg)%mg%iy-1)
        if ( iy_shift<=ny ) then
        do ix=1,shgat(isg)%mg%iLx
        ix_shift = ix+(shgat(isg)%mg%ix-1)
        if ( ix_shift<=nx ) then
          vel_x(iy) = vel_x(iy) + sqrt(buoy(ix_shift,iy_shift,iz_shift)*cmpr(ix_shift,iy_shift,iz_shift))
        end if
        end do
        end if
        end do
        end if
        end do
        vel_x = vel_x/real(shgat(isg)%mg%iLx*shgat(isg)%mg%iLz)
        
!         vel_z_mean = c_mean
        vel_z_mean = sum( vel_x )/real(shgat(isg)%mg%iLy)
        
        call apply_smooth_op_small( shgat(isg), mat_aux_mg, freq_inv, vel_z_mean, vel_x )
        
        !----------------------------------------------------------!
        !write(*,*) 'xxx 6'
        ! Some gradient statistics at first iter (no normalization).
        first_normalization = .true.
        if ( first_normalization ) then
          g_max = 0._sp
          cont  = 0
          mean  = 0._sp
          mean2 = 0._sp
          do iz=1,shgat(isg)%mg%iLz
          iz_shift = iz+(shgat(isg)%mg%iz-1)
          if ( iz_shift<=nz ) then
          do iy=1,shgat(isg)%mg%iLy
          iy_shift = iy+(shgat(isg)%mg%iy-1)
          if ( iy_shift<=ny ) then
          do ix=1,shgat(isg)%mg%iLx
          ix_shift = ix+(shgat(isg)%mg%ix-1)
          if ( ix_shift<=nx ) then
            raux = mat_aux_mg(ix,iy,iz)
!             if ( .not. isnan(real(raux)) ) then
              mean  = mean  + abs(raux)
              mean2 = mean2 + raux*raux
              g_max = max( g_max, abs(raux) )
!             end if
            cont = cont + 1
          end if
          end do
          end if
          end do
          end if
          end do
          
          !----------------------------------------------------------!
          if ( g_max == 0._sp ) then
            write(*,*) " ***** WARNING solv3Dahv_dir_adj_all: g_max == 0._sp ***** ", shgat(isg)%id
!             stop " ***** ERROR solv3Dahv_dir_adj_all: g_max == 0._sp ***** "
            g_max = 1.d9
          end if
          
          if ( isnan(real(g_max)) ) then
            write(*,*) " ***** WARNING solv3Dahv_dir_adj_all: isnan(real(g_max)) ***** "
            stop " ***** ERROR solv3Dahv_dir_adj_all: isnan(real(g_max)) ***** "
          end if
          
          mean  = mean/real(cont,sp)
          mean2 = mean2/real(cont,sp)
          mean_grad0_sg(shgat(isg)%ord) = mean
          sig_grad0_sg(shgat(isg)%ord)  = sqrt(mean2-mean*mean)
          
!         open(unit=1,file='grad_test.txt',status='unknown',action='write')
!         do iy=1,shgat(isg)%mg%iLy
!           write(1,'(100000(es13.6E2,1x))') (real(mat_aux_mg(ix,iy,1)),ix=1,shgat(isg)%mg%iLx)
!         end do
!         close(1)
          
          ! Saving the first mean of the grad to normalized it at the first iter.
!           mean_grad0_it0_sg(shgat(isg)%ord) = 1._sp/g_max
          if ( error > 0. ) then
            mean_grad0_it0_sg(shgat(isg)%ord) = sqrt(abs(error))/g_max
          else
            mean_grad0_it0_sg(shgat(isg)%ord) = sqrt(abs(1./error))/g_max
          end if
        end if
        
        !----------------------------------------------------------!
        ! Normalize in function of the mean of the first gradient.
        norm = mean_grad0_it0_sg(shgat(isg)%ord)
        
        !----------------------------------------------------------!
!         mf_sg_tmp(shgat(isg)%ord) = norm*error
        mf_sg_tmp(shgat(isg)%ord) = error
        
        !----------------------------------------------------------!
!         error_sum = error_sum + norm*error
        error_sum = error_sum + error
        
!         !----------------------------------------------------------!
!         ! 
! !         if ( wei_grad ) then
! ! !           coef_sea = (mf_sg_a - mf_sg_b)/mf_sg_a
! ! !           raux1 = sum(coef_sea*coef_sea)
! ! !           coef_sea = (mf_sg_b - mf_sg_tmp/mf_i0_f0)/mf_sg_b
! ! !           raux2 = sum(coef_sea*coef_sea)
! ! !           coef_sea = (mf_sg_tmp-mf_sg_b) + (mf_sg_b-mf_sg_a)*(raux2/raux1)
! !           coef_sea(shgat(isg)%ord) = (mf_sg_b(shgat(isg)%ord) - mf_sg_tmp(shgat(isg)%ord)/mf_i0_f0)/mf_sg_b(shgat(isg)%ord) + &
! !             (mf_sg_a(shgat(isg)%ord) - mf_sg_b(shgat(isg)%ord))/mf_sg_a(shgat(isg)%ord)
! !           coef_sea(shgat(isg)%ord) = max( -0.5_sp, coef_sea(shgat(isg)%ord) )
! !           coef_sea(shgat(isg)%ord) = 1._sp + 1._sp*coef_sea(shgat(isg)%ord)
! !         else
!           coef_sea(shgat(isg)%ord) = 1._sp
! !         end if
        
        !----------------------------------------------------------!
        !write(*,*) 'xxx 7'
        ! Gradient statistics normnalized.
        cont  = 0
        mean  = 0._sp
        mean2 = 0._sp
        do iz=1,shgat(isg)%mg%iLz
        iz_shift = iz+(shgat(isg)%mg%iz-1)
        if ( iz_shift<=nz ) then
        do iy=1,shgat(isg)%mg%iLy
        iy_shift = iy+(shgat(isg)%mg%iy-1)
        if ( iy_shift<=ny ) then
        do ix=1,shgat(isg)%mg%iLx
        ix_shift = ix+(shgat(isg)%mg%ix-1)
        if ( ix_shift<=nx ) then
          raux = norm*mat_aux_mg(ix,iy,iz)
!           if ( .not. isnan(real(raux)) ) then
            grad_sg(isg,ix,iy,iz) = raux
            grad_sum(ix_shift,iy_shift,iz_shift) = grad_sum(ix_shift,iy_shift,iz_shift) + coef_sea(shgat(isg)%ord)*raux
            mean  = mean  + abs(raux)
            mean2 = mean2 + raux*raux
!           end if
          cont = cont + 1
        end if
        end do
        end if
        end do
        end if
        end do
        
        mean  = mean/real(cont,sp)
        mean2 = mean2/real(cont,sp)
        mean_grad0_sg(shgat(isg)%ord) = mean
        sig_grad0_sg(shgat(isg)%ord)  = sqrt(mean2-mean*mean)
        
        !----------------------------------------------------------!
        if ( dyw_Wd .and. dyw_read_Wd ) then
          deallocate( voffset, stat=stal ); if ( stal/=0 ) stop 'dAE0 solv3Dahv_dir_adj_all'
        end if
        
        deallocate( vel_x, stat=stal ); if ( stal/=0 ) stop 'dAE solv3Dahv_dir_adj_all 5'
        
      end do
      
      !**********************************************************!
      !**********************************************************!
      !**********************************************************!
      ! 
      dyw_read = .false.
      dyw_read_pos = .false.
      dyw_read_J = .false.
      dyw_reread_tr = .false.
      
      first_normalization = .false.
      
      !**********************************************************!
      ! 
      ! Sum gradient subsums for each process.
! #ifdef usempi
!       mat_aux = grad_sum
!       if ( sp == 4 ) then
!         call MPI_ALLREDUCE( mat_aux, grad_sum, nx*ny*nz, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr )
!       else if ( sp == 8 ) then
!         call MPI_ALLREDUCE( mat_aux, grad_sum, nx*ny*nz, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr )
!       else
!         stop '***** ERROR - solv3Dahv_dir_adj_all: MPI_ALLREDUCE grad: sp wrong *****'
!       end if
! #endif

      !**********************************************************!
      ! 
      ! Sum gradient subsums for each process.
        !write(*,*) 'xxx 8'
! #ifdef usempi
!       do iy=1,ny
!       do iz=1,nz
!         v_aux = grad_sum(:,iy,iz)
!         call MPI_ALLREDUCE( v_aux, grad_sum(:,iy,iz), nx, mpi_fwi, MPI_SUM, MPI_COMM_WORLD, ierr )
!       end do
!       end do
!       
!       ! Sum gradient subsums for each process.
!       do iy=1,ny
!       do iz=1,nz
!         v_aux = illu_sum(:,iy,iz)
!         call MPI_ALLREDUCE( v_aux, illu_sum(:,iy,iz), nx, mpi_fwi, MPI_SUM, MPI_COMM_WORLD, ierr )
!       end do
!       end do
! #endif

#ifdef usempi
      ! Max value of all grads.
      raux = maxval(abs(grad_sum))
      call MPI_ALLREDUCE( raux, r_max_mpi, 1, mpi_fwi, MPI_MAX, MPI_COMM_WORLD, ierr )
      
      ! Sum gradient subsums for each process + recover norm.
      do iy=1,ny
      do iz=1,nz
        ! Reduce precision.
        sr_aux = grad_sum(:,iy,iz)/r_max_mpi
        ! Send + sum.
        call MPI_ALLREDUCE( sr_aux, sr_aux2, nx, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr )
        ! Recover precision.
        grad_sum(:,iy,iz) = r_max_mpi*dble(sr_aux2)
      end do
      end do
      
      !----------------------------------------------------------!
      ! Max value of all illus.
      raux = maxval(abs(illu_sum))
      call MPI_ALLREDUCE( raux, r_max_mpi, 1, mpi_fwi, MPI_MAX, MPI_COMM_WORLD, ierr );
      
      ! Sum illumination subsums for each process.
      do iy=1,ny
      do iz=1,nz
        ! Reduce precision.
        sr_aux = illu_sum(:,iy,iz)/r_max_mpi
        ! Send + sum.
        call MPI_ALLREDUCE( sr_aux, sr_aux2, nx, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr )
        ! Recover precision.
        illu_sum(:,iy,iz) = r_max_mpi*dble(sr_aux2)
      end do
      end do
#endif
      
        !write(*,*) 'xxx 9'
      !*********************************************************************/      
      if ( dyw_Wd .and. dyw_read_Wd ) then
        deallocate( vp_0, stat=stal ); if ( stal/=0 ) stop 'dAE5 solv3Dahv_dir_adj_all'
        deallocate( tt_0, stat=stal ); if ( stal/=0 ) stop 'dAE7 solv3Dahv_dir_adj_all'
        deallocate( nvp_0, stat=stal ); if ( stal/=0 ) stop 'dAE9 solv3Dahv_dir_adj_all'
      end if
      
      !**********************************************************!
      ! Deallocate.
      !write(*,*) 'xxx 10'
      deallocate( v_aux, sr_aux, sr_aux2, stat=stal ); if ( stal/=0 ) stop 'dAE solv3Dahv_dir_adj_all 1'
      deallocate( mat_aux_mg, stat=stal ); if ( stal/=0 ) stop 'dAE solv3Dahv_dir_adj_all 1b'
      
    !*********************************************************************/
    end subroutine solv3Dahv_dir_adj_all
  
  !***********************************************************************/
  end module m_solv3Dahv_dir_adj_all    




