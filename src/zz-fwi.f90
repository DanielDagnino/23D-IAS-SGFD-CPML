  !*********************************************************************/
  !** This code has been done in the Barcelona Center for Subsurface 
  !** Imaging (BCSI).
  !** Authors: Daniel Dagnino.
  !*********************************************************************/
  !** Short resume of this code/file:
  !** Main FWI program. Here, we performe the FWI by calling the necessary
  !** subroutines/functions. The algorithm is basically (without writes)
  !** the following:
  !** 
  !** Loop over freq. + loop over constant freq.:
  !**   Source inversion.
  !**   Solve wave eq. forward+adjoint: gradient, misfit and trace calculation.
  !**   Calculate the search direction: Preconditioning + Smoothing + Optimization.
  !**   Calculate initial step length.
  !**   Minimization: Line search.
  !**   Check confident creterium.
  !** Repeat for the loops.
  !** 
  !*********************************************************************/
  module m_fwi
    implicit none
    
    !*********************************************************************/
    contains
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine fwi( n_sg, shgat, n_freq, freq_0_k, freq_k, iter_k )
      use m_mpi
      use m_dyw
      use m_geo
      use m_wrt_tr
      use m_sg_type
      use m_sg_data, only: folder_output, folder_strat
      use m_sg_path_file
      use m_model
      use m_type_inversion
      use m_filter
      use m_variable
      use m_solv3Dahv_dir_adj_all
      use m_change_var
      use m_copy_mat_2D
      use m_copy_mat_4D
      use m_copy_mat_3D
      use m_relat_change
      use m_alp_cota
      use m_normalize_grad
      use m_search_min
      use m_minimize
      use m_mod_prec
      use m_coef_norm_grad
      use m_control_grad
      use m_solv3Dahv_deallo
      use m_Hessian
      use m_poly_fit_min
      use m_source_inversion
      use m_smooth_op
      use m_inv_data
      use m_work
      use m_solv3Dahv_test_all
      use m_Hessian_red
      use m_param_control
#ifdef usempi
      use mpi
#endif
      implicit none
      
      !**********************************************************!
      ! The variables that are passed to the function.
      integer, allocatable, intent(in)              :: n_sg(:)
      type(shot_gather), allocatable, intent(inout) :: shgat(:)
      integer, intent(in)               :: n_freq
      real(sp), allocatable, intent(in) :: freq_0_k(:), freq_k(:)
      integer, allocatable, intent(in)  :: iter_k(:)
      
      !**********************************************************!
      ! The variables that are generated inside the function. 
      logical  :: wol_cur, wol_cur_str, wol_arm
      logical  :: inv_QC, pol_min_QC, succes_red
      
      integer :: ny0
      integer :: k, ix, iy, iz, ifr, isg, iwtr
      
      real(sp) :: raux
      real(sp) :: mf_0, mf_1, mf_2, mf_fit
      real(sp) :: mf_min, mf_0_old
      real(sp) :: perc_max_alp_cota, perc_min_alp_cota
      real(sp) :: mod_change
      real(sp) :: alp_1, alp_2, alp_fit, alp_min, alp_cota, alp_scaled
      real(sp) :: prod_grad_search, prod_grad_search0
      real(sp) :: norm_grad
      real(sp) :: freq_0_inv
      real(sp) :: mean, sigma, z_mean
      real(sp) :: dx0, dx01, dy0, dy1, dz0, dz01
      
      integer :: n_saved, ind_saved, n_max_saved
      real(sp), allocatable :: grad(:,:,:,:), var(:,:,:,:)
      real(sp), allocatable :: grad_aux(:,:,:), mat_aux(:,:,:), illu(:,:,:)
      real(sp), allocatable :: v_sg_aux(:)
      real(sp), allocatable :: grad_sg(:,:,:,:)
      logical :: start_grad_sg
      
      integer :: stal
      character(200) :: path_file     
      character(50)  :: str_tmp
      
      logical :: wei_grad
      
      !**********************************************************!
      ! Some initial variables.
      
      ! Number maximum of gradients and models to be save for the l-BFGS.
      n_max_saved = 5
      
      ! Change at the first step.
!       change_cota = 0.005_sp   ! High-freq.
!       change_cota = 0.10_sp   ! High-freq.
!       change_cota = 0.050_sp   ! Background.
!       change_cota = 0.10_sp   ! Background: CC_W
      perc_max_alp_cota = 0.40_sp
      perc_min_alp_cota = 0.05_sp
      
      ! Border of the gradient that it will be reduced.
      dx0  = 0.100_sp   ! x0 left and rigth.
      dx01 = 0.000_sp   ! Increment at depth.
      
      dz0  = 0.000_sp   ! z0 back and front.
      dz01  = 0.000_sp   ! Increment at depth.
      
!       dy0  = 0.090_sp   ! y up.
      dy0  = 0.000_sp   ! y up.
      dy1  = 0.000_sp   ! y down.
      
      !**********************************************************!
      wol_cur     = .false.
      wol_cur_str = .false.
      wol_arm     = .false.
      
      prod_grad_search0 = 0._sp
      prod_grad_search  = 0._sp
      
      !**********************************************************!
      ! Files.
      if ( rank == 0 ) then
        write(str_tmp,*)  strat_inv
        path_file = trim(adjustl( folder_output )) // '/strat_inv=' // trim(adjustl(str_tmp)) // '.txt'
        open(unit=1,file=path_file,status='unknown',action='write')
          write(1,*) '************************************************'
          write(1,*) ' ----- Type of inversion ----- '
          write(1,'(a,i3)') ' lay_frz  = ', lay_frz
          write(1,'(a,i3)') ' meth_inv = ', meth_inv
          write(1,*) 'source_inver = ', source_inver
          write(1,*) 'fast_search  = ', fast_search
          write(1,*) ' ----- Frequencies ----- '
          write(1,*) 'ifr  iter_k   freq_0_k    freq_k'
          do ifr=1,n_freq
            write(1,'(i4,4x,i4,1x,f10.2,f10.2)') ifr, iter_k(ifr), freq_0_k(ifr), freq_k(ifr)
          end do
          write(1,*) '************************************************'
          write(1,*) ' ----- Acquisition time for the inversion ----- '
          write(1,'(a,f10.2,a)') ' shgat(1)%t1_inv (s) = ', shgat(1)%t1_inv
          write(1,*) '************************************************'
          write(1,*) ' ----- Space ----- '
          write(1,'(a,f10.2,a)') ' dx (m) = ', 1000._sp*dx
          write(1,'(a,i8)') ' nx = ', nx
          write(1,'(a,i8)') ' ny = ', ny
          write(1,'(a,i8)') ' nz = ', nz
          write(1,*) '************************************************'
          write(1,*) ' ----- Filter ----- '
          write(1,'(a,i3)') ' zero_pad_time  = ', zero_pad_time
          write(1,'(a,i3)') ' zero_pad_space = ', zero_pad_space
          write(1,'(a,i3)') ' order_filt     = ', order_filt
          write(1,*) '************************************************'
          write(1,*) ' ----- Writing options ----- '
          write(1,*) 'write_model   = ', write_model
          write(1,*) 'write_trace   = ', write_trace
          write(1,*) 'write_sou     = ', write_sou
          write(1,*) 'write_adj_sou = ', write_adj_sou
          write(1,*) 'write_grad    = ', write_grad
          write(1,*) 'write_sear    = ', write_sear
          write(1,*) '************************************************'
        close(1)
      end if
      
      !**********************************************************!
      ! Files.
      if ( rank == 0 ) then
        write(str_tmp,*) strat_inv
        path_file = trim(adjustl( folder_output )) // '/resume/strat_inv=' // trim(adjustl(str_tmp)) // '.txt'
        open(unit=101,file=path_file,status='unknown',action='write')
      end if
      
      !**********************************************************!
      ! Initialize variables.
      norm_grad0 = 0._sp
      
      mf_0   = 0._sp
      mf_1   = 0._sp
      mf_2   = 0._sp
      mf_fit = 0._sp
      mf_min = 0._sp
      
      alp_cota  = 0._sp
      alp_1     = 0._sp
      alp_2     = 0._sp
      alp_fit   = 0._sp
      alp_min   = 0._sp
      alp_scaled = 0._sp
      
      !**********************************************************!
      ! Allocate.
      allocate( m_0(nx,ny,nz), m_1(nx,ny,nz), m_2(nx,ny,nz), m_fit(nx,ny,nz), m_min(nx,ny,nz), stat=stal ); if ( stal/=0 ) stop 'AE fwi2'
      m_0   = 0._sp
      m_1   = 0._sp
      m_2   = 0._sp
      m_fit = 0._sp
      m_min = 0._sp
      allocate( grad_aux(nx,ny,nz), search(nx,ny,nz), stat=stal ); if ( stal/=0 ) stop 'AE fwi6'
      grad_aux = 0._sp
      search   = 0._sp
      allocate( mat_aux(nx,ny,nz),  stat=stal ); if ( stal/=0 ) stop 'AE fwi6c'
      mat_aux  = 0._sp
      allocate( illu(nx,ny,nz), prec_mod(nx,ny,nz), stat=stal ); if ( stal/=0 ) stop 'AE fwi6d'
      illu     = 0._sp
      prec_mod = 0._sp
      allocate( mf_sg(n_shot), mf_sg_1(n_shot), mf_sg_2(n_shot), mf_sg_fit(n_shot), mf_sg_ck(n_shot), mf_sg_min(n_shot), stat=stal ); if ( stal/=0 ) stop 'AE fwi6e'
      allocate( mf_sg_a(n_shot), mf_sg_b(n_shot), stat=stal ); if ( stal/=0 ) stop 'AE fwi6e'
      allocate( es_sg(n_shot), mean_grad0_it0_sg(n_shot), mean_grad0_sg(n_shot), sig_grad0_sg(n_shot), stat=stal ); if ( stal/=0 ) stop 'AE fwi6e'
      allocate( coef_sea(n_shot), stat=stal ); if ( stal/=0 ) stop 'AE fwi6f'
      mean_grad0_it0_sg = 0._sp
      mf_sg         = 0._sp
      mf_sg_a       = 0._sp
      mf_sg_b       = 0._sp
      mf_sg_1       = 0._sp
      mf_sg_2       = 0._sp
      mf_sg_fit     = 0._sp
      mf_sg_ck      = 0._sp
      mf_sg_min     = 0._sp
      coef_sea      = 0._sp
      es_sg         = 0._sp
      mean_grad0_sg = 0._sp
      sig_grad0_sg  = 0._sp
      allocate( v_sg_aux(n_shot), stat=stal ); if ( stal/=0 ) stop 'AE fwi6e'
      v_sg_aux = 0._sp
      
      !**********************************************************************!
      !**********************************************************************!
      !*** Loop for the inversion of different frequencies (multi-scale). ***!
      !**********************************************************************!
      !**********************************************************************!
      if ( rank == 0 ) write(*,*) 'Loop over all frequencies:'
      
      ! Initial variables.
      dyw_read_pos = .true.
      iter = 0
      ifr  = 1
      
      start_grad_sg = .true.
      
      !**********************************************************************!
      ! Recover the model in the non-canonical variables.
      call change_var_2( cmpr, buoy , m_0 )
      
      !**********************************************************************!
      ! Loop for the inversion of different frequencies (multi-scale).
      do ifr=1,n_freq
        !----------------------------------------------------------!
        ! Frequency
        freq_0_inv = freq_0_k(ifr)
        freq_inv = freq_k(ifr)
        dyw_read = .true.
        dyw_reread_tr = .true.
        first_normalization = .true.
        
        !----------------------------------------------------------!
        ! Number of grad and model necessary to save according with the selected minimization model.
        if ( meth_inv == 0 ) then
          n_saved = 1
        else if ( meth_inv == 1 ) then
          n_saved = 2
        else if ( meth_inv == 2 ) then
          n_saved = 0
          do k=1,n_freq
            n_saved = max( iter_k(k), n_saved )
          end do
          n_saved = min( n_max_saved, n_saved ) 
        else
          stop ' ********** ERROR: meth_inv is not a correct value. ********** '
        end if
        
        !----------------------------------------------------------!
        ! Allocate.
        if ( .not. allocated( var ) ) then
          allocate( var(n_saved,nx,ny,nz), stat=stal ); if ( stal/=0 ) stop 'AE fwi5'
          var = 0._sp
        end if
        
        if ( .not. allocated( grad ) ) then
          allocate( grad(n_saved,nx,ny,nz), stat=stal ); if ( stal/=0 ) stop 'AE fwi5'
          grad = 0._sp
        end if
        
        ! 
!         if ( dyw_Wd .and. dyw_read_Wd ) then
!           call get_Wd( folder_strat )
!           dyw_read_Wd = .false.
!         end if
        
        !----------------------------------------------------------!
        ! 
        if ( rank == 0 ) then
          raux = sqrt(minval(cmpr*buoy))
          write(*,*) '******************************************************************************'
          write(*,'(a,f10.2,a)') ' dx       = ', 1000._sp*dx, ' m'
          write(*,'(a,f9.2,a)')  ' freq_inv = ', freq_inv, ' Hz'
          write(*,'(a,f6.2,a)')  ' Advise (relaxed) dx < (min(vP)/(1.0*freq_inv))/5 = ', 1000.*(raux/(1.0*freq_inv))/5., ' m'
          write(*,'(a,f6.2,a)')  ' Advise (middle)  dx < (min(vP)/(2.0*freq_inv))/5 = ', 1000.*(raux/(2.0*freq_inv))/5., ' m'
          write(*,'(a,f6.2,a)')  ' Advise (strong)  dx < (min(vP)/(2.5*freq_inv))/5 = ', 1000.*(raux/(2.5*freq_inv))/5., ' m'
          write(*,*) '******************************************************************************'
          
          write(101,*) '******************************************************************************'
          write(101,'(a,f10.2,a)') ' dx       = ', 1000._sp*dx, ' m'
          write(101,'(a,f9.2,a)')  ' freq_inv = ', freq_inv, ' Hz'
          write(101,'(a,f6.2,a)')  ' Advise (relaxed) dx < (min(vP)/(1.0*freq_inv))/5 = ', 1000.*(raux/(1.0*freq_inv))/5., ' m'
          write(101,'(a,f6.2,a)')  ' Advise (middle)  dx < (min(vP)/(2.0*freq_inv))/5 = ', 1000.*(raux/(2.0*freq_inv))/5., ' m'
          write(101,'(a,f6.2,a)')  ' Advise (strong)  dx < (min(vP)/(2.5*freq_inv))/5 = ', 1000.*(raux/(2.5*freq_inv))/5., ' m'
          write(101,*) '******************************************************************************'
        end if
        
        !********************************************************!
        !********************************************************!
        !*** Loop: Inversion for a frequency band. **************!
        !********************************************************!
        !********************************************************!
        inv_QC    = .true.   ! Used in the first iteration.
        ind_saved = 0
        iter_freq = 0
        
        do while ( iter_freq + 1 <= iter_k(ifr) )
          !----------------------------------------------------------!
          ! Counting the iterations, total and by each frequency.
          iter      = iter + 1
          iter_freq = iter_freq + 1
          
          !----------------------------------------------------------!
          ! Iteration number and grad saved.
          
          ! For source inversion and first iter ...
          if ( iter_freq <= 1 ) then
            
            ind_saved = 1
            inv_QC    = .true.
            
          ! Then ...
          else
            
            ! If previous iter ok we make another.
            if ( inv_QC ) then
              ind_saved = min( n_saved, ind_saved + 1 )
            ! Otherwise we finish here the loop and remove the last iter.
            else
              iter      = iter - 1
              iter_freq = iter_freq - 1
              exit
            end if
            
          end if
          
          !----------------------------------------------------------!
          ! Write.
          if ( rank == 0 ) then
            write(101,*) '******************************************************************************'
            write(101,*) 'iter  fr0_inv  fr1_inv  it_fr  ind_sav  it_k(ifr)'
            write(101,'(i5,3x,f6.2,3x,f6.2,3x,i4,5x,i4,7x,i4)') iter, freq_0_inv, freq_inv, iter_freq, ind_saved, iter_k(ifr)
            
            write(*,*) '******************************************************************************'
            write(*,*) 'iter  fr0_inv  fr1_inv  it_fr  ind_sav  it_k(ifr)'
            write(*,'(i5,3x,f6.2,3x,f6.2,3x,i4,5x,i4,7x,i4)') iter, freq_0_inv, freq_inv, iter_freq, ind_saved, iter_k(ifr)
            write(*,*) '******************************************************************************'
          end if
          
          !----------------------------------------------------------!          
          ! Source inversion.
          if ( source_inver ) then
            
            if ( rank == 0 ) write(*,*) 'source inver'
            
            ! Source inversion.
            es_sg = 0._sp
            shgat(:)%alg = 0
            call source_inversion( n_sg, shgat, freq_0_inv, freq_inv, iter_freq, source_fix )
            
            ! Write starting point.
            if ( iter_freq == 1 ) then
              ! 
              if ( rank == rank_wrt ) then
                
                ! Write initial source.
                if ( write_sou ) then
                  call copy_mat_2D_unformatted( strat_inv, freq_inv, 0, '/source/sou', shgat(isg_wrt)%sou, shgat(isg_wrt)%t1_sou/real(shgat(isg_wrt)%nt_sou-1,sp) )
                end if
                
!                 ! Write initial tr (calculate with the initial source).
!                 if ( write_trace ) then
!                   call copy_mat_2D_unformatted( strat_inv, freq_inv, 0, '/trace/tr_s', shgat(isg_wrt)%tr_syn, shgat(isg_wrt)%t1_syn/real(shgat(isg_wrt)%nt_tr_syn-1,sp) )
!                   call copy_mat_2D_unformatted( strat_inv, freq_inv, 0, '/trace/tr', shgat(isg_wrt)%tr_obs, shgat(isg_wrt)%t1_obs/real(shgat(isg_wrt)%nt_tr_obs-1,sp) )
!                   call copy_mat_2D_unformatted( strat_inv, freq_inv, 0, '/adj_sou/prec_data', shgat(isg_wrt)%pr_dat, shgat(isg_wrt)%t1_obs/real(shgat(isg_wrt)%nt_tr_obs-1,sp) )
!                 end if
                
              end if
              
              ! 
              ! Write initial tr (calculate with the initial source).
              if ( write_trace ) then
                do iwtr=1,n_id_wrt_tr
                do isg=1,n_sg(rank+1)
                  if ( shgat(isg)%id == id_wrt_tr(iwtr) ) then
                    call copy_mat_2D_B_unformatted( strat_inv, freq_inv, 0, '/trace/tr_s_end', shgat(isg)%id, shgat(isg)%tr_syn, shgat(isg)%t1_obs/real(shgat(isg)%nt_tr_obs-1,sp) )
                    call copy_mat_2D_B_unformatted( strat_inv, freq_inv, 0, '/trace/tr_s_end', shgat(isg)%id, shgat(isg)%tr_obs, shgat(isg)%t1_obs/real(shgat(isg)%nt_tr_obs-1,sp) )
                    call copy_mat_2D_B_unformatted( strat_inv, freq_inv, 0, '/trace/prec_data_end', shgat(isg)%id, shgat(isg)%pr_dat, shgat(isg)%t1_obs/real(shgat(isg)%nt_tr_obs-1,sp) )
                  end if
                end do
                end do
              end if
              
            end if
            
          end if
          
          !**********************************************************!
          !* 0er step - Gradient. ***********************************!
          !**********************************************************!
          
          !***************************************************!
          ! Hessian.
          if ( use_iH .and. iter_freq == 1 ) then
            if ( .not. full_H ) then
              call get_iH_diag( shgat(1)%mg%iLx, shgat(1)%mg%iLy, shgat(1)%mg%iLz, folder_strat )
!             else
!               call get_iH( folder_strat )
            end if
          end if
          
          !***************************************************!
          ! Gradient calculation.
          if ( rank == 0 ) write(*,*) 'solv3Dahv_dir_adj_all'
          
          shgat(:)%alg = 1
          if ( iter_freq >= 3 ) then
            wei_grad = .true.
          else
            wei_grad = .false.
          end if
          call solv3Dahv_dir_adj_all( n_sg, shgat, buoy, cmpr, freq_0_inv, freq_inv, &
                                      mf_0, mf_sg, grad_aux, illu, start_grad_sg, grad_sg, mf_sg_a, mf_sg_b, wei_grad )
          
          !***************************************************!
          ! Model preconditioner.
          if ( iter_freq == 1 ) then
            ! get_m_prec
            call get_mod_prec( prec_m, prec_mod )
            ! apply_border_red
            if ( border_red ) call apply_border_red( prec_mod, dx0, dx01, dy0, dy1, dz0, dz01 )
          end if
          
          ! Model preconditioner.
          grad_aux = grad_aux*prec_mod
          
          !***************************************************!
          ! Smooth operator.          
          call apply_smooth_op( grad_aux, dy0, freq_inv )
          
          !***************************************************!
          ! Clean water + Layer frezzing.
          
          do iz=1,nz
          do ix=1,nx
            ny0 = grid%ind_bath(ix,iz)-1
            do iy=1,ny0+lay_frz
              prec_mod(ix,iy,iz) = 1._rl
              illu(ix,iy,iz)  = 0._rl
            end do
          end do
          end do
          
          do iz=1,nz
          do ix=1,nx
            ny0 = grid%ind_bath(ix,iz)-1
            do iy=1,ny0+lay_frz
              grad_aux(ix,iy,iz)  = 0._rl
            end do
          end do
          end do
          
          !***************************************************!
          ! Grad.
          grad(ind_saved,:,:,:) = grad_aux
          
          ! Common mf.
#ifdef usempi
          raux = mf_0
          call MPI_ALLREDUCE( raux, mf_0, 1, mpi_fwi, MPI_SUM, MPI_COMM_WORLD, ierr )
#endif
          
          ! Normalization mf at first iter per freq.
          if ( iter_freq == 1 ) mf_i0_f0 = abs(mf_0)
          
          mf_0 = mf_0/mf_i0_f0
          mf_sg = mf_sg/mf_i0_f0
          
          if ( rank == 0 ) write(*,*) ' mf_0   = ', mf_0,   ' mf_i0_f0   = ', mf_i0_f0
          
          !***************************************************!
          ! Mf and energy per each shot.
#ifdef usempi
            v_sg_aux = mf_sg
            call MPI_ALLREDUCE( v_sg_aux, mf_sg, n_shot, mpi_fwi, MPI_SUM, MPI_COMM_WORLD, ierr )
            
            v_sg_aux = mean_grad0_it0_sg
            call MPI_ALLREDUCE( v_sg_aux, mean_grad0_it0_sg, n_shot, mpi_fwi, MPI_SUM, MPI_COMM_WORLD, ierr )
            
            v_sg_aux = es_sg
            call MPI_ALLREDUCE( v_sg_aux, es_sg, n_shot, mpi_fwi, MPI_SUM, MPI_COMM_WORLD, ierr )
            
            v_sg_aux = mean_grad0_sg
            call MPI_ALLREDUCE( v_sg_aux, mean_grad0_sg, n_shot, mpi_fwi, MPI_SUM, MPI_COMM_WORLD, ierr )
            
            v_sg_aux = sig_grad0_sg
            call MPI_ALLREDUCE( v_sg_aux, sig_grad0_sg, n_shot, mpi_fwi, MPI_SUM, MPI_COMM_WORLD, ierr )
#endif
          
          !***************************************************!
          ! Write misfit and energy of the source per each shot.
          if ( rank == 0 ) then
            
            path_file = trim(adjustl( folder_output )) // '/resume/mf1'
            write(str_tmp,*) strat_inv
            path_file = trim( path_file ) // '_strat_'// trim(adjustl(str_tmp))
            write(str_tmp,fmt='(f5.2)') freq_inv
            path_file = trim( path_file ) // '_freq_'// trim(adjustl(str_tmp))
            write(str_tmp,*) iter_freq
            path_file = trim( path_file ) // '_iter_freq_'// trim(adjustl(str_tmp))
            path_file = trim( path_file ) // '.txt'
            
            open(unit=102,file=path_file,status='unknown',action='write')
            do k=1,n_shot
              write(102,'(6(e12.5,1x))') real(mf_sg(k)), real(coef_sea(k)), real(es_sg(k)), real(mean_grad0_sg(k)), real(sig_grad0_sg(k)), real(mean_grad0_it0_sg(k))
            end do
            close(102)
            
          end if
          
          !***************************************************!
          ! Write results at the current iteration.
          if ( rank == rank_wrt ) then
            ! Info.
            write(*,*) 'info ---> shgat(isg_wrt)%id = ', shgat(isg_wrt)%id
            write(*,*) 'all: " %pos_sou(1)%x = ', dx*real(shgat(isg_wrt)%mg%ix-1,sp) + shgat(isg_wrt)%pos_sou(1)%x
            write(*,*) 'all: " %pos_sou(1)%y = ', dx*real(shgat(isg_wrt)%mg%iy-1,sp) + shgat(isg_wrt)%pos_sou(1)%y
            write(*,*) 'all: " %pos_sou(1)%z = ', dx*real(shgat(isg_wrt)%mg%iz-1,sp) + shgat(isg_wrt)%pos_sou(1)%z
            write(*,*) 'all: " %pos_rec(1)%x = ', dx*real(shgat(isg_wrt)%mg%ix-1,sp) + shgat(isg_wrt)%pos_rec(1)%x
            write(*,*) 'all: " %pos_rec(1)%y = ', dx*real(shgat(isg_wrt)%mg%iy-1,sp) + shgat(isg_wrt)%pos_rec(1)%y
            write(*,*) 'all: " %pos_rec(1)%z = ', dx*real(shgat(isg_wrt)%mg%iz-1,sp) + shgat(isg_wrt)%pos_rec(1)%z
            write(*,*) 'all: " %pos_rec(shgat%n_rec)%x = ', dx*real(shgat(isg_wrt)%mg%ix-1,sp) + shgat(isg_wrt)%pos_rec(shgat(isg_wrt)%n_rec)%x
            write(*,*) 'all: " %pos_rec(shgat%n_rec)%y = ', dx*real(shgat(isg_wrt)%mg%iy-1,sp) + shgat(isg_wrt)%pos_rec(shgat(isg_wrt)%n_rec)%y
            write(*,*) 'all: " %pos_rec(shgat%n_rec)%z = ', dx*real(shgat(isg_wrt)%mg%iz-1,sp) + shgat(isg_wrt)%pos_rec(shgat(isg_wrt)%n_rec)%z
            ! Grad.
!             if ( write_grad ) call copy_mat_4D( strat_inv, freq_inv, iter_freq, '/grad/grad', ind_saved, grad )
            if ( write_grad ) call copy_mat_4D_unformatted( strat_inv, freq_inv, iter_freq, '/grad/grad', ind_saved, grad, dx )
            ! Illumination.
!             if ( write_illu ) call copy_mat_3D( strat_inv, freq_inv, iter_freq, '/illu/illu', illu )
!             if ( write_illu ) call copy_mat_3D( strat_inv, freq_inv, iter_freq, '/illu/prec_mod', prec_mod )
            if ( write_illu ) call copy_mat_3D_unformatted( strat_inv, freq_inv, iter_freq, '/illu/illu', illu, dx )
            if ( write_illu ) call copy_mat_3D_unformatted( strat_inv, freq_inv, iter_freq, '/illu/prec_mod', prec_mod, dx )
!             ! Trace.
!             if ( write_trace ) then
!               call copy_mat_2D_unformatted( strat_inv, freq_inv, iter_freq, '/trace/tr_s', shgat(isg_wrt)%tr_syn, shgat(isg_wrt)%t1_syn/real(shgat(isg_wrt)%nt_tr_syn-1,sp) )
!               call copy_mat_2D_unformatted( strat_inv, freq_inv, iter_freq, '/trace/tr', shgat(isg_wrt)%tr_obs, shgat(isg_wrt)%t1_obs/real(shgat(isg_wrt)%nt_tr_obs-1,sp) )
!               call copy_mat_2D_unformatted( strat_inv, freq_inv, iter_freq, '/adj_sou/prec_data', shgat(isg_wrt)%pr_dat, shgat(isg_wrt)%t1_obs/real(shgat(isg_wrt)%nt_tr_obs-1,sp) )
!             end if
            ! Source.
            if ( write_sou ) then
              call copy_mat_2D_unformatted( strat_inv, freq_inv, iter_freq, '/source/sou', shgat(isg_wrt)%sou, shgat(isg_wrt)%t1_sou/real(shgat(isg_wrt)%nt_sou-1,sp) )
            end if
!             ! Adj Source.
!             if ( write_adj_sou ) then
!               call copy_mat_2D_unformatted( strat_inv, freq_inv, iter_freq, '/adj_sou/adj_sou', shgat(isg_wrt)%adj_sou, shgat(isg_wrt)%t1_adj/real(shgat(isg_wrt)%nt_adj-1,sp) )
!             end if
          end if
          
          ! 
          do iwtr=1,n_id_wrt_tr
          do isg=1,n_sg(rank+1)
            if ( shgat(isg)%id == id_wrt_tr(iwtr) ) then
              if ( write_trace ) then
                call copy_mat_2D_B_unformatted( strat_inv, freq_inv, iter_freq, '/trace/tr_s', shgat(isg)%id, shgat(isg)%tr_syn, shgat(isg)%t1_obs/real(shgat(isg)%nt_tr_obs-1,sp) )
                call copy_mat_2D_B_unformatted( strat_inv, freq_inv, iter_freq, '/trace/tr', shgat(isg)%id, shgat(isg)%tr_obs, shgat(isg)%t1_obs/real(shgat(isg)%nt_tr_obs-1,sp) )
                call copy_mat_2D_B_unformatted( strat_inv, freq_inv, iter_freq, '/adj_sou/prec_data', shgat(isg)%id, shgat(isg)%pr_dat, shgat(isg)%t1_obs/real(shgat(isg)%nt_tr_obs-1,sp) )
              end if
              if ( write_adj_sou ) then
                call copy_mat_2D_B_unformatted( strat_inv, freq_inv, iter_freq, '/adj_sou/adj_sou', shgat(isg)%id, shgat(isg)%adj_sou, shgat(isg)%t1_adj/real(shgat(isg)%nt_adj-1,sp) )
              end if
            end if
          end do
          end do
                 
          !***************************************************!
          ! Control gradient: Calculation of some statistical values.
          if ( iter_freq == 1 ) then
            
            call control_grad( m_0, grad_aux, mean, sigma, z_mean )
            
            norm_grad0 = 1./mean
            
            if ( rank == 0 ) then
              write(*,'(a)') '   mean                    sigma/mean              z_mean'
              write(*,*) mean, sigma, 1000._sp*z_mean
            end if
            
          end if
          
          ! Avoiding overloading: Grad normalization.
          call normalize_grad( ind_saved, norm_grad0, grad )
          
          !**********************************************************!
          ! Calculate gradient norm.
          norm_grad = norm_gradient( ind_saved, grad )
          
          ! Check.
! #ifndef valgrind
!           if ( isnan(real(norm_grad)) ) stop '***** ERROR - fwi: isnan(real(norm_grad)) = T *****'
! #else
          if ( isnan(real(norm_grad)) ) then
            if ( rank == 0 ) write(*,*) '***** WARNING - fwi: isnan(real(norm_grad)) = T *****'
            grad = 1.
            norm_grad = 1.
            mean_grad0_it0_sg = 1.
          end if
! #endif
          
          !**********************************************************!
          ! Wolfe conditions.
          if ( ind_saved > 1 ) then
            wol_cur_str = wolfe_curvature_strong( ind_saved-1, grad, search )
            wol_cur     = wolfe_curvature( ind_saved-1, grad, search )
          end if
          
          !**********************************************************!
          ! Inversion QC.
          ! First iter is true then:
          if ( iter_freq > 1 ) then
            
            if ( meth_min_QC == 0 ) then
              inv_QC = .true.
            else if ( meth_min_QC == 1 ) then
              inv_QC = wol_arm
            else if ( meth_min_QC == 2 ) then
              inv_QC = wol_cur
            else if ( meth_min_QC == 3 ) then
              inv_QC = wol_cur_str
            end if
            
          end if
          
          !**********************************************************!
          !**********************************************************!
          !**********************************************************!
          ! Define QC=true
          ! Step 1: grad, QC?, var=m_0, min, update m_0=min
          ! QC true
          ! Step X: grad, QC?, var=m_0, min, update m_0=min
          ! QC flase
          ! Step X: grad, QC?, m_0=var old
          !**********************************************************!          
          ! If inv fails we recovered the last saved model.
          if ( .not. inv_QC ) then
            
            ! Recover the last.
            m_0 = var(max(1,ind_saved-1),:,:,:)
            
!             ! Recover the model in the canonical variables.
!             call change_var_inv_2( m_0, cmpr, buoy  )
            
            ! Impose constraints.
            if ( vel_cut_off ) call param_control_3( m_0 )
            ! Gardner.
            if ( gardner ) call gardner_density_2( m_0, buoy )
            ! From parametrization to cmpr and dnst.
            call change_var_inv_3( m_0, buoy, cmpr  )
            
          !----------------------------------------------------------!
          ! Otherwise we continue.
          else
            
            ! If not wol_cur_str with l-BFGS, we start by applying the SD.
            if ( ( meth_inv == 2 ) .and. ( .not. wol_cur_str ) ) then
              grad(1,:,:,:) = grad(ind_saved,:,:,:)
              ind_saved = 1
            end if
            
            ! Save previous model.
            var(ind_saved,:,:,:) = m_0
            
            ! Search direction.
            call search_min( meth_inv, ind_saved, var, grad, search )
            
            ! alp_cota
            alp_cota = cal_alp_cota( m_0, search, change_cota )
            
            ! prod_grad_search
            prod_grad_search = sum(sum(sum( grad(ind_saved,:,:,:)*search ,3),2),1)
            
            ! Synthetic with all OK.
            ! First step.
            if ( ( meth_inv == 0 .or. meth_inv == 1 ) .and. iter_freq == 1 ) then
              alp_1 = perc_max_alp_cota*alp_cota
            else if ( meth_inv == 2 .and. ind_saved == 1 ) then
              alp_1 = perc_max_alp_cota*alp_cota
            ! Following steps.
            else
              ! Initial step length.
              alp_scaled = alp_min*( prod_grad_search0/prod_grad_search )
              ! If prod_grad_search0/prod_grad_search is negative we reject alp_scaled since is wrong.
              if ( alp_scaled > 0._sp ) then
                alp_1 = max( perc_min_alp_cota*alp_cota, min( alp_scaled, perc_max_alp_cota*alp_cota ) )
              else
                alp_1 = perc_max_alp_cota*alp_cota
              end if
            end if
            
            !**********************************************************!
            ! Write.
            if ( rank == 0 .and. iter_freq > 1 ) then
              write(*,*) ' prod_grad_search0/prod_grad_search = ', real(prod_grad_search0/prod_grad_search)
              write(*,'(a,e12.5)') ' alp_1/alp_cota      = ', real(alp_1/alp_cota)
              write(*,'(a,e12.5)') ' alp_scaled/alp_cota = ', real(alp_scaled/alp_cota)
              write(*,'(a,e12.5)') ' alp_min/alp_cota    = ', real(alp_min/alp_cota)
            end if
            
            ! Write.
            if ( rank == rank_wrt .and. write_sear ) then
              mat_aux = ( alp_cota*search )/m_0
!               call copy_mat_3D( strat_inv, freq_inv, iter_freq, '/grad/search', mat_aux )
              call copy_mat_3D_unformatted( strat_inv, freq_inv, iter_freq, '/grad/search', mat_aux, dx )
            end if
            
            !**********************************************************!
            !*** Minimization. ****************************************!
            !**********************************************************!
            ! Old mf.
            mf_0_old = mf_0
            
            ! Minimize.
            shgat(:)%alg = 2
            call minimize( n_sg, shgat, freq_0_inv, freq_inv, mf_0, mf_i0_f0, alp_1, &
                           alp_cota, mf_1, alp_2, mf_2, alp_fit, mf_fit, alp_min, mf_min, &
                           ind_saved, grad, grad_sg, prod_grad_search, wol_arm, mod_change, &
                           pol_min_QC, succes_red, mat_aux )
            
            !***************************************************!
            ! Write the evolution of the FWI.
            if ( rank == 0 ) then
              write(*,*) 'pol_min_QC = ', pol_min_QC
              write(*,*) 'succes_red = ', succes_red
              write(*,'(a,f7.3,a)') ' red. mf = ', 100._sp*(1._sp-mf_min/mf_0_old), '%'
              write(*,'(a,f7.3,a)') ' change  = ', 100._sp*mod_change, '%'
              write(*,*) '/////////////////////////////////////////////////'
              write(*,*)   'norm_grad      = ', norm_grad
              write(*,*)   'wol_arm        = ', wol_arm
              if ( ind_saved > 1 ) then
                write(*,*) 'wolfe_curv     = ', wol_cur
                write(*,*) 'wolfe_curv_str = ', wol_cur_str
              end if
              write(*,*) '/////////////////////////////////////////////////'
              
              write(101,*) 'pol_min_QC = ', pol_min_QC
              write(101,*) 'succes_red = ', succes_red
              write(101,'(a,f9.5,a)') ' red. mf = ', 100._sp*(1._sp-mf_min/mf_0_old), '%'
              write(101,'(a,f9.5,a)') ' change  = ', 100._sp*mod_change, '%'
              if ( iter_freq == 1 ) write(101,*) 'mf_i0_f0   = ', mf_i0_f0
              write(101,'(1(a,f9.5))') ' mf_0   = ', mf_0
              write(101,'(3(a,f9.5))') ' mf_1   = ', mf_1,   ' a_1/a_ct   = ', alp_1/alp_cota,   ' a_1   = ', alp_1
              write(101,'(3(a,f9.5))') ' mf_2   = ', mf_2,   ' a_2/a_ct   = ', alp_2/alp_cota,   ' a_2   = ', alp_2
              write(101,'(3(a,f9.5))') ' mf_fit = ', mf_fit, ' a_fit/a_ct = ', alp_fit/alp_cota, ' a_fit = ', alp_fit
              write(101,'(3(a,f9.5))') ' mf_min = ', mf_min, ' a_fit/a_ct = ', alp_min/alp_cota, ' a_fit = ', alp_min
              if ( iter_freq > 1 ) write(101,*) 'prod_grad_search0/prod_grad_search = ', prod_grad_search0/prod_grad_search
              write(101,*) '//////////////////////////////////////////////////////////////////////////////'
              write(101,*)            'norm_grad      = ', norm_grad
              write(101,*)            'wolfe_arm      = ', wol_arm
              if ( ind_saved > 1 ) then
                write(101,*) ' wolfe_curv     = ', wol_cur
                write(101,*) ' wolfe_curv_str = ', wol_cur_str
              end if
            end if
            
            !***************************************************!
            ! Update.
            prod_grad_search0 = prod_grad_search
            
            !***************************************************!
            ! Sum of the mf coming from all the cores.
#ifdef usempi
            v_sg_aux = mf_sg
            call MPI_ALLREDUCE( v_sg_aux, mf_sg, n_shot, mpi_fwi, MPI_SUM, MPI_COMM_WORLD, ierr )
#endif
            
            mf_sg_a = mf_sg_b
            mf_sg_b = mf_sg
            
            !***************************************************!
            ! Write mf per shot.
            if ( rank == 0 ) then
              
              path_file = trim(adjustl( folder_output )) // '/resume/mf2'
              write(str_tmp,*) strat_inv
              path_file = trim( path_file ) // '_strat_'// trim(adjustl(str_tmp))
              write(str_tmp,fmt='(f5.2)') freq_inv
              path_file = trim( path_file ) // '_freq_'// trim(adjustl(str_tmp))
              write(str_tmp,*) iter_freq
              path_file = trim( path_file ) // '_iter_freq_'// trim(adjustl(str_tmp))
              path_file = trim( path_file ) // '.txt'
              
              open(unit=102,file=path_file,status='unknown',action='write')
              do k=1,n_shot
                write(102,'(e12.5)') real(mf_sg(k))
              end do
              close(102)
              
            end if
            
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !             !***************************************************!
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !             ! Saving and reorganising the grad and model.
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !             if ( ind_saved == n_saved ) then
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !               do k=2,n_saved
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !                 grad(k-1,:,:,:) = grad(k,:,:,:)
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !                 var(k-1,:,:,:)  = var(k,:,:,:)
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !               end do
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !             end if
            
            !**********************************************************!
!             ! Recover the model in the canonical variables.
!             call change_var_inv_2( m_0, cmpr, buoy  )
            
            ! Impose constraints.
            if ( vel_cut_off ) call param_control_3( m_0 )
            ! Gardner.
            if ( gardner ) call gardner_density_2( m_0, buoy )
            ! From parametrization to cmpr and dnst.
            call change_var_inv_3( m_0, buoy, cmpr  )
            
            !***************************************************!
            ! Writing the updated common model.
            if ( write_model .and. rank == rank_wrt ) then
              mat_aux = 1000._sp*sqrt(cmpr*buoy)
              if ( any(isnan(mat_aux)) ) write(*,*) ' any(isnan(mat_aux)) '
!               call copy_mat_3D( strat_inv, freq_inv, iter_freq, '/model/vP', mat_aux )
              call copy_mat_3D_unformatted( strat_inv, freq_inv, iter_freq, '/model/vP', mat_aux, dx )
            end if
            
          end if
          
          !**********************************************************!          
          ! Final results at the last frequency iteration.
          if ( rank == rank_wrt ) then
            ! model
            if ( .not. write_model ) then
              mat_aux = 1000._sp*sqrt(cmpr*buoy)
              call copy_mat_3D_unformatted( strat_inv, 0._sp, 0, '/results/vP_last', mat_aux, dx )
            end if
            ! Grad.
            if ( .not. write_sear ) then
              call copy_mat_3D_unformatted( strat_inv, 0._sp, 0, '/results/search_last', search, dx )
            end if
            ! Illumination.
            if ( .not. write_illu ) then
              call copy_mat_3D_unformatted( strat_inv, 0._sp, 0, '/results/illu_last', illu, dx )
              call copy_mat_3D_unformatted( strat_inv, 0._sp, 0, '/results/prec_mod_last', prec_mod, dx )
            end if
            ! Trace.
!             if ( .not. write_trace ) then
!               call copy_mat_2D_unformatted( strat_inv, 0._sp, 0, '/results/tr_s_last', shgat(isg_wrt)%tr_syn, shgat(isg_wrt)%t1_obs/real(shgat(isg_wrt)%nt_tr_obs-1,sp) )
!               call copy_mat_2D_unformatted( strat_inv, 0._sp, 0, '/results/tr_last', shgat(isg_wrt)%tr_obs, shgat(isg_wrt)%t1_obs/real(shgat(isg_wrt)%nt_tr_obs-1,sp) )
!               call copy_mat_2D_unformatted( strat_inv, 0._sp, 0, '/results/prec_data_last', shgat(isg_wrt)%pr_dat, shgat(isg_wrt)%t1_obs/real(shgat(isg_wrt)%nt_tr_obs-1,sp) )
!             end if
          end if
          
          ! 
          do iwtr=1,n_id_wrt_tr
          do isg=1,n_sg(rank+1)
            if ( shgat(isg)%id == id_wrt_tr(iwtr) ) then
              if ( .not. write_trace ) then
                call copy_mat_2D_B_unformatted( strat_inv, 0._sp, 0, '/trace/tr_s', shgat(isg)%id, shgat(isg)%tr_syn, shgat(isg)%t1_obs/real(shgat(isg)%nt_tr_obs-1,sp) )
                call copy_mat_2D_B_unformatted( strat_inv, 0._sp, 0, '/trace/tr', shgat(isg)%id, shgat(isg)%tr_obs, shgat(isg)%t1_obs/real(shgat(isg)%nt_tr_obs-1,sp) )
                call copy_mat_2D_B_unformatted( strat_inv, 0._sp, 0, '/adj_sou/prec_data', shgat(isg)%id, shgat(isg)%pr_dat, shgat(isg)%t1_obs/real(shgat(isg)%nt_tr_obs-1,sp) )
              end if
              if ( .not. write_adj_sou ) then
                call copy_mat_2D_B_unformatted( strat_inv, 0._sp, 0, '/adj_sou/adj_sou', shgat(isg)%id, shgat(isg)%adj_sou, shgat(isg)%t1_adj/real(shgat(isg)%nt_adj-1,sp) )
              end if
            end if
          end do
          end do
          
          !**********************************************************!          
          ! Break loop and pass to the next freq?
          if ( .not. inv_QC ) exit   ! We have recovered the last ok iteration, now we pass to next iter.
          if ( .not. succes_red ) exit   ! Nothing to be recover, we just pass to next iter.
          if ( abs(1._sp-mf_min/mf_0_old) < chng_min_QC ) exit   ! The mf reduction is too small. Abs to also maxmize instead of minimize.
          
          !**********************************************************!          
        end do     ! Loop single freq.
      end do     ! Loop group of freq.
      
      
      
      !**********************************************************!
      !**********************************************************!
      !**********************************************************!
      !*** We have finish the FWI.   ****************************!
      !**********************************************************!
      ! Writting final results.
      if ( rank == rank_wrt ) then
        ! model
        if ( write_model ) then
          mat_aux = 1000.*sqrt(cmpr*buoy)
          call copy_mat_3D_unformatted( strat_inv, 0._sp, 0, '/results/vP_end', mat_aux, dx )
        end if
        ! Grad.
        if ( write_model ) then
          call copy_mat_4D_unformatted( strat_inv, 0._sp, 0, '/results/grad_end', ind_saved, grad, dx )
        end if
        ! Illumination.
        if ( write_model ) then
          call copy_mat_3D_unformatted( strat_inv, 0._sp, 0, '/results/illu_end', illu, dx )
          call copy_mat_3D_unformatted( strat_inv, 0._sp, 0, '/results/prec_m_end', prec_mod, dx )
        end if
!         ! Trace.
!         if ( write_trace ) then
!           call copy_mat_2D_unformatted( strat_inv, 0._sp, 0, '/results/tr_s_end', shgat(isg_wrt)%tr_syn, shgat(isg_wrt)%t1_obs/real(shgat(isg_wrt)%nt_tr_obs-1,sp) )
!           call copy_mat_2D_unformatted( strat_inv, 0._sp, 0, '/results/tr_end', shgat(isg_wrt)%tr_obs, shgat(isg_wrt)%t1_obs/real(shgat(isg_wrt)%nt_tr_obs-1,sp) )
!           call copy_mat_2D_unformatted( strat_inv, 0._sp, 0, '/results/prec_data_end', shgat(isg_wrt)%pr_dat, shgat(isg_wrt)%t1_obs/real(shgat(isg_wrt)%nt_tr_obs-1,sp) )
!         end if
        ! Source.
        if ( write_sou ) then
          call copy_mat_2D_unformatted( strat_inv, 0._sp, 0, '/results/sou_end', shgat(isg_wrt)%sou, shgat(isg_wrt)%t1_sou/real(shgat(isg_wrt)%nt_sou-1,sp) )
        end if
!         ! Adj Source.
!         if ( write_adj_sou ) then
!           call copy_mat_2D_unformatted( strat_inv, 0._sp, 0, '/results/adj_sou_end', shgat(isg_wrt)%adj_sou, shgat(isg_wrt)%t1_obs/real(shgat(isg_wrt)%nt_tr_obs-1,sp) )
!         end if
        
      end if
      
      ! 
      do iwtr=1,n_id_wrt_tr
      do isg=1,n_sg(rank+1)
        if ( shgat(isg)%id == id_wrt_tr(iwtr) ) then
          if ( write_trace ) then
            call copy_mat_2D_B_unformatted( strat_inv, 0._sp, 0, '/results/tr_s_end', shgat(isg)%id, shgat(isg)%tr_syn, shgat(isg)%t1_obs/real(shgat(isg)%nt_tr_obs-1,sp) )
            call copy_mat_2D_B_unformatted( strat_inv, 0._sp, 0, '/results/tr_s_end', shgat(isg)%id, shgat(isg)%tr_obs, shgat(isg)%t1_obs/real(shgat(isg)%nt_tr_obs-1,sp) )
            call copy_mat_2D_B_unformatted( strat_inv, 0._sp, 0, '/results/prec_data_end', shgat(isg)%id, shgat(isg)%pr_dat, shgat(isg)%t1_obs/real(shgat(isg)%nt_tr_obs-1,sp) )
          end if
          if ( write_adj_sou ) then
            call copy_mat_2D_B_unformatted( strat_inv, 0._sp, 0, '/results/adj_sou', shgat(isg)%id, shgat(isg)%adj_sou, shgat(isg)%t1_adj/real(shgat(isg)%nt_adj-1,sp) )
          end if
        end if
      end do
      end do
      
      !**********************************************************!
      ! Deallocate shgat.
      call solv3Dahv_deallo
      
      !**********************************************************!
      ! Deallocate.
      deallocate( grad_aux, stat=stal ); if ( stal/=0 ) stop 'AE fwi14c'
      deallocate( mat_aux, stat=stal ); if ( stal/=0 ) stop 'AE fwi14c'
      deallocate( grad, stat=stal ); if ( stal/=0 ) stop 'AE fwi18'
      deallocate( var, stat=stal ); if ( stal/=0 ) stop 'AE fwi18b'
      deallocate( search, stat=stal ); if ( stal/=0 ) stop 'AE fwi9'
      deallocate( m_0, m_1, m_2, m_fit, m_min, stat=stal ); if ( stal/=0 ) stop 'AE fwi20'
      deallocate( illu, stat=stal ); if ( stal/=0 ) stop 'AE fwi23'
      deallocate( prec_mod, stat=stal ); if ( stal/=0 ) stop 'AE fwi24'
      deallocate( v_sg_aux, stat=stal ); if ( stal/=0 ) stop 'AE fwi25'
      deallocate( mean_grad0_it0_sg, mf_sg, mf_sg_1, mf_sg_2, mf_sg_fit, mf_sg_ck, mf_sg_min, es_sg, mean_grad0_sg, sig_grad0_sg, stat=stal ); if ( stal/=0 ) stop 'AE fwi26'
      deallocate( mf_sg_a, mf_sg_b, stat=stal ); if ( stal/=0 ) stop 'AE fwi26'
      deallocate( coef_sea, stat=stal ); if ( stal/=0 ) stop 'AE fwi26b'
      
      deallocate( grad_sg, stat=stal ); if ( stal/=0 ) stop 'AE fwi30'
      
      ! Deallocate H things.
      if ( use_iH ) then
        if ( full_H ) then
          deallocate( inv_H, stat=stal ); if ( stal/=0 ) stop 'dAE fwi27'
        else
          deallocate( inv_H_diag, stat=stal ); if ( stal/=0 ) stop 'dAE fwi28'
        end if
      end if
      
      ! Deallocate Wd.
      if ( dyw_Wd ) then
        deallocate( Wd, stat=stal ); if ( stal/=0 ) stop 'dAE fwi29'
        deallocate( iW, stat=stal ); if ( stal/=0 ) stop 'dAE fwi30'
        deallocate( jW, stat=stal ); if ( stal/=0 ) stop 'dAE fwi31'
      end if
      
      !**********************************************************!
      ! Close unit of master core.
      if ( rank == 0 ) close(101)
      
      !**********************************************************!
      ! Ending.
      write(*,*) 'End of the inversion ---> rank = ',rank
      
    !*********************************************************************/
    end subroutine fwi
  
  !***********************************************************************/
  end module m_fwi




