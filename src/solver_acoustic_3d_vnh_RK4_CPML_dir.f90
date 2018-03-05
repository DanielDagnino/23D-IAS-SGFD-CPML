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
  !** Short resume of this code/file:
  !** 2D/3D acoustic adjoint wave propagator
  !** 
  !** 
  !** 
  !*********************************************************************/
  module m_solv3Dahv_dir
    implicit none
    
    !*********************************************************************/
    contains
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine solv3Dahv_dir( time_all, shgat, buoy_in, cmpr_in, freq_inv, i_save, &
                              ind_save, vf_dir, illu_sum_mg )
      use m_mpi
      use m_dyw, only: dyw_allo, dyw_allo_b, dyw_deallo, dyw_deallo_b
      use m_sg_type
      use m_sg_data, only: folder_tmp
      use m_sg_tr
      use m_sg_sou
      use m_solv3Dahvc_solver
      use m_solv3Dahvc_PML_building
      use m_solv3Dahvc_solver_sources
      use m_solv3Dahvc_RK_loop
      use m_solv3Dahvc_interpol_sg
      use m_geo
      use m_dag_time
      use m_solv3Dahvc_write_field
      use m_solv3Dahvc_RK_save_boundary
      use m_solv3Dahvc_RK_loop_back_aux
      use m_solv3Dahvc_deriv
      use m_solv3Dahvc_deriv_o
      implicit none
      ! The variables which are passed to the function.
      real(8), intent(out)                :: time_all
      type(shot_gather), intent(inout)    :: shgat
      real(rl), allocatable, intent(in)   :: buoy_in(:,:,:), cmpr_in(:,:,:)
      real(rl), intent(in)                :: freq_inv
      type(field), intent(out)            :: vf_dir   !last_vf
      integer, intent(inout)              :: i_save
      integer, allocatable, intent(inout) :: ind_save(:)
      real(rl), allocatable, intent(inout):: illu_sum_mg(:,:,:)
      ! The variables which are generated inside the function.
!       type(field), target :: vf_dir
      integer :: ix, iy, iz
      integer :: ix_QC, iy_QC, iz_QC
      integer :: cont, rest
      integer :: ir
      integer :: it
      integer :: stal
      integer :: cont_save, ind_bd
      real(rl) :: vp_max
      real(rl) :: c_simpson
#ifdef test_plot
      integer  :: n_plot_plane
#endif

!       !*********************************************************************/
!       ! 
! #ifdef testing
!       write(*,*) '*****     Testing     *****'
! #else
!       write(*,*) '*****     No Testing     *****'
! #endif
      
      !*********************************************************************/
      ! Time of the CPU.
      time_all = rtc()      
      
!       !*********************************************************************/
!       ! Checking things.
!       if ( shgat%ts0 /= shgat%t0_obs ) then
!         write(*,*) shgat%ts0, shgat%t0_obs
!         stop '***** ERROR - source time does not match time shot: shgat%ts0 /= shgat%t0_obs *****'
! !       else if ( shgat%ts1 < shgat%t1_obs ) then
! !         write(*,*) shgat%ts1, shgat%t1_obs
! !         stop '***** ERROR - source time does not match time shot: shgat%ts1 < shgat%t1_obs *****'
!       end if
      
      !*********************************************************************/
      !----------------------------------------------------------!
      ! Definition of numercial constants.
!       write(*,*) 'init_solver_coef_space'
      call init_solver_coef_space
      
      ! V_p max.
      vp_max = sqrt(maxval(cmpr_in*buoy_in))
      ! Check.
      if ( vp_max > 10._rl ) then
        write(*,*) '***** WARNING: solv3Dahv_test: vp_max >  10._rl *****'
        write(*,*) 'vp_max          = ', vp_max, ' km/s'
        write(*,*) 'maxval(cmpr_in) = ', maxval(cmpr_in)
        write(*,*) 'maxval(buoy_in) = ', maxval(buoy_in)
      else if ( vp_max < 0.5_rl ) then
        write(*,*) '***** WARNING: solv3Dahv_test: vp_max < 0.5_rl *****'
        write(*,*) 'vp_max          = ', vp_max, ' km/s'
        write(*,*) 'maxval(cmpr_in) = ', maxval(cmpr_in)
        write(*,*) 'maxval(buoy_in) = ', maxval(buoy_in)
      end if
      
!       write(*,*) 'init_solver_coef_time'
      call init_solver_coef_time( shgat%t1_inv, real(vp_max,rl) )
!       if ( rank == 0 ) write(*,'(a,1x,i4)') ' nt_solver (dir) = ', nt_solver
      
!       write(*,*) 'adapt_one_sg_sou'
      call adapt_one_sg_sou( shgat, nt_solver, shgat%t1_inv )
      
!       write(*,*) 'allo_one_sg_tr'
      call allo_one_sg_tr_s( shgat, nt_store( nt_solver, 1 ) )
      shgat%t0_syn = 0._sp
      shgat%t1_syn = shgat%t1_inv
      
      !----------------------------------------------------------!
      ! Allocation.
!       write(*,*) 'init_solver_allo'
      call allo_field( vf_dir, shgat%mg )
      
      if ( dyw_allo ) then
        call init_solver_allo( shgat%mg )
        call init_PML_allo( shgat%mg )
        dyw_allo = .false. 
      end if
      
      if ( dyw_allo_b ) then
        call init_b_allo( shgat%mg )
        dyw_allo_b = .false. 
      end if
      
!       write(*,*) 'source allocate'
      allocate( int_sou(shgat%n_sou,nt_solver), stat=stal ); if ( stal/=0 ) stop 'AE m_solv3Dahvc_test 2'
      int_sou = 0._rl
      
      !----------------------------------------------------------!
      ! Copy model for transformation.
!$omp parallel do default(none) shared(buoy,buoy_in,cmpr,cmpr_in,tau,tau_in) firstprivate(shgat%mg%iLx,shgat%mg%iLy,shgat%mg%iLz) schedule(dynamic,1000) collapse(3)
      do iz=1,shgat%mg%iLz
      iz_QC = min( nz, (shgat%mg%iz-1)+iz )
      do iy=1,shgat%mg%iLy
      iy_QC = min( ny, (shgat%mg%iy-1)+iy )
      do ix=1,shgat%mg%iLx
      ix_QC = min( nx, (shgat%mg%ix-1)+ix )
        buoy(ix,iy,iz) = buoy_in(ix_QC,iy_QC,iz_QC)
        cmpr(ix,iy,iz) = cmpr_in(ix_QC,iy_QC,iz_QC)
      end do
      end do
      end do
!$omp end parallel do
      
      !----------------------------------------------------------!
      ! Staggered grid model requiered.
!       write(*,*) 'init_staggered_model'
      call init_staggered_model( shgat%mg )
      
      !----------------------------------------------------------!
      ! Definition of PML constants and model parameters.
!       write(*,*) 'PML_definitions_model'
      call PML_definitions_model( shgat%mg )
!       write(*,*) 'PML_definitions_param'
      call PML_definitions_param( freq_inv, shgat%mg )
      
      !----------------------------------------------------------!
      ! Source.
!       write(*,*) 'init_source'
      call init_source( shgat, int_sou )
      
      !*********************************************************************/
      !*********************************************************************/
      !* Solver - Elastic wave equation. ***********************************/
      !*********************************************************************/
      !*********************************************************************/
!       if ( nPML_B == 0 ) then
!         apply_FS = .true.
!       else
!         apply_FS = .false.
!       end if
      
      if ( nPML_B == 0 ) then
        if ( .not. apply_FS ) stop ' ********** ERROR: solv3Dahv_dir: nPML_B == 0 .and. .not. apply_FS ********** '
      else
        if ( apply_FS ) stop ' ********** ERROR: solv3Dahv_dir: nPML_B == 0 .and. apply_FS ********** '
      end if
      
!       write(*,*) 'loop'
      
      !*********************************************************************/
      ! 
!       if ( nt_solver<save_every ) then
!         stop '***** ERROR solv3Dahv_dir - nt_solver<save_every *****'
!       end if
      
      ! Index value for saving boundary.
      if ( nt_solver>save_every ) then
        
        ! Save field every ??? iteration.
        rest = mod(nt_solver,save_every)
        if ( rest==0 ) then
          i_save = nt_solver/save_every
        else
          i_save = 1 + floor(dble(nt_solver)/dble(save_every))
        end if
        
        ! Index of the saved fields.
        allocate( ind_save(i_save), stat=stal ); if ( stal/=0 ) stop 'AE solv3Dahv_dir'
        
        if ( rest == 0 ) then
          do cont_save=1,i_save
            ind_save(cont_save) = save_every*cont_save -1
          end do
        else
          do cont_save=1,i_save
            ind_save(cont_save) = save_every*(cont_save-1) + (rest-1)
          end do
        end if
        
      else
        
        i_save = 1
        allocate( ind_save(1), stat=stal ); if ( stal/=0 ) stop 'AE solv3Dahv_dir'
        ind_save(1) = nt_solver - 1
      
      end if
      
      !*********************************************************************/
      ! Wave equation integration.
      ! NOTE: Only nt_solver-1 steps are computed. The first one is t=0 and it is only saved.
      cont      = 0   ! Trace saving.
      ind_bd    = 0
      cont_save = 0   ! Boundary saving.
      
      c_simpson   = 1._rl
      illu_sum_mg = 0._rl
      
      !----------------------------------------------------------!
      !----------------------------------------------------------!
      do it=1,nt_solver
        
#ifdef test_plot
        ! Iteration number.
        if ( mod(it-1,100)==0 ) write(*,*) 'it , nt_solver = ',  it, nt_solver
        ! For plotting.
        n_plot_plane = 1 + floor(0.5*dble(shgat%mg%iLy))
        include 'plot_include/prop_shot_plot_xz.f90'
        n_plot_plane = 1 + floor(0.5*dble(shgat%mg%iLx))
        include 'plot_include/prop_shot_plot_yz.f90'
        n_plot_plane = 1 + floor(0.5*dble(shgat%mg%iLz))
        include 'plot_include/prop_shot_plot_xy.f90'
#endif
        
        !----------------------------------------------------------!
        ! 
        if ( it==1 ) then
          c_simpson = 7._rl
        else if ( mod(it-1,4)==0 ) then
          c_simpson = 14._rl
        else if ( mod(it-1,2)==0 ) then
          c_simpson = 12._rl
        else
          c_simpson = 32._rl
        end if
        
        !----------------------------------------------------------!
        ! 
        do iz=1,shgat%mg%iLz
        do iy=1,shgat%mg%iLy
        do ix=1,shgat%mg%iLx
          illu_sum_mg(ix,iy,iz) = illu_sum_mg(ix,iy,iz) + dt_solver*c_simpson*(vf_dir%u(ix,iy,iz))**2
        end do
        end do
        end do
        
        !----------------------------------------------------------!
        ! Seismogram.
!         if ( mod(it-1,store)==0 ) then
          
          cont = cont + 1
          if ( cont<=shgat%nt_tr_syn ) then ! En caso de no ser el último.
            
            do ir=1,shgat%n_rec
              shgat%tr_syn(ir,cont) = interpol_1o_3d_ug( shgat%pos_rec(ir), vf_dir%u, dx )
            end do
            
          else
            
            stop '***** ERROR: solv3Dahv_dir: cont > shgat%nt_tr_syn *****'
            
          end if
          
!         end if
        
        !----------------------------------------------------------!
        ! Runge-Kutta.        
        ! RK step.
        if ( it/=nt_solver ) then
          !----------------------------------------------------------!
#ifndef test_no_RK_loop
          call RK_loop( it, shgat, vf_dir, ktmp, k1, k2, k3, k4 )
#endif
          
          !----------------------------------------------------------!
          ! Save field in RAM.
          ! Iteration number of the boundary saving.
          ind_bd = ind_bd + 1
!         write(*,*) 'it ind_bd cont_save = ', it, ind_bd, cont_save
          
          dt_solver = -dt_solver
          dt_o_6 = dt_solver/6._rl
          dt_o_3 = dt_solver/3._rl
          dt_o_2 = dt_solver/2._rl
          
            ! RK back propagation to calculate the RK steps backward.
#ifndef test_no_RK_loop
            call RK_loop_back_aux( shgat, vf_dir, ktmp, k1, k2, k3, k4 )
#endif
!             write(*,*) 'it, ind_bd, cont_save, save_every = ', it, ind_bd, cont_save, save_every
!             write(*,*) 'allocated(bd1(ind_bd)) = ', allocated(bd1(ind_bd)%L%u)
            call save_boundary( bd1(ind_bd), k1, shgat%mg )
            call save_boundary( bd2(ind_bd), k2, shgat%mg )
            call save_boundary( bd3(ind_bd), k3, shgat%mg )
            call save_boundary( bd4(ind_bd), k4, shgat%mg )
            
          dt_solver = -dt_solver
          dt_o_6 = dt_solver/6._rl
          dt_o_3 = dt_solver/3._rl
          dt_o_2 = dt_solver/2._rl
          
          !----------------------------------------------------------!
          ! Save the boundary in ROM.
!           if ( 0 == cont_save ) then
            
            if ( it==ind_save(cont_save+1) .and. it/=nt_solver-1 ) then
              cont_save = cont_save + 1
  !             write(*,*) 'write - cont_save = ', cont_save
              call write_field( bd1, 1, shgat%id, cont_save, folder_tmp, save_every, shgat%mg )
              call write_field( bd2, 2, shgat%id, cont_save, folder_tmp, save_every, shgat%mg )
              call write_field( bd3, 3, shgat%id, cont_save, folder_tmp, save_every, shgat%mg )
              call write_field( bd4, 4, shgat%id, cont_save, folder_tmp, save_every, shgat%mg )
              ind_bd = 0
            end if
            
!           else if ( i_save > cont_save+1 ) then
!             
!             if ( it==ind_save(cont_save+1) ) then
!               cont_save = cont_save + 1
!   !             write(*,*) 'write - cont_save = ', cont_save
!               call write_field( bd1, 1, shgat%id, cont_save, folder_tmp, save_every )
!               call write_field( bd2, 2, shgat%id, cont_save, folder_tmp, save_every )
!               call write_field( bd3, 3, shgat%id, cont_save, folder_tmp, save_every )
!               call write_field( bd4, 4, shgat%id, cont_save, folder_tmp, save_every )
!               ind_bd = 0
!             end if
!             
!           end if
        
        end if
        
        !----------------------------------------------------------!
        if ( it>nt_solver ) stop 'ERROR solv3Dahv_dir it>nt_solver'
        
      end do
      
#ifdef test_no_RK_loop
      shgat%tr_syn = 1._rl
#endif

      !*********************************************************************/
      ! Dealloaction.
!       call deallo_field( vf_dir )

      if ( dyw_deallo ) then
        call init_solver_deallo
        call init_PML_deallo
      end if
      
      if ( dyw_deallo_b ) then
        call init_b_deallo
      end if
      
      deallocate( int_sou, stat=stal ); if ( stal/=0 ) stop 'dAE m_solv3Dahvc_test 2'
      
      !*********************************************************************/      
#ifdef testing
      write(*,*) 'cont        = ', cont
      write(*,*) 'shgat%nt_tr_syn  = ', shgat%nt_tr_syn
      write(*,*) 'nt_solver   = ', nt_solver
#endif
      time_all = rtc()-time_all

    end subroutine solv3Dahv_dir
  
  !***********************************************************************/
  end module m_solv3Dahv_dir
  
  


