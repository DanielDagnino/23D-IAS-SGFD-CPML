  !*********************************************************************/
  !** This code has been done in the Barcelona Center for Subsurface 
  !** Imaging (BCSI).
  !** Authors: Daniel Dagnino.
  !*********************************************************************/
  !** Short resume of this code/file:
  !** 
  !** 
  !** 
  !*********************************************************************/
  module m_correct_search
    implicit none
    
    !*********************************************************************/
    contains
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    logical function correct_search( n_sg, shgat, mf_sg, mf_sg_new, search, grad_sg, coef_sea, mat_aux, perc_wrong, &
                                     alp_1_global, alp_2_global, alp_min_global )
      use m_mpi
      use m_data_kind
      use m_sg_type
      use m_sg_type
      use m_geo
      use m_mygeo
      use m_smooth_op
      use m_inv_data
      use m_type_inversion, only: lay_frz
      use m_work, only: mf_sg_1, mf_sg_2, n_shot
      use m_poly_fit_min, only: poly_fit_min
#ifdef usempi
      use mpi
#endif
      implicit none
      ! The variables which are passed to the function.
      integer, allocatable, intent(in) :: n_sg(:)
      type(shot_gather), allocatable, intent(in) :: shgat(:)
      real(rl), allocatable, intent(in) :: mf_sg(:), mf_sg_new(:)
      real(rl), intent(in) :: alp_1_global, alp_2_global, alp_min_global
      real(rl), allocatable, intent(inout) :: coef_sea(:)
      real(rl), allocatable, intent(inout) :: search(:,:,:), mat_aux(:,:,:)
      real(rl), allocatable, intent(inout) :: grad_sg(:,:,:,:)
      real(rl), intent(out) :: perc_wrong
      ! The variables which are generated inside the function.    
      integer :: isg, ix, iy, iz
      integer :: ix_shift, iy_shift, iz_shift
      integer :: ny0
      integer :: stal
      real(rl) :: perc_ok, raux, r_max_mpi
      real(4), allocatable :: sr_aux(:), sr_aux2(:)
      logical :: laux, minimum
      
      real(rl) :: coef, alp_fit, mf_fit_inter
      
!dir$ assume_aligned mf_sg(1):64,mf_sg_new(1):64
!dir$ assume_aligned search(1,1,1):64,mat_aux(1,1,1):64,grad_sg(1,1,1,1):64
      
      !**********************************************************!
      ! 
      correct_search = .false.
      perc_ok = 0._rl
      perc_wrong = 0._rl
      
      allocate( sr_aux(nx), sr_aux2(nx), stat=stal ); if ( stal/=0 ) stop 'AE correct_search 1'
      sr_aux  = 0._rl
      sr_aux2 = 0._rl
      
      !**********************************************************!
      !**********************************************************!
      !**********************************************************!
      ! Shot.
      do isg=1,n_sg(rank+1)
        
        ! 
        if ( mf_sg_new(shgat(isg)%ord) - mf_sg(shgat(isg)%ord) <  0._rl ) then
!           coef_sea(shgat(isg)%ord) = 1._rl
          perc_ok = perc_ok + (mf_sg_new(shgat(isg)%ord)-mf_sg(shgat(isg)%ord))**2
        else
          coef_sea(shgat(isg)%ord) = 0._rl
          perc_wrong = perc_wrong + (mf_sg_new(shgat(isg)%ord)-mf_sg(shgat(isg)%ord))**2
          correct_search = .true.
        end if
        
      end do
      
      !**********************************************************!
#ifdef usempi
        ! Max value of all grads.
        laux = correct_search
        call MPI_ALLREDUCE( laux, correct_search, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr )
        
        ! 
        raux = perc_wrong
        call MPI_ALLREDUCE( raux, perc_wrong, 1, mpi_fwi, MPI_SUM, MPI_COMM_WORLD, ierr )
        
        ! 
        raux = perc_ok
        call MPI_ALLREDUCE( raux, perc_ok, 1, mpi_fwi, MPI_SUM, MPI_COMM_WORLD, ierr )
#endif
      
      perc_wrong = sqrt(perc_wrong/(perc_wrong+perc_ok))
      
!       !**********************************************************!
!       ! 
!       if ( correct_search ) then
!       do isg=1,n_sg(rank+1)
!         
!         ! 
!         call poly_fit_min( minimum, alp_1_global, alp_2_global, mf_sg(shgat(isg)%ord), mf_sg_1(shgat(isg)%ord), &
!                            mf_sg_2(shgat(isg)%ord), alp_fit, mf_fit_inter )
!         
!         ! 
!         if ( minimum ) then
!           
!           coef = max( 0.05_sp, min( 2._sp, alp_fit/alp_min_global ) )
!           coef_sea(shgat(isg)%ord) = coef*coef_sea(shgat(isg)%ord)
!           write(*,'(i4,8(1x,f7.3))') shgat(isg)%id, coef, alp_1_global, alp_2_global, 1000.*mf_sg(shgat(isg)%ord), &
!             1000.*mf_sg_1(shgat(isg)%ord),1000.*mf_sg_2(shgat(isg)%ord), alp_fit, 1000.*mf_fit_inter
!           
!         else
!           
!           if ( mf_sg_2(shgat(isg)%ord) < mf_sg(shgat(isg)%ord) ) then
!             
!             coef = max( 0.05_sp, min( 2._sp, alp_fit/alp_min_global ) )
!             coef_sea(shgat(isg)%ord) = coef*coef_sea(shgat(isg)%ord)
!             write(*,'(i4,8(1x,f7.3))') shgat(isg)%id, coef, alp_1_global, alp_2_global, 1000.*mf_sg(shgat(isg)%ord), &
!               1000.*mf_sg_1(shgat(isg)%ord),1000.*mf_sg_2(shgat(isg)%ord), alp_fit, 1000.*mf_fit_inter
!             
!           else
!             
!             coef = 0.05_sp
!             coef_sea(shgat(isg)%ord) = coef*coef_sea(shgat(isg)%ord)
!             write(*,'(i4,8(1x,f7.3))') shgat(isg)%id, coef, alp_1_global, alp_2_global, 1000.*mf_sg(shgat(isg)%ord), &
!               1000.*mf_sg_1(shgat(isg)%ord),1000.*mf_sg_2(shgat(isg)%ord), alp_fit, 1000.*mf_fit_inter
!             
!           end if
!           
!         end if
!         
!       end do
!       
!       end if
      
      !**********************************************************!
      !**********************************************************!
      !**********************************************************!
      ! 
      if ( correct_search ) then
        
        !----------------------------------------------------------!
        ! 
        mat_aux = 0._rl
        do isg=1,n_sg(rank+1)
          
          !----------------------------------------------------------!
          ! 
          do iz=1,shgat(isg)%mg%iLz
          iz_shift = iz+(shgat(isg)%mg%iz-1)
          if ( iz_shift<=nz ) then
          do iy=1,shgat(isg)%mg%iLy
          iy_shift = iy+(shgat(isg)%mg%iy-1)
          if ( iy_shift<=ny ) then
          do ix=1,shgat(isg)%mg%iLx
          ix_shift = ix+(shgat(isg)%mg%ix-1)
          if ( ix_shift<=nx ) then
            mat_aux(ix_shift,iy_shift,iz_shift) = mat_aux(ix_shift,iy_shift,iz_shift) + &
                                                  coef_sea(shgat(isg)%ord)*grad_sg(isg,ix,iy,iz)
          end if
          end do
          end if
          end do
          end if
          end do
          
        end do
        
#ifdef usempi
        !----------------------------------------------------------!
        ! Max value of all grads.
        raux = maxval(abs(mat_aux))
        call MPI_ALLREDUCE( raux, r_max_mpi, 1, mpi_fwi, MPI_MAX, MPI_COMM_WORLD, ierr )
        
        ! Sum gradient subsums for each process + recover norm.
        do iy=1,ny
        do iz=1,nz
          ! Reduce precision.
          sr_aux = mat_aux(:,iy,iz)/r_max_mpi
          ! Send + sum.
          call MPI_ALLREDUCE( sr_aux, sr_aux2, nx, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr )
          ! Recover precision.
          mat_aux(:,iy,iz) = r_max_mpi*dble(sr_aux2)
        end do
        end do
#endif
        
!         !----------------------------------------------------------!
!         ! Model preconditioner.
!         mat_aux = mat_aux*prec_mod
        
        !----------------------------------------------------------!
        ! Smooth operator.
        call apply_smooth_op( mat_aux, 0.000_sp, freq_inv )
        
        !----------------------------------------------------------!
        ! Clean water + Layer frezzing.
        do iz=1,nz
        do ix=1,nx
          ny0 = grid%ind_bath(ix,iz)-1
          do iy=1,ny0+lay_frz
            mat_aux(ix,iy,iz) = 0._rl
          end do
        end do
        end do
        
        !----------------------------------------------------------!
        ! 
        search = (sum(sum(sum(search*mat_aux,3),2),1)/sum(sum(sum(mat_aux*mat_aux,3),2),1)) * mat_aux
        
        !----------------------------------------------------------!
        ! 
#ifdef test_no_RK_loop
        search = 0.
#else
        if ( any(isnan(real(search))) ) then
          write(*,*) " ***** WARNING correct_search: any(isnan(real(search))) ***** "
          stop " ***** ERROR correct_search: any(isnan(real(search))) ***** "
        end if
#endif

      end if
      
      !**********************************************************!
      ! Deallocate.
      deallocate( sr_aux, sr_aux2, stat=stal ); if ( stal/=0 ) stop 'dAE correct_search 1'
      
    !*********************************************************************/
    end function correct_search
  
  !***********************************************************************/
  end module m_correct_search    




