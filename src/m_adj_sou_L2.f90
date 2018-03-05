  !*********************************************************************/
  !*****     Author: D. Dagnino     ************************************/
  !*********************************************************************/
  module m_adj_sou_L2
    implicit none
    
    !**********************************************************!
    !**********************************************************!
    contains
    
    !**********************************************************!
    real(sp) function misfit_L2( n_sg, shgat )
      use m_data_kind
      use m_mpi
      use m_sg_type
#ifdef usempi
      use mpi
#endif
      implicit none
      ! The variables which are passed to the function.
      integer, allocatable, intent(in) :: n_sg(:)
      type(shot_gather), allocatable, intent(inout) :: shgat(:)
      ! The variables which are generated inside the function.
      integer :: isg
      
      !**********************************************************!
      ! L2
      misfit_L2 = 0._sp
      
      ! Shot.
      do isg=1,n_sg(rank+1)
        ! error
        misfit_L2 = misfit_L2 + misfit_one_L2( shgat(isg), isg )
      end do
      
      !----------------------------------------------------------!
    end function misfit_L2
    
    !**********************************************************!
    real(sp) function misfit_one_L2( shgat, isg )
      use m_data_kind
      use m_mpi
      use m_sg_type
      use m_sg_tr
      use m_sg_pr
      use m_calib
      use m_type_inversion, only: calib, rank_wrt, isg_wrt, prec_d, offset_pr_min, offset_pr_max
      use m_inv_data
      implicit none
      ! The variables which are passed to the function.
      type(shot_gather), intent(inout) :: shgat
      integer, intent(in) :: isg
      ! The variables which are generated inside the function.
      integer :: ir
      real(sp) :: dt, weigth
      character(200) :: file_name     
      character(50)  :: str_tmp
      
      !----------------------------------------------------------!
      ! Adapt traces.
      call adapt_one_sg_tr_both( shgat, shgat%nt_sou, shgat%t1_obs )
      call adapt_prec( shgat, shgat%nt_sou )
      
      ! New time grid.
      dt = shgat%t1_obs/real(shgat%nt_sou-1,sp)
      
      !----------------------------------------------------------!
      ! Calibrate.
      if ( calib ) then
        call get_one_sg_pr( shgat, 31, freq_inv, offset_pr_min, offset_pr_max )
        call coef_calib2( shgat )
        call get_one_sg_pr( shgat, prec_d, freq_inv, offset_pr_min, offset_pr_max )
      else
        shgat%calib = 1._sp
      end if
      
      !----------------------------------------------------------!
      ! 
      misfit_one_L2 = 0._sp
      weigth = 0._sp
      do ir=1,shgat%n_rec
        
        ! L2        
        misfit_one_L2 = misfit_one_L2 + dt*sum( shgat%pr_dat(ir,:)*( ( shgat%calib(ir)*shgat%tr_obs(ir,:) - shgat%tr_syn(ir,:) )**2 ) )
        ! weigth
        weigth = weigth + dt*sum( shgat%pr_dat(ir,:)*( (shgat%calib(ir)*shgat%tr_obs(ir,:))**2 ) )
        
      end do
      
!       misfit_one_L2 = 100._sp*misfit_one_L2
!       misfit_one_L2 = 100._sp*misfit_one_L2/weigth
      misfit_one_L2 = 100._sp*(misfit_one_L2/weigth)*dble(shgat%n_rec)
      
      !----------------------------------------------------------!
      ! Write shgat%calib.
      if ( rank == rank_wrt .and. isg == isg_wrt ) then
        ! File.
        file_name = 'data_output/calib/a_calib_mf'
        write(str_tmp,*) strat_inv
        file_name = trim( file_name ) // '_strat_'// trim(adjustl(str_tmp))
        write(str_tmp,fmt='(f5.2)') freq_inv
        file_name = trim( file_name ) // '_freq_'// trim(adjustl(str_tmp))
        write(str_tmp,*) iter_freq
        file_name = trim( file_name ) // '_iter_freq_'// trim(adjustl(str_tmp))
        file_name = trim( file_name ) // '.txt'
        
        ! Write.
        open(unit=1,file=file_name,status='unknown',action='write')
          do ir=1,shgat%n_rec
            write(1,'(i5,1x,es13.6E2)') ir, real(shgat%calib(ir))
          end do
        close(1)
      end if
      
      !----------------------------------------------------------!
    end function misfit_one_L2
    
    
    
    !**********************************************************!
    subroutine adj_sou_L2( shgat, error, isg )
      use m_data_kind
      use m_mpi
      use m_sg_type
      use m_sg_tr
      use m_sg_pr
      use m_sg_adj_sou
      use m_calib
      use m_type_inversion, only: calib, rank_wrt, isg_wrt, prec_d, offset_pr_min, offset_pr_max
      use m_inv_data
      implicit none
      ! The variables which are passed to the function.
      type(shot_gather), intent(inout) :: shgat
      real(sp), intent(out) :: error
      integer, intent(in) :: isg
      ! The variables which are generated inside the function.
      real(sp) :: raux
      integer :: ir, it
      real(sp) :: dt, weigth
      character(200) :: file_name     
      character(50)  :: str_tmp
      
      !----------------------------------------------------------!
      ! Adapt traces.
      call adapt_one_sg_tr_both( shgat, shgat%nt_sou, shgat%t1_obs )
      call adapt_prec( shgat, shgat%nt_sou )
      
      ! Adapt adj.
      call allo_one_sg_adj_sou( shgat, shgat%nt_sou )
      shgat%t0_adj = shgat%t0_obs
      shgat%t1_adj = shgat%t1_obs
      
      ! New time grid.
      dt = shgat%t1_adj/real(shgat%nt_adj-1,sp)
      
      !----------------------------------------------------------!
      ! Calibrate.
      if ( calib ) then
        call get_one_sg_pr( shgat, 31, freq_inv, offset_pr_min, offset_pr_max )
        call coef_calib2( shgat )
        call get_one_sg_pr( shgat, prec_d, freq_inv, offset_pr_min, offset_pr_max )
      else
        shgat%calib = 1._sp
      end if
      
      !----------------------------------------------------------!
      ! 
      error = 0._sp
      weigth = 0._sp
      do it=1,shgat%nt_adj
      do ir=1,shgat%n_rec
        
        ! L2
        raux = -( shgat%calib(ir)*shgat%tr_obs(ir,it) - shgat%tr_syn(ir,it) )
        ! weigth
        weigth = weigth + dt*shgat%pr_dat(ir,it)*( (shgat%calib(ir)*shgat%tr_obs(ir,it))**2 )
        ! error
        error = error + dt*shgat%pr_dat(ir,it)*(raux**2)
        ! adj sou
        shgat%adj_sou(ir,it) = shgat%pr_dat(ir,it)*raux
        
      end do
      end do
      
!       error = 100._sp*error
      error = 100._sp*(error/weigth)*dble(shgat%n_rec)
      shgat%adj_sou = (shgat%adj_sou/weigth)*dble(shgat%n_rec)
!       error = 100._sp*error/weigth
!       shgat%adj_sou = shgat%adj_sou/weigth
      
      !----------------------------------------------------------!
      ! Write a_calib.
      if ( rank == rank_wrt .and. isg == isg_wrt ) then
        ! File.
        file_name = 'data_output/calib/a_calib_adj'
        write(str_tmp,*) strat_inv
        file_name = trim( file_name ) // '_strat_'// trim(adjustl(str_tmp))
        write(str_tmp,fmt='(f5.2)') freq_inv
        file_name = trim( file_name ) // '_freq_'// trim(adjustl(str_tmp))
        write(str_tmp,*) iter_freq
        file_name = trim( file_name ) // '_iter_freq_'// trim(adjustl(str_tmp))
        file_name = trim( file_name ) // '.txt'
        
        ! Write.
        open(unit=1,file=file_name,status='unknown',action='write')
          do ir=1,shgat%n_rec
            write(1,'(i5,1x,es13.6E2)') ir, real(shgat%calib(ir))
          end do
        close(1)
      end if
      
      !----------------------------------------------------------!
    end subroutine adj_sou_L2
    
    !***********************************************************************/
  end module m_adj_sou_L2
  
  
  
  
  