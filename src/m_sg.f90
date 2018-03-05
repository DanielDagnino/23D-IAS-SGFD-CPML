  !*********************************************************************/
  !** This code has been done in the Barcelona Center for Subsurface 
  !** Imaging (BCSI).
  !** Authors: Daniel Dagnino.
  !*********************************************************************/
  !** Short resume of this code/file:
  !** 
  !** 
  !** 
  !** 
  !** 
  !** 
  !*********************************************************************/
  module m_sg
    use m_sg_data
    use m_sg_type
    use m_sg_func
    use m_sg_pos
    use m_sg_tr
    use m_sg_sou
    implicit none
    
    !**********************************************************!
    !**********************************************************!
    !**********************************************************!
    
    contains
    
    !**********************************************************!
    subroutine rmv_sg( shgat, n_shot )
      implicit none
      ! The variables which are passed to the function.
      type(shot_gather), allocatable, intent(inout) :: shgat(:)
      integer, intent(in) :: n_shot
      ! The variables which are generated inside the function.
      integer :: stal
      integer :: ish
      
      !----------------------------------------------------------!    
      do ish=1,n_shot
        call rmv_one_sg( shgat(ish) )
      end do
      
      if ( allocated( shgat ) ) deallocate( shgat, stat=stal ); if ( stal/=0 ) stop "dAE rmv_sg"
      
    end subroutine rmv_sg
    
    !**********************************************************!
    subroutine rmv_one_sg( shgat )
      implicit none
      ! The variables which are passed to the function.
      type(shot_gather), intent(inout) :: shgat
      ! The variables which are generated inside the function.
      integer :: stal
      
      !----------------------------------------------------------!
      if ( allocated( shgat%pos_rec ) ) then
        deallocate( shgat%pos_rec, stat=stal ); if ( stal/=0 ) stop "dAE rmv_one_shot_gather 1"
      end if
      
      if ( allocated( shgat%pos_sou ) ) then
        deallocate( shgat%pos_sou, stat=stal ); if ( stal/=0 ) stop "dAE rmv_one_shot_gather 2"
      end if
      
      if ( allocated( shgat%tr_obs ) ) then
        deallocate( shgat%tr_obs, stat=stal ); if ( stal/=0 ) stop "dAE rmv_one_shot_gather 3"
      end if
      
      if ( allocated( shgat%tr_syn ) ) then
        deallocate( shgat%tr_syn, stat=stal ); if ( stal/=0 ) stop "dAE rmv_one_shot_gather 5"
      end if
      
      if ( allocated( shgat%sou ) ) then
        deallocate( shgat%sou, stat=stal ); if ( stal/=0 ) stop "dAE rmv_one_shot_gather 6"
      end if
      
      if ( allocated( shgat%adj_sou ) ) then
        deallocate( shgat%adj_sou, stat=stal ); if ( stal/=0 ) stop "dAE rmv_one_shot_gather 7"
      end if
      
      if ( allocated( shgat%pr_dat ) ) then
        deallocate( shgat%pr_dat, stat=stal ); if ( stal/=0 ) stop "dAE rmv_one_shot_gather 8"
      end if
      
      if ( allocated( shgat%kill_tr ) ) then
        deallocate( shgat%kill_tr, stat=stal ); if ( stal/=0 ) stop "dAE rmv_one_shot_gather 9"
      end if
      
      if ( allocated( shgat%calib ) ) then
        deallocate( shgat%calib, stat=stal ); if ( stal/=0 ) stop "dAE rmv_one_shot_gather 13"
      end if
      
      if ( allocated( shgat%pick_o ) ) then
        deallocate( shgat%pick_o, stat=stal ); if ( stal/=0 ) stop "dAE rmv_one_shot_gather 10"
      end if
      
      if ( allocated( shgat%pick_s ) ) then
        deallocate( shgat%pick_s, stat=stal ); if ( stal/=0 ) stop "dAE rmv_one_shot_gather 11"
      end if
      
      if ( allocated( shgat%dly ) ) then
        deallocate( shgat%dly, stat=stal ); if ( stal/=0 ) stop "dAE rmv_one_shot_gather 12"
      end if
      
      if ( allocated( shgat%kill_rec_pos_out ) ) then
        deallocate( shgat%kill_rec_pos_out, stat=stal ); if ( stal/=0 ) stop "dAE rmv_one_shot_gather 13"
      end if
      
    end subroutine rmv_one_sg
    
    !**********************************************************!
  end module m_sg
  
  
  
  
  