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
  !*********************************************************************/
  module m_sg_pos
    implicit none
    
    !*********************************************************************/
    contains
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine get_sg_pos( folder, shgat, n_shot, offset_min, offset_max )
      use m_data_kind
      use m_sg_type
      use m_geo
      implicit none
      
      ! The variables which are passed to the function.
      type(shot_gather), allocatable, intent(inout) :: shgat(:)
      integer, intent(inout) :: n_shot
      character(200), intent(in) :: folder
      real(sp), intent(in) :: offset_min, offset_max
      ! The variables which are generated inside the function.
      integer :: ish
      
      !----------------------------------------------------------!    
      ! Read all shots.
      do ish=1,n_shot
        call get_one_sg_pos( folder, shgat(ish), offset_min, offset_max )
      end do
      
      !----------------------------------------------------------!    
    end subroutine get_sg_pos
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine get_sg_bath( shgat, n_shot )
      use m_data_kind
      use m_sg_type
      use m_geo
      implicit none
      
      ! The variables which are passed to the function.
      type(shot_gather), allocatable, intent(inout) :: shgat(:)
      integer, intent(inout) :: n_shot
      ! The variables which are generated inside the function.
      integer :: ish
      
      !----------------------------------------------------------!    
      ! Read all shots.
      do ish=1,n_shot
        call get_one_sg_bath( shgat(ish) )
      end do
      
      !----------------------------------------------------------!    
    end subroutine get_sg_bath
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine get_one_sg_pos( folder, shgat, offset_min, offset_max )
      use m_data_kind
      use m_sg_type
      use m_sg_func
      use m_dyw, only: dyw_3d
      use m_sg_path_file
      use m_geo
      implicit none
      ! The variables which are passed to the function.
      type(shot_gather), intent(inout) :: shgat
      character(200), intent(in) :: folder      
      real(sp), intent(in) :: offset_min, offset_max
      ! The variables which are generated inside the function.
      character(100) :: str_aux
      integer :: is, ir
      real(sp) :: x, y, z
      character(200) :: path_file
      integer :: stal
      
      !----------------------------------------------------------!    
      ! Read all shots.
      ! Open file.
      path_file = shot_gather_path_file( folder, 'sou_rec/ps', shgat%id )
      open(unit=1,file=path_file,status='old',action='read')
      
      !----------------------------------------------------------!
      ! Reading source position.
      read(1,'(a)') str_aux
      read(1,*) shgat%n_sou
      
      ! 
#ifdef testing
      if ( allocated( shgat%pos_sou ) ) stop "ERROR get_one_sg_pos: shgat%pos_sou already allocated"
#endif
      
      allocate( shgat%pos_sou(shgat%n_sou), stat=stal ); if ( stal/=0 ) stop "AE get_one_sg_pos"
      
      ! 
      do is=1,shgat%n_sou
        read(1,*) x, y, z
        x = x/1000._sp
        y = y/1000._sp
        z = z/1000._sp
        if ( dyw_3d ) z = 0._sp
        shgat%pos_sou(is)%x = x
        shgat%pos_sou(is)%y = y
        shgat%pos_sou(is)%z = z
        shgat%pos_sou(is)%ix = 1+anint(x/dx)
        shgat%pos_sou(is)%iy = 1+anint(y/dx)
        shgat%pos_sou(is)%iz = 1+anint(z/dx)
        
        if ( dx*real(nx-1,sp)<x .or. 0._sp>x ) stop "***** ERROR get_one_sg_pos: pos_sou x out of the model. *****"
        if ( dx*real(ny-1,sp)<y .or. 0._sp>y ) stop "***** ERROR get_one_sg_pos: pos_sou y out of the model. *****"
        if ( dx*real(nz-1,sp)<z .or. 0._sp>z ) stop "***** ERROR get_one_sg_pos: pos_sou z out of the model. *****"
      end do
      
      !----------------------------------------------------------!  
      ! Reading receiver positions.
      read(1,'(a)') str_aux
      read(1,'(a)') str_aux
      read(1,*) shgat%n_rec_file_pos
      
#ifdef testing
      if ( allocated( shgat%pos_rec ) ) stop "ERROR get_sg_pos: shgat%pos_rec already allocated"
#endif
      
      allocate( shgat%pos_rec(shgat%n_rec_file_pos), stat=stal ); if ( stal/=0 ) stop "AE get_sg_pos"
      
      ! 
      do ir=1,shgat%n_rec_file_pos
        read(1,*) x, y, z
        x = x/1000._sp
        y = y/1000._sp
        z = z/1000._sp
        if ( dyw_3d ) z = 0._sp
        
        shgat%pos_rec(ir)%x = x
        shgat%pos_rec(ir)%y = y
        shgat%pos_rec(ir)%z = z
        shgat%pos_rec(ir)%ix = 1+floor(x/dx)
        shgat%pos_rec(ir)%iy = 1+floor(y/dx)
        shgat%pos_rec(ir)%iz = 1+floor(z/dx)
      end do
      
      ! 
      call fix_pos_offset_sou_rec_inside( shgat, offset_min, offset_max )
      
      ! 
      do ir=1,shgat%n_rec
        x = shgat%pos_rec(ir)%x
        y = shgat%pos_rec(ir)%y
        z = shgat%pos_rec(ir)%z
        if ( dx*real(nx-1,sp)<x .or. 0._sp>x ) stop "***** ERROR get_one_sg_pos: pos_rec x out of the model. *****"
        if ( dx*real(ny-1,sp)<y .or. 0._sp>y ) stop "***** ERROR get_one_sg_pos: pos_rec y out of the model. *****"
        if ( dx*real(nz-1,sp)<z .or. 0._sp>z ) stop "***** ERROR get_one_sg_pos: pos_rec z out of the model. *****"
        if ( shgat%pos_rec(ir)%ix<0 ) stop "***** ERROR get_one_sg_pos: shgat%pos_rec(ir)%ix<0 *****"
        if ( shgat%pos_rec(ir)%iy<0 ) stop "***** ERROR get_one_sg_pos: shgat%pos_rec(ir)%iy<0 *****"
        if ( shgat%pos_rec(ir)%iz<0 ) stop "***** ERROR get_one_sg_pos: shgat%pos_rec(ir)%iz<0 *****"
      end do
      
!       !----------------------------------------------------------!    
!       ! 
!       write(*,*) 'sou x', minval(shgat%pos_sou(:)%ix), maxval(shgat%pos_sou(:)%ix)
!       write(*,*) 'sou y', minval(shgat%pos_sou(:)%iy), maxval(shgat%pos_sou(:)%iy)
!       write(*,*) 'sou z', minval(shgat%pos_sou(:)%iz), maxval(shgat%pos_sou(:)%iz)
!       
!       write(*,*) 'rec x', minval(shgat%pos_rec(:)%ix), maxval(shgat%pos_rec(:)%ix)
!       write(*,*) 'rec y', minval(shgat%pos_rec(:)%iy), maxval(shgat%pos_rec(:)%iy)
!       write(*,*) 'rec z', minval(shgat%pos_rec(:)%iz), maxval(shgat%pos_rec(:)%iz)
      
      !----------------------------------------------------------!
      close(1)
      
    end subroutine get_one_sg_pos
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine get_one_sg_bath( shgat )
      use m_data_kind
      use m_sg_type
      use m_geo
      use m_solv3Dahvc_grid_bath
      implicit none
      ! The variables which are passed to the function.
      type(shot_gather), intent(inout) :: shgat
      ! The variables which are generated inside the function.
      integer :: ix, iz
      real(sp) :: x, z, x0, z0, dist, dist_min
      
      !----------------------------------------------------------!    
      ! Surface position in the global model.
      x0 = shgat%pos_sou(1)%x
      z0 = shgat%pos_sou(1)%z
      
      x0 = x0 + dx*dble(shgat%mg%ix-1)
      z0 = z0 + dx*dble(shgat%mg%iz-1)
      
      ! Searc min distance.
      dist_min = huge(1.)
      do iz=1,nz
      z = dx*dble(iz-1)
        do ix=1,nx
        x = dx*dble(ix-1)
          
          ! Surface distance: point surface-source.
          dist = (x-x0)**2 + (z-z0)**2
          
          ! Min distance.
          if ( dist<dist_min ) then
            shgat%shot_bath = grid%bath(ix,iz)
            dist_min = dist
          end if
          
        end do
      end do
      
!         dist_min = huge(1.)
!         do iz=1,shgat%mg%iLz
!         iz_shift = iz+(shgat%mg%iz-1)
!         z  = dx*dble(iz_shift-1)
!         if ( iz_shift<=nz ) then
!           do ix=1,shgat%mg%iLx
!           ix_shift = ix+(shgat%mg%ix-1)
!           x  = dx*dble(ix_shift-1)
!           if ( ix_shift<=nx ) then
!             
!             ! Surface distance: point surface-source.
!             dist = (x-x0)**2 + (z-z0)**2
!             
!             ! Min dist.
!             if ( dist<dist_min ) then
!               shgat%shot_bath = grid%bath(ix_shift,iz_shift)
!               dist_min = dist
!             end if
!             
!           end if
!           end do
!         end if
!         end do
      
    end subroutine get_one_sg_bath
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine allo_one_sg_calib( shgat )
      use m_data_kind
      use m_sg_type
      implicit none
      ! The variables which are passed to the function.
      type(shot_gather), intent(inout) :: shgat
      ! The variables which are generated inside the function.
      integer :: stal
      
      if ( .not. allocated( shgat%calib ) ) then
        allocate( shgat%calib(shgat%n_rec), stat=stal ); if ( stal/=0 ) stop "AE allo_one_sg_calib"
      end if
      
      shgat%calib = 0._sp
      
    end subroutine allo_one_sg_calib
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine allo_one_sg_kill( shgat )
      use m_data_kind
      use m_sg_type
      implicit none
      ! The variables which are passed to the function.
      type(shot_gather), intent(inout) :: shgat
      ! The variables which are generated inside the function.
      integer :: stal
      
      if ( .not. allocated( shgat%kill_tr ) ) then
        allocate( shgat%kill_tr(shgat%n_rec), stat=stal ); if ( stal/=0 ) stop "AE allo_one_sg_kill"
      end if
      
      shgat%kill_tr = .false.
      
    end subroutine allo_one_sg_kill
    
    
    
  !**********************************************************!
  end module m_sg_pos
  
  
  
  
  
