  !*********************************************************************/
  !*****     Author: D. Dagnino     ************************************/
  !*********************************************************************/
  module m_Hessian
    use m_data_kind
    use m_geo
    use m_type_inversion, only: use_iH, full_H
    implicit none
    real(sp) :: alpha_step
    real(sp), allocatable :: inv_H(:,:)
    real(sp), allocatable :: inv_H_diag(:)
    real(sp), allocatable :: Wd(:)
    integer, allocatable :: iW(:), jW(:)
    integer :: nxH, nyH, nzH
    integer :: dimWd
    
    !**********************************************************!
    !**********************************************************!
    !**********************************************************!
    contains
    
    
    
!     !**********************************************************!
!     !**********************************************************!
!     subroutine get_iH( folder )
!       use m_data_kind
!       implicit none
!       ! The variables which are passed to the function.
!       character(200), intent(in) :: folder
!       ! The variables which are generated inside the function.
!       integer :: ix_H, iy_H, jx_H, jy_H, i_H, j_H, k
!       integer :: stal
!       real(sp), allocatable :: ar_aux(:)
!       character(200) :: path_file
!       
!       !----------------------------------------------------------!    
!       ! Open file.
!       path_file = trim( folder ) // '/inv_Hij.txt'
!       
!       open(unit=1,file=path_file,status='old',action='read')
!       
!       ! Reading time
!       read(1,*) nxH, nyH
!       nzH = 1
!       
!       !----------------------------------------------------------!
! !     if ( use_iH ) then
! !       if ( .not. allocated( inv_H ) ) then
!         if ( full_H ) then
!           allocate( inv_H(nxH*nyH*nzH,nxH*nyH*nzH), stat=stal ); if ( stal/=0 ) stop 'AE get_iH'
!           inv_H = 0.
!         else
!           allocate( inv_H(nx,ny,nz), stat=stal ); if ( stal/=0 ) stop 'AE get_iH'
!           inv_H = 0.
!         end if
! !       end if
! !     end if
!       
!       allocate( ar_aux(nxH*nyH), stat=stal ); if ( stal/=0 ) stop 'AE get_iH 2'
!       
!       !----------------------------------------------------------!
!       ! 
!       do iy_H=1,nyH
!       do ix_H=1,nxH
!         
!         read(1,*) ( ar_aux(k), k=1,nxH*nyH )
!         
!         i_H = ix_H+nxH*(iy_H-1)
!         
!         do jy_H=1,nyH
!         do jx_H=1,nxH
!           
!           j_H = jx_H+nxH*(jy_H-1)
!           
!           inv_H(i_H,j_H) = ar_aux(j_H)
!           
!         end do
!         end do
!         
!       end do
!       end do
!       
!       !----------------------------------------------------------!    
!       ! 
! !       write(*,*) ' nxH, nyH = ', nxH, nyH
! !       write(*,*) ' inv_H(,) = ', inv_H(1,1), inv_H(1,2), inv_H(1,3), inv_H(1,nxH), inv_H(1,nxH*nyH)
! !       write(*,*) ' inv_H(,) = ', inv_H(nxH,nxH), inv_H(nxH,2), inv_H(nxH,3), inv_H(nxH,nxH), inv_H(nxH,nxH*nyH)
! !       write(*,*) ' inv_H(,) = ', inv_H(5*nxH,5*nxH), inv_H(5*nxH,2), inv_H(5*nxH,3), inv_H(5*nxH,nxH), inv_H(5*nxH,nxH*nyH)
! !       write(*,*) ' inv_H(,) = ', inv_H(10*nxH,10*nxH)
! !       write(*,*) ' inv_H(,) = ', inv_H(nxH*nyH,nxH*nyH)
!       
!       do ix_H=1,nxH
!       do iy_H=1,nyH
!         
!         i_H = ix_H+nxH*(iy_H-1)
!         
!         if ( inv_H(i_H,i_H) <= 0._sp ) then
!           write(*,*) ' i_H = ', i_H
!           stop ' ********** ERROR: inv_H(i_H,i_H) <= 0._sp ********** '
!         end if
!         
!       end do
!       end do
!       
!       !----------------------------------------------------------!    
!       deallocate( ar_aux, stat=stal ); if ( stal/=0 ) stop 'dAE get_iH 2'
!       
!       !----------------------------------------------------------!
!       close(1)
!       
!     end subroutine get_iH
    
    
    
    !**********************************************************!
    !**********************************************************!
    subroutine get_iH_diag( iLx, iLy, iLz, folder )
      use m_data_kind
      implicit none
      ! The variables which are passed to the function.
      integer, intent(in) :: iLx, iLy, iLz
      character(200), intent(in) :: folder
      ! The variables which are generated inside the function.
      integer :: i, j, k, dimH
      integer :: stal
      character(200) :: path_file
      real(sp) :: raux
      
      !----------------------------------------------------------!    
      ! Open file.
      path_file = trim( folder ) // '/inv_Hessian_diag.txt'
      
      open(unit=1,file=path_file,status='old',action='read')
      
      !       
!       if ( nx/=nxH .or. ny/=nyH ) then
!         write(*,*) ' ERROR - get_iH_diag: nx, nxH, ny, nyH = ', nx, nxH, ny, nyH
!         stop ' ********** ERROR - get_iH_diag: nx/=nxH .or. ny/=nyH ********** '
!       end if
      
      !----------------------------------------------------------!
      dimH = iLx*iLy
      if ( .not. allocated( inv_H_diag ) ) then
        allocate( inv_H_diag(dimH), stat=stal ); if ( stal/=0 ) stop 'AE get_iH 1'
      end if
      
      inv_H_diag = 0._sp
      
      !----------------------------------------------------------!
      ! 
      do k=1,dimH
        read(1,*) i, j, raux
        if ( i/=j ) stop ' ********** ERROR: get_iH_diag i/=j ********** '
        inv_H_diag(i) = 1._sp/raux
      end do
      
      !----------------------------------------------------------!    
      ! 
      do i=1,dimH
        
        if ( inv_H_diag(i) <= 0._sp ) then
          write(*,*) ' i,j = ', i,j
          stop ' ********** ERROR: inv_H_diag(i,j) <= 0._sp ********** '
        end if
        
      end do
      
      !----------------------------------------------------------!
      close(1)
      
    end subroutine get_iH_diag
    
    
    
    !**********************************************************!
    !**********************************************************!
    subroutine get_Wd( folder )
      use m_data_kind
      implicit none
      ! The variables which are passed to the function.
      character(200), intent(in) :: folder
      ! The variables which are generated inside the function.
      integer :: k
      integer :: stal
      character(200) :: path_file
      
      !----------------------------------------------------------!    
      ! Open file.
      path_file = trim( folder ) // '/Wd.txt'
      
      open(unit=1,file=path_file,status='old',action='read')
      
      !----------------------------------------------------------!
      read(1,*) dimWd
      
      if ( .not. allocated( Wd ) ) then
        allocate( Wd(dimWd), stat=stal ); if ( stal/=0 ) stop 'AE get_iH 1'
        allocate( iW(dimWd), stat=stal ); if ( stal/=0 ) stop 'AE get_iH 2'
        allocate( jW(dimWd), stat=stal ); if ( stal/=0 ) stop 'AE get_iH 3'
      end if
      
      !----------------------------------------------------------!
      ! 
      do k=1,dimWd
        read(1,*) iW(k), jW(k), Wd(k)
      end do
      
      !----------------------------------------------------------!
      close(1)
      
    end subroutine get_Wd    
    
    
    
    
    
!     !**********************************************************!
!     !**********************************************************!
!     subroutine get_cal_Wd( offset, t_nmo, v_nmo, ny, tf, n_rec )
!       use m_data_kind
!       use m_mpi
!       implicit none
!       ! The variables which are passed to the function.
!       integer, intent(in) :: ny, n_rec
!       real(rl), intent(in) :: tf
!       real(4), allocatable, intent(in) :: t_nmo(:), v_nmo(:), offset(:)
!       ! The variables which are generated inside the function.
!       integer :: k, k0, nt_Wd, n_rec_aux
!       integer :: stal
!       integer :: ir1, ir2, cont, nr
!       integer :: it, it0, it1, it2, it_last, its0, it1_aux, it2_aux
!       real(4) :: t0, t1, t2, t0a, e1, e2, dt, raux
!       integer :: ie1, ie2
!       real(4), allocatable :: iv_nmo(:)
!       logical :: succes
!       
!       !----------------------------------------------------------!    
!       ! Half correlation offset in receiver points (ir-nr,ir+nr).
!       nr = 10
!       
!       ! Fix shot-gather.
!       nt_Wd = 1021
!       n_rec_aux = 192
!       
!       !----------------------------------------------------------!
!       ! Allocate.
! !       write(*,*) 'A - hola 1'
!       allocate( iv_nmo(nt_Wd), stat=stal ); if ( stal/=0 ) stop 'AE get_cal_Wd 1'
!       
!       ! Time step.
!       dt = real(tf,4)/real(nt_Wd-1,4)
!       
!       ! 
!       its0 = nint(0.03/dt)
!       
!       !----------------------------------------------------------!    
!       ! Interpolation of a monotonus increasing function in a irregular grid.
! !       write(*,*) 'A - hola 2'
!       iv_nmo = 0.
!       k0 = 2
!       do it=1,nt_Wd
!         t0 = dt*real(it-1,4)
!         
!         ! Loop to look the middle point.
!         do k=k0,ny          
!           ie1 = nint( (t_nmo(k-1)-t0)/dt )
!           ie2 = nint( (t_nmo(k)-t0)/dt )
!           ! Are we in the middle?
!           if ( ie1*ie2 <= 0 ) then
!             e1 = t_nmo(k-1)-t0
!             e2 = t_nmo(k)-t0
!             raux = -e1 + e2
!             iv_nmo(it) = (1.+(e1/raux))*v_nmo(k-1) + (1.-(e2/raux))*v_nmo(k)
!             k0 = k   ! Since it is monot. we start the new search in this last point.
!             it_last = it   ! Last point that we interpolate.
!             exit
!           end if
!         end do
!         
!       end do
!       
!       ! Repeat last point to end the interpolation (extrapolation: nearest-neighbor).
!       do it=it_last+1,nt_Wd
!         iv_nmo(it) = iv_nmo(it_last)
!       end do
!       
! !       write(*,*) 'iv_nmo(1) = ', iv_nmo(1)
! !       write(*,*) 'iv_nmo(100) = ', iv_nmo(100)
! !       write(*,*) 'iv_nmo(nt_Wd-100) = ', iv_nmo(nt_Wd-100)
! !       write(*,*) 'iv_nmo(nt_Wd) = ', iv_nmo(nt_Wd)
!       
!       !----------------------------------------------------------!    
!       ! Inverse.
! !       write(*,*) 'A - hola 3'
!       do it=1,nt_Wd
!         iv_nmo(it) = 1./iv_nmo(it)
!       end do
!       
!       !----------------------------------------------------------!
!       ! 
! !       write(*,*) 'A - hola 4'
!       cont = 0
!       ! Point in the shot-gather.
!       do ir1=1,n_rec
!       it0 = 2+its0
!       do it1=1+its0,nt_Wd
!         ! Time.
!         t1 = dt*real(it1-1,4)
!         
!         ! Looking for the point in another receiver with the same reflexion curve.
!         ! Loop over zero-offset time where we start the curve that pass throw (ir1,it1). We solve the eq:
!         ! t0 = sqrt( t1**2 - (offset(ir1)/vel(t0))**2 )
!         succes = .false.
! !         do it=it0,it1+its0
!         do it=min(it1+its0+1,nt_Wd),it0,-1
! !           t0a = dt*real(it-2,4)
! !           t0  = dt*real(it-1,4)
! !           ie1 = nint( (t0a - sqrt( t1**2 - (offset(ir1)*iv_nmo(it-1))**2 ) )/dt )
! !           ie2 = nint( (t0  - sqrt( t1**2 - (offset(ir1)*iv_nmo(it  ))**2 ) )/dt )
!           ie1 = (it-2) - nint(sqrt( real((it1-1-its0)**2,4) - (offset(ir1)*iv_nmo(it-1)/dt)**2 ))
!           ie2 = (it-1) - nint(sqrt( real((it1-1-its0)**2,4) - (offset(ir1)*iv_nmo(it  )/dt)**2 ))
!           ! Min?
!           if ( ie1*ie2 <= 0 ) then
!             if ( abs(ie1) < abs(ie2) ) then
!               it0 = it-1
!             else
!               it0 = it
!             end if
! !             t0 = dt*real(it0-1,4)
!             succes = .true.
!             exit
!           end if
!         end do
!         
!         ! Connect point in the curve trought Wd: (ir1,it1)--->(ir2,it2).
!         if ( ( it0 > nint(1.62/dt) ) .and.( it0 < nint(1.74/dt) ) ) then
!         if ( succes ) then
! !         do ir2=max(1,ir1-nr),min(n_rec,ir1+nr)
!         do ir2=min(n_rec,ir1+nr),min(n_rec,ir1+nr)
!           ! t2
! !           t2 = sqrt( t0**2 + (offset(ir2)*iv_nmo(it0))**2 )
! !           it2 = nint(t2/dt)
!           it2 = nint(sqrt( real((it0-1)**2,4) + (offset(ir2)*iv_nmo(it0)/dt)**2 ))
!           ! Wd and index.
!           it1_aux = it1
!           it2_aux = its0 + it2
!           if ( (it1_aux<=nt_Wd) .and. (it1_aux>0) ) then
!           if ( (it2_aux<=nt_Wd) .and. (it2_aux>0) ) then
!             cont = cont + 1
!           end if
!           end if
!           
!         end do
!         end if
!         end if
!         
!       end do
!       end do
!       
!       !----------------------------------------------------------!    
!       ! Dimension Wd.
!       dimWd = cont
!       
!       if ( rank == 0 ) write(*,*) 'dimWd = ', dimWd
!       
!       ! Allocate.
!       if ( .not. allocated( Wd ) ) then
!         allocate( Wd(dimWd), stat=stal ); if ( stal/=0 ) stop 'AE get_cal_Wd 1'
!         allocate( iW(dimWd), stat=stal ); if ( stal/=0 ) stop 'AE get_cal_Wd 2'
!         allocate( jW(dimWd), stat=stal ); if ( stal/=0 ) stop 'AE get_cal_Wd 3'
!       else
!         deallocate( Wd, stat=stal ); if ( stal/=0 ) stop 'dAE get_cal_Wd 1'
!         deallocate( iW, stat=stal ); if ( stal/=0 ) stop 'dAE get_cal_Wd 2'
!         deallocate( jW, stat=stal ); if ( stal/=0 ) stop 'dAE get_cal_Wd 3'
!         allocate( Wd(dimWd), stat=stal ); if ( stal/=0 ) stop 'AE get_cal_Wd 1b'
!         allocate( iW(dimWd), stat=stal ); if ( stal/=0 ) stop 'AE get_cal_Wd 2b'
!         allocate( jW(dimWd), stat=stal ); if ( stal/=0 ) stop 'AE get_cal_Wd 3b'
!       end if
!       
!       iW = 0
!       jW = 0
!       Wd = 0.
!       
!       !----------------------------------------------------------!
!       ! 
! !       write(*,*) 'A - hola 5'
!       cont = 0
!       ! Point in the shot-gather.
!       do ir1=1,n_rec
!       it0 = 2+its0
!       do it1=1+its0,nt_Wd
!         ! Time.
!         t1 = dt*real(it1-1,4)
!         
!         ! Looking for the point in another receiver with the same reflexion curve.
!         ! Loop over zero-offset time where we start the curve that pass throw (ir1,it1). We solve the eq:
!         ! t0 = sqrt( t1**2 - (offset(ir1)/vel(t0))**2 )
!         succes = .false.
! !         do it=it0,it1+its0
!         do it=min(it1+its0+1,nt_Wd),it0,-1
! !           t0a = dt*real(it-2,4)
! !           t0  = dt*real(it-1,4)
! !           ie1 = nint( (t0a - sqrt( t1**2 - (offset(ir1)*iv_nmo(it-1))**2 ) )/dt )
! !           ie2 = nint( (t0  - sqrt( t1**2 - (offset(ir1)*iv_nmo(it  ))**2 ) )/dt )
!           ie1 = (it-2) - nint(sqrt( real((it1-1-its0)**2,4) - (offset(ir1)*iv_nmo(it-1)/dt)**2 ))
!           ie2 = (it-1) - nint(sqrt( real((it1-1-its0)**2,4) - (offset(ir1)*iv_nmo(it  )/dt)**2 ))
!           ! Min?
!           if ( ie1*ie2 <= 0 ) then
!             if ( abs(ie1) < abs(ie2) ) then
!               it0 = it-1
!             else
!               it0 = it
!             end if
! !             t0 = dt*real(it0-1,4)
!             succes = .true.
!             exit
!           end if
!         end do
!         
!         ! Connect point in the curve trought Wd: (ir1,it1)--->(ir2,it2).
!         if ( ( it0 > nint(1.62/dt) ) .and.( it0 < nint(1.74/dt) ) ) then
!         if ( succes ) then
! !         do ir2=max(1,ir1-nr),min(n_rec,ir1+nr)
!         do ir2=min(n_rec,ir1+nr),min(n_rec,ir1+nr)
!           ! t2
! !           t2 = sqrt( t0**2 + (offset(ir2)*iv_nmo(it0))**2 )
! !           it2 = nint(t2/dt)
!           it2 = nint(sqrt( real((it0-1)**2,4) + (offset(ir2)*iv_nmo(it0)/dt)**2 ))
!           ! Wd and index.
!           it1_aux = it1
!           it2_aux = its0 + it2
!           if ( (it1_aux<=nt_Wd) .and. (it1_aux>0) ) then
!           if ( (it2_aux<=nt_Wd) .and. (it2_aux>0) ) then
!             cont = cont + 1
!             iW(cont) = ir1+n_rec_aux*(it1_aux-1)
!             jW(cont) = ir2+n_rec_aux*(it2_aux-1)
!             Wd(cont) = 1._sp
!           end if
!           end if
!           
!         end do
!         end if
!         end if
!         
!       end do
!       end do
!       
! !       write(*,*) 'cont = ', cont
!       
!       if ( cont /= dimWd ) stop 'ERROR - get_cal_Wd: cont /= dimWd'
!       
!       !----------------------------------------------------------!
! !       write(*,*) 'A - hola 6'
!       deallocate( iv_nmo, stat=stal ); if ( stal/=0 ) stop 'dAE get_cal_Wd 1'
! !       write(*,*) 'A - hola 7'
!       
!     end subroutine get_cal_Wd    
    
    
    
    
    
    !**********************************************************!
    !**********************************************************!
    subroutine get_cal_Wd( offset, t_nmo, v_nmo, ny, tf, n_rec, isg )
      use m_data_kind
      use m_mpi
      use m_type_inversion, only: rank_wrt, isg_wrt
      implicit none
      ! The variables which are passed to the function.
      integer, intent(in) :: ny, n_rec, isg
      real(rl), intent(in) :: tf
      real(4), allocatable, intent(in) :: t_nmo(:), v_nmo(:), offset(:)
      ! The variables which are generated inside the function.
      integer :: k, k0, nt_Wd, n_rec_aux
      integer :: stal
      integer :: ir1, ir2, cont, nr
      integer :: it, it0, it1, it2, its0, it_last
      real(4) :: t0, e1, e2, dt, raux, e1_min, rit0
      integer :: ie1, ie2
      real(4), allocatable :: iv_nmo(:)
      logical :: succes
      
      !----------------------------------------------------------!    
      ! Half correlation offset in receiver points (ir-nr,ir+nr).
      nr = 10
      
      ! Fix shot-gather.
      nt_Wd = 1021
      n_rec_aux = 192
      
      !----------------------------------------------------------!
      ! Allocate.
      allocate( iv_nmo(nt_Wd), stat=stal ); if ( stal/=0 ) stop 'AE get_cal_Wd 1'
      
      ! Time step.
      dt = real(tf,4)/real(nt_Wd-1,4)
      
      ! 
      its0 = nint(0.03/dt)
      
      !----------------------------------------------------------!    
      ! Interpolation of a monotonus increasing function in a irregular grid.
!       write(*,*) 'A - hola 2'
      iv_nmo = 0.
      k0 = 2
      do it=1,nt_Wd
        t0 = dt*real(it-1,4)
        
        ! Loop to look the middle point.
        do k=k0,ny          
          ie1 = nint( (t_nmo(k-1)-t0)/dt )
          ie2 = nint( (t_nmo(k)-t0)/dt )
          ! Are we in the middle?
          if ( ie1*ie2 <= 0 ) then
            e1 = t_nmo(k-1)-t0
            e2 = t_nmo(k)-t0
            raux = -e1 + e2
            iv_nmo(it) = (1.+(e1/raux))*v_nmo(k-1) + (1.-(e2/raux))*v_nmo(k)
            k0 = k   ! Since it is monot. we start the new search in this last point.
            it_last = it   ! Last point that we interpolate.
            exit
          end if
        end do
        
      end do
      
      ! Repeat last point to end the interpolation (extrapolation: nearest-neighbor).
      do it=it_last+1,nt_Wd
        iv_nmo(it) = iv_nmo(it_last)
      end do
      
      !----------------------------------------------------------!    
      ! Inverse.
!       write(*,*) offset(1)
!       write(*,*) offset(n_rec)
!       write(*,*) n_rec
      do it=1,nt_Wd
!         write(*,*) iv_nmo(it)
        iv_nmo(it) = 1./iv_nmo(it)
      end do
      
      !----------------------------------------------------------!
      ! 
      cont = 0
      ! Point in the shot-gather.
      do it1=1+its0,nt_Wd
      do ir1=1,n_rec
        ! Looking for the point in another receiver with the same reflexion curve.
        ! Loop over zero-offset time where we start the curve that pass throw (ir1,it1). We solve the eq:
        ! t00 = sqrt( (t1-ts0)**2 - (offset(ir1)/vel(t0))**2 )
        e1_min = huge(1.)
        succes = .false.
        do it=1+its0,min(it1,nt_Wd)
          ! Error.
          e1 = abs( real(it,4) - sqrt( real((it1-1-its0)**2,4) - (offset(ir1)*iv_nmo(it)/dt)**2 ) )
          ! Min?
          if ( e1 < e1_min ) then
            it0 = it
            rit0 = sqrt( real((it1-1-its0)**2,4) - (offset(ir1)*iv_nmo(it)/dt)**2 )
            e1_min = e1
            succes = .true.
          end if
        end do
        
        ! Connect point in the curve trought Wd: (ir1,it1)--->(ir2,it2).
        if ( succes ) then
        do ir2=min(n_rec,ir1-nr),min(n_rec,ir1+nr)
          ! t2
          it2 = its0 + nint(sqrt( rit0**2 + (offset(ir2)*iv_nmo(it0)/dt)**2 ))
          ! Wd and index.
          if ( (it2<=nt_Wd) .and. (it2>0) ) then
            cont = cont + 1
          end if
        end do
        end if
        
      end do
      end do
      
      !----------------------------------------------------------!    
      ! Dimension Wd.
      dimWd = cont
      
      if ( ( rank == rank_wrt ) .and. ( isg == isg_wrt ) ) write(*,*) 'dimWd = ', dimWd
      
      ! Allocate.
      if ( .not. allocated( Wd ) ) then
        allocate( Wd(dimWd), stat=stal ); if ( stal/=0 ) stop 'AE get_cal_Wd 1'
        allocate( iW(dimWd), stat=stal ); if ( stal/=0 ) stop 'AE get_cal_Wd 2'
        allocate( jW(dimWd), stat=stal ); if ( stal/=0 ) stop 'AE get_cal_Wd 3'
      else
        deallocate( Wd, stat=stal ); if ( stal/=0 ) stop 'dAE get_cal_Wd 1'
        deallocate( iW, stat=stal ); if ( stal/=0 ) stop 'dAE get_cal_Wd 2'
        deallocate( jW, stat=stal ); if ( stal/=0 ) stop 'dAE get_cal_Wd 3'
        allocate( Wd(dimWd), stat=stal ); if ( stal/=0 ) stop 'AE get_cal_Wd 1b'
        allocate( iW(dimWd), stat=stal ); if ( stal/=0 ) stop 'AE get_cal_Wd 2b'
        allocate( jW(dimWd), stat=stal ); if ( stal/=0 ) stop 'AE get_cal_Wd 3b'
      end if
      
      iW = 0
      jW = 0
      Wd = 0.
      
      !----------------------------------------------------------!
      ! 
      cont = 0
      ! Point in the shot-gather.
!       do it1=1+its0,nt_Wd
      do it1=nt_Wd,1+its0,-1
      do ir1=1,n_rec
      if ( cont < dimWd ) then
        ! Looking for the point in another receiver with the same reflexion curve.
        ! Loop over zero-offset time where we start the curve that pass throw (ir1,it1). We solve the eq:
        ! t00 = sqrt( (t1-ts0)**2 - (offset(ir1)/vel(t0))**2 )
        e1_min = huge(1.)
        succes = .false.
        do it=1+its0,min(it1,nt_Wd)
          ! Error.
          e1 = abs( real(it,4) - sqrt( real((it1-1-its0)**2,4) - (offset(ir1)*iv_nmo(it)/dt)**2 ) )
          ! Min?
          if ( e1 < e1_min ) then
            it0 = it
            rit0 = sqrt( real((it1-1-its0)**2,4) - (offset(ir1)*iv_nmo(it)/dt)**2 )
            e1_min = e1
            succes = .true.
          end if
        end do
        
        ! Connect point in the curve trought Wd: (ir1,it1)--->(ir2,it2).
        if ( succes ) then
        do ir2=min(n_rec,ir1-nr),min(n_rec,ir1+nr)
          ! t2
          it2 = its0 + nint(sqrt( rit0**2 + (offset(ir2)*iv_nmo(it0)/dt)**2 ))
          ! Wd and index.
          if ( (it2<=nt_Wd) .and. (it2>0) ) then
            cont = cont + 1
            iW(cont) = ir1+n_rec_aux*(it1-1)
            jW(cont) = ir2+n_rec_aux*(it2-1)
            Wd(cont) = 1._sp
!             write(*,*) 'p1<p2 = ', it0, it1, it2, ir1, ir2
          end if
        end do
        end if
        
      end if
      end do
      end do
      
      dimWd = cont
!       if ( cont /= dimWd ) stop 'ERROR - get_cal_Wd: cont /= dimWd'
      
      !----------------------------------------------------------!
      deallocate( iv_nmo, stat=stal ); if ( stal/=0 ) stop 'dAE get_cal_Wd 1'
      
    end subroutine get_cal_Wd   
    
    
    
    
    
!     !**********************************************************!
!     !**********************************************************!
!     subroutine comp_Wd( vp0 )
!       use m_data_kind
!       implicit none
!       ! The variables which are passed to the function.
!       character(200), intent(in) :: folder
!       ! The variables which are generated inside the function.
!       integer :: k
!       integer :: stal
!       character(200) :: path_file
!       
! !       !----------------------------------------------------------!
! !       dimWd = 0
! !       
! !         read(1,*) iW(k), jW(k), Wd(k)
! ! 
! !       
! !       
! !       !----------------------------------------------------------!
! !       if ( .not. allocated( Wd ) ) then
! !         allocate( Wd(dimWd), stat=stal ); if ( stal/=0 ) stop 'AE comp_Wd 1'
! !         allocate( iW(dimWd), stat=stal ); if ( stal/=0 ) stop 'AE comp_Wd 2'
! !         allocate( jW(dimWd), stat=stal ); if ( stal/=0 ) stop 'AE comp_Wd 3'
! !       end if
!       
!       !----------------------------------------------------------!
!       
!     end subroutine comp_Wd  
    
    
    
!     !**********************************************************!
!     !**********************************************************!
!     subroutine mult_inv_H_v( v )
!       use m_data_kind
!       use m_geo
!       implicit none
!       ! The variables which are passed to the function.
!       real(sp), allocatable, intent(inout) :: v(:,:,:)
!       ! The variables which are generated inside the function.
!       integer :: ixg, iyg, izg, jxg, jyg, ix_H, iy_H, jx_H, jy_H, i_H, j_H
!       integer :: stal
!       real(sp) :: dxH, dyH, ux, uy, uux, uuy
!       real(sp) :: i_H_x1, i_H_x2, i_H_y1, i_H_y2, i_H_x1y1, i_H_x2y1, i_H_x1y2, i_H_x2y2
!       real(sp) :: j_H_x1, j_H_x2, j_H_y1, j_H_y2, j_H_x1y1, j_H_x2y1, j_H_x1y2, j_H_x2y2 
!       real(sp) :: d_x0, d_x1, d_x2,  d_y0, d_y1, d_y2
!       real(sp) :: d_x0y0, d_x1y0, d_x2y0, d_x0y1, d_x1y1, d_x2y1, d_x0y2, d_x1y2, d_x2y2, d_sum
!       real(sp) :: inv_Ha_interp, inv_Hb_interp, inv_H_interp
!       real(sp), allocatable :: w(:,:,:)
!       
!       !----------------------------------------------------------!
!       ! 
!       allocate( w(nx,ny,nz), stat=stal ); if ( stal/=0 ) stop 'AE mult_inv_H_v 1'
!       
!       w = v
!       v = 0._sp
!       
!       !----------------------------------------------------------!
!       ! 
! !       if ( dx/=16._sp/1000._sp ) stop 'ERROR: mult_inv_H_v'
!       
!       dxH = dx*20._sp
!       dyH = dx*4._sp
!       ux = dx/dxH
!       uy = dx/dyH
!       uux = 1._sp/ux
!       uuy = 1._sp/uy
!       
!       do izg=1,nz
!         
!         do ixg=1,nx
!         ix_H = 1+floor(dble(ixg-1)*ux)
!         do iyg=1,ny
!         iy_H = 1+floor(dble(iyg-1)*uy)
!         i_H = ix_H+nxH*(iy_H-1)
!           
!           do jxg=1,nx
!           jx_H = 1+floor(dble(jxg-1)*ux)
!           do jyg=1,ny
!           jy_H = 1+floor(dble(jyg-1)*uy)
!           j_H = jx_H+nxH*(jy_H-1)
!           
! !             v(ixg,iyg,izg) = v(ixg,iyg,izg) + inv_H(i_H,j_H)*w(jxg,jyg,izg)
! 
! !           H(p)_ii = d²X/dp² = dX/dp*dp²/dm² + d²X/dp²*dm/dp = grad(m)*dm/dp*dp²/dm² + H(m)*dm²/dp²
! !           H^p_ij = ...
! !           p = (cmpr/prec_mod)^(1-n)
! !           dp²/dm² =  d/dm abs(1-n) (cmpr/prec_mod)^(-n) = abs(n*(-1+n))*(cmpr/prec_mod)^(-n-1) = 
! !                   = abs(n*(-1+n))*p^(-n-1)/(1-n) = abs(n*(-1+n))*p^((n+1)/(n-1))
!             
!             i_H_x1 = ix_H+1+nxH*(iy_H-1)
!             i_H_x2 = ix_H+2+nxH*(iy_H-1)
!             j_H_x1 = jx_H+1+nxH*(jy_H-1)
!             j_H_x2 = jx_H+2+nxH*(jy_H-1)
!             
!             i_H_y1 = ix_H+nxH*(iy_H+1-1)
!             i_H_y2 = ix_H+nxH*(iy_H+2-1)
!             j_H_y1 = jx_H+nxH*(jy_H+1-1)
!             j_H_y2 = jx_H+nxH*(jy_H+2-1)
!             
!             i_H_x1y1 = ix_H+1+nxH*(iy_H+1-1)
!             i_H_x2y1 = ix_H+2+nxH*(iy_H+1-1)
!             j_H_x1y1 = jx_H+1+nxH*(jy_H+1-1)
!             j_H_x2y1 = jx_H+2+nxH*(jy_H+1-1)
!             
!             i_H_x2y1 = ix_H+2+nxH*(iy_H+1-1)
!             i_H_x2y2 = ix_H+2+nxH*(iy_H+2-1)
!             j_H_x2y1 = jx_H+2+nxH*(jy_H+1-1)
!             j_H_x2y2 = jx_H+2+nxH*(jy_H+2-1)
!             
!             ! i
!             d_x0 = ( uux*dble(ix_H  -1) - dble(ixg-1) )
!             d_x1 = ( uux*dble(ix_H+1-1) - dble(ixg-1) )
!             d_x2 = ( uux*dble(ix_H+2-1) - dble(ixg-1) )
!             
!             d_y0 = ( uuy*dble(iy_H  -1) - dble(iyg-1) )
!             d_y1 = ( uuy*dble(iy_H+1-1) - dble(iyg-1) )
!             d_y2 = ( uuy*dble(iy_H+2-1) - dble(iyg-1) )
!             
!             d_x0y0 = sqrt( d_x0**2 + d_y0**2 )
!             d_x1y0 = sqrt( d_x1**2 + d_y0**2 )
!             d_x2y0 = sqrt( d_x2**2 + d_y0**2 )
!             
!             d_x0y1 = sqrt( d_x0**2 + d_y1**2 )
!             d_x1y1 = sqrt( d_x1**2 + d_y1**2 )
!             d_x2y1 = sqrt( d_x2**2 + d_y1**2 )
!             
!             d_x0y2 = sqrt( d_x0**2 + d_y2**2 )
!             d_x1y2 = sqrt( d_x1**2 + d_y2**2 )
!             d_x2y2 = sqrt( d_x2**2 + d_y2**2 )
!             
!             d_sum = 0.5_sp/(d_x0y0 + d_x1y0 + d_x2y0 + d_x0y1 + d_x1y1 + d_x2y1 + d_x0y2 + d_x1y2 + d_x2y2)
!             
!             inv_Ha_interp = d_sum*( &
!                          d_x0y0*inv_H(i_H   ,j_H) + d_x1y0*inv_H(i_H_x1  ,j_H) + d_x2y0*inv_H(i_H_x2  ,j_H) + &
!                          d_x0y1*inv_H(i_H_y1,j_H) + d_x1y1*inv_H(i_H_x1y1,j_H) + d_x2y1*inv_H(i_H_x2y1,j_H) + &
!                          d_x0y2*inv_H(i_H_y2,j_H) + d_x1y2*inv_H(i_H_x1y2,j_H) + d_x2y2*inv_H(i_H_x2y2,j_H)   &
!                          )
!             
!             ! j
!             d_x0 = ( uux*dble(jx_H  -1) - dble(jxg-1) )
!             d_x1 = ( uux*dble(jx_H+1-1) - dble(jxg-1) )
!             d_x2 = ( uux*dble(jx_H+2-1) - dble(jxg-1) )
!             
!             d_y0 = ( uuy*dble(jy_H  -1) - dble(jyg-1) )
!             d_y1 = ( uuy*dble(jy_H+1-1) - dble(jyg-1) )
!             d_y2 = ( uuy*dble(jy_H+2-1) - dble(jyg-1) )
!             
!             d_x0y0 = sqrt( d_x0**2 + d_y0**2 )
!             d_x1y0 = sqrt( d_x1**2 + d_y0**2 )
!             d_x2y0 = sqrt( d_x2**2 + d_y0**2 )
!             
!             d_x0y1 = sqrt( d_x0**2 + d_y1**2 )
!             d_x1y1 = sqrt( d_x1**2 + d_y1**2 )
!             d_x2y1 = sqrt( d_x2**2 + d_y1**2 )
!             
!             d_x0y2 = sqrt( d_x0**2 + d_y2**2 )
!             d_x1y2 = sqrt( d_x1**2 + d_y2**2 )
!             d_x2y2 = sqrt( d_x2**2 + d_y2**2 )
!             
!             d_sum = 0.5_sp/(d_x0y0 + d_x1y0 + d_x2y0 + d_x0y1 + d_x1y1 + d_x2y1 + d_x0y2 + d_x1y2 + d_x2y2)
!             
!             inv_Hb_interp = d_sum*( &
!                          d_x0y0*inv_H(i_H,j_H   ) + d_x1y0*inv_H(i_H,j_H_x1  ) + d_x2y0*inv_H(i_H,j_H_x2  ) + &
!                          d_x0y1*inv_H(i_H,j_H_y1) + d_x1y1*inv_H(i_H,j_H_x1y1) + d_x2y1*inv_H(i_H,j_H_x2y1) + &
!                          d_x0y2*inv_H(i_H,j_H_y2) + d_x1y2*inv_H(i_H,j_H_x1y2) + d_x2y2*inv_H(i_H,j_H_x2y2)   &
!                          )
!             
!             ! 
!             inv_H_interp = inv_Ha_interp + inv_Hb_interp
!             
! !             v(ixg,iyg,izg) = v(ixg,iyg,izg) + inv_H(i_H,j_H)*w(jxg,jyg,izg)
!             v(ixg,iyg,izg) = v(ixg,iyg,izg) + inv_H_interp*w(jxg,jyg,izg)
!             
!           end do
!           end do
!           
!         end do
!         end do
!         
!       end do
!       
!       !----------------------------------------------------------!
!       ! 
!       deallocate( w, stat=stal ); if ( stal/=0 ) stop 'dAE mult_inv_H_v 1'
!       
!       !----------------------------------------------------------!    
!       
!     end subroutine mult_inv_H_v
    
    
    
    !**********************************************************!
    !**********************************************************!
    subroutine mult_inv_H_v_diag( iLx, iLy, iLz, v )
      use m_data_kind
      use m_geo
      implicit none
      ! The variables which are passed to the function.
      integer, intent(in) :: iLx, iLy, iLz
      real(sp), allocatable, intent(inout) :: v(:,:,:)
      ! The variables which are generated inside the function.
      integer :: ix, iy, iz, i
      
      !----------------------------------------------------------!
      !       
      do iz=1,iLz
      do iy=1,iLy
      do ix=1,iLx
        i = ix + iLx*(iy-1)
        v(ix,iy,iz) = inv_H_diag(i)*v(ix,iy,iz)
      end do
      end do
      end do
      
      !----------------------------------------------------------!    
      
    end subroutine mult_inv_H_v_diag
    
    
!     !**********************************************************!
!     !**********************************************************!
!     subroutine mult_inv_H_v( v )
!       use m_data_kind
!       use m_geo
!       implicit none
!       ! The variables which are passed to the function.
!       real(sp), allocatable, intent(inout) :: v(:,:,:)
!       ! The variables which are generated inside the function.
!       integer :: ixg, iyg, izg, jxg, jyg, ix_H, iy_H, jx_H, jy_H, i_H, j_H
!       integer :: stal
!       real(sp) :: dxH, dyH, ux, uy
!       integer :: i_H_x1, i_H_y1, i_H_x1y1
!       integer :: j_H_x1, j_H_y1, j_H_x1y1
!       real(sp) :: d_x0, d_x1,  d_y0, d_y1
!       real(sp) :: d_x0y0, d_x1y0, d_x0y1, d_x1y1
!       real(sp) :: inv_Ha_interp, inv_Hb_interp, inv_H_interp
!       real(sp), allocatable :: w(:,:,:)
!       
!       !----------------------------------------------------------!
!       ! 
!       allocate( w(nx,ny,nz), stat=stal ); if ( stal/=0 ) stop 'AE mult_inv_H_v 1'
!       
!       w = v
!       v = 0._sp
!       
!       !----------------------------------------------------------!
!       !       
!       dxH = dx*20._sp
!       dyH = dx*3._sp
!       ux = dx/dxH
!       uy = dx/dyH
!       
!       do izg=1,nz
!         
!         do ixg=1,nx
!         ix_H = 1+floor(dble(ixg-1)*ux)
!         do iyg=1,ny
!         iy_H = 1+floor(dble(iyg-1)*uy)
!         i_H = ix_H+nxH*(iy_H-1)
!           
!           do jxg=1,nx
!           jx_H = 1+floor(dble(jxg-1)*ux)
!           do jyg=1,ny
!           jy_H = 1+floor(dble(jyg-1)*uy)
!           j_H = jx_H+nxH*(jy_H-1)
!             
!             i_H_x1 = ix_H+1+nxH*(iy_H-1)
!             j_H_x1 = jx_H+1+nxH*(jy_H-1)
!             
!             i_H_y1 = ix_H+nxH*(iy_H+1-1)
!             j_H_y1 = jx_H+nxH*(jy_H+1-1)
!             
!             i_H_x1y1 = ix_H+1+nxH*(iy_H+1-1)
!             j_H_x1y1 = jx_H+1+nxH*(jy_H+1-1)
!             
!             ! i
!             d_x0 = ( dble(ix_H  -1) - ux*dble(ixg-1) )
!             d_x1 = ( dble(ix_H+1-1) - ux*dble(ixg-1) )
!             d_y0 = ( dble(iy_H  -1) - uy*dble(iyg-1) )
!             d_y1 = ( dble(iy_H+1-1) - uy*dble(iyg-1) )
!             
!             d_x0y0 = (1._sp-d_x0)*(1._sp-d_y0)
!             d_x1y0 = d_x1*(1._sp-d_y0)
!             d_x0y1 = (1._sp-d_x0)*d_y1
!             d_x1y1 = d_x1*d_y1
!             
!             inv_Ha_interp = ( &
!                          d_x0y0*inv_H(i_H   ,j_H) + d_x1y0*inv_H(i_H_x1  ,j_H) + &
!                          d_x0y1*inv_H(i_H_y1,j_H) + d_x1y1*inv_H(i_H_x1y1,j_H)   &
!                          )
!             
!             ! j
!             d_x0 = ( dble(jx_H  -1) - ux*dble(jxg-1) )
!             d_x1 = ( dble(jx_H+1-1) - ux*dble(jxg-1) )
!             d_y0 = ( dble(jy_H  -1) - uy*dble(jyg-1) )
!             d_y1 = ( dble(jy_H+1-1) - uy*dble(jyg-1) )
!             
!             d_x0y0 = (1._sp-d_x0)*(1._sp-d_y0)
!             d_x1y0 = d_x1*(1._sp-d_y0)
!             d_x0y1 = (1._sp-d_x0)*d_y1
!             d_x1y1 = d_x1*d_y1
!             
!             inv_Hb_interp = ( &
!                          d_x0y0*inv_H(i_H,j_H   ) + d_x1y0*inv_H(i_H,j_H_x1  ) + &
!                          d_x0y1*inv_H(i_H,j_H_y1) + d_x1y1*inv_H(i_H,j_H_x1y1)   &
!                          )
!             
!             ! 
!             inv_H_interp = inv_Ha_interp + inv_Hb_interp
!             
! !             v(ixg,iyg,izg) = v(ixg,iyg,izg) + inv_H(i_H,j_H)*w(jxg,jyg,izg)
!             v(ixg,iyg,izg) = v(ixg,iyg,izg) + inv_H_interp*w(jxg,jyg,izg)
!             
!           end do
!           end do
!           
!         end do
!         end do
!         
!       end do
!       
!       !----------------------------------------------------------!
!       ! 
!       deallocate( w, stat=stal ); if ( stal/=0 ) stop 'dAE mult_inv_H_v 1'
!       
!       !----------------------------------------------------------!    
!       
!     end subroutine mult_inv_H_v
    
    
    
  !**********************************************************!
  end module m_Hessian








