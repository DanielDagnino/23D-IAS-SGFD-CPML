  !*********************************************************************/
  !** This code has been done in the Barcelona Center for Subsurface 
  !** Imaging (BCSI).
  !** Authors: Daniel Dagnino.
  !*********************************************************************/
  !** Short resume of this code/file:
  !** 
  !*********************************************************************/
  module m_Hessian_red
    use m_data_kind
    implicit none
    real(sp), allocatable :: base(:,:)
    integer :: dim_H, nt_J
    
    !*********************************************************************/
    contains
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine get_J_red( folder, shgat, nt_J )
      use m_sg_type
      use m_mpi
      implicit none
      ! The variables which are passed to the function.
      integer, intent(out) :: nt_J
      character(200), intent(in) :: folder
      type(shot_gather), intent(inout) :: shgat
      ! The variables which are generated inside the function.
      integer :: imp, it, ir, iaux
      integer :: stal
      character(200) :: path_file
      character(50)  :: str_tmp, str_tmp2
      
      !----------------------------------------------------------!
      ! Read Jacobian.
      do imp=1,dim_H
        
        ! Open file.
        write(str_tmp,*)  shgat%id
        write(str_tmp2,*)  imp
        path_file = trim(adjustl( folder )) // '/Jacobian/tr_imp' // trim(adjustl(str_tmp2)) // '_' // trim(adjustl(str_tmp)) // '.txt'
        open(unit=1,file=path_file,status='old',action='read')
        
        read(1,*) iaux, nt_J
        
        if ( imp == 1 ) then
!           if ( .not. allocated( shgat%J ) ) then
            allocate( shgat%J(dim_H,nt_J,shgat%n_rec), stat=stal ); if ( stal/=0 ) stop 'AE get_J_red 1'
!           end if
          shgat%J = 0._sp
        end if
        
        do it=1,nt_J
          read(1,*) ( shgat%J(imp,it,ir), ir=1,shgat%n_rec )
        end do
        
        ! Close file.
        close(1)
        
      end do
      
    end subroutine get_J_red
    
    
    
    !**********************************************************!
    !**********************************************************!
    subroutine get_iH_red( folder, shgat )
      use mod_mat_inv
      use m_sg
      use m_sg_type
      use m_sg_pr
      use m_mpi
      use m_Hessian
      use m_type_inversion, only: dyw_Wd
      implicit none
      ! The variables which are passed to the function.
      character(200), intent(in) :: folder
      type(shot_gather), intent(inout) :: shgat
      ! The variables which are generated inside the function.
      integer :: imp, jmp, nt_end, it, ir
      integer :: stal
!       character(200) :: path_file
      real(sp) :: damp
!       real(sp) :: error_inv, mat_ij
!       real(sp) :: mat(dim_H,dim_H)
!       character(50)  :: str_tmp
!       real(sp), allocatable :: J_aux(:,:), J_aux2(:,:)
      integer :: nt_Wd
      
      nt_Wd = 1021
      
!       !----------------------------------------------------------!    
!       ! Open file.
! !       path_file = trim( folder ) // '../z-new_version_H/data_output/hessians'
!       path_file = '../z-new_version_H/data_output/hessians/'
!       write(str_tmp,*)  shgat%id
!       path_file = trim(adjustl( path_file )) // 'H_' // trim(adjustl(str_tmp)) // '.txt'
!       
!       open(unit=1,file=path_file,status='old',action='read')
      
      !----------------------------------------------------------!
      if ( .not. allocated( shgat%iH ) ) then
        allocate( shgat%iH(dim_H,dim_H), stat=stal ); if ( stal/=0 ) stop 'AE get_iH 1'
      end if
      
      shgat%iH = 0._sp
      
!       !----------------------------------------------------------!
!       ! 
!       do imp=1,dim_H
!         read(1,*) ( shgat%iH(imp,jmp), jmp=1,dim_H )
!       end do
      
      !----------------------------------------------------------!
      nt_end = shgat%nt_sou
      
      call adapt_prec( shgat, nt_J )
      
!       if ( dyw_Wd ) then
!         allocate( J_aux(shgat%n_rec,nt_Wd), stat=stal ); if ( stal/=0 ) stop 'AE get_iH_red'
!         allocate( J_aux2(shgat%n_rec,nt_Wd), stat=stal ); if ( stal/=0 ) stop 'AE get_iH_red'
!         J_aux = 0.
!         J_aux2 = 0.
!       end if
      
      do imp=1,dim_H
      
!         if ( dyw_Wd ) then
!           J_aux2 = 0.
!           do k=1,dimWd
!             it = 1+floor(dble(iW(k)-1)/dble(192))
!             ir = iW(k)-192*(it-1)
!             jt = 1+floor(dble(jW(k)-1)/dble(192))
!             jr = jW(k)-192*(jt-1)
!             if ( ( ir <= shgat%n_rec ) .and. ( jr <= shgat%n_rec ) ) then
!               J_aux2(it,ir) = J_aux2(it,ir) + sqrt(shgat%pr_d(jr,jt)) * Wd(k) * shgat%J(imp,jt,jr)
!             end if
!           end do
!         end if
      
      do jmp=1,imp
        
!         if ( dyw_Wd ) then
!           
!           J_aux = 0.
!           do k=1,dimWd
!             it = 1+floor(dble(iW(k)-1)/dble(192))
!             ir = iW(k)-192*(it-1)
!             jt = 1+floor(dble(jW(k)-1)/dble(192))
!             jr = jW(k)-192*(jt-1)
!             if ( ( ir <= shgat%n_rec ) .and. ( jr <= shgat%n_rec ) ) then
!               J_aux(it,ir) = J_aux(it,ir) + sqrt(shgat%pr_d(jr,jt)) * Wd(k) * shgat%J(jmp,jt,jr)
!             end if
!           end do
!           
!           do ir=1,shgat%n_rec
!           do it=1,nt_J
!             shgat%iH(imp,jmp) = shgat%iH(imp,jmp) + J_aux2(it,ir)*J_aux(it,ir)
!           end do
!           end do
!           
!         else
        
          do ir=1,shgat%n_rec
          do it=1,nt_J
            shgat%iH(imp,jmp) = shgat%iH(imp,jmp) + shgat%J(imp,it,ir)*shgat%J(jmp,it,ir)
          end do
          end do
          
!         end if
        
        shgat%iH(jmp,imp) = shgat%iH(imp,jmp)
        
      end do
      end do
      
      call adapt_prec( shgat, nt_end )
      
      !----------------------------------------------------------!
      !----------------------------------------------------------!
      !----------------------------------------------------------!
      ! 
      damp = 1.e-3*maxval(abs(shgat%iH))
!       damp = 1.e-1*maxval(abs(shgat%iH))
      do imp=1,dim_H
        shgat%iH(imp,imp) = shgat%iH(imp,imp) + damp
      end do
      
!       write(*,*) 'maxval(shgat%iH), rank a = ', maxval(shgat%iH), rank    
      
      !----------------------------------------------------------!
      ! 
!       mat = shgat%iH
      call mat_inv( dim_H, shgat%iH )
      
!       !----------------------------------------------------------!
!       ! 
!       error_inv = 0.
!       do imp=1,dim_H
!       do jmp=jmp+1,dim_H
!         mat_ij = 0.
!         do kmp=1,dim_H
!           mat_ij = mat_ij + mat(imp,kmp)*shgat%iH(kmp,jmp)
!         end do
!         error_inv = max( error_inv, abs(mat_ij) )
!       end do
!       end do
!       
!       do imp=1,dim_H
!         mat_ij = 0.
!         do kmp=1,dim_H
!           mat_ij = mat_ij + mat(imp,kmp)*shgat%iH(kmp,imp)
!         end do
!         error_inv = max( error_inv, abs(mat_ij-1.) )
!       end do
!       
!       write(*,*) 'error_inv = ', error_inv, rank   

      !----------------------------------------------------------!
      ! Error check.
!       write(*,*) 'maxval(shgat%iH), rank b = ', maxval(shgat%iH), rank       
      
!       !----------------------------------------------------------!
!       !----------------------------------------------------------!
!       !----------------------------------------------------------!
!       close(1)

      !----------------------------------------------------------!
!       if ( dyw_Wd ) then
!         deallocate( J_aux, stat=stal ); if ( stal/=0 ) stop 'dAE get_iH_red'
!         deallocate( J_aux2, stat=stal ); if ( stal/=0 ) stop 'dAE get_iH_red'
!       end if
      
    end subroutine get_iH_red
    
    
    
    !**********************************************************!
    !**********************************************************!
    subroutine get_base( folder, shgat )
      use m_data_kind
      use m_sg_type
      use m_geo, only: ny, dx
      use m_mpi
      use m_splines
      implicit none
      ! The variables which are passed to the function.
      character(200), intent(in) :: folder
      type(shot_gather), intent(inout) :: shgat
      ! The variables which are generated inside the function.
      integer :: iy, imp
!       real(sp) :: k_num
      real(sp) :: y, dbase
      integer :: stal
      integer :: stype(2)
      real(8) :: svalue(2)
      real(rl), allocatable :: yi(:), basei(:)
      real(rl), allocatable :: c(:,:)
      
      !----------------------------------------------------------!
      if ( .not. allocated( base ) ) then
        allocate( base(dim_H,ny), stat=stal ); if ( stal/=0 ) stop 'AE get_base 1'
      end if
      
      !----------------------------------------------------------!
      allocate( yi(dim_H+4), stat=stal ); if ( stal/=0 ) stop 'AE get_base 1'
      allocate( basei(dim_H+4), stat=stal ); if ( stal/=0 ) stop 'AE get_base 2'
      allocate( c(0:4,(dim_H+4)-1), stat=stal ); if ( stal/=0 ) stop 'AE get_base 3'
      
      !----------------------------------------------------------!
      do imp=1,dim_H
        
!         k_num = -floor( dble(imp)/2. )
!         
!         if ( mod(imp,2) == 0 ) then
!           
!           do iy=1,ny
!             base(imp,iy) = sin( 2.*3.14159*k_num*(dble(iy-1)/dble(ny-1)) )
!           end do
!           
!         else
!           
!           do iy=1,ny
!             base(imp,iy) = cos( 2.*3.14159*k_num*(dble(iy-1)/dble(ny-1)) )
!           end do
!           
!         end if





!         !----------------------------------------------------------!
!         if ( imp == 1 ) then
!           y0 = 350./2000.
!         else if ( imp == 2 ) then
!           y0 = 700./2000.
!         else if ( imp == 3 ) then
!           y0 = 1050./2000.
!         else if ( imp == 4 ) then
!           y0 = 1400./2000.
!         else if ( imp == 5 ) then
!           y0 = 1600./2000.
!         end if
! 
!         do iy=1,ny
!           y = dble(iy-1)/dble(ny-1)
!           base(imp,iy) = exp( -((y-y0)/0.065)**2 )
!         end do
        
        
        
        
        
        
        !----------------------------------------------------------!
        ! 
        basei = 0.
        if ( imp == 1 ) then
          basei(imp+2) = 1.
        else if ( imp == 2 ) then
          basei(imp+2) = 1.
        else if ( imp == 3 ) then
          basei(imp+2) = 1.
        else if ( imp == 4 ) then
          basei(imp+2) = 1.
        else if ( imp == 5 ) then
          basei(imp+2) = 1.
          basei(imp+3) = 1.
          basei(imp+4) = 1.
        end if
        
        yi(1) =    0.
        yi(2) =   50./1000.
        yi(3) =  300./1000.   ! 1st element of the base.
        yi(4) =  600./1000.   ! 2on element of the base.
        yi(5) =  900./1000.   ! 3th element of the base.
        yi(6) = 1200./1000.   ! 4th element of the base.
        yi(7) = 1500./1000.   ! 5th element of the base.
        yi(8) = 1800./1000.
        yi(9) = 2000./1000.
                
        stype(1) = 1
        stype(2) = 1
        svalue(1) = 0.
        svalue(2) = 0.
        call spline3pars( yi, basei, stype, svalue, c )
        
        do iy=1,ny
          y = dx*dble(iy-1)
          call spline3valder( y, yi, c, base(imp,iy), dbase )
        end do
        
        base(imp,:) = base(imp,:)/maxval(abs(base(imp,:))) 
        
        
        
      end do
      
      deallocate( c, stat=stal ); if ( stal/=0 ) stop 'dAE get_base 2'
      deallocate( yi, stat=stal ); if ( stal/=0 ) stop 'dAE get_base 3'
      deallocate( basei, stat=stal ); if ( stal/=0 ) stop 'dAE get_base 4'
      
    end subroutine get_base
    
    
    
    !**********************************************************!
    !**********************************************************!
    subroutine mult_inv_H_v_red( iLx, iLy, iLz, grad_red, shgat )
      use m_data_kind
      use m_sg_type
!       use m_geo
      implicit none
      ! The variables which are passed to the function.
      integer, intent(in) :: iLx, iLy, iLz
!       real(sp), allocatable, intent(inout) :: grad_mg(:,:,:)
      real(sp), allocatable, intent(inout) :: grad_red(:)
      type(shot_gather), intent(in) :: shgat
      ! The variables which are generated inside the function.
      integer :: imp, jmp
      real(sp) :: grad_aux(dim_H)
      
      !----------------------------------------------------------!
      ! 
      grad_aux = 0.
      do imp=1,dim_H
      do jmp=1,dim_H
        grad_aux(imp) = grad_aux(imp) + shgat%iH(imp,jmp)*grad_red(jmp)
      end do
      end do
      
      grad_red = grad_aux
      
    end subroutine mult_inv_H_v_red
    
    
    
    
    
  !**********************************************************!
  end module m_Hessian_red








