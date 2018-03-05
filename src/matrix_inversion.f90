  !*********************************************************************/
  !** This code has been done in the Barcelona Center for Subsurface 
  !** Imaging (BCSI).
  !** Authors: Daniel Dagnino.
  !*********************************************************************/
  !** Short resume of this code/file:
  !** Matrix inversion using the LAPACK libraries.
  !*********************************************************************/
  module mod_mat_inv
    
    !*********************************************************************/
    contains
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine mat_inv( dim_H, H )
      use m_data_kind
      implicit none
      ! The variables which are passed to the function.
      integer, intent(in) :: dim_H
      real(rl), allocatable, intent(inout) :: H(:,:)
      ! The variables which are generated inside the function.
      real(rl), allocatable :: ap(:)
      integer :: i, j
      integer :: stal
      ! Variables for LAPACK
      real(rl), allocatable :: work(:), ipiv(:)
      integer :: info1, info2
!dir$ assume_aligned H(1,1):64
      
      !**********************************************************!
      ! 
      allocate( work(dim_H), ipiv(dim_H), stat=stal ); if ( stal/=0 ) stop 'AE mat_inv'
      allocate( ap( (dim_H*(dim_H+1))/2 ), stat=stal ); if ( stal/=0 ) stop 'AE mat_inv'
      
      !**********************************************************!
      ! 
      ! (N*(N+1)/2)
      ! AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
      do j=1,dim_H
        do i=1,j
          ap(i+((j-1)*j)/2) = H(i,j)
  !         write(*,'(f3.0,1x)',advance='no') H(i,j)
        end do
  !       write(*,*)
      end do
      
#ifndef valgrind
      !**********************************************************!
      ! Inversion.
      if ( rl==4 ) then
        
  !       call spotri( 'U', dim_H, H, dim_H, info )
        
        call ssptrf( 'U', dim_H, ap, ipiv, info1 )
        call ssptri( 'U', dim_H, ap, ipiv, work, info2 )
        
      else if ( rl==8 ) then
        
  !       call dpotri( 'U', dim_H, H, dim_H, info )

        call dsptrf( 'U', dim_H, ap, ipiv, info1 )
        call dsptri( 'U', dim_H, ap, ipiv, work, info2 )
        
      else
        
        stop ' ********** ERROR: Invalid value for rl ( = 4 or 8 ) ********** '
        
      end if
#endif
      
      ! Succesful inversion ?
      if ( info1/=0 .or. info2/=0 ) then
        write(*,*) ' ********** ERROR: matrix_inversion - LAPACK - info/=0 ********** '
        write(*,*) ' info1 = ', info1
        write(*,*) ' info2 = ', info2
        stop
      end if
      
      !**********************************************************!
      ! 
      do j=1,dim_H
      do i=1,j
        H(i,j) = ap(i+((j-1)*j)/2)
        H(j,i) = H(i,j)
      end do
      end do
      
      !**********************************************************!
      ! 
      deallocate( work, ipiv, stat=stal ); if ( stal/=0 ) stop 'dAE mat_inv'
      deallocate( ap, stat=stal ); if ( stal/=0 ) stop 'dAE mat_inv'
      
    end subroutine mat_inv
  
  !**********************************************************!
  end module mod_mat_inv




