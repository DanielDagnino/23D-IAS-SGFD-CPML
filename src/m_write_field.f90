  !*********************************************************************/
  !** This code has been done in the Barcelona Center for Subsurface 
  !** Imaging (BCSI).
  !** Authors: Daniel Dagnino.
  !*********************************************************************/
  !** Short resume of this code/file:
  !** Write boundary fields for the back-propagation field reconstraction.
  !*********************************************************************/
  module m_solv3Dahvc_write_field
    implicit none
    
    !*********************************************************************/
    contains
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine write_field( bd, ii, ind_shot, cont_save, folder_field, max_cont_save, mg )
      use m_solv3Dahvc_data_kind
      use m_sg_type
      use m_dyw, only: dyw_3d
      implicit none
#ifdef __INTEL_COMPILER
      include 'for_iosdef.for'
#endif
      
      ! The variables which are passed to the function.
      type(boundary), intent(in) :: bd(:)
      integer, intent(in)        :: ii, ind_shot, cont_save
      integer, intent(in)        :: max_cont_save
      character(200), intent(in) :: folder_field
      type(my_geo), intent(in) :: mg
      
      ! The variables which are generated inside the function.
      character(100) :: str_tmp1, str_tmp2, str_tmp0
      character(200) :: path_file
      integer :: k
      integer :: ierr
      
!       LOGICAL :: I_OPENED, L_EXISTS
!       character(50)   :: I_ACTION
!       CHARACTER (200) :: I_NAME, C_DIRSPEC
      
      !----------------------------------------------------------!
      ! File name.
      write(str_tmp0,*) ind_shot
      write(str_tmp1,*) ii
      write(str_tmp2,*) cont_save
      
      path_file = trim(adjustl( folder_field)) // '/shot_'
      path_file = trim( path_file) // adjustl(str_tmp0)
      path_file = trim( path_file) // '_k'
      path_file = trim( path_file) // adjustl(str_tmp1)
      path_file = trim( path_file) // '_step_'
      path_file = trim( path_file) // adjustl(str_tmp2)
      path_file = trim( path_file) // '.txt'
      
      !----------------------------------------------------------!
      ! Open.
      open(unit=1,file=path_file,status='replace',action='write',&
           form='unformatted',access='sequential',IOSTAT=ierr,ERR=100)
!       open(unit=1,file=path_file,status='unknown',action='write',form='unformatted')
!       open(unit=1,file=path_file,status='unknown',action='write',form='unformatted',BUFFERED="NO")
!       open(unit=1,file=path_file,status='unknown',action='write',ACCESS="STREAM")
      
      !----------------------------------------------------------!
      ! Write.
        write(1) max_cont_save
!         write(*,*) 'max_cont_save = ', max_cont_save
        
        do k=1,max_cont_save
          
          if ( dyw_3d ) then
            write(1) bd(k)%L%u(       -1:3       ,       -1:mg%iLy+2,       -1:mg%iLz+2 )
            write(1) bd(k)%R%u( mg%iLx-2:mg%iLx+2,       -1:mg%iLy+2,       -1:mg%iLz+2 )
            write(1) bd(k)%B%u(       -1:mg%iLx+2,       -1:3       ,       -1:mg%iLz+2 )
            write(1) bd(k)%F%u(       -1:mg%iLx+2, mg%iLy-2:mg%iLy+2,       -1:mg%iLz+2 )
            write(1) bd(k)%U%u(       -1:mg%iLx+2,       -1:mg%iLy+2,       -1:3        )
            write(1) bd(k)%D%u(       -1:mg%iLx+2,       -1:mg%iLy+2, mg%iLz-2:mg%iLz+2 )
          else
            write(1) bd(k)%L%u(       -1:3       ,       -1:mg%iLy+2,        1:mg%iLz   )
            write(1) bd(k)%R%u( mg%iLx-2:mg%iLx+2,       -1:mg%iLy+2,        1:mg%iLz   )
            write(1) bd(k)%B%u(       -1:mg%iLx+2,       -1:3       ,        1:mg%iLz   )
            write(1) bd(k)%F%u(       -1:mg%iLx+2, mg%iLy-2:mg%iLy+2,        1:mg%iLz   )
          end if
          
        end do
        
      !----------------------------------------------------------!
      ! Close.
      close(1)
      
      !----------------------------------------------------------!
      ! Error handling.
#ifdef __INTEL_COMPILER
      IF ( ierr .ne. FOR$IOS_SUCCESS ) THEN
100   IF ( ierr .eq. FOR$IOS_FILNOTFOU ) THEN
        write(*,*) ' ***** ERROR - write_field: File ', path_file, ' does not exist  ***** '
        
!         open(unit=1,file='./a/shot_1_k1_step_1.txt')
!         INQUIRE ( 1, OPENED=I_OPENED, NAME=I_NAME, ACTION=I_ACTION )
!         write(*,*) ' 1 - I_OPENED, I_NAME, I_ACTION = ', I_OPENED, I_NAME, I_ACTION
!         close(1)
!         
!         open(unit=1,file='./a/shot_1_k1_step_1.txt',status='replace',action='write',disp='keep')
!         INQUIRE ( 1, OPENED=I_OPENED, NAME=I_NAME, ACTION=I_ACTION )
!         write(*,*) ' 2 - I_OPENED, I_NAME, I_ACTION = ', I_OPENED, I_NAME, I_ACTION
!         close(1)
!         
!         open(unit=1,file='./a/shot_1_k1_step_1.txt',status='replace',action='write',form='unformatted',disp='keep')
!         INQUIRE ( 1, OPENED=I_OPENED, NAME=I_NAME, ACTION=I_ACTION )
!         write(*,*) ' 3 - I_OPENED, I_NAME, I_ACTION = ', I_OPENED, I_NAME, I_ACTION
!         close(1)
!         
!         open(unit=1,file='./a/shot_1_k1_step_1.txt',status='replace',action='write',form='unformatted',disp='keep',access='sequential')
!         INQUIRE ( 1, OPENED=I_OPENED, NAME=I_NAME, ACTION=I_ACTION )
!         write(*,*) ' 5 - I_OPENED, I_NAME, I_ACTION = ', I_OPENED, I_NAME, I_ACTION
!         close(1)
!         
!         open(unit=1,file=path_file,status='replace',action='write',form='unformatted',disp='keep',access='sequential')
!         INQUIRE ( 1, OPENED=I_OPENED, NAME=I_NAME, ACTION=I_ACTION )
!         write(*,*) ' 6 - I_OPENED, I_NAME, I_ACTION = ', I_OPENED, I_NAME, I_ACTION
!         close(1)
        
        stop
      ELSE IF ( ierr .eq. FOR$IOS_FILNAMSPE ) THEN
        write(*,*) ' ***** ERROR - write_field: File ', path_file, ' was bad, enter new file name ***** '
        stop
      ELSE
        write(*,*) ' ***** ERROR - write_field: Unrecoverable error, code = ', ierr, ' ***** '
        stop
      END IF
      END IF
#else
100   if ( ierr .ne. 0 ) then
        write(*,*) ' ***** ERROR - write_field: Error code: iostat says ', ierr, ' ***** '
        stop
      end if
#endif

      !----------------------------------------------------------!    
    end subroutine write_field
  
  !***********************************************************************/
  end module m_solv3Dahvc_write_field
  
  
  
  
  