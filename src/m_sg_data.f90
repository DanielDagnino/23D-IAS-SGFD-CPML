  !*********************************************************************/
  !** This code has been done in the Barcelona Center for Subsurface 
  !** Imaging (BCSI).
  !** Authors: Daniel Dagnino.
  !*********************************************************************/
  !** Short resume of this code/file:
  !** Definition of input/output folder/files + check existence.
  !*********************************************************************/
  module m_sg_data
    implicit none
    
    ! Inputs.
    character(200) :: folder_strat     ! < Basic + Freq. strategy + Wd + H.
    character(200) :: folder_model     ! < Initial model.
    character(200) :: folder_sg        ! < Traces.
    character(200) :: folder_pos_wv    ! < Source/Receiver positons + Source wavelet.
    ! Outputs.
    character(200) :: folder_output    ! > Common outputs.
    character(200) :: folder_tmp       ! > Temporal fields.
    character(200) :: folder_test      ! > In testing mode, field at various time steps.
    
    !*********************************************************************/
    contains
    
    
    
    !*********************************************************************/
    !*********************************************************************/
    !*********************************************************************/
    subroutine check_exit_folders_files()
      implicit none
#ifdef __INTEL_COMPILER
      include 'for_iosdef.for'
#endif
      integer :: ierr
      logical :: dir_exists
      character(100) :: path_file
      character(10) :: path_extra
      
      !----------------------------------------------------------!
      ! Check main input files:
      path_file = 'data_input/initial_data.txt'
      open(unit=1,file=path_file,status='old',action='read',iostat=ierr,err=100)
      close(1)
      
      path_file = trim(adjustl( folder_strat )) // '/basic.txt'
      open(unit=1,file=path_file,status='old',action='read',iostat=ierr,err=100)
      close(1)
      
      path_file = trim(adjustl( folder_strat )) // '/freq_inv.txt'
      open(unit=1,file=path_file,status='old',action='read',iostat=ierr,err=100)
      close(1)
      
      path_file = trim(adjustl( folder_strat )) // '/strategy.txt'
      open(unit=1,file=path_file,status='old',action='read',iostat=ierr,err=100)
      close(1)
      
      ! Check main input files for the model and geometry.
      path_file = trim(adjustl( folder_model )) // '/bathymetry.txt'
      open(unit=1,file=path_file,status='old',action='read',iostat=ierr,err=100)
      close(1)
      
      path_file = trim(adjustl( folder_model )) // '/dnst'
      open(unit=1,file=path_file,status='old',action='read',iostat=ierr,err=100)
      close(1)
      
      path_file = trim(adjustl( folder_model )) // '/geometry.txt'
      open(unit=1,file=path_file,status='old',action='read',iostat=ierr,err=100)
      close(1)
      
      path_file = trim(adjustl( folder_model )) // '/vP'
      open(unit=1,file=path_file,status='old',action='read',iostat=ierr,err=100)
      close(1)
      
      ! Error file output.
#ifdef __INTEL_COMPILER
      if ( ierr .ne. FOR$IOS_SUCCESS ) then
100    if ( ierr .eq. FOR$IOS_FILNOTFOU ) then
          write(*,*) ' ***** ERROR - check_exit_folders_files: File ', path_file, ' does not exist  ***** '
          stop
        else if ( ierr .eq. FOR$IOS_FILNAMSPE ) then
          write(*,*) ' ***** ERROR - check_exit_folders_files: File ', path_file, ' was bad, enter new file name ***** '
          stop
        else
          write(*,*) ' ***** ERROR - check_exit_folders_files: Unrecoverable error, code = ', ierr, ' ***** '
          stop
        end if
      end if
#else
100   if ( ierr .ne. 0 ) then
        write(*,*) ' ***** ERROR - check_exit_folders_files: File ', path_file, ' ***** '
        write(*,*) ' ***** ERROR - check_exit_folders_files: Error code: iostat says ', ierr, ' ***** '
        stop
      end if
#endif
      
      !----------------------------------------------------------!
      ! Check folders.
      
      ! Determine compiler to use the correct inquire form for folders.
#ifdef __INTEL_COMPILER
#define INQUIREFOLDER directory
#define COMPILEROK 1
      write(*,*) 'check_exit_folders_files: Check using IFORT '
      path_extra = ''
#elif __GFORTRAN__
#define INQUIREFOLDER file
#define COMPILEROK 1
      write(*,*) 'check_exit_folders_files: Check using GFORTRAN '
      path_extra = '/.'
#else
#define COMPILEROK 0
      write(*,*) 'check_exit_folders_files: Using an unknown compiler => NO CHECKING if files and folders exist!!! '
#endif
      
#if COMPILEROK == 1
      ! Check folders.
      path_file = trim(adjustl( trim(adjustl( folder_strat )) // path_extra ))
      inquire( INQUIREFOLDER=path_file, exist=dir_exists )
      if ( .not. dir_exists ) then
        write(*,*) ' ***** ERROR - check_exit_folders_files: Directory '
        write(*,*) path_file
        write(*,*) ' does not exist  ***** '
        stop
      end if
      
      path_file = trim(adjustl( trim(adjustl( folder_model )) // path_extra ))
      inquire( INQUIREFOLDER=path_file, exist=dir_exists )
      if ( .not. dir_exists ) then
        write(*,*) ' ***** ERROR - check_exit_folders_files: Directory '
        write(*,*) path_file
        write(*,*) ' does not exist  ***** '
        stop
      end if
      
      path_file = trim(adjustl( trim(adjustl( folder_sg )) // path_extra ))
      inquire( INQUIREFOLDER=path_file, exist=dir_exists )
      if ( .not. dir_exists ) then
        write(*,*) ' ***** ERROR - check_exit_folders_files: Directory '
        write(*,*) path_file
        write(*,*) ' does not exist  ***** '
        stop
      end if
      
      path_file = trim(adjustl( trim(adjustl( folder_pos_wv )) // path_extra ))
      inquire( INQUIREFOLDER=path_file, exist=dir_exists )
      if ( .not. dir_exists ) then
        write(*,*) ' ***** ERROR - check_exit_folders_files: Directory '
        write(*,*) path_file
        write(*,*) ' does not exist  ***** '
        stop
      end if
      
      path_file = trim(adjustl( trim(adjustl( folder_output )) // path_extra ))
      inquire( INQUIREFOLDER=path_file, exist=dir_exists )
      if ( .not. dir_exists ) then
        write(*,*) ' ***** ERROR - check_exit_folders_files: Directory '
        write(*,*) path_file
        write(*,*) ' does not exist  ***** '
        stop
      end if
      
      path_file = trim(adjustl( trim(adjustl( folder_tmp )) // path_extra ))
      inquire( INQUIREFOLDER=path_file, exist=dir_exists )
      if ( .not. dir_exists ) then
        write(*,*) ' ***** ERROR - check_exit_folders_files: Directory '
        write(*,*) path_file
        write(*,*) ' does not exist  ***** '
        stop
      end if
      
      path_file = trim(adjustl( trim(adjustl( folder_test )) // path_extra ))
      inquire( INQUIREFOLDER=path_file, exist=dir_exists )
      if ( .not. dir_exists ) then
        write(*,*) ' ***** ERROR - check_exit_folders_files: Directory '
        write(*,*) path_file
        write(*,*) ' does not exist  ***** '
        stop
      end if
      
      !----------------------------------------------------------!
      ! Check default folders in data_output.
      path_file = trim(adjustl( trim(adjustl( folder_output )) // '/adj_sou' // path_extra ))
      inquire( INQUIREFOLDER=path_file, exist=dir_exists )
      if ( .not. dir_exists ) then
        write(*,*) ' ***** ERROR - check_exit_folders_files: Directory '
        write(*,*) path_file
        write(*,*) ' does not exist  ***** '
        stop
      end if
      
      path_file = trim(adjustl( trim(adjustl( folder_output )) // '/calib' // path_extra ))
      inquire( INQUIREFOLDER=path_file, exist=dir_exists )
      if ( .not. dir_exists ) then
        write(*,*) ' ***** ERROR - check_exit_folders_files: Directory '
        write(*,*) path_file
        write(*,*) ' does not exist  ***** '
        stop
      end if
      
#ifdef testing
      path_file = trim(adjustl( trim(adjustl( folder_output )) // '/field' // path_extra ))
      inquire( INQUIREFOLDER=path_file, exist=dir_exists )
      if ( .not. dir_exists ) then
        write(*,*) ' ***** ERROR - check_exit_folders_files: Directory '
        write(*,*) path_file
        write(*,*) ' does not exist  ***** '
        stop
      end if
#endif
      
      path_file = trim(adjustl( trim(adjustl( folder_output )) // '/grad' // path_extra ))
      inquire( INQUIREFOLDER=path_file, exist=dir_exists )
      if ( .not. dir_exists ) then
        write(*,*) ' ***** ERROR - check_exit_folders_files: Directory '
        write(*,*) path_file
        write(*,*) ' does not exist  ***** '
        stop
      end if
      
      path_file = trim(adjustl( trim(adjustl( folder_output )) // '/illu' // path_extra ))
      inquire( INQUIREFOLDER=path_file, exist=dir_exists )
      if ( .not. dir_exists ) then
        write(*,*) ' ***** ERROR - check_exit_folders_files: Directory '
        write(*,*) path_file
        write(*,*) ' does not exist  ***** '
        stop
      end if
      
      path_file = trim(adjustl( trim(adjustl( folder_output )) // '/model' // path_extra ))
      inquire( INQUIREFOLDER=path_file, exist=dir_exists )
      if ( .not. dir_exists ) then
        write(*,*) ' ***** ERROR - check_exit_folders_files: Directory '
        write(*,*) path_file
        write(*,*) ' does not exist  ***** '
        stop
      end if
      
      path_file = trim(adjustl( trim(adjustl( folder_output )) // '/results' // path_extra ))
      inquire( INQUIREFOLDER=path_file, exist=dir_exists )
      if ( .not. dir_exists ) then
        write(*,*) ' ***** ERROR - check_exit_folders_files: Directory '
        write(*,*) path_file
        write(*,*) ' does not exist  ***** '
        stop
      end if
      
      path_file = trim(adjustl( trim(adjustl( folder_output )) // '/resume' // path_extra ))
      inquire( INQUIREFOLDER=path_file, exist=dir_exists )
      if ( .not. dir_exists ) then
        write(*,*) ' ***** ERROR - check_exit_folders_files: Directory '
        write(*,*) path_file
        write(*,*) ' does not exist  ***** '
        stop
      end if
      
      path_file = trim(adjustl( trim(adjustl( folder_output )) // '/source' // path_extra ))
      inquire( INQUIREFOLDER=path_file, exist=dir_exists )
      if ( .not. dir_exists ) then
        write(*,*) ' ***** ERROR - check_exit_folders_files: Directory '
        write(*,*) path_file
        write(*,*) ' does not exist  ***** '
        stop
      end if
      
      path_file = trim(adjustl( trim(adjustl( folder_output )) // '/tmp' // path_extra ))
      inquire( INQUIREFOLDER=path_file, exist=dir_exists )
      if ( .not. dir_exists ) then
        write(*,*) ' ***** ERROR - check_exit_folders_files: Directory '
        write(*,*) path_file
        write(*,*) ' does not exist  ***** '
        stop
      end if
      
      path_file = trim(adjustl( trim(adjustl( folder_output )) // '/trace' // path_extra ))
      inquire( INQUIREFOLDER=path_file, exist=dir_exists )
      if ( .not. dir_exists ) then
        write(*,*) ' ***** ERROR - check_exit_folders_files: Directory '
        write(*,*) path_file
        write(*,*) ' does not exist  ***** '
        stop
      end if
#endif
      
      !----------------------------------------------------------!
    end subroutine check_exit_folders_files
    
    !*********************************************************************/
  end module m_sg_data
  
  
  
  
  