  !*********************************************************************/
  !*****     Author: D. Dagnino     ************************************/
  !*********************************************************************/
  module m_geo
    use m_data_kind
    implicit none
    
    !----------------------------------------------------------!
    type model_grid
      ! Bathymetry.
      real(sp), allocatable :: bath(:,:)
      integer, allocatable :: ind_bath(:,:)
    end type model_grid
    
    !----------------------------------------------------------!
    integer :: nx, ny, nz
    real(sp) :: dx
    type(model_grid) :: grid
    
    !**********************************************************!
    !**********************************************************!
    !**********************************************************!
    contains
    
    !**********************************************************!
    subroutine get_geo_size( folder_model, nx_model, ny_model, nz_model, dx_model )      
      ! The variables which are passed to the function.
      integer, intent(out) :: nx_model, ny_model, nz_model
      real(sp), intent(out) :: dx_model
      character(200), intent(in) :: folder_model
      ! The variables which are generated inside the function.
      character(200) :: path_file, str_aux
      
      !----------------------------------------------------------!    
      ! Open file.
      path_file = trim(adjustl( folder_model )) // '/'
      path_file = trim( path_file ) // 'geometry.txt'
      
      open(unit=1,file=path_file,status='old',action='read')
        
        ! Read number of shot.
        read(1,'(a)') str_aux
        read(1,*) nx_model
        
        read(1,'(a)') str_aux
        read(1,*) ny_model
        
        read(1,'(a)') str_aux
        read(1,*) nz_model
        
        read(1,'(a)') str_aux
        read(1,*) dx_model
        
        dx_model = dx_model/1000._sp
        
      close(1)
      
    end subroutine get_geo_size
    
    !***********************************************************************/
  end module m_geo
  
  
  
  
  