  !*********************************************************************/
  !*****     Author: D. Dagnino     ************************************/
  !*********************************************************************/
  module m_alp_cota
    implicit none
    
    !**********************************************************!
    !**********************************************************!
    contains
    
    !**********************************************************!
    function cal_alp_cota( model, search, change_cota )
      use m_data_kind
      implicit none
      ! The variables which are passed to the function.
      real(sp)             :: cal_alp_cota
      real(sp), intent(in) :: change_cota
      real(sp), allocatable, intent(in)    :: model(:,:,:), search(:,:,:)
!dir$ assume_aligned model(1,1,1):64,search(1,1,1):64
      
      !**********************************************************!
      ! 
      cal_alp_cota = change_cota*minval(abs( model/search ))
      
    !**********************************************************!
    end function cal_alp_cota
    
  !**********************************************************!
  end module m_alp_cota




