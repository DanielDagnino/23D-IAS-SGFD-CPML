
  !**********************************************************!
  module m_mpi
    implicit none
    integer :: n_procs
    integer :: rank
    integer :: ierr
    integer :: mpi_fwi
    
   !**********************************************************!
   contains
   
   !**********************************************************!
#ifdef usempi
   
   subroutine prec_mpi_fwi( prec )
      use m_data_kind
      use mpi
      implicit none
      ! The variables which are passed to the function.
      integer, intent(in) :: prec
      
      !----------------------------------------------------------!
      if ( prec == 4 ) then
        mpi_fwi = MPI_REAL
      else if ( prec == 8 ) then
!         mpi_fwi = MPI_DOUBLE
        mpi_fwi = MPI_DOUBLE_PRECISION
      else
        stop '***** ERROR - prec_mpi_fwi: mpi_fwi /= 4 .nor. 8 *****'
      end if
      
   end subroutine prec_mpi_fwi
   
#endif
   
   !**********************************************************!
    
  end module m_mpi




