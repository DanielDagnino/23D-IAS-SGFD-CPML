  !*********************************************************************/
  !** This code has been done in the Barcelona Center for Subsurface 
  !** Imaging (BCSI).
  !** Authors: Daniel Dagnino.
  !*********************************************************************/
  !** Short resume of this code/file:
  !** Here, we define several modules for common variables.
  !*********************************************************************/
  
  
  
  !**********************************************************!
  module m_dyw
    implicit none
    logical, parameter :: dyw_3d = .false.
    logical :: dyw_allo     = .true.
    logical :: dyw_allo_b   = .true.
    logical :: dyw_deallo   = .false.
    logical :: dyw_deallo_b = .false.
    logical :: dyw_read, dyw_reread_tr
    logical :: dyw_read_pos
    logical :: dyw_reread_pr = .false.
    logical :: ck_mem        = .true.
    logical :: dyw_read_iH   = .false.
    logical :: dyw_read_J   = .false.
    logical :: dyw_get_iH_J = .false.
  end module m_dyw
  
  !**********************************************************!
  module m_type_inversion
    use m_data_kind
    implicit none
    
    real(sp), allocatable :: prec_mod(:,:,:)
    
    logical :: norm_in_tr_1p1
    logical :: different_time
    
    integer :: geometry_type
    integer :: meth_OF
    real(8) :: pow_OF = 2.
    integer :: meth_min_QC
    real(sp) :: chng_min_QC
    real(sp) :: change_cota
    
    logical :: use_iH_red
    logical :: use_H0
    real(sp) :: power_H0
    logical :: use_iH
    logical :: full_H
    
    integer :: lay_frz
    integer :: meth_inv
    logical :: source_inver, source_fix
    logical :: fast_search
    
    logical :: calib, calib_sou, corr_dly
    logical :: calib_sr, calib_sr_sou
    
    integer :: prec_d
    integer :: prec_d_w
    logical :: dyw_Wd
    integer :: prec_m
    logical :: border_red
    
    logical :: vel_cut_off, dnst_cut_off
    real(sp) :: vel_max, vel_min, dnst_max, dnst_min
    
    real(sp) :: offset_cal_min, offset_cal_max
    real(sp) :: offset_sou_max, offset_sou_min
    real(sp) :: offset_pr_min, offset_pr_max
    real(sp) :: offset_allowed_grad, long_streamer
    
    integer  :: filt_grad
    real(sp) :: c_mean
    real(sp) :: pl_x_filt, pl_y_filt, l_x_filt, l_y_filt
    
    real(sp) :: p0_cut, p1_cut
    
    logical :: write_model, write_trace, write_illu, write_sou, write_adj_sou, write_grad, write_sear
    
    integer :: rank_wrt, isg_wrt
    
    logical, parameter :: gardner = .true.
    logical, parameter :: mult_tr_sqrt_t = .true.
    
    ! NOTE: To use the code to a forward modeling use only_source_inver=.true. and sou_inv_water=.false.
    logical, parameter :: only_source_inver = .false.
    logical :: sou_inv_water = .true.
    
  end module m_type_inversion
  
  !**********************************************************!
  module m_inv_data
    use m_data_kind
    implicit none
    integer :: strat_inv, iter, iter_freq
    real(sp) :: freq_inv
  end module m_inv_data
  
  !**********************************************************!
  module m_work
    use m_data_kind
    implicit none
    integer :: n_shot
    real(sp) :: mf_i0_f0, norm_grad0
    real(sp), allocatable :: mf_sg(:), mf_sg_1(:), mf_sg_2(:), mf_sg_fit(:), mf_sg_ck(:), mf_sg_min(:), coef_sea(:)
    real(sp), allocatable :: mf_sg_a(:), mf_sg_b(:)
    real(sp), allocatable :: es_sg(:)
    real(sp), allocatable :: mean_grad0_sg(:), sig_grad0_sg(:), mean_grad0_it0_sg(:)
    logical :: first_normalization
  end module m_work
  
  !**********************************************************!
  module m_model
    use m_data_kind
    implicit none
    real(sp), allocatable :: cmpr(:,:,:), buoy(:,:,:)
!     real(sp), allocatable :: cmpr0(:,:,:), cmpr(:,:,:), buoy(:,:,:)
  end module m_model
  
  !**********************************************************!
  module m_variable
    use m_data_kind
    implicit none
    real(sp), allocatable :: m_0(:,:,:), m_1(:,:,:), m_2(:,:,:), m_fit(:,:,:), m_min(:,:,:)
    real(sp), allocatable :: search(:,:,:)
  end module m_variable
  
  !**********************************************************!
  module m_filter
    use m_data_kind
    implicit none
    integer :: filter
    
    integer, parameter  :: zero_pad_time = 3   ! 3 (recommended)
    real(sp), parameter :: pt_time = 0.025     ! seconds, 0.25*period at 10Hz
    
    integer, parameter :: zero_pad_space = 1      ! 2-5 (recommended)
    ! Marmousi size.
    real(sp), parameter :: px_space = 101.  ! 100. meters (recommended)
    real(sp), parameter :: py_space = 101.  ! 100. meters (recommended)
!     ! Nicaragua size.
!     real(sp), parameter :: px_space = 2001.  ! 100. meters (recommended)
!     real(sp), parameter :: py_space = 2001.  ! 100. meters (recommended)
    
    integer, parameter :: order_filt = 6      ! 5-10 (recommended)
    integer, parameter :: o_spc_filt = 6      ! 5-10 (recommended)
  end module m_filter
  
  !**********************************************************!
  module m_phys
    use m_data_kind
    implicit none
    real(sp), parameter :: c_water = 1.520_sp   ! 1.500_sp 
  end module m_phys
  
  !**********************************************************!
  module m_wrt_tr
    implicit none
    integer :: n_id_wrt_tr
    integer, allocatable :: id_wrt_tr(:)
  end module m_wrt_tr


