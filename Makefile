#!/bin/bash
include $(CWPROOT)/src/Makefile.config
D = $L/libcwp.a $L/libpar.a $L/libsu.a

FC = mpif90
FFLAGS = -g -fpp -Dusempi \
         -O3 -ipo -unroll-aggressive\
         -m64 -xSSE4.2 \
         -heap-arrays -mcmodel=large -shared-intel

#FC = gfortran
#FFLAGS = -g -cpp -Dusempi \
         -O3 -m64 -msse4.2
         
LIBS = -L./1-lib_base_v3/lib -lbase \
       -L../../../fftpack5.1/lib -lfftpack \
       -L../../../dfftpack5.1/lib -ldfftpack \
HEAD = -I./1-lib_base_v4/src

EXE = fwi

SRC = 00-m_data_kind.f90 \
	 m_mpi.f90 \
	 00-module.f90 \
	 m_solv3Dahvc_data_kind.f90 \
	 01-correctness_criterion.f90 \
	 01-m_mat.f90 \
	 01-mod_definit.f90 \
	 m_support_funct.f90 \
	 m_sg_path_file.f90 \
	 m_sg_data.f90 \
	 m_geo.f90 \
	 m_sg_type.f90 \
	 m_grid_bath.f90 \
	 m_sg_func.f90 \
	 m_sg_get_pos.f90 \
	 m_sg_get_sou.f90 \
	 m_sg_get_tr.f90 \
	 m_sg_filt.f90 \
	 m_sg.f90 \
	 m_get_freq_inv.f90 \
	 m_mygeo.f90 \
	 m_solver_source.f90 \
	 m_interpol_sg.f90 \
	 m_solver.f90 \
	 m_deriv_b.f90 \
	 m_deriv.f90 \
	 m_PML_definitions_param_face.f90 \
	 m_PML_building.f90 \
	 m_coll_allo.f90 \
	 m_change_var.f90 \
	 m_clear_div_sou.f90 \
	 m_copy_mat_2d.f90 \
	 m_copy_mat_3D.f90 \
	 m_copy_mat_4.f90 \
	 m_get_id_shots.f90 \
	 m_get_initial_data.f90 \
	 m_get_strategy.f90 \
	 m_get_model_and_QC.f90 \
	 m_source_ricker.f90 \
	 RK_der_edge_BD.f90 \
	 RK_der_edge_BU.f90 \
	 RK_der_edge_FD.f90 \
	 RK_der_edge_FU.f90 \
	 RK_der_edge_LB.f90 \
	 RK_der_edge_LD.f90 \
	 RK_der_edge_LF.f90 \
	 RK_der_edge_LU.f90 \
	 RK_der_edge_RB.f90 \
	 RK_der_edge_RD.f90 \
	 RK_der_edge_RF.f90 \
	 RK_der_edge_RU.f90 \
	 RK_der.f90 \
	 RK_der_face_B.f90 \
	 RK_der_face_D.f90 \
	 RK_der_face_F.f90 \
	 RK_der_face_L.f90 \
	 RK_der_face_R.f90 \
	 RK_der_face_U.f90 \
	 RK_der_vertex_LBD.f90 \
	 RK_der_vertex_LBU.f90 \
	 RK_der_vertex_LFD.f90 \
	 RK_der_vertex_LFU.f90 \
	 RK_der_vertex_RBD.f90 \
	 RK_der_vertex_RBU.f90 \
	 RK_der_vertex_RFD.f90 \
	 RK_der_vertex_RFU.f90 \
	 RK_dif_op_FS.f90 \
	 RK_dif_op_FS_back.f90 \
	 RK_dif_op_split_back.f90 \
	 RK_dif_op_split_bd.f90 \
	 RK_dif_op_split.f90 \
	 RK_new_step_edge.f90 \
	 RK_new_step.f90 \
	 RK_new_step_face.f90 \
	 RK_new_step_FS.f90 \
	 RK_new_step_vertex.f90 \
	 RK_new_step_split.f90 \
	 RK_sou_adj.f90 \
	 RK_sou.f90 \
	 RK_step_edge.f90 \
	 RK_step.f90 \
	 RK_step_face.f90 \
	 RK_step_FS.f90 \
	 RK_step_vertex.f90 \
	 RK_step_split_bd.f90 \
	 RK_step_split.f90 \
	 RK_loop.f90 \
	 solver_acoustic_3d_vnh_RK4_CPML_deallo.f90 \
	 solver_acoustic_3d_vnh_RK4_CPML_test.f90 \
	 solv3Dahv_test_all.f90 \
	 zz-main.f90
	 
# No need to edit below this line
.SUFFIXES: .f90 .o

OBJ = $(SRC:.f90=.o)

.f90.o:
	$(FC) $(FFLAGS) $(HEAD) -c $<

all: $(EXE)

$(EXE): $(OBJ)
	$(FC) $(OBJ) $(HEAD) $(LIBS) -o $@ 

$(OBJ): $(MF)

tar:
	tar cvf $(EXE).tar $(MF) $(SRC)

clean:
	rm -f $(OBJ) $(EXE) core *.mod *.MOD
