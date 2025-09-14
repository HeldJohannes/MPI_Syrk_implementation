//
// Created by Johannes Held on 10.07.24.
//

#ifndef MPI_SYRK_IMPLEMENTATION_TWO_D_SYRK_H
#define MPI_SYRK_IMPLEMENTATION_TWO_D_SYRK_H

#ifdef __cplusplus
extern "C" {
#endif

#include "MPI_Syrk_implementation.h"


void distribute_input_matrix_2D(run_config *s, int rank, floatArray input_array, floatMatrix rank_input, MPI_Comm communicator);

int calculate_Q_i(intArray q_i, int m_bloc_i, int c);

void copy_to_2D(float *B, float **A, int k, int block_height, int block_length, int index);

void copy_to_d(double *A_i, float *A, int c, int block_height, int block_length, int index, _Bool trans);

void copy_to_f(float *A_i, float *A, int c, int block_height, int block_length, int index);

void copy_to_1D(float *dest, float *source, int r_pos, int w_pos, int block_height, int block_length, int row_length, int w_offset);


int cal_block_size(run_config *s);

void accumulate_B_into_A(run_config *s, int k, floatArray A, floatArray B, intArray R_k, intArray Q_i, floatMatrix input);

void two_d_syrk(run_config *s, int k, floatArray rank_result, floatMatrix input, MPI_Comm communicator);

#ifdef __cplusplus
}
#endif

#endif //MPI_SYRK_IMPLEMENTATION_TWO_D_SYRK_H
