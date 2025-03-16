//
// Created by Johannes Held on 10.07.24.
//

#ifndef MPI_SYRK_IMPLEMENTATION_TWO_D_SYRK_H
#define MPI_SYRK_IMPLEMENTATION_TWO_D_SYRK_H

#ifdef __cplusplus
extern "C" {
#endif

#include "MPI_Syrk_implementation.h"


void distribute_input_matrix_2D(run_config *s, int rank, float *input_array, float **rank_input);

int calculate_Q_i(int *q_i, int m_bloc_i, int c);

void copy_to_2D(float *B, float **A, int k, int block_height, int block_length, int index);

int cal_block_size(run_config *s);

void two_d_syrk(run_config *s, int k, float *rank_result, float **input);

#ifdef __cplusplus
}
#endif

#endif //MPI_SYRK_IMPLEMENTATION_TWO_D_SYRK_H
