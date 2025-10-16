#ifndef MPI_SYRK_IMPLEMENTATION_THREE_D_SYRK_H
#define MPI_SYRK_IMPLEMENTATION_THREE_D_SYRK_H

#ifdef __cplusplus
extern "C" {
#endif

#include "MPI_Syrk_implementation.h"
#include "two_d_syrk.h"

void copy_array_slice(floatArray A_i, floatArray A, int block_height, int block_length, int shift);

void three_d_syrk(run_config *s, int rank, int comm_rank, floatArray rank_result, floatMatrix input, MPI_Comm communicator);

void distribute_input_matrix_3D(run_config *s, int rank, int *comm_rank, floatArray input_array, floatMatrix rank_input, MPI_Comm *pMpiCommunicators);

#ifdef __cplusplus
}
#endif

#endif //MPI_SYRK_IMPLEMENTATION_THREE_D_SYRK_H