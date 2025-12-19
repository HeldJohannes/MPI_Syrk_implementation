#include "one_d_syrk.h"

#ifdef USE_CBLAS_64
    #define CBLAS_SYRK cblas_ssyrk64_
#else
    #define CBLAS_SYRK cblas_ssyrk
#endif

/**
 * @brief This function computes the SYRK operation using a triple nested for-loop.
 * 
 * @param s A pointer to a run_config structure containing the configuration values.
 * @param rank The rank of the current processor.
 * @param index_arr_rank The size of the input slice for the current processor.
 * @param rank_input A pointer to a 2D float array containing the input values.
 * @param rank_input_t A pointer to a 2D float array containing the transposed input values.
 * @param rank_result A pointer to a float array where the result will be stored.
 */
void syrkIterative(run_config *s, int rank, int index_arr_rank, floatMatrix rank_input, floatMatrix rank_input_t,
                   floatArray rank_result) {
    log_trace("[rank %d] syrkIterative()", rank);
    // for each result row:
    for (int row = 0; row < s->m; ++row) {
        // for each result column
        for (int col = 0; col < s->m; ++col) {
            // run for slice of the input:
            for (int c = 0; c < index_arr_rank; ++c) {
                rank_result.data[row * s->m + col] += rank_input.data[row][c] * rank_input_t.data[c][col];
            }
        }
    }
}

/**
 * @brief This function computes the SYRK operation using an improved triple nested for-loop. 
 *  The improvement comes from reducing the number of iterations in the inner loop by taking advantage of the symmetry of the result matrix.
 *  In this version, the inner loop iterates only over the upper triangular part of the result matrix,
 *  and the result is stored in a 1D array.
 * 
 * @param s A pointer to a run_config structure containing the configuration values.
 * @param rank The rank of the current processor.
 * @param index_arr_rank The size of the input slice for the current processor.
 * @param rank_input A pointer to a 2D float array containing the input values.
 * @param rank_input_t A pointer to a 2D float array containing the transposed input values.
 * @param rank_result A pointer to a float array where the result will be stored.
 */
void improved_syrkIterative(run_config *s, int rank, const int index_arr_rank, float **rank_input, float **rank_input_t,
                            float *rank_result) {
    log_trace("[rank %d] improved_syrkIterative()", rank);
    for (int row = 0; row < s->m; ++row) {
        for (int col = row; col < s->m; ++col) {
            for (int c = 0; c < index_arr_rank; ++c) {
                rank_result[row * s->m + col] += rank_input[row][c] * rank_input_t[c][col];
            }
        }
    }
}

/**
 * @brief This function computes the SYRK operation using the OpenBLAS library.
 * 
 * @param s A pointer to a run_config structure containing the configuration values.
 * @param rank The rank of the current processor.
 * @param index_arr_rank The size of the input slice for the current processor.
 * @param rank_input A pointer to a 2D float array containing the input values.
 * @param rank_result A pointer to a float array where the result will be stored.
 */
void syrk_withOpenBLAS(run_config *config, int rank, int index_arr_rank, float **rank_input, float *rank_result) {
    log_trace("[rank %d] syrk_withOpenBLAS()", rank);
    // transform 2d array to 1d:
    float *A = (float *) calloc(config->m * index_arr_rank, sizeof(float));
    if (A == NULL) {
        log_error("Calloc failed");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    for (int i = 0; i < config->m; ++i) {
        for (int j = 0; j < index_arr_rank; ++j) {
            A[i * index_arr_rank + j] = rank_input[i][j];
        }
    }
    
    // compute syrk:
    CBLAS_SYRK(
            CblasRowMajor,
            CblasUpper,
            CblasConjNoTrans,
            config->m,
            index_arr_rank,
            1.0f,
            A,
            index_arr_rank,
            0.0f,
            rank_result,
            config->m);
    
    // cleanup
    free(A);
}
