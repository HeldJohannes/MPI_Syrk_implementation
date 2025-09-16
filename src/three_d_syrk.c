#include "three_d_syrk.h"

#define TEST_RANK -1

/**
 * @brief This function is used to compute the result of the SYRK operation for a given rank.
 * 
 * @param s A pointer to a run_config structure containing the configuration values.
 * @param rank The rank of the current processor.
 * @param rank_result A pointer to a float array where the result will be stored.
 * @param input A pointer to a 2D float array containing the input values.
 * 
 * @see run_config
 * @see two_d_syrk
 * @see MPI_Syrk_implementation
 */
void three_d_syrk(run_config *s, int rank, floatArray rank_result, floatMatrix input, MPI_Comm communicator) {

    //TODO: remove 
    //fprintf(stderr, "TEST: three_d_syrk with rank %d \n", rank);

    //TODO remove:
    // print input matrix for a specific processor to test if correct
    if (rank == TEST_RANK) {
        FILE *fp;
        fp = fopen("log_input", "w");
        printMatrix(input, fp);
        fclose(fp);
    }

    // requires the number of processors to be |π| = p1 * p2
    // where p1 = c * (c + 1) and c is a prime number

    /* ****************************************
    STEP 1: Split the rank into two parts
    ******************************************/ 

    // in order to work the rank has to be split into 2 parts:
    // the first part k which is given by the modulo of the rank
    // k (0 <= k < p1) and p1 = c * (c + 1)
    int p1 = s->c * (s->c + 1);
    int k = rank % p1;

    // and into the second part l which is given by the division of the rank
    // l (0 <= l < p2)
    int l = rank / p1;


    /* ****************************************
    STEP 2: Compute intermediate results C_kl using the two_d_syrk function on the input slice A_*l and π_*l
    ******************************************/ 

    // input Blocks have size (m / c^2) x (n / p2)

    // change n:
    run_config *run_config_copy = malloc(sizeof(run_config));
    if (!run_config_copy) {
        log_fatal("Memory allocation failed for copy of run_config");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    run_config_copy->algo = s->algo;         
    run_config_copy->world_size = s->world_size;     
    run_config_copy->m = s->m;              
    run_config_copy->n = s->n;              
    run_config_copy->c = s->c;              
    run_config_copy->P2 = s->P2;             
      
    two_d_syrk(run_config_copy, rank, rank_result, input, communicator);

    free(run_config_copy);

    /* ****************************************
    STEP 3: Compute the final result C_kl by summing up the intermediate results C_kl 
    unsing REDUCE-SCATTER on C_kl and π_k*
    ******************************************/


    log_trace("rank_result %p", rank_result);
}


void copy_array_slice(floatArray A_i, floatArray A, int block_height, int block_length, int shift) {

    assert(A.data != NULL);
    assert(A_i.data != NULL);
    assert(A.row_length > 0);
    assert(A_i.row_length > 0);


    for (int i = 0; i < block_height; ++i) {
        for (int j = 0; j < block_length; ++j) {

            int index_A_i = i * block_length + j;
            int index_A = i * A.row_length + j + shift;

            assert__(index_A_i < A_i.length) {
                fprintf(stderr, "assertion will fail\n"); 
            }
            assert(index_A < A.length) ;

            A_i.data[index_A_i] = A.data[index_A];
        }
    }
}
