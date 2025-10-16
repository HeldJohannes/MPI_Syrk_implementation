#include "three_d_syrk.h"

#define TEST_RANK -1
#define PRINT_DEBUG false

/**
 * @brief This function is used to compute the result of the SYRK operation for a given rank.
 * 
 * @param s A pointer to a run_config structure containing the configuration values.
 * @param rank The rank of the current processor.
 * @param comm_rank The rank within the communicator.
 * @param rank_result A pointer to a float array where the result will be stored.
 * @param input A pointer to a 2D float array containing the input values.
 * 
 * @see run_config
 * @see two_d_syrk
 * @see MPI_Syrk_implementation
 */
void three_d_syrk(run_config *s, int rank, int comm_rank, floatArray rank_result, floatMatrix input, MPI_Comm communicator) {

    //TODO: remove 
    //fprintf(stderr, "TEST: three_d_syrk with rank %d \n", rank);

    //TODO remove:
    // print input matrix for a specific processor to test if correct
    if (PRINT_DEBUG) {
        FILE *fp;
        char filename[256];
        snprintf(filename, sizeof(filename), "log_input_3D_%d", rank);
        fp = fopen(filename, "w");
        printMatrix(input, fp);
        fclose(fp);
    }

    // requires the number of processors to be |π| = p1 * p2
    // where p1 = c * (c + 1) and c is a prime number

    /* ****************************************
    STEP 1: Split the comm_rank into two parts
    ******************************************/ 

    // in order to work the rank has to be split into 2 parts:
    // the first part k which is given by the modulo of the rank
    // k (0 <= k < p1) and p1 = c * (c + 1)
    // int p1 = s->c * (s->c + 1);
    // int k = rank % p1;

    // and into the second part l which is given by the division of the rank
    // l (0 <= l < p2)
    // int l = rank / p1;


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
    run_config_copy->world_size = s->world_size / s->P2;     
    run_config_copy->m = s->m;              
    run_config_copy->n = s->n / s->P2;              
    run_config_copy->c = s->c;              
    run_config_copy->P2 = s->P2;             
      
    two_d_syrk(run_config_copy, comm_rank, rank_result, input, communicator);

    free(run_config_copy);

    /* ****************************************
    * STEP 3: Compute the final result C_kl by summing up the intermediate results C_kl 
    * unsing REDUCE-SCATTER on C_kl and π_k*
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

void distribute_input_matrix_3D(run_config *s, int rank, int *comm_rank, floatArray input_array, floatMatrix rank_input, MPI_Comm *pMpiCommunicators) {

    MPI_Comm communicator = MPI_COMM_WORLD;
            
    // Create P2 groups:
    MPI_Group main_group;
    MPI_Comm_group(communicator, &main_group);

    // create P2 mpi - groups:
    MPI_Group *pMpiGroups = (MPI_Group *) malloc(s->P2 * sizeof (MPI_Group ));
    if (!pMpiGroups) {
        log_fatal("Memory allocation failed for pMpiGroups", 0);
        MPI_Abort(communicator, EXIT_FAILURE);
    }

    assert(s->world_size % s->P2 == 0);

    int P1 = s->c * (s->c +1);

    intMatrix processor_Ranks;
    allocate_int_matrix(&processor_Ranks, s->P2, P1);

    for (int i = 0; i < s->world_size; i++){
        //TODO find out why l = i % s->P2 works and the other not (???)
        //int l = i / P1; // [0 .. P2] because Π = P1 * P2
        int l = i % s->P2;
        //int k = i % P1; // [0 .. P1]
        int k = i / s->P2;
        processor_Ranks.data[l][k] = i;
    }

    if (false) {
        FILE *fp;
        char filename[256];
        snprintf(filename, sizeof(filename), "log_processor_Ranks");
        fp = fopen(filename, "w");
        printIntMatrix(processor_Ranks, fp);
        fclose(fp);
    }


    // create the new groups and communicators
    for (int i = 0; i < s->P2; ++i) {

        // calculate the Q_i values
        // create the new group
        int ret = MPI_Group_incl(
            main_group,                // parent group
            P1,                        // number of nodes in new group
            processor_Ranks.data[i],   // ranks of the processes in the new group
            &pMpiGroups[i]             // new group
        );
        // check if the MPI_Group_incl was successful
        assert(pMpiGroups[i] != NULL);
        if (ret != MPI_SUCCESS) {
            log_fatal("MPI_Group_incl failed for group %d", i);
            MPI_Abort(communicator, ret);
        }
        // create the new communicator
        int err = MPI_Comm_create(
            communicator,           // parent communicator
            pMpiGroups[i],          // group
            &pMpiCommunicators[i]   // new communicator
        );
        // check if the MPI_Comm_create was successful
        int l = rank % s->P2;
        if (l == i && err == MPI_SUCCESS) {
            assert(pMpiCommunicators[i] != MPI_COMM_NULL);
        } else {
            assert(pMpiCommunicators[i] == MPI_COMM_NULL);
        }
    }

    /** ****************************************
     * split the data of A into P2 slices (0 <= l < P2) [A_*l] in the first processors of each MPI Π_*l group
     **************************************** */


    //TODO find a better solution:
    // broadcast the input matrix to all processors
    MPI_Bcast(
        input_array.data,       // buffer to broadcast
        input_array.length,     // number of elements in the buffer
        MPI_FLOAT,              // data type of the buffer
        0,                      // root process
        communicator            // communicator
    );

    floatArray input_slice;
    allocate_float_array(&input_slice, s->m, (s->n / s->P2));
    

    // split the input array in such way that the first Π / P2 processors get A_*1 and so on.

    int block_length = s->n / s->P2;
    int shift = (rank % s->P2) * block_length;
    copy_array_slice(input_slice, input_array, s->m, block_length, shift);

    // call distribute_input_matrix_2D for all P2 A_*l input matrices with the corresponding communicator

    

    // validate that the conmmunicator size is as expected:
    int comm_size;
    MPI_Comm_size(pMpiCommunicators[rank % s->P2], &comm_size);
    assert(comm_size == s->c * (s->c+1));

    // Compute the rank of this process in the new communicator group:
    MPI_Comm_rank(pMpiCommunicators[rank % s->P2], comm_rank);
    assert(*comm_rank < s->c * (s->c+1));

    //TODO: remove
    // to test if distribute works correctly 
    if (PRINT_DEBUG) {
        FILE *fp;
        char filename[256]; 
        snprintf(filename, sizeof(filename), "log_input_slice_%d", rank);
        fp = fopen(filename, "w");
        printArray(input_slice, s->m, (s->n / s->P2), fp);
        fclose(fp);
    }

    run_config *run_config_copy = malloc(sizeof(run_config));
    if (!run_config_copy) {
        log_fatal("Memory allocation failed for copy of run_config");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    run_config_copy->algo = s->algo;         
    run_config_copy->world_size = s->world_size;     
    run_config_copy->m = s->m;              
    run_config_copy->n = s->n / s->P2;              
    run_config_copy->c = s->c;              
    run_config_copy->P2 = s->P2;

    // distribute the input matrix slice A_*l acros c^2 processors.
    distribute_input_matrix_2D(run_config_copy, *comm_rank, input_slice, rank_input, pMpiCommunicators[rank % s->P2]);

    free(run_config_copy);

    if (PRINT_DEBUG) {
        FILE *fp;
        char filename[256];
        snprintf(filename, sizeof(filename), "log_rank_input_%d", rank);
        fp = fopen(filename, "w");
        printMatrix(rank_input, fp);
        fclose(fp);
    }
}

