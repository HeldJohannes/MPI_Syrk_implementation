#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>
#include <float.h>
#include "log.h"
#include "MPI_Syrk_implementation.h"
#include "one_d_syrk.h"
#include "two_d_syrk.h"
#include "three_d_syrk.h"

#define ROOT 0

/**
 *
 *
 *
 * @param argc number of input parameters
 * @param argv  available options are -n and -m; where -m specifies the number of input rows and -n the number of input columns
 * @return 0 if successful
 */
int main(int argc, char *argv[]) {

    /** ************************************************************************************************
     * STEP 1: Initialize the MPI environment and parse the input parameters
     ************************************************************************************************ */
    log_set_level(LOG_DEBUG);
    static run_config config;
    config.fileName = NULL;

    int world_size, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    config.world_size = world_size;

    // parse and check the input parameters  
    // this is done on all processors  
    parseInput(&config, argc, argv, rank);

    log_debug("Size = %d (SIZE_MAX = %zu) => %zu", config.m * config.n, SIZE_MAX, config.m * config.n * sizeof(float));

    /** ************************************************************************************************
     * STEP 2: Allocate memory for the input matrix and the index array and distribute the input matrix
     ************************************************************************************************ */
    
    // allocate memory for the input matrix (of size m * n) 
    // this is done on all processors
    floatArray input_array;
    allocate_float_array(&input_array, config.m * config.n);

    // TODO: check if this is necessary
    // allocate memory for the ... (of size m)
    float **input = (float **) calloc(config.m, sizeof(float *));
    for (int i = 0; i < config.m; ++i) {
        input[i] = &(input_array.data[config.n * i]);
    }


    // variable to store the index array value at the current rank
    int index_arr_rank;
    // array to store the index array values for all processors
    intArray index_arr;
    allocate_int_array(&index_arr, world_size);
    
    // array to store the cumulate index array values for all processors
    intArray cumulate_index_arr;
    allocate_int_array(&cumulate_index_arr, world_size);

    
    if (rank == ROOT) {

        // if a file is given read the input from the file
        // otherwise generate the input
        if (config.fileName != NULL) {
            read_input_file(rank, &config, input);
        } else {
            generate_input(&config, input);
        }

        // calculate how many cols each node gets and save it in index_arr
        index_calculation(index_arr, config.n, world_size);

        // TODO check if this is necessary:
        int rank_count = 0;
        for (int i = 0; i < world_size; ++i) {
            assert(cumulate_index_arr.data + i != NULL);
            cumulate_index_arr.data[i] = rank_count;
            assert(index_arr.length > i);
            rank_count += index_arr.data[i];
        }
    }

    // send index_arr[rank] to all processors using MPI_Scatter 
    // and continune work with index_arr_rank
    MPI_Scatter(
        index_arr.data,     // send buffer
        1,                  // number of elements to send to each processor
        MPI_INT,            // senddata type
        &index_arr_rank,    // receive buffer
        1,                  // number of elements to receive
        MPI_INT,            // recive data type
        0,                  // root process
        MPI_COMM_WORLD      // communicator (in this case the default communicator)
    );
    if(rank == ROOT) log_debug("successfully scattered the index_arr to all processors");

    // send cumulate_index_arr[rank] to all processors and work with cumulate_index_arr_rank
    //MPI_Scatter(cumulate_index_arr, 1, MPI_INT, &cumulate_index_arr_rank, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //input matrix for each node:
    floatMatrix rank_input;

    // transposed input matrix for each node:
    // needed only for 1D-Algo
    floatMatrix rank_input_t;

    if (config.algo != 3) {

        // if the algo is not 3 (2D-Algo) we can split the input matrix into rows and distribute them to the processors
        // this is done by using MPI_Scatterv to split the input matrix



        /** ************************************************************************************************
         * STEP 2.1:
         * allocate memory for the node input matrix of size m * index_arr[rank] (= index_arr_rank)
         ************************************************************************************************ */
        allocate_float_matrix(&rank_input, config.m, index_arr_rank);

        /** ************************************************************************************************
         * STEP 2.2:
         * Distribute the input matrix to all processors using MPI_Scatterv 
         * which allows to scatter the input matrix with different sizes 
         ************************************************************************************************ */
        for (int i = 0; i < config.m; ++i) {
            // for each row split and distribute across all processors
            MPI_Scatterv(
                input[i],                   // row to be scatterd
                index_arr.data,             // array holding number of elements to send to each processor
                cumulate_index_arr.data,    // array holding the displacement to apply to the message sent by each processor
                MPI_FLOAT,                  // data type of the buffer
                rank_input.data[i],         // receive buffer of size index_arr_rank
                index_arr_rank,             // number of elements to receive
                MPI_FLOAT,                  // data type of the recive buffer  
                ROOT,                       // root process
                MPI_COMM_WORLD              // communicator
            );
        }

        //transposed input matrix for each node:
        allocate_float_matrix(&rank_input_t, index_arr_rank, config.m);
        
        // compute the input matrix, and it's transpose,
        // which consists of the columns and all rows in that column:
        //computeInputAndTransposed(&config, rank, index_arr_rank, cumulate_index_arr_rank, input, rank_input, rank_input_t);

        transposeMatrix(config.m, index_arr_rank, rank_input.data, rank_input_t.data);

    } else {
        assert(input_array.data != NULL);

        // Allocate the arrays:
        // each node gets c block parts of A (A_i)
        // each block has m / (c * c) rows and n / (c+1) columns
        int row_block_height = config.m / (config.c * config.c);
        int row_block_length = config.n / (config.c + 1);

        // allocate memory for the node input matrix of size row_block_height * n
        allocate_float_matrix(&rank_input, config.c * row_block_height, row_block_length);
        
        distribute_input_matrix_2D(&config, rank, input_array.data, rank_input.data);
        assert(rank_input.data != NULL);
    }

    // free index_arr and cumulate_index_arr because they are no longer needed
    free(index_arr.data);
    free(cumulate_index_arr.data);

    // TODO: remove
    // to test if distribute_input_matrix_2D works correctly 
    // I print the input matrix and the rank_input matrix for a specific processor:
    /*if (rank == 11) {
        log_info("input array:");
        for (int i = 0; i < config.m; ++i) {
            for (int j = 0; j < config.n; ++j) {
                fprintf(stderr, "%0.0f ", input_array[i* config.n + j]);
            }
            fprintf(stderr, "\n");
        }

        log_info("rank_input:");
        for (int i = 0; i < config.m / config.c; i++)
        {
            for (int j = 0; j < config.n / (config.c + 1); j++)
            {
                fprintf(stderr, "%0.0f ", rank_input[i][j]);
            }
            fprintf(stderr, "\n"); 
        }
    }*/

    // input no longer needed
    free_float_array(input_array);

    // SYRK:
    // Compute the result matrix for each node which gets
    //  summed up by the MPI_Reduce_scatter() methode 
    //  and store the result in rank_syrk_result
    float *rank_syrk_result = (float *) calloc((long) config.m * config.m, sizeof(float));
    if (!rank_syrk_result) {
        log_fatal("[processor %d] Memory allocation failed for input with errno", rank);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }


    // Synchronize before starting time
    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();

    switch (config.algo) {
        case 0:
            syrkIterative(&config, rank, index_arr_rank, rank_input.data, rank_input_t.data, rank_syrk_result);
            break;
        case 1:
            improved_syrkIterative(&config, rank, index_arr_rank, rank_input.data, rank_input_t.data, rank_syrk_result);
            break;
        case 2:
            // 1D SYRK - OpenBLAS
            //In the case that m ≤ n and P is not too large
            syrk_withOpenBLAS(&config, rank, index_arr_rank, rank_input.data, rank_syrk_result);
            break;
        case 3:
            // 2D SYRK - OpenBLAS
            // In the case that m > n and P is not too large, a 2D algorithm is optimal
            assert(rank_input.data != NULL);
            assert(rank_syrk_result != NULL);
            two_d_syrk(&config, rank, rank_syrk_result, rank_input.data);
            break;
        case 4:
            // 3D SYRK - OpenBLAS
            three_d_syrk(&config, rank, rank_syrk_result, rank_input.data);
            break;
        default:
            log_fatal("no SYRK operator selected --> error ALOG %d not in [0..2]", config.algo);
            error_exit(rank, argv[0], "no SYRK operator selected");
    }
    // Synchronize again before obtaining the time
    //MPI_Barrier(MPI_COMM_WORLD);
    //log_info("Syrk algo(%d) took %f sec", ALGO, MPI_Wtime() - start);
    log_debug("Successfully freed the buffer -> cumulate_index_arr");
    free_float_matrix(&rank_input);
    log_debug("Successfully freed the buffer -> rank_input");
    free_float_matrix(&rank_input_t);
    log_debug("Successfully freed the buffer -> rank_input_t");


    intArray counts;
    allocate_int_array(&counts, world_size);
    
    index_calculation(counts, (long) config.m * config.m, world_size);
    log_debug("Successfully allocated the buffer -> counts");

    floatArray reduction_result;
    allocate_float_array(&reduction_result, counts.data[rank]);

    // Save the time before the MPI_Reduce_scatter (Start the colock)
    double start_mpi_reduce_scatter = MPI_Wtime();

    // Reduce the results of all processors and scatter the result to all processors
    int result = MPI_Reduce_scatter(
        rank_syrk_result,       // send buffer
        reduction_result.data,  // receive buffer
        counts.data,            // number of elements to receive from each processor
        MPI_FLOAT,              // data type of the buffer
        MPI_SUM,                // operation to perform
        MPI_COMM_WORLD          // communicator
    );
    // check if the MPI_Reduce_scatter was successful
    if (result != MPI_SUCCESS) {
        log_error("MPI_Reduce_scatter returned exit coed %d", result);
    }

    // Save the time after the MPI_Reduce_scatter (Stop the colock)
    double runtime_mpi_reduce_scatter = MPI_Wtime() - start_mpi_reduce_scatter;
    // log the time needed for the MPI_Reduce_scatter
    log_debug("[rank %d]: MPI_Reduce took %f sec", rank, runtime_mpi_reduce_scatter);

    if (rank == ROOT) {
        double runtime = MPI_Wtime() - start;
        printf("The process took %f seconds to run.\n", runtime);

        log_debug("m = %d", config.m);
        log_debug("[int] config.m * config.m = %d", config.m * config.m);
        log_debug("[long] config.m * config.m = %ld", config.m * config.m);

        floatArray buffer;
        allocate_float_array(&buffer, config.m * config.m);

        intArray displacements;
        allocate_int_array(&displacements, world_size);

        displacements.data[0] = 0;
        for (int i = 1; i < world_size; ++i) {
            displacements.data[i] = displacements.data[i - 1] + counts.data[i - 1];
        }
//        printf("counts:\n");
//        printResult(rank, world_size, counts);
//        printf("displacements:\n");
//        printResult(rank, world_size, displacements);
        int status = MPI_Gatherv(
            reduction_result.data,      // send buffer
            counts.data[rank],          // number of elements to send
            MPI_FLOAT,                  // data type of the
            buffer.data,                // receive buffer
            counts.data,                // number of elements to receive from each processor
            displacements.data,         // An array containing the displacement to apply to the message received by each process
            MPI_FLOAT,                  // data type of the buffer
            0,                          // root process
            MPI_COMM_WORLD              // communicator
        );

        if (status != MPI_SUCCESS) {
            log_error("MPI_Gatherv returned %d", status);
        }

        //Print the result:
        log_info("Values gathered in the buffer on process %d\n", rank);

        if (config.print_result) {
            // No synchronization needed because only processor 0 operates here
            double start_print_results = MPI_Wtime();

            printResult(&config, config.m, buffer.data);

            double runtime_print_results = MPI_Wtime() - start_print_results;
            log_info("runtime_print_results = %f", runtime_print_results);
        }

        free_float_array(buffer);
        free_int_array(displacements);
    } else {
        // all processes in the communicator must invoke Gatherv 
        // the reciver buffer, size and displacements are ignored and therefore NULL
        int status = MPI_Gatherv(
            reduction_result.data,  // send buffer
            counts.data[rank],      // number of elements to send
            MPI_FLOAT, 
            NULL,
            NULL, 
            NULL, 
            MPI_FLOAT, 
            0, 
            MPI_COMM_WORLD
        );
        if (status != MPI_SUCCESS) {
            log_error("MPI_Gatherv returned %d", status);
        }
    }

    free_float_array(reduction_result);

    free(rank_syrk_result);
    log_debug("Successfully freed the buffer -> rank_syrk_result");

    log_info("[rank %d] finished the program --> MPI_Finalize()", rank);
    
    // finish the program:
    int status;
    if ((status = MPI_Finalize()) != MPI_SUCCESS) {
        log_error("MPI Failed with %d", status);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

