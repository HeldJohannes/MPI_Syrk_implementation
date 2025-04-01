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
#define TEST_RANK 11

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
    int ret = parseInput(&config, argc, argv, rank);
    if (ret != 0) {
        log_error("Error while parsing the input parameters");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    log_debug("Size = %d (SIZE_MAX = %zu) => %zu", config.m * config.n, SIZE_MAX, config.m * config.n * sizeof(float));

    /** ************************************************************************************************
     * STEP 2: Allocate memory for the input matrix and the index array and distribute the input matrix
     ************************************************************************************************ */
    
    // allocate memory for the input matrix (of size m * n) 
    // this is done on all processors
    floatArray input_array;
    allocate_float_array(&input_array, config.m, config.n);

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

    int comm_rank = -1;
    MPI_Comm *pMpiCommunicators = NULL;

    switch (config.algo) {
        case 0:
        case 1:
        case 2: 
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
            break;
        case 3: 
            assert(input_array.data != NULL);

            // Allocate the arrays:
            // each node gets c block parts of A (A_i)
            // each block has m / (c * c) rows and n / (c+1) columns
            int row_block_height = config.m / (config.c * config.c);
            int row_block_length = config.n / (config.c + 1);

            // allocate memory for the node input matrix of size row_block_height * n
            allocate_float_matrix(&rank_input, config.c * row_block_height, row_block_length);
            
            distribute_input_matrix_2D(&config, rank, input_array, rank_input, MPI_COMM_WORLD);
            assert(rank_input.data != NULL);
            break;
        case 4:
            assert(input_array.data != NULL);

            /** ************************************
             * STEP ...
             * 
             * compute the groupe Π_*l
             ************************************ */

            MPI_Comm communicator = MPI_COMM_WORLD;

            // Create P2 groups:
            MPI_Group main_group;
            MPI_Comm_group(communicator, &main_group);

            // create P2 mpi - groups:
            MPI_Group *pMpiGroups = (MPI_Group *) malloc(config.P2 * sizeof (MPI_Group ));
            if (!pMpiGroups) {
                log_fatal("Memory allocation failed for pMpiGroups", 0);
                MPI_Abort(communicator, EXIT_FAILURE);
            }

            // create P2 mpi - communicators:
            pMpiCommunicators = (MPI_Comm *) malloc(config.P2 * sizeof (MPI_Comm ));
            if (!pMpiCommunicators) {
                log_fatal("[rank %d] Memory allocation failed for pMpiCommunicators", rank);
                MPI_Abort(communicator, EXIT_FAILURE);
            }

            assert(config.world_size % config.P2 == 0);

            int P1 = config.c * (config.c +1);

            intMatrix processor_Ranks;
            allocate_int_matrix(&processor_Ranks, config.P2, P1);

            for (int i = 0; i < config.world_size; i++)
            {
                //TODO find out why l = i % config.P2 works and the other not (???)
                //int l = i / P1; // [0 .. P2] because Π = P1 * P2
                int l = i % config.P2;
                //int k = i % P1; // [0 .. P1]
                int k = i / config.P2;
                processor_Ranks.data[l][k] = i;
            }

            if (rank == TEST_RANK) {
                FILE *fp;
                fp = fopen("log_processor_Ranks", "w");
                printIntMatrix(processor_Ranks, fp);
                fclose(fp);
            }


            // create the new groups and communicators
            for (int i = 0; i < config.P2; ++i) {

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
                int l = rank % config.P2;
                if (l == i && err == MPI_SUCCESS) {
                    assert(pMpiCommunicators[i] != MPI_COMM_NULL);
                } else {
                    assert(pMpiCommunicators[i] == MPI_COMM_NULL);
                }
            }

            /**
             * STEP .2
             * 
             * split the data of A into P2 slices (0 <= l < P2) [A_*l] in the first processors of each MPI Π_*l group
             */
            //

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
            allocate_float_array(&input_slice, config.m, (config.n / config.P2));
            

            // split the input array in such way that the first Π / P2 processors get A_*1 and so on.

            int block_length = config.n / config.P2;
            int shift = (rank % config.P2) * block_length;
            copy_array_slice(input_slice, input_array, config.m, block_length, shift);

            // call distribute_input_matrix_2D for all P2 A_*l input matrices with the corresponding communicator

            // Each processor has after the distribution c blocks of A_il of size (m/c^2)×(n/p2)
            int block_height = (config.m / (config.c * config.c));
            int block_length_ri =  config.n / config.P2;
            allocate_float_matrix(&rank_input, config.c * block_height, block_length_ri / (config.c + 1));

            // validate that the conmmunicator size is as expected:
            int comm_size;
            MPI_Comm_size(pMpiCommunicators[rank % config.P2], &comm_size);
            assert(comm_size == config.c * (config.c+1));

            // Compute the rank of this process in the new communicator group:
            MPI_Comm_rank(pMpiCommunicators[rank % config.P2], &comm_rank);
            assert(comm_rank < config.c * (config.c+1));

            //TODO: remove
            // to test if distribute works correctly 
            if (rank == TEST_RANK) {
                log_info("printing input_slice:");
                FILE *fp;
                fp = fopen("log_input_slice", "w");
                printArray(input_slice, config.m, (config.n / config.P2), fp);
                fclose(fp);
            }

            run_config *run_config_copy = malloc(sizeof(run_config));
            if (!run_config_copy) {
                log_fatal("Memory allocation failed for copy of run_config");
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }

            run_config_copy->algo = config.algo;         
            run_config_copy->world_size = config.world_size;     
            run_config_copy->m = config.m;              
            run_config_copy->n = config.n / config.P2;              
            run_config_copy->c = config.c;              
            run_config_copy->P2 = config.P2;

            // distribute the input matrix slice A_*l acros c^2 processors.
            distribute_input_matrix_2D(run_config_copy, comm_rank, input_slice, rank_input, pMpiCommunicators[rank % config.P2]);

            free(run_config_copy);

            if (rank == TEST_RANK) {
                log_info("printing rank_input:");
                FILE *fp;
                fp = fopen("log_rank_input", "w");
                printMatrix(rank_input, fp);
                fclose(fp);
            }
            
            break;
        default:
            log_fatal("Algorithm %d doesn't exist --> Abort", config.algo);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    // free index_arr and cumulate_index_arr because they are no longer needed
    free_int_array(index_arr);
    free_int_array(cumulate_index_arr);

    //TODO: remove
    // to test if distribute works correctly 
    if (rank == TEST_RANK) {
        log_info("printing input_array:");
        FILE *fp;
        fp = fopen("log_input_array", "w");
        printArray(input_array, config.m, config.n, fp);
        fclose(fp);
    }

    // input no longer needed
    free_float_array(input_array);

    // SYRK:
    // Compute the result matrix for each node which gets
    //  summed up by the MPI_Reduce_scatter() methode 
    //  and store the result in rank_syrk_result
    floatArray rank_syrk_result;
    allocate_float_array(&rank_syrk_result, config.m, config.m);

    // Synchronize before starting time
    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();

    switch (config.algo) {
        case 0:
            syrkIterative(&config, rank, index_arr_rank, rank_input, rank_input_t, rank_syrk_result);
            break;
        case 1:
            improved_syrkIterative(&config, rank, index_arr_rank, rank_input.data, rank_input_t.data, rank_syrk_result.data);
            break;
        case 2:
            // 1D SYRK - OpenBLAS
            //In the case that m ≤ n and P is not too large
            syrk_withOpenBLAS(&config, rank, index_arr_rank, rank_input.data, rank_syrk_result.data);
            break;
        case 3:
            // 2D SYRK - OpenBLAS
            // In the case that m > n and P is not too large, a 2D algorithm is optimal
            assert(rank_input.data[0] != 0);
            assert(rank_syrk_result.data != NULL);
            two_d_syrk(&config, rank, rank_syrk_result, rank_input, MPI_COMM_WORLD);
            break;
        case 4:
            // 3D SYRK - OpenBLAS
            assert(comm_rank != -1);
            assert(pMpiCommunicators != NULL);
            assert(pMpiCommunicators[rank % config.P2] != MPI_COMM_NULL);
            three_d_syrk(&config, comm_rank, rank_syrk_result, rank_input, pMpiCommunicators[rank % config.P2]);
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
    allocate_float_array(&reduction_result, counts.data[rank], 1);

    // Save the time before the MPI_Reduce_scatter (Start the colock)
    double start_mpi_reduce_scatter = MPI_Wtime();

    // TODO: remove
    if(rank == TEST_RANK) {
        log_info("printing rank [%d] rank_syrk_result:", TEST_RANK);
        FILE *fp;
        fp = fopen("log_rank_syrk_result_2", "w");
        printArray(rank_syrk_result, config.m, config.m, fp);
        fclose(fp);
    }
    // TODO: remove
    // if(rank == 3) {
    //     log_info("printing rank [3] rank_syrk_result:");
    //     FILE *fp;
    //     fp = fopen("log_rank_syrk_result_1", "w");
    //     printArray(rank_syrk_result, config.m, config.m, fp);
    //     fclose(fp);
    // }
    // TODO: remove
    // if(rank == 5) {
    //     log_info("printing rank [5] rank_syrk_result:");
    //     FILE *fp;
    //     fp = fopen("log_rank_syrk_result_3", "w");
    //     printArray(rank_syrk_result, config.m, config.m, fp);
    //     fclose(fp);
    // }


    // Reduce the results of all processors and scatter the result to all processors
    int result = MPI_Reduce_scatter(
        rank_syrk_result.data,      // send buffer
        reduction_result.data,      // receive buffer
        counts.data,                // number of elements to receive from each processor
        MPI_FLOAT,                  // data type of the buffer
        MPI_SUM,                    // operation to perform
        MPI_COMM_WORLD              // communicator
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
        allocate_float_array(&buffer, config.m, config.m);

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

            printResult(&config, config.m, buffer);

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

    free_float_array(rank_syrk_result);
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

