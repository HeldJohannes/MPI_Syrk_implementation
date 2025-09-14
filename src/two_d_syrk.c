#include "two_d_syrk.h"
#define TEST_RANK -1
#define ROOT 0

#ifdef USE_CBLAS_64
    #define CBLAS_SYRK cblas_ssyrk64_
    #define CBLAS_DGEMM cblas_dgemm64_
#else
    #define CBLAS_SYRK cblas_ssyrk
    #define CBLAS_DGEMM cblas_dgemm
#endif

/**
 * computes (⌊k/c⌋(u − 1) + k) mod c + cu
 * @param k rank of the current processor (0 ≤ k < P)
 * @param u result array position (0 ≤ u < c)
 * @param c prime number such that P=c(c+1)
 */
int func_f(int k, int u, int c) {
    return (k / c * (u - 1) + k) % c + c * u;
}

/**
 * computes h_i(q) = (i − (⌊i/c⌋ − 1) q) mod c + cq
 * @param i (0 ≤ i < c^2)
 * @param q (0 ≤ q < c)
 * @param c prime number such that P=c(c+1)
 * @return the processor assigned bloc
 */
int func_h(int i, int q, int c) {
    return (i - (i / c - 1) * q) % c + c * q;
}

/**
 * specify the set of row block indices that defines the triangle block for a particular processor k
 * @param r_k row block indices array of size c
 * @param k rank of the current processor (0 ≤ k < P)
 */
void calculate_R_k(int *r_k, int k, int c) {
    if (k < c * c) {
        assert(r_k + 0 != NULL);
        r_k[0] = k / c;
        for (int i = 1; i < c; ++i) {
            assert(r_k + i != NULL);
            r_k[i] = func_f(k, i, c);
        }
    } else if (k < (c * c + c)) {
        for (int i = 0; i < c; ++i) {
            assert(r_k + i != NULL);
            r_k[i] = (k - c * c) * c + i;
        }
    }
}

/**
 * d_k defines the diagonal block owned by processor k (|d_k| ≤ 1)
 * @return the index of the diagonal block owned by processor k or -1 if there is none
 */
int calculate_D_k(int k, int c) {
    if (k < c) {
        return -1;
    } else if (k < c * c && k % c == 0) {
        //log_info("[rank == %d] k / c == %d", k, k/c);
        return k / c;
    } else if (k < c * c && k % c != 0) {
        //log_info("[rank == %d] func_f(k, (k / c), c) == %d", k, func_f(k, (k / c), c));
        return func_f(k, (k / c), c);
    }
    return func_f(c * (k - c * c), k - c * c, c);
}

/**
 * calculate the set of processors Q_i that are assigned to a particular row block i
 *
 * @param q_i array of size c+1
 * @param m_bloc_i (0 ≤ i < c^2)
 * @param c prime number such that P=c(c+1)
 * 
 * @return -1 if m_bloc_i is not in the range [0, c^2) or 0 otherwise
 */
int calculate_Q_i(intArray q_i, int m_bloc_i, int c) {

    // Überprüfen, ob 0 ≤ m_bloc_i < c^2 erfüllt ist
    if (!(m_bloc_i >= 0 && m_bloc_i < c * c)) {
        return -1; 
    }

    if (m_bloc_i < c) {
        for (int i = 0; i < c; ++i) {
            assert(q_i.data + i != NULL);
            q_i.data[i] = c * m_bloc_i + i;
        }
        assert(q_i.data + c != NULL);
        q_i.data[c] = c * c;
    } else {
        for (int i = 0; i < c; ++i) {
            assert(q_i.data + i != NULL);
            q_i.data[i] = func_h(m_bloc_i, i, c);
        }
        assert(q_i.data + c != NULL);
        q_i.data[c] = c * c + m_bloc_i / c;
    }
    return 0;
}

/**
 * calculate the block size of the matrix A
 * @param s run_config struct
 * @return the block size of the matrix A [m*n / c^2(c+1)]
 */
int cal_block_size(run_config *s) {
    return (s->m * s->n) / (s->c * s->c * (s->c + 1));
}

void copy_array(const double *src, double *des, int h, int l, int offset_src, int offset_des) {
    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < l; ++j) {
            des[(i * l) + j + offset_des] = src[(i * l) + j + offset_src];
        }
    }
}

/**
 * copies a block of values from a 2D array to a 1D array
 * @param B 1D destination array
 * @param A 2D source array 
 * @param k index of the block in B 
 * @param block_height height of the block
 * @param block_length length of the block
 * @param index block index in the 2D array A
 */
void copy_to_2D(float *B, float **A, int k, int block_height, int block_length, int index) {

    assert(B != NULL);
    assert(A != NULL);
    assert(A[index] != NULL);

    int block_size = (block_height * block_length);
    for (int i = 0; i < block_height; ++i) { //row 
        for (int j = 0; j < block_length; ++j) { //col
            
            int B_index = k * block_size + (i * block_length) + j;
            int A_index = i + (block_height * index);

            B[B_index] = A[A_index][j];
        }
    }
}

/**
 * copy the values of a 1D array to a 1D array at a specific position
 * @param dest destination array
 * @param source source array
 * @param r_pos read position in the 2D array source
 * @param w_pos write position in the 1D array dest
 * @param block_height height of the block; needs to be (m / (c * (c+1)))
 * @param block_length length of the block; needs to be n
 * @param row_length length of the rows in the 2D array source
 * @param w_offset offset in the write position
 */
void copy_to_1D(float *dest, float *source, int r_pos, int w_pos, int block_height, int block_length, int row_length, int w_offset) {
   
    // Pointer checks:
    assert(dest != NULL);
    assert(source != NULL);

    int block_size = block_height * block_length;

    assert(source + r_pos * block_size != NULL);
    assert(dest + w_pos * block_size != NULL);
    
    for (int i = 0; i < block_height; ++i) {
        for (int j = 0; j < block_length; ++j) {

            int A_i_index = w_pos * block_size + w_offset + (i * row_length) + j;
            int B_index = r_pos * block_size + (i * block_length) + j;

            assert(dest + A_i_index != NULL);
            assert(source + B_index != NULL);

            dest[A_i_index] = source[B_index];
        }
    }
}

/**
 * copy the values of a 1D array to a 1D array and cast values to double
 * @param A_i destination array
 * @param A source array
 * @param c prime number such that P=c(c+1)
 * @param block_height height of the block; needs to be (m / (c * (c+1)))
 * @param block_length length of the block; needs to be n
 * @param index block index in the 2D array A; needs to be in [0, c]
 * @param trans transpose the array
 */
void copy_to_d(double *A_i, float *A, int c, int block_height, int block_length, int index, _Bool trans) {

    assert(A_i != NULL);
    assert(A != NULL);

    int block_size = (block_height * block_length);
    int shift = index * block_size * (c + 1);
    if (!trans) {
        for (int i = 0; i < block_height; ++i) {
            for (int j = 0; j < c + 1; ++j) {
                for (int k = 0; k < block_length; ++k) {

                    int index_A_i = k + j * block_length + i * block_length * (c + 1);
                    int index_A = k + j * block_size + i * block_length + shift;

                    A_i[index_A_i] = A[index_A];
                }
            }
        }
    } else {
        for (int j = 0; j < c + 1; ++j) {
            for (int k = 0; k < block_length; ++k) {
                for (int i = 0; i < block_height; ++i) {

                    int index_A_i = i + k * block_height + j * block_size;
                    int index_A = i * block_length + k + j * block_size + shift;

                    A_i[index_A_i] = A[index_A];
                }
            }
        }
    }
}

/**
 * copy the values of a 1D array to a 1D array
 * @param A_i 1D destination array
 * @param A 1D source array
 * @param c prime number such that P=c(c+1)
 * @param block_height height of the block; needs to be (m / (c * (c+1)))
 * @param block_length length of the block; needs to be n
 * @param index block index in the 2D array A; needs to be in [0, c]
 * @param trans transpose the array
 */
void copy_to_f(float *A_i, float *A, int c, int block_height, int block_length, int index) {

    assert(A_i != NULL);
    assert(A != NULL);

    int block_size = (block_height * block_length);
    int shift = index  * (c + 1) * block_size;
    
    for (int i = 0; i < block_height; ++i) {            //rows
        for (int j = 0; j < c + 1; ++j) {               //cols
            for (int k = 0; k < block_length; ++k) {    //elements in the block

                int index_A_i = k + j * block_length + i * block_length * (c + 1);
                int index_A = k + j * block_size + i * block_length + shift;

                A_i[index_A_i] = A[index_A];
            }
        }
    }

}

/**
 * cast double array to float array
 * @param pDouble double array
 * @param pFloat float array
 * @param size size of the arrays
 */
void cast_d_to_f(double *pDouble, float *pFloat, int size) {

    assert(pDouble != NULL);
    assert(pFloat != NULL);

    for (int j = 0; j < size; ++j) {
        pFloat[j] = (float) pDouble[j];
    }
}

void accumulate_B_into_A(run_config *s, int k, floatArray A, floatArray B, intArray R_k, intArray Q_i, floatMatrix input) {
    
    int block_height = s->m / (s->c * s->c);
    //fprintf(stderr, "block_height = %d\n", block_height);
    int block_length = s->n / (s->c + 1);
    //fprintf(stderr, "block_length = %d\n", block_length);
    
    for (int i = 0; i < s->c; ++i) {
        // for each i ∈ R_k calculate Q_i:
        int ret = calculate_Q_i(Q_i, R_k.data[i], s->c);
        if (ret == -1) {
            log_fatal("calculate_Q_i failed for i == %d", i);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        for (int j = 0; j < s->c + 1; ++j) {
            // for each k' ∈ Q_i do:
            // Accumulate B_k' into A
            if (Q_i.data[j] != k) {

                int read_position = Q_i.data[j];
                int write_position = j + ((s->c + 1) * i);

                copy_to_1D(A.data, B.data, read_position, write_position, block_height, block_length, block_length, 0);
            } else {

                int write_position = j + ((s->c + 1) * i);

                copy_to_2D(A.data, input.data, write_position, block_height, block_length, i);
            }
        }
    }

    /*
    for (int i = 0; i < s->c; ++i) {
        // for each i ∈ R_k calculate Q_i:
        int ret = calculate_Q_i(Q_i, R_k.data[i], s->c);
        if (ret == -1) {
            log_fatal("calculate_Q_i failed for i == %d", i);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        for (int j = 0; j < s->c + 1; ++j) {
            // for each k' ∈ Q_i do:
            // Accumulate B_k' into A
            if (Q_i.data[j] != k) {
                int read_position = Q_i.data[j];
                int write_position = j + ((s->c + 1) * i);
                
                assert(read_position < s->world_size);
                assert(write_position < s->world_size);
                copy_to_1D(A.data, B.data, read_position, write_position, block_height, block_length);
            } else {
                int write_position = j + ((s->c + 1) * i);
                assert(write_position < s->world_size);
                copy_to_2D(A.data, input.data, write_position, block_height, block_length, i);
            }
        }
    }
    */
}

/**
 * <p> main function of this class </p>
 *
 * <p> <b>Require:</b> |Π|=P=c(c+1) for prime c </p>
 * <p> <b>Require:</b> A is evenly subdivided into c^2 row blocks, and each row
 *   block A_i is evenly divided across a set of c + 1 processors Q_i </p>
 *
 * @param s run_config struct containing the configuration of the run
 * @param k rank of the processor [0 ≤ k < world_size]
 * @param rank_result result array of size ...
 * @param input input matrix A_i 
 */
void two_d_syrk(run_config *s, int k, floatArray rank_result, floatMatrix input, MPI_Comm communicator) {
    
    assert(input.data != NULL);
    assert(rank_result.data != NULL);

    fprintf(stderr, "[rank %d] s->m = %d, s->n = %d, s->c = %d, s->world_size = %d\n", k, s->m, s->n, s->c, s->world_size);

    // block size is the number of elements in a block A_i^(k)
    int block_size = cal_block_size(s);
    // block height is the number of rows in a block A_i^(k)
    int block_height = s->m / (s->c * s->c);
    // block length is the number of cols in a block A_i^(k)
    int block_length = s->n / (s->c + 1);

    fprintf(stderr ,"block_size = %d, block_height = %d, block_length = %d\n", block_size, block_height, block_length);

    assert(block_size != 0);
    assert(block_height != 0);
    assert(block_length != 0);
    assert(block_size == block_height * block_length);

    //TODO remove:
    // print input matrix for a specific processor to test if correct
    //if (k == RANK) {
    //    FILE *fp;
    //    fp = fopen("log_input", "w");
    //    printMatrix(input, fp);
    //    fclose(fp);
    //}

    // R_k defines the row block indices that defines the triangle block for a particular processor k
    intArray R_k;
    allocate_int_array(&R_k, s->c);
    
    calculate_R_k(R_k.data, k, s->c);

    // Q_i defines the set of processors that are assigned to a particular row block i (A_i)
    intArray Q_i;
    allocate_int_array(&Q_i, s->c + 1);

    /** ************************************************************************************************
     * STEP 1: 
     * Gather c row blocks in row block set 
     ************************************************************************************************ */

    // Step 1.1:
    // allocate array B of P blocks, each of size block_size [m*n / c^2 (c+1)]
    floatArray B;
    allocate_float_array(&B, s->world_size, block_size);

    /** ************************************************************************************************
     * Step 2:
     * copy the input matrix block A^(k) to position B_k\'
     ************************************************************************************************ */
    for (int i = 0; i < s->c; ++i) {
        // for each i ∈ R_k do:
        calculate_Q_i(Q_i, R_k.data[i], s->c);
        for (int j = 0; j < s->c + 1; ++j) {
            // for each j ∈ Q_i \{k} do:
            // (ignore the own block)
            if (Q_i.data[j] != k) {
                // copy the blocs A_i^(k) to B_k'
                int copy_to = Q_i.data[j];// needs to be in [0, c * (c+1)]
                int copy_from = i; // needs to be in [0, c]
                copy_to_2D(B.data, input.data, copy_to, block_height, block_length, copy_from);
            }
        }
    }

    // TODO: remove 
    // print the matrix B for a specific processor after copy to test if correct
    if (k == TEST_RANK) {
        FILE *fp;
        fp = fopen("log_B", "w");
        printArray(B, s->world_size, block_size, fp);
        fclose(fp);
    }

    assert(B.data != NULL);
    assert(B.data + block_size * (s->world_size +1) -1 != NULL);

    /** ************************************************************************************************
     * Step 3:
     * Communicate B ALL-TO-ALL
     ************************************************************************************************ */
    MPI_Alltoall(
        B.data,         // send buffer
        block_size,     // send count
        MPI_FLOAT,      // send datatype
        B.data,         // receive buffer
        block_size,     // receive count
        MPI_FLOAT,      // receive datatype
        communicator   // communicator
    );

    // TODO remove:
    // print the matrix B after ALLtoALL:
    if (k == TEST_RANK) {
        FILE *fp;
        fp = fopen("log_B_after_ATA", "w");
        printArray(B, s->world_size, block_size, fp);
        fclose(fp);
    }

    floatArray A;
    allocate_float_array(&A, s->world_size, block_size);

    /** ************************************************************************************************
     * STEP 4:
     * Accumulate B_k′ into A
     ****************************************** */
    // R_k contains the row block indices that defines the triangle block for a particular processor k 
    // and has size c

    accumulate_B_into_A(s, k, A, B, R_k, Q_i, input);

    //TODO remove:
    // print A after the accumulation
    // Attention after the accumulation the blocks are row wise in A_i
    if (k == TEST_RANK) {
        FILE *fp;
        fp = fopen("log_A_i", "w");
        printArray(A, s->world_size, block_size, fp);
        fclose(fp);
    }
    

    /** ************************************************************************************************
     * STEP 5: 
     * Compute c(c − 1)/2 off-diagonal blocks
     * 
     * To do this, we need to:
     * compute C_ij= Local-GEMM(Ai , Aj^T) for each (i, j) ∈ Rk with i > j
     * C_ij is a block of size block_height x block_height
     ************************************************************************************************* */

    // Temporary matrix to hold the result of DGEMM which is in double
    doubleArray result;
    allocate_double_array(&result, block_height * block_height);
    
    // Temporary matrix to hold the result of DGEMM as float
    // this is needed to work with the rest of the code 
    floatArray result_f;
    allocate_float_array(&result_f, block_height, block_height);
    
    // Temporary matrix A_i of size block_height x n 
    // to store the values of A_i as double
    doubleArray A_i;
    allocate_double_array(&A_i, block_height * s->n);
    
    // Temporary matrix A_j of size n x block_height 
    // to store the values of A_j as double
    doubleArray A_j;
    allocate_double_array(&A_j, block_height * s->n);

    for (int i = 0; i < s->c; ++i) {
        for (int j = 0; j < s->c; ++j) {
            if (R_k.data[i] > R_k.data[j]) {

                /** ************************************************************************************************
                 * STEP 5.1:
                 * copy  the required index values of A to A_i and A_j 
                 * and covert them from float to double
                 ************************************************************************************************  */ 
                copy_to_d(A_i.data, A.data, s->c, block_height, block_length, i, false);
                if (k == TEST_RANK) {
                    FILE *fp;
                    fp = fopen("log_A_i_tmp", "w");
                    printDoubleArray(A_i, block_height, s->n, fp);
                    fclose(fp);
                }
                copy_to_d(A_j.data, A.data, s->c, block_height, block_length, j, true);
                if (k == TEST_RANK) {
                    FILE *fp;
                    fp = fopen("log_A_j_tmp", "w");
                    printDoubleArray(A_j, block_height, s->n, fp);
                    fclose(fp);
                }

                /** ************************************************************************************************
                 * STEP 5.2:
                 * Compute C_ij= Local-GEMM(Ai , Aj^T) for each (i, j) ∈ Rk with i > j
                 * 
                 * DGEMM := C = αAB + βC
                 ************************************************************************************************ */
                
                 CBLAS_DGEMM(
                        CblasRowMajor,          // matrix layout, Specifies row-major (C) or column-major (Fortran) data ordering.
                        CblasNoTrans,           // transA, Specifies whether to transpose matrix A.
                        CblasNoTrans,           // transB, Specifies whether to transpose matrix B
                        block_height,           // m, Number of rows in matrices A and C
                        block_height,           // n, Number of columns in matrices B and C
                        s->n,                   // k, Number of columns in matrix A; number of rows in matrix B
                        1.0,                    // alpha, Scaling factor for the product of matrices A and B
                        A_i.data,               // Matrix A
                        s->n,                   // lda, The size of the first dimension of matrix A; if you are passing a matrix A[m][n], the value should be ...
                        A_j.data,               // Matrix B
                        block_height,           // ldb, The size of the first dimension of matrix B; Number of Elements between rows; if you are passing a matrix B[m][n], the value should be n
                        0.0,                    // beta, Scaling factor for matrix C
                        result.data,            // Result matrix C
                        block_height            // ldc, The size of the first dimension of matrix C; if you are passing a matrix C[m][n], the value should be ...
                );
                
                /**
                 * STEP 5.3:
                 * write result to C_ij
                 */

                //TODO: remove
                // print result
                if (k == TEST_RANK) {
                    FILE *fp;
                    fp = fopen("log_result_tmp", "w");
                    printDoubleArray(result, block_height, block_height, fp);
                    fclose(fp);
                }
                
                // cast the result to float
                cast_d_to_f(result.data, result_f.data, result.length);
                // copy the result to the correct position in the rank_result array

                copy_to_1D(
                    rank_result.data, 
                    result_f.data, 
                    0, 
                    R_k.data[i] * s->c * s->c, 
                    block_height,
                    block_height, 
                    s->m, 
                    R_k.data[j] * block_height
                );
            }
        }
    }

    //TODO: remove
    // print result
    // if (k == RANK) {
    //     FILE *fp;
    //     fp = fopen("log_result", "w");
    //     printArray(rank_result, s->m, s->m, fp);
    //     fclose(fp);
    // }

    // the tmporary arrays are no longer needed 
    // because the values are copied to the rank_result array 
    // and can be freed
    free_double_array(result);
    
    free_float_array(result_f);
    
    free_double_array(A_i);
    
    free_double_array(A_j);

    /** ************************************************************************************************
     * STEP 6:
     * Compute diagonal block if assigned
     * Not all processors have a diagonal block assigned to them
     ************************************************************************************************ */

    // for each i ∈ D_k do
    int d_k = calculate_D_k(k, s->c);

    if (d_k != -1) {

        // Temporary array to store the result of SYRK for the diagonal block
        floatArray result_D_k;
        allocate_float_array(&result_D_k, block_height, block_height);
        
        // Temporary array to store the values of A_i for the diagonal block as float
        floatArray A_i_D_k;
        allocate_float_array(&A_i_D_k, block_height, s->n);

        for (int i = 0; i < s->c; ++i) {
            if (R_k.data[i] == d_k) {
                //log_info("[rank == %d] i == %d", k, i);
                copy_to_f(A_i_D_k.data, A.data, s->c, block_height, block_length, i);
                
                //C_ii = Local-SYRK(A_i)
                CBLAS_SYRK(
                    CblasRowMajor,
                    CblasLower,
                    CblasNoTrans,
                    block_height,
                    s->n,
                    1.0f,
                    A_i_D_k.data,
                    s->n,
                    0.0f,
                    result_D_k.data,
                    block_height
                );
                copy_to_1D(
                    rank_result.data, 
                    result_D_k.data, 
                    0, 
                    d_k * s->c * s->c, 
                    block_height,
                    block_height, 
                    s->m, 
                    d_k * block_height
                );
            }
        }

        // TODO: remove
        if (k == TEST_RANK) {
            FILE *fp;
            fp = fopen("log_A_i_D_k", "w");
            printArray(A_i_D_k, block_height, s->n, fp);
            fclose(fp);
        }


        // TODO: remove
        // print rank_result
        if (k == TEST_RANK) {
            FILE *fp;
            fp = fopen("log_result_D_k", "w");
            printArray(result_D_k, block_height, block_height, fp);
            fclose(fp);
        }

        // free the no longer neede temporary arrays 
        // because the result is copied to the rank_result array
        free_float_array(result_D_k);
    
        free_float_array(A_i_D_k);
        
    }

    // cleanup: free the no longer needed arrays
    // the computation is done and the results are stored in rank_result
    free_int_array(R_k);

    free_int_array(Q_i);
    
    free_float_array(B);
    
    free_float_array(A);
}

bool includes(intArray array, int value) {
    for (int i = 0; i < array.length; ++i) {
        if (array.data[i] == value) {
            return true;
        }
    }
    return false;
}

/**
 * Methode to distribute the data of an input matrix.
 * 
 * @param s run_config struct containing the configuration of the run
 * @param rank rank of the processor
 * @param input_array input to be distributed
 * @param rank_input input for each rank after the distribution
 * @param communicator MPI_Comm to distribut the data over
 */
void distribute_input_matrix_2D(run_config *s, int rank, floatArray input_array, floatMatrix rank_input, MPI_Comm communicator) {
    if (rank == ROOT) log_debug("Starting 2D-Algo array distribution");

    // If algo == 3 (2D-Algo), we need to split the input matrix A into c^2 blocks (A_i) 
    // and distribute them to the processors Q_i
    // This is done by creating new communicators based on Q_i and using MPI_Scatter to split A_i as A_i^(k)

    //TODO find a better solution:
    // broadcast the input matrix to all processors
    MPI_Bcast(
        input_array.data,       // buffer to broadcast
        input_array.length,     // number of elements in the buffer
        MPI_FLOAT,              // data type of the buffer
        0,                      // root process
        communicator          // communicator
    );

    if (rank == TEST_RANK) {
        FILE *fp;
        fp = fopen("log_input_array", "w");
        printArray(input_array, s->m, s->n, fp);
        fclose(fp);
    }
    

    if (rank == ROOT) log_debug("Successfully broadcasted the input matrix");

    /** ************************************************************************************************
     * STEP 2.1:
     * Create new Communicators Based on Q_i to take advantage of MPI_SCATTER to split A_i as A_i^(k).
     * note: this step might not be necessary but I couldn't find a better current solution
     ************************************************************************************************ */
    
    // create a main group for the MPI_COMM_WORLD
    MPI_Group main_group;
    MPI_Comm_group(communicator, &main_group);

    // create c^2 mpi - groups:
    MPI_Group *pMpiGroups = (MPI_Group *) malloc(s->c * s->c * sizeof (MPI_Group ));
    if (!pMpiGroups) {
        log_fatal("Memory allocation failed for pMpiGroups", 0);
        MPI_Abort(communicator, EXIT_FAILURE);
    }

    // create c^2 mpi - communicators:
    MPI_Comm *pMpiCommunicators = (MPI_Comm *) malloc(s->c * s->c * sizeof (MPI_Comm ));
    if (!pMpiCommunicators) {
        log_fatal("Memory allocation failed for pMpiCommunicators", 0);
        MPI_Abort(communicator, EXIT_FAILURE);
    }

    // create an array to store the Q_i values which are used to create the new communicators
    intArray Q_i;
    allocate_int_array(&Q_i, s->c + 1);
        
            
    // create the new groups and communicators
    for (int i = 0; i < s->c * s->c; ++i) {
        // calculate the Q_i values
        calculate_Q_i(Q_i, i, s->c);
        // create the new group
        int ret = MPI_Group_incl(
            main_group,         // parent group
            s->c + 1,           // number of nodes in new group
            Q_i.data,           // ranks of the processes in the new group
            &pMpiGroups[i]      // new group
        );
        // check if the MPI_Group_incl was successful
        assert(pMpiGroups[i] != NULL);
        if (ret != MPI_SUCCESS) {
            log_fatal("MPI_Group_incl failed for group %d", i);
            MPI_Abort(communicator, ret);
        }
        // create the new communicator
        int err = MPI_Comm_create(
            communicator,         // parent communicator
            pMpiGroups[i],          // group
            &pMpiCommunicators[i]   // new communicator
        );
        // check if the MPI_Comm_create was successful
        if (includes(Q_i, rank) && err == MPI_SUCCESS) {
            assert(pMpiCommunicators[i] != MPI_COMM_NULL);
        } else {
            assert(pMpiCommunicators[i] == MPI_COMM_NULL);
        }
        
    }

    /** ************************************************************************************************
     * STEP 2.2:
     *
     * Spilt the input:
     ************************************************************************************************ */

    assert(input_array.data != NULL);
    assert(rank_input.data != NULL);

    int row_block_height = s->m / (s->c * s->c);
    int pos = 0;
    for (int i = 0; i < s->c * s->c; ++i) {
        calculate_Q_i(Q_i, i, s->c);
        // if Q_i includes the current rank otherwise continue
        bool is_included = includes(Q_i, rank);
        if (is_included) {
            //TODO comunicate the input matrix to the processors
            for (int j = 0; j < row_block_height; j++)
            {
                int index = i * row_block_height * s->n + j * s->n;
                int count = s->n / (s->c + 1);
                MPI_Scatter(
                    input_array.data +index,    // send buffer
                    count,                      // send count
                    MPI_FLOAT,                  // send datatype
                    rank_input.data[pos*row_block_height + j],          // receive buffer
                    count,                      // receive count
                    MPI_FLOAT,                  // receive datatype
                    0,                          // root process
                    pMpiCommunicators[i]        // communicator
                );
            }
            pos++;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    /*int row_block_length = s->n / (s->c + 1);
    

    int pos = 0;
    // Distribute the rows of A along the Q_i processors
    for (int i = 0; i < s->c * s->c; ++i) {
        calculate_Q_i(Q_i, i, s->c);
        // if Q_i includes the current rank otherwise continue
        bool is_included = false;
        for (int i = 0; i < s->c + 1; i++)
        {
            if (Q_i[i] == rank)
            {
                is_included = true;
                break;
            }
        }
        if (is_included) {
            for (int j = 0; j < row_block_height; ++j) {
                for (int k = 0; k < s->n; k++)
                {
                    // calculate the index of the input array
                    // where j [0, row_block_height-1] is the row index of the rank_input array 
                    // and pos [0, c] is the position of the block
                    int rank_input_row = j + pos * row_block_height;
                    assert(pos < s->c);
                    assert(rank_input_row < s->m);

                    int input_array_row_index = j*s->n + i*row_block_height * s->n;
                    int input_array_index = input_array_row_index + k;
                    assert(input_array_index < s->m * s->n);


                    rank_input[j + pos*row_block_height][k] = input_array[input_array_index];
                }
            }
            pos++;
        }
    }
        */

    free_int_array(Q_i);

    // if (rank == RANK) {
    //     FILE *fp;
    //     fp = fopen("log_rank_input", "w");
    //     printMatrix(rank_input, fp);
    //     fclose(fp);
    // }
    
    if (rank == ROOT) log_info("Successfully distributed the input matrix");
    MPI_Barrier(communicator);
}