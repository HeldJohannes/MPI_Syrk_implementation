//
// Created by Johannes Held on 04.10.23.
//

#ifndef MPI_SYRK_IMPLEMENTATION_MPI_SYRK_IMPLEMENTATION_H
#define MPI_SYRK_IMPLEMENTATION_MPI_SYRK_IMPLEMENTATION_H

#ifdef __cplusplus
extern "C" {
#endif

#include <assert.h>
#include <mpi.h>
#include <stdlib.h>
#include <cblas.h>
#include "log.h"
#include "utils.h"

/**
 * Structure to store the configuration values for the program.
 */
typedef struct {
    int algo;           // 0 for tripple for-loop, 1 for improved tripple for-loop, 2 for 1D, 3 for 2D and 4 for 3D
    int world_size;     // total number of processors
    int m;              // total number of rows
    int n;              // total number of columns
    int c;              // c is a prime number (e.g. 3, so that P = c(c+1) = 12) and is required for the 2D and 3D algorithm
    int P2;             // p2 is required for the 3D algorithm |Π|= P * P2
    char *fileName;     // name of the input file if not provided random input will be generated
    char *result_File;  // name of the output file
    _Bool print_result; // flag to print the result
} run_config;

/**
 * read_input function parses command-line arguments and populates a run_config structure.
 *
 * @param s A pointer to a run_config structure where the parsed values will be stored.
 * @param argc The number of command-line arguments.
 * @param argv An array of strings containing the command-line arguments.
 * @return
 *  EXIT_SUCCESS: if the parsing is successful, otherwise an error code.
 */
int read_input(run_config *s, int argc, char* argv[]);

/**
 * read_input_file function reads the input file and populates the array A with its contents.
 *
 * @param s A pointer to A run_config structure containing the file name.
 * @param A A pointer to a float array where the input values will be stored.
 * @param rank The rank of the current processor
 * @return EXIT_SUCCESS if the file is read successfully, otherwise an error code.
 */
int read_input_file(const int rank, run_config *s, float **A);

/**
 * This function is used to print error messages and exit the program.
 *
 * @param rank Rank of the caller
 * @param name Name of the program or function where the error occurred
 * @param msg Format string for the error message
 * @param ... Additional arguments for the error message (variable argument list)
 */
void error_exit(int rank, char *name, const char *msg, ...);

int parseInput(run_config *s, int argc, char **argv, int rank);

void printResult(run_config *s, int cols, floatArray array);

/**
 * This function prints an array to the specified file.
 *
 * @param array The array to be printed
 * @param row The number of rows in the array
 * @param cols The number of columns in the array
 * @param file The file where the array will be printed
 */
void printArray(floatArray array, int row, int cols, FILE *file);

void printDoubleArray(doubleArray array, int row, int cols, FILE *file);

void printMatrix(floatMatrix matrix, FILE *file);

/**
 * This function calculates the block size for the input array.
 */
void index_calculation(intArray arr, long n, int p);

/**
 * This function reads the input file and populates the input array.
 *
 * @param input A pointer to an integer array where the input values will be stored.
 * @param rank The rank of the current processor
 * @param argv An array of strings containing the command-line arguments.
 */
void readInputFile(int *input, int rank, char **argv);

void computeInputAndTransposed(run_config *s, int rank, int index_arr_rank, int cum_index_arr_rank, float **input, float **rank_input, float **rank_input_t);

/**
 * Transposes a matrix.
 * @param m Number of rows in the matrix
 * @param n Number of columns in the matrix
 * @param matrix The matrix to be transposed
 * @param result The transposed matrix
 */
void transposeMatrix(long m, long n, float** matrix, float** result);

/**
 * Generates random input values for the matrix A.
 *
 * @param s A pointer to a run_config structure containing the matrix dimensions.
 * @param A A pointer to a float array where the input values will be stored.
 */
void generate_input(run_config *s, float **A);

/**
 * Prints usage information for the program.
 *
 * @param prog_name Name of the program
 */
void print_usage(char *prog_name);

/**
 * This function is used to print error messages and exit the program.
 *
 * @param rank Rank of the caller
 * @param name Name of the program or function where the error occurred
 * @param msg Format string for the error message
 * @param ... Additional arguments for the error message (variable argument list)
 */
void error_exit(int rank, char *name, const char *msg, ...);

#ifdef __cplusplus
}
#endif

#endif //MPI_SYRK_IMPLEMENTATION_MPI_SYRK_IMPLEMENTATION_H
