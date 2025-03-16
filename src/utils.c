#include "utils.h"
#include "log.h"
#include <stdlib.h>
#include <mpi.h>

void allocate_float_array(floatArray *array, int length) {
    array->length = length;
    array->data = (float *) calloc(array->length, sizeof(float));
    if (array->data == NULL) {
        log_fatal("Memory allocation failed for input");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
}

void free_float_array(floatArray array) {
    log_debug("Freeing array %p", array.data);
    free(array.data);
    log_debug("Successfully freed -> float array...");
}

void allocate_int_array(intArray *array, int length) {
    array->length = length;
    array->data = (int *) calloc(array->length, sizeof(int));
    if (array->data == NULL) {
        log_fatal("Memory allocation failed for input");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
}

void free_int_array(intArray array) {
    log_debug("Freeing array %p", array.data);
    free(array.data);
    log_debug("Successfully freed -> int array...");
}

void allocate_double_array(doubleArray *array, int length) {
    array->length = length;
    array->data = (double *) calloc(array->length, sizeof(double));
    if (array->data == NULL) {
        log_fatal("Memory allocation failed for input");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
}

void free_double_array(doubleArray array) {
    log_debug("Freeing array %p", array.data);
    free(array.data);
}

void allocate_float_matrix(floatMatrix *matrix, int rows, int cols) {

    matrix->length = rows * cols;
    matrix->rows = rows;
    matrix->cols = cols;
    matrix->data = (float **) calloc(matrix->rows, sizeof(float *));
    if (matrix->data == NULL) {
        log_fatal("Memory allocation failed for input");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    for (int i = 0; i < matrix->rows; ++i) {
        matrix->data[i] = (float *) calloc(matrix->cols, sizeof(float));
        if (matrix->data[i] == NULL) {
            log_fatal("Memory allocation failed for input");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
    }
}

void free_float_matrix(floatMatrix *matrix) {
    log_debug("Freeing matrix %p", matrix);
    for (int i = 0; i < matrix->rows; ++i) {
        free(matrix->data[i]);
    }
    free(matrix->data);
}