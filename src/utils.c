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