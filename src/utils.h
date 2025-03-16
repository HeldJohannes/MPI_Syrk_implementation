#ifndef UTILS_H
#define UTILS_H

#ifdef __cplusplus
extern "C" {
#endif


typedef struct {
    float* data;
    int length;  
} floatArray;

/**
 * Allocates memory for a float array.
 * 
 * @param array The floatArray struct to be allocated.
 * @param length The length of the array.
 */
void allocate_float_array(floatArray *array, int length);

void free_float_array(floatArray array);

typedef struct {
    double* data;
    int length;
} doubleArray;

void allocate_double_array(doubleArray *array, int length);

void free_double_array(doubleArray array);

typedef struct {
    int* data;
    int length;
} intArray;

void allocate_int_array(intArray *array, int length);

void free_int_array(intArray array);

typedef struct {
    float** data;
    int length;
    int rows;
    int cols;
} floatMatrix;

void allocate_float_matrix(floatMatrix *matrix, int rows, int cols);

void free_float_matrix(floatMatrix *matrix);

typedef struct {
    double** data;
    int length;
    int rows;
    int cols;
} doubleMatrix;

#ifdef __cplusplus
}
#endif

#endif //UTILS_H