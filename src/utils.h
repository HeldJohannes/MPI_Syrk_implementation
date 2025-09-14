#ifndef UTILS_H
#define UTILS_H

#ifdef __cplusplus
extern "C" {
#endif

#define assert__(x) for ( ; !(x) ; assert(x) )

typedef struct {
    float* data;        // pointer to the float data
    int length;         // total length of the array
    int row_length;     // length of each row (for 2D representation)
} floatArray;

/**
 * Allocates memory for a float array. Total length = rows * row_length
 * 
 * @param array The floatArray struct to be allocated.
 * @param rows Number of rows in the array.
 * @param row_length Length of each row in the array.
 */
void allocate_float_array(floatArray *array, int rows, int row_length);

void free_float_array(floatArray array);

typedef struct {
    double* data;   // pointer to the double data
    int length;     // total length of the array
} doubleArray;

/**
 * Allocates memory for a double array. Total length = length
 * 
 * @param array The doubleArray struct to be allocated.
 * @param length Total length of the array.
 */
void allocate_double_array(doubleArray *array, int length);

void free_double_array(doubleArray array);

typedef struct {
    int* data;      // pointer to the int data
    int length;     // total length of the array
} intArray;

/**
 * Allocates memory for an int array. Total length = length
 * 
 * @param array The intArray struct to be allocated.
 * @param length Total length of the array.
 */
void allocate_int_array(intArray *array, int length);

void free_int_array(intArray array);

typedef struct {
    float** data;   // pointer to the float data
    int length;     // total length of the matrix (1D representation)
    int rows;       // number of rows
    int cols;       // number of columns
} floatMatrix;

/**
 * Allocates memory for a float matrix with given rows and columns.
 * 
 * @param matrix The floatMatrix struct to be allocated.
 * @param rows Number of rows in the matrix.
 * @param cols Number of columns in the matrix.
 */
void allocate_float_matrix(floatMatrix *matrix, int rows, int cols);

void free_float_matrix(floatMatrix *matrix);

typedef struct {
    double** data;  // pointer to the double data
    int length;     // total length of the matrix (1D representation)
    int rows;       // number of rows
    int cols;       // number of columns
} doubleMatrix;

typedef struct {
    int** data;     // pointer to the int data
    int length;     // total length of the matrix (1D representation)
    int rows;       // number of rows
    int cols;       // number of columns
} intMatrix;

/**
 * Allocates memory for an int matrix with given rows and columns.
 * 
 * @param matrix The intMatrix struct to be allocated.
 * @param rows Number of rows in the matrix.
 * @param cols Number of columns in the matrix.
 */
void allocate_int_matrix(intMatrix *matrix, int rows, int cols);

void free_int_matrix(intMatrix *matrix);

#ifdef __cplusplus
}
#endif

#endif //UTILS_H