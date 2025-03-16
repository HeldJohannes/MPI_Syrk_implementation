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


typedef struct {
    double* data;
    int length;
} doubleArray;

typedef struct {
    int* data;
    int length;
} intArray;

typedef struct {
    float** data;
    int length;
    int rows;
    int cols;
} floatMatrix;

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