#ifndef MPI_SYRK_IMPLEMENTATION_THREE_D_SYRK_H
#define MPI_SYRK_IMPLEMENTATION_THREE_D_SYRK_H

#include "MPI_Syrk_implementation.h"

void three_d_syrk(run_config *s, int rank, float *rank_result, float **input);

#endif //MPI_SYRK_IMPLEMENTATION_THREE_D_SYRK_H