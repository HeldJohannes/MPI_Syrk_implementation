#include "three_d_syrk.h"


/**
 * @brief This function is used to compute the result of the SYRK operation for a given rank.
 * 
 * @param s A pointer to a run_config structure containing the configuration values.
 * @param rank The rank of the current processor.
 * @param rank_result A pointer to a float array where the result will be stored.
 * @param input A pointer to a 2D float array containing the input values.
 * 
 * @see run_config
 * @see two_d_syrk
 * @see MPI_Syrk_implementation
 */
void three_d_syrk(run_config *s, int rank, float *rank_result, float **input) {

    // requires the number of processors to be |π| = p1 * p2
    // where p1 = c * (c + 1) and c is a prime number

    /* ****************************************
    STEP 1: Split the rank into two parts
    ******************************************/ 

    // in order to work the rank has to be split into 2 parts:
    // the first part k which is given by the modulo of the rank
    // k (0 <= k < p1) and p1 = c * (c + 1)
    int p1 = s->c * (s->c + 1);
    int k = rank % p1;

    // and into the second part l which is given by the division of the rank
    // l (0 <= l < p2)
    int l = rank / p1;


    /* ****************************************
    STEP 2: Compute intermediate results C_kl using the two_d_syrk function on the input slice A_*l and π_*l
    ******************************************/ 

    // input Blocks have size (m / c^2) x (n / p2)

    // TODO: implement the three_d_syrk function

    /* ****************************************
    STEP 3: Compute the final result C_kl by summing up the intermediate results C_kl 
    unsing REDUCE-SCATTER on C_kl and π_k*
    ******************************************/


    log_trace("rank_result %p", rank_result);
}
