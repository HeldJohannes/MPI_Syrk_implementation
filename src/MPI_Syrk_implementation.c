#include "MPI_Syrk_implementation.h"
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>
#include <float.h>
#include "log.h"
#include <getopt.h>

#define ROOT 0

int parseInput(run_config *s, int argc, char **argv, int rank) {

    if (argc <= 5) {
        error_exit(rank, argv[0], "To many or not enough input variables!");
    }

    log_trace("Enter parseInput");
    //total_row_number
    s->m = -1;
    //total_col_number
    s->n = -1;

    // define the long options
    static struct option long_options[] = {
        // option, has_arg, flag, val
        {"algorithm", required_argument, 0, 'a'},
        {"rows", required_argument, 0, 'm'},
        {"columns", required_argument, 0, 'n'},
        {"output", required_argument, 0, 'o'},
        {"c", required_argument, 0, 'c'},
        {"print-result", no_argument, 0, 'p'},
        {"prozessor-2", required_argument, 0, 'i'}
        {0, 0, 0, 0}
    };

    int opt;
    char *end;
    while ((opt = getopt_long(argc, argv, "a:m:n:o:c:i:", long_options, NULL)) != -1) {
        switch (opt) {
            case 'a':
                s->algo = (int) strtol(optarg, &end, 10);
                break;
            case 'm':
                s->m = (int) strtol(optarg, &end, 10);
                break;
            case 'n':
                s->n = (int) strtol(optarg, &end, 10);
                break;
            case 'o':
                s->result_File = optarg;
                break;
            case 'c':
                s->c = (int) strtol(optarg, &end, 10);
                break;
            case 'p':
                s->print_result = true;
                break;
            case 'i':
                s->P2 = (int) strtol(optarg, &end, 10);
                break;
            default:
            case '?':
                if (rank == ROOT) {
                    fprintf(stderr, "wrong usage: option %c doesn't exist", optopt);
                    fprintf(stderr, "Usage: %s -m <rows> -n <columns> <input_file>\n", argv[0]);
                }
                //print_usage(argv[0]);
                return EXIT_FAILURE; 
        }
    }

    log_trace("m = %d; n = %d", s->m, s->n);

    if (s->m == -1)
    {
        if (rank == ROOT) {
            fprintf(stderr, "missing parameter m\n");
            fprintf(stderr, "Usage: %s -m <rows> -n <columns> <input_file>\n", argv[0]);
        }
        return EXIT_FAILURE;
    }
    if (s->n == -1)
    {
        if (rank == ROOT) {
            fprintf(stderr, "missing parameter n\n");
            fprintf(stderr, "Usage: %s -m <rows> -n <columns> <input_file>\n", argv[0]);
        }
        return EXIT_FAILURE;
    }
    if (s->m <= 0 || s->n <= 0) {
        if (rank == ROOT) {
            fprintf(stderr, "parameters m (= %d) and n (= %d) have to be bigger than 0\n", s->m, s->n);
            fprintf(stderr, "Usage: %s -m <rows> -n <columns> <input_file>\n", argv[0]);
        }
        return EXIT_FAILURE;
    }

    log_trace("optind = %d, argc = %d", optind, argc);
    if (optind < argc) {
        s->fileName = argv[optind];
    } else {
        if (rank == ROOT) {
            fprintf(stderr, "missing input file name --> Generate random input\n");
        }
    }
    log_trace("Exit parseInput");
    return EXIT_SUCCESS;
}


void printResult(run_config *s, int cols, floatArray array) {
    FILE *file;
    if (s->result_File != NULL) {
        log_info("Printing result to set file = %s", s->result_File);
        file = fopen(s->result_File, "w");
    } else {
        log_info("Printing result to default file = result.csv");
        file = fopen("syrk_result.csv", "w");
    }

    printArray(array, cols, cols, file);

    log_debug("Finished printResults()");
}

void printArray(floatArray array, int row, int cols, FILE *file) {
    log_info("Printing float array to file");
    assert(array.data != NULL);
    assert(array.length == row * cols);
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < cols; ++j) {
            //log_debug("array[%d][%d] = %f", i, j, array[i * cols + j]);
            fprintf(file, j == cols - 1 ? "%0.0f" : "%0.0f; ", array.data[i * cols + j]);
        }
        fprintf(file, "\n");
    }
}

void printDoubleArray(doubleArray array, int row, int cols, FILE *file) {
    log_info("Printing float array to file");
    assert(array.data != NULL);
    assert(array.length == row * cols);
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < cols; ++j) {
            fprintf(file, j == cols - 1 ? "%0.0lf" : "%0.0lf; ", array.data[i * cols + j]);
        }
        fprintf(file, "\n");
    }
}

void printMatrix(floatMatrix matrix, FILE *file) {
    log_info("Printing matrix to file");
    assert(matrix.data != NULL);
    for (int i = 0; i < matrix.rows; ++i) {
        for (int j = 0; j < matrix.cols; ++j) {
            if (j == matrix.cols - 1) {
                fprintf(file, "%0.0f", matrix.data[i][j]);
            } else {
                fprintf(file, "%0.0f; ", matrix.data[i][j]);
            }

        }
        fprintf(file, "\n");
    }
}

void index_calculation(intArray arr, long n, int p) {
    long input_size = n / p;

    long rest = n % p;
    log_trace("rest = %d", rest);

    for (int i = 0; i < p; ++i) {
        assert(arr.data + i != NULL);
        arr.data[i] = (int) input_size;
        if (i < rest) {
            arr.data[i] += 1;
        }
    }
}

int read_input_file(const int rank, run_config *s, float **A) {
    //    MPI_File mpiFile;
    //    if (MPI_File_open(MPI_COMM_WORLD, s->fileName, MPI_MODE_RDONLY, MPI_INFO_NULL, &mpiFile)) {
    //        printf("[MPI process %d] Failure in opening the file.\n", rank);
    //        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    //    }
    //
    //    MPI_Offset filesize;
    //    MPI_File_get_size(mpiFile, &filesize);
    //    printf("[MPI process %d] File size == %lld\n", rank, filesize);
    //
    //    MPI_File_read(mpiFile, )
    //
    //
    //    MPI_File_close(&mpiFile);
    
        FILE *stream = fopen(s->fileName, "r");
    
        if (stream == NULL) {
            log_fatal("[MPI process %d] Failure in opening the file.", rank);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
    
        char *line = NULL;
        size_t len = 0;
        int row = 0;
        int col = 0;
        while (1) {
            ssize_t n = getline(&line, &len, stream);
            log_debug("ssize_t n = %d", n);
            if (n == -1) break;
            char *subtoken = strtok(line, ";");
            while (subtoken) {
                char *pEnd;
                float res = strtof(subtoken, &pEnd);
                if (res > FLT_MIN || res < FLT_MAX) {
                    A[row][col] = res;
                } else {
                    log_error("Input is not an float --> Abort");
                    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
                }
                //log_trace("[MPI process %d] input[%d] = %d", rank, row, A[row]);
                col++;
                subtoken = strtok(NULL, ";");
            }
            col = 0;
            row++;
        }
    
        fclose(stream);
        return EXIT_SUCCESS;
    }
    

void computeInputAndTransposed(run_config *s, int rank, int index_arr_rank, int cum_index_arr_rank, float **input,
    float **rank_input,
    float **rank_input_t) {
    log_trace("[rank %d] computeInputAndTransposed()", rank);
    for (int row = 0; row < s->m; row++) {
        for (int col = 0; col < index_arr_rank; ++col) {
            float tmp = input[row][cum_index_arr_rank + col];
            rank_input[row][col] = tmp;
            log_debug("&input[(rank_count) + row_count * s->n] = %p => %f", &tmp, tmp);
            log_debug("rank_input[row_count] = %p => %f", rank_input[row], *rank_input[row]);
        }
    }
    log_debug("for loop success");

    transposeMatrix(s->m, index_arr_rank, rank_input, rank_input_t);
}

void transposeMatrix(long m, long n, float **matrix, float **result) {
    for (long i = 0; i < m; i++) {
        for (long j = 0; j < n; j++) {
            assert(matrix[i] + j != NULL);
            assert(result[j] + i != NULL);
            result[j][i] = matrix[i][j];
        }
    }
}

void generate_input(run_config *s, float **A) {
    static int seeded = 0;
    if (!seeded) {
        srandom((unsigned int) time(NULL));
        seeded = 1;
    }

    for (int i = 0; i < s->m; i++) {
        for (int j = 0; j < s->n; j++) {
            A[i][j] = ((float) random() / RAND_MAX) * 10.0;
        }
    }
}

void print_usage(char *prog_name) {
    fprintf(stderr, "Usage: \"mpiexec -np <CORES> %s -o <output_file> -m <ROWS> -n <COLS> <input_file> \"\n", prog_name);
}

void error_exit(int rank, char *name, const char *msg, ...) {
    if (rank == 0) {

        char buf[16];
        time_t t = time(NULL);
        buf[strftime(buf, sizeof(buf), "%H:%M:%S", localtime(&t))] = '\0';
        fprintf(stderr, "[%s] %s %-5s : ", name, buf, "ERROR");

        va_list ap;
        va_start(ap, msg);
        vfprintf(stderr, msg, ap);
        va_end(ap);
        fprintf(stderr, "\n");
    }
    MPI_Finalize();
    exit(EXIT_FAILURE);
}