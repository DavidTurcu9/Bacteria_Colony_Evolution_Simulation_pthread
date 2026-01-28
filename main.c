#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>

#define FILENAME_SIZE 1024

typedef unsigned char cell_t; // 0 if cell is empty

typedef struct {
    int rows;
    int cols;
    cell_t *buf; // 1D array representing 2D array of bacteria, uses dynamic allocation
} grid_t;

grid_t* grid_create(int rows, int cols) {
    grid_t *grid = NULL;
    grid = (grid_t *) malloc(sizeof(grid_t));
    if (grid == NULL) {
        perror("malloc");
        exit(EXIT_FAILURE);
    }
    grid->rows = rows;
    grid->cols = cols;
    // initialize bacteria grid with all cells as 0
    grid->buf = (cell_t *) calloc(rows * cols, sizeof(cell_t));
    if (grid->buf == NULL) {
        perror("calloc");
        exit(EXIT_FAILURE);
    }
    return grid;
}

void grid_free(grid_t *g) {
    if (!g) {
        return;
    }
    free(g->buf);
    free(g);
}

void print_grid10(grid_t* grid) {
    for (int row = 0; row < grid->rows; row++) {
        for (int col = 0; col < grid->cols; col++) {
            printf("%d ", grid->buf[row * grid->cols + col]);
        }
        printf("\n");
    }
}

void print_gridX(grid_t* grid) {
    for (int row = 0; row < grid->rows; row++) {
        for (int col = 0; col < grid->cols; col++) {
            if (grid->buf[row * grid->cols + col] != 0) {
                printf("X ");
            } else {
                printf(". ");
            }

        }
        printf("\n");
    }
}

void read_file_and_init_grid(char* file_name, grid_t** grid) {
    FILE* fin = fopen(file_name, "r");
    if (fin == NULL) {
        perror("fopen");
        exit(EXIT_FAILURE);
    }

    // read rows and cols from top of file
    int rows, cols;
    if (fscanf(fin, "%d %d\n", &rows, &cols) != 2) {
        perror("fscanf");
        fclose(fin);
        exit(EXIT_FAILURE);
    }

    grid_t* grid_temp = grid_create(rows, cols); // creates grid with rows rows, cols columns and array of rows * columns cells initialized with 0

    // read grid from file with . and X format
    char* buffer = (char *) malloc(cols + 2);
    if (buffer == NULL) {
        perror("malloc");
        exit(EXIT_FAILURE);
    }
    for (int row = 0; row < rows; row++) {
        if (!fgets(buffer, cols + 2, fin)) {
            perror("fgets");
            free(buffer);
            exit(EXIT_FAILURE);
        }
        buffer[strcspn(buffer, "\r\n")] = 0; // remove newline from buffer
        for (int col = 0; col < cols; col++) {
            if (buffer[col] == 'X') {
                grid_temp->buf[row * cols + col] = 1;
            } else {
                grid_temp->buf[row * cols + col] = 0;
            }
        }
    }
    free(buffer);
    fclose(fin);
    *grid = grid_temp;
}

// counts how many bacteria are in the 8 surrounding cells
int count_neighbors(grid_t* grid, int row, int col) {
    int count = 0;
    // checks the surrounding 8 cells
    for (int dr = -1; dr <= 1; dr++) {
        for (int dc = -1; dc <= 1; dc++) {
            // skip current cell
            if (dr == 0 && dc == 0) {
                continue;
            }
            // get indices of the cell it's checking
            int nr = row + dr;
            int nc = col + dc;
            // if grid->buf[nr][nc] is not outside the grid and a bacteria, count++
            if (nr >= 0 && nr <grid->rows && nc >= 0 && nc < grid->cols) {
                count += grid->buf[nr * grid->cols + nc];
            }
        }
    }
    return count;
}

// for testing count_neighbors function
void print_grid_nneighbors(grid_t* grid) {
    for (int row = 0; row < grid->rows; row++) {
        for (int col = 0; col < grid->cols; col++) {
            printf("%d ", count_neighbors(grid, row, col));
        }
        printf("\n");
    }
}

void step_serial(grid_t* currentgen, grid_t* nextgen) {
    for (int r = 0; r < currentgen->rows; r++) {
        for (int c = 0; c < currentgen->cols; c++) {
            int nr_neighbors = count_neighbors(currentgen, r, c);
            cell_t new_value = 0;
            if (currentgen->buf[r * currentgen->cols + c] == 1) { // if current cell is alive
                if (nr_neighbors < 2 || nr_neighbors > 3) {
                    new_value = 0; // dies if isolated or suffocated
                } else {
                    new_value = 1;
                }
            } else { // if current cell is dead
                if (nr_neighbors == 3) {
                    new_value = 1; // multiplies if there are exactly 3 neighbors
                } else {
                    new_value = 0;
                }
            }
            nextgen->buf[r * currentgen->cols + c] = new_value;
        }
    }
}

double run_serial(grid_t* currentgen, grid_t* nextgen, int generations, bool debug) {
    struct timespec start_time, end_time;
    clock_gettime(CLOCK_MONOTONIC, &start_time);

    // run serial algorithm for generations generations
    for (int g = 0; g < generations; g++) {
        step_serial(currentgen, nextgen);

        // swap pointers(currentgen with nextgen)
        cell_t* tmp = currentgen->buf;
        currentgen->buf = nextgen->buf;
        nextgen->buf = tmp;

        if (debug && g != generations - 1) {
            printf("Serial grid after generation %d:\n", g + 1);
            print_gridX(currentgen);
            printf("\n");
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &end_time);
    double elapsed_time = (end_time.tv_sec - start_time.tv_sec) + (end_time.tv_nsec - start_time.tv_nsec) / 1e9; // 1e9 is 10^9 (1000000000.0)
    return elapsed_time;
}

typedef struct {
    grid_t *cur;
    grid_t *next;
    int row_start;
    int row_end;
    int gens;
    pthread_barrier_t* barrier;
} pthread_args_t;

void* thread_func(void* varg) {
    pthread_args_t* args = (pthread_args_t *) varg;
    int gens = args->gens;
    int row_start = args->row_start;
    int row_end = args->row_end;
    for (int g = 0; g < gens; g++) {
        grid_t* cur = args->cur;
        grid_t* next = args->next;
        // this for would be step_parallel
        for (int row = row_start; row < row_end; row++) {
            for (int col = 0; col < cur->cols; col++) {
                int nr_neighbors = count_neighbors(cur, row, col);
                cell_t new_value = 0;
                if (cur->buf[row * cur->cols + col] == 1) { // if current cell is alive
                    if (nr_neighbors < 2 || nr_neighbors > 3) {
                        new_value = 0; // dies if isolated or suffocated
                    } else {
                        new_value = 1;
                    }
                } else { // if current cell is dead
                    if (nr_neighbors == 3) {
                        new_value = 1; // multiplies if there are exactly 3 neighbors
                    } else {
                        new_value = 0;
                    }
                }
                next->buf[row * cur->cols + col] = new_value;
            }
        }
        // wait for others to finish computation
        pthread_barrier_wait(args->barrier);
        // wait for master to swap cur/next pointers in main thread
        pthread_barrier_wait(args->barrier);
    }
    return NULL;
}

double run_parallel(grid_t* currentgen, grid_t* nextgen, int ngenerations, int nthreads, bool debug) {
    if (nthreads < 1) {
        nthreads = 1;
    }
    pthread_t* threads = (pthread_t *) malloc(sizeof(pthread_t) * nthreads);
    pthread_args_t* pthread_args = (pthread_args_t *) malloc(sizeof(pthread_args_t) * nthreads);
    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier, NULL, nthreads + 1); // +1 for master thread

    // compute workload for each thread and arguments for pthread_create
    grid_t *cur = currentgen;
    grid_t *next = nextgen;
    int R = currentgen->rows;
    int base = R / nthreads;
    int rem = R % nthreads;
    int r0 = 0;
    for (int thr = 0; thr < nthreads; thr++) {
        int rcount = 0;
        if (thr < rem) {
            rcount = base + 1;
        } else {
            rcount = base;
        }
        pthread_args[thr].row_start = r0;
        pthread_args[thr].row_end = r0 + rcount;
        pthread_args[thr].barrier = &barrier;
        pthread_args[thr].gens = ngenerations;
        pthread_args[thr].cur = cur;
        pthread_args[thr].next = next;
        r0 += rcount;
    }

    // create threads
    for (int t = 0; t < nthreads; t++) {
        if (pthread_create(&threads[t], NULL, thread_func, &pthread_args[t]) != 0) {
            perror("pthread_create");
            exit(EXIT_FAILURE);
        }
    }

    // start timing
    struct timespec start_time, end_time;
    clock_gettime(CLOCK_MONOTONIC, &start_time);

    for (int g = 0; g < ngenerations; g++) {
        pthread_barrier_wait(&barrier);
        grid_t *tmp = cur;
        cur = next;
        next = tmp;
        for (int t = 0; t < nthreads; t++) {
            pthread_args[t].cur = cur;
            pthread_args[t].next = next;
        }
        if (debug && g != ngenerations - 1) {
            printf("Parallel grid after generation %d:\n", g + 1);
            print_gridX(cur);
            printf("\n");
        }
        pthread_barrier_wait(&barrier);
    }

    clock_gettime(CLOCK_MONOTONIC, &end_time);
    double elapsed = (end_time.tv_sec - start_time.tv_sec) + (end_time.tv_nsec - start_time.tv_nsec) / 1e9;

    // join threads
    for (int t = 0; t < nthreads; t++) {
        pthread_join(threads[t], NULL);
    }

    if (cur != currentgen) {
        memcpy(currentgen->buf, cur->buf, (size_t)currentgen->rows * currentgen->cols * sizeof(cell_t));
    }

    pthread_barrier_destroy(&barrier);
    free(threads);
    free(pthread_args);
    return elapsed;
}

bool grids_equal(const grid_t *a, const grid_t *b) {
    if (a->rows != b->rows || a->cols != b->cols) {
        return false;
    }
    return memcmp(a->buf, b->buf, (size_t)a->rows * a->cols * sizeof(cell_t)) == 0;
}

void write_grid_to_file(char filename[], grid_t* grid) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("fopen");
        exit(EXIT_FAILURE);
    }
    fprintf(file, "%d %d\n", grid->rows, grid->cols);
    for (int i = 0; i < grid->rows; i++) {
        for (int j = 0; j < grid->cols; j++) {
            if (grid->buf[i * grid->cols + j] == 1) {
                fprintf(file, "X ");
            } else {
                fprintf(file, ". ");
            }
        }
        fprintf(file, "\n");
    }
    fclose(file);
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <input-file> <generations> <num-threads> [DEBUG]\n", argv[0]);
        return 1;
    }

    char* input_file = argv[1];
    int ngenerations = atoi(argv[2]);
    int nthreads = atoi(argv[3]);
    bool debug;
    if (argc >= 5 && strcmp(argv[4], "DEBUG") == 0) {
        debug = true;
    } else {
        debug = false;
    }

    // give initial configuration of bacteria (loaded from file, also nr of rows and columns loaded from file)
    grid_t* initial_grid = NULL;
    read_file_and_init_grid(input_file, &initial_grid);

    // here the code makes different grids for the current and next generations of bacteria
    // 2 for the serial version and 2 for the parallel version of the algorithm.
    // 4 grids are needed for the automatic verification method to compare that the serial and parallel versions produce the same result.
    grid_t* grid_serial_currentgen = grid_create(initial_grid->rows, initial_grid->cols);
    grid_t* grid_serial_nextgen = grid_create(initial_grid->rows, initial_grid->cols);
    grid_t* grid_parallel_currentgen = grid_create(initial_grid->rows, initial_grid->cols);
    grid_t* grid_parallel_nextgen = grid_create(initial_grid->rows, initial_grid->cols);
    // error handle malloc
    if (!grid_serial_currentgen || !grid_serial_nextgen || !grid_parallel_currentgen || !grid_parallel_nextgen) {
        perror("malloc");
        exit(EXIT_FAILURE);
    }

    // copy buf from initial_grid to the currentgen grids
    memcpy(grid_serial_currentgen->buf, initial_grid->buf, (size_t)initial_grid->cols * (size_t)initial_grid->rows * sizeof(cell_t));
    memcpy(grid_parallel_currentgen->buf, initial_grid->buf, (size_t)initial_grid->cols * (size_t)initial_grid->rows * sizeof(cell_t));

    if (debug) {
        printf("Initial grid:\n");
        print_gridX(initial_grid);
    }

    // start algorithm
    double serial_time = run_serial(grid_serial_currentgen, grid_serial_nextgen, ngenerations, debug);
    if (debug) {
        printf("Serial final after %d generations:\n", ngenerations);
        print_gridX(grid_serial_currentgen);
    }

    double parallel_time = run_parallel(grid_parallel_currentgen, grid_parallel_nextgen, ngenerations, nthreads, debug);
    if (debug) {
        printf("Parallel final after %d generations:\n", ngenerations);
        print_gridX(grid_parallel_currentgen);
    }

    if (grids_equal(grid_serial_currentgen, grid_parallel_currentgen)) {
        printf("Serial and parallel results match!\n");
    } else {
        fprintf(stderr, "Mismatch between serial and parallel results!\n");
        exit(EXIT_FAILURE);
    }

    // Output final results to files
    char fout_serial[FILENAME_SIZE], fout_parallel[FILENAME_SIZE];
    snprintf(fout_serial, FILENAME_SIZE, "%s_serial_out.txt", input_file);
    snprintf(fout_parallel, FILENAME_SIZE, "%s_parallel_out.txt", input_file);
    write_grid_to_file(fout_serial, grid_serial_currentgen);
    write_grid_to_file(fout_parallel, grid_parallel_currentgen);

    printf("Serial time:   %.6f s\n", serial_time);
    printf("Parallel time: %.6f s (threads = %d)\n", parallel_time, nthreads);
    if (parallel_time > 0.0) {
        printf("Speedup (serial/parallel): %.4f\n", serial_time / parallel_time);
    }

    grid_free(initial_grid);
    grid_free(grid_serial_currentgen);
    grid_free(grid_serial_nextgen);
    grid_free(grid_parallel_currentgen);
    grid_free(grid_parallel_nextgen);

    return 0;
}