#include <omp.h>
#include <stdio.h>

void main() {
    int num_threads = omp_get_max_threads();
    printf("Number of threads: %d\n", num_threads);

    #pragma omp parallel
    {
        #pragma omp single
        {
            num_threads = omp_get_num_threads();
            printf("Number of threads in parallel region: %d\n", num_threads);
        }

        // Parallelized code here
    }
}
