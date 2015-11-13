#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>

double*** threadPartialBofZ;
double** preX;
double*** threadPartialOutMM;
int targetDimension = 3;
int blockSize = 64;

void matrixMultiply(double** A, double** B, int aHeight, int bWidth, int comm, int bz, double** C) {

    int aHeightBlocks = aHeight / bz; // size = Height of A
    int aLastBlockHeight = aHeight - (aHeightBlocks * bz);
    if (aLastBlockHeight > 0) {
        aHeightBlocks++;
    }

    int bWidthBlocks = bWidth / bz; // size = Width of B
    int bLastBlockWidth = bWidth - (bWidthBlocks * bz);
    if (bLastBlockWidth > 0) {
        bWidthBlocks++;
    }

    int commnBlocks = comm / bz; // size = Width of A or Height of B
    int commLastBlockWidth = comm - (commnBlocks * bz);
    if (commLastBlockWidth > 0) {
        commnBlocks++;
    }

    int aBlockHeight = bz;
    int bBlockWidth = bz;
    int commBlockWidth = bz;

    for (int ib = 0; ib < aHeightBlocks; ib++) {
        if (aLastBlockHeight > 0 && ib == (aHeightBlocks - 1)) {
            aBlockHeight = aLastBlockHeight;
        }
        bBlockWidth = bz;
        commBlockWidth = bz;
        for (int jb = 0; jb < bWidthBlocks; jb++) {
            if (bLastBlockWidth > 0 && jb == (bWidthBlocks - 1)) {
                bBlockWidth = bLastBlockWidth;
            }
            commBlockWidth = bz;
            for (int kb = 0; kb < commnBlocks; kb++) {
                if (commLastBlockWidth > 0 && kb == (commnBlocks - 1)) {
                    commBlockWidth = commLastBlockWidth;
                }

                for (int i = ib * bz; i < (ib * bz) + aBlockHeight; i++) {
                    for (int j = jb * bz; j < (jb * bz) + bBlockWidth;
                         j++) {
                        for (int k = kb * bz;
                             k < (kb * bz) + commBlockWidth; k++) {
                            if (A[i][k] != 0 && B[k][j] != 0) {
                                C[i][j] += A[i][k] * B[k][j];
                            }
                        }
                    }
                }
            }
        }
    }
}

void bcReplica(int threadCount, int iterations, int globalColCount, int rowCountPerUnit) {
    preX = (double**) malloc(sizeof(double*)*globalColCount);
    int i;
    for (i = 0; i < globalColCount; ++i ){
        preX[i] = (double*) malloc(sizeof(double)*targetDimension);
    }

    threadPartialBofZ = (double***) malloc(sizeof(double**)*threadCount);
    threadPartialOutMM = (double***) malloc(sizeof(double**)*threadCount);
    int j;
    for (i = 0; i < threadCount; ++i){
        threadPartialBofZ[i] = (double**) malloc(sizeof(double*)*rowCountPerUnit);
        threadPartialOutMM[i] = (double**) malloc(sizeof(double*)*rowCountPerUnit);
        for (j = 0; j < rowCountPerUnit; ++j){
            threadPartialBofZ[i][j] = (double*) malloc(sizeof(double)*globalColCount);
            threadPartialOutMM[i][j] = (double*) malloc(sizeof(double)*targetDimension);
        }
    }

    int itr;
    for (itr = 0; itr < iterations; ++itr){

    }
}


int main() {
    printf("Hello World");
    const int total_threads = omp_get_max_threads();
    printf("There are %d available threads.\n", total_threads); fflush(stdout);

    /*int num_threads = omp_get_max_threads();*/
    int num_points = 1000;
    long num_itr = 1000000000;

    /*int array[100];*/

    time_t t;
    srand((unsigned)time(&t));
    //parallelize this part
    #pragma omp parallel
    {
        double result =0.0;
        const int thread_id = omp_get_thread_num();

        long i;
        for (i = 0; i < num_itr; ++i) {
            result += sqrt(rand() % num_points);
        }
        printf("Hello world from thread %d result %lf iterations %ld\n", thread_id, result, num_itr);
    }

    return 0;
}

