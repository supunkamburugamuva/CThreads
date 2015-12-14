#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

double*** threadPartialBofZ;
double** preX;
double*** threadPartialOutMM;
int targetDimension = 3;
int blockSize = 64;

double currentTimeInSeconds(void)
{

    int flag;
    clockid_t cid = CLOCK_REALTIME; // CLOCK_MONOTONE might be better
    struct timespec tp;
    double timing;

    flag = clock_gettime(cid, &tp);
    if (flag == 0) timing = tp.tv_sec + 1.0e-9*tp.tv_nsec;
    else           timing = -17.0;         // If timer failed, return non-valid time

    return(timing);
};


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

    int ib, jb, kb, i, j, k;
    for (ib = 0; ib < aHeightBlocks; ib++) {
        if (aLastBlockHeight > 0 && ib == (aHeightBlocks - 1)) {
            aBlockHeight = aLastBlockHeight;
        }
        bBlockWidth = bz;
        commBlockWidth = bz;
        for (jb = 0; jb < bWidthBlocks; jb++) {
            if (bLastBlockWidth > 0 && jb == (bWidthBlocks - 1)) {
                bBlockWidth = bLastBlockWidth;
            }
            commBlockWidth = bz;
            for (kb = 0; kb < commnBlocks; kb++) {
                if (commLastBlockWidth > 0 && kb == (commnBlocks - 1)) {
                    commBlockWidth = commLastBlockWidth;
                }

                for (i = ib * bz; i < (ib * bz) + aBlockHeight; i++) {
                    for (j = jb * bz; j < (jb * bz) + bBlockWidth;
                         j++) {
                        for (k = kb * bz;
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

    double totalTime = 0.0;
    double times[threadCount];

    int itr;
    int k;
    for (itr = 0; itr < iterations; ++itr){
        for (i = 0; i < globalColCount; ++i){
            for (j = 0; j < targetDimension; ++j){
                preX[i][j] = (double)rand() / (double)RAND_MAX;
            }
        }

        for (k = 0; k < threadCount; ++k){
            for (i = 0; i < rowCountPerUnit; ++i){
                for (j = 0; j < globalColCount; ++j){
                    threadPartialBofZ[k][i][j] = (double)rand() / (double)RAND_MAX;
                }
            }
        }

        for (k = 0; k < threadCount; ++k){
            for (i = 0; i < rowCountPerUnit; ++i){
                for (j = 0; j < targetDimension; ++j){
                    threadPartialOutMM[k][i][j] = (double)0.0;
                }
            }
        }


        {
 //           int num_t = omp_get_num_threads();
 //           if (num_t != threadCount){
 //               printf("Error num_t %d expected %d", num_t, threadCount);
 //           }

            double t1, t2;
            //const int threadIdx = omp_get_thread_num();
            const int threadIdx = 0;
            t1 = currentTimeInSeconds();
            //t1 = omp_get_wtime();
            matrixMultiply(threadPartialBofZ[threadIdx], preX, rowCountPerUnit, targetDimension, globalColCount, blockSize, threadPartialOutMM[threadIdx]);
            t2 = currentTimeInSeconds() - t1;
            //times[threadIdx] += t2;
	    totalTime += t2;
        }
/*
        double max = -1.0;
        for (i = 0; i < threadCount; ++i){
            if (times[i] > max){
                max = times[i];
            }
        }
        totalTime += max;*/
    }

    printf("%d,%d,%d,%d,%lf s\n", rowCountPerUnit, globalColCount, iterations, threadCount, totalTime);
}


int main(int argc, char **args) {
    if (argc < 5) {
        printf("We need 4 arguments");
        exit(1);
    }

    //double t1 = currentTimeInSeconds();
    //sleep(4);
    //double t2 = currentTimeInSeconds();
    //printf("%lf s\n",(t2-t1));

    //const int total_threads = omp_get_max_threads();
    //printf("There are %d total available threads.\n", total_threads); fflush(stdout);

    /* Take these as command line args
     * 1. num threads -- t
     * 2. iterations -- i
     * 3. rowcount -- r
     * 4. colcount  -- c*/
    int t = 0;
    int i = 0;
    int c = 0;
    int r = 0;

    t = atoi(args[1]);
    i = atoi(args[2]);
    r = atoi(args[3]);
    c = atoi(args[4]);

//    omp_set_num_threads(t);


    bcReplica(t, i, c, r);

    return 0;
}

