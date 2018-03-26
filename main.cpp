#include <iostream>
#include <sstream>
#include <fstream>
#include <chrono>
#include <omp.h>

// Macro used for printing to standard output
#define PRINT std::cout

/*Limit value that is compared to actual size of multiplicated polynomials. When their size is equal or bellow
 * this value, the polynomials are multiplicated using the naive method*/
int naive_limit;

/*If defined, Karatsuba algorithm is called. Otherwise, the naive algorithm computes all multiplications*/
#define KARATSUBA

#define START_LIMIT 1
#define END_LIMIT 1
#define STEP 1

/*Debug the result of multiplication - if defined, the result is printed to standard output*/
//#define DEBUG

// read polynomial from input file
long* readPolynomial(const std::string &file, int *size){
    int degree;
    long *result = nullptr;

    std::ifstream inputFile(file);
    if (inputFile) {
        inputFile >> degree;
        *size = degree;
        result = new long[degree];
        for(int i=0; i<degree; i++){
            inputFile>>result[i];
        }
    }
    return result;
}

// print polynomial
void print(long *A, int size){
    for (int i = 0; i < size ; i++) {
        PRINT<<A[i]<<" ";
    }
    PRINT<<std::endl;
}

// naive multiplication
long* naive(long *A, long *B, int size_A, int size_B){
    auto *result = new long[size_A + size_B - 1];

    // initialization
    for(int iterator=0; iterator<size_A + size_B - 1; iterator++)
        result[iterator] = 0;

    // naive multiplication
    for (int i = 0; i < size_A; i++) {
        // VECTORIZED
        for (int j = 0; j < size_B; ++j) {
            result[i+j] += A[i] * B[j];
        }
    }

    return result;
}

// recursive Karatsuba algorithm
long* karatsuba(long *A, long *B, int size) {
    long *lowA, *highA, *lowB, *highB, *midA, *midB;

    // when the size of polynomial is bellow the limit, naive algorithm is called
    if (size <= ::naive_limit)
        return naive(A, B, size, size);

    // compute the half for splitting the polynomial
    int half = size / 2;

    // if size is odd number
    if(size % 2 == 1)
        half++;

    // prepare arrays for splitted parts
    lowA = new long[half];  lowB = new long[half];
    midA = new long[half];  midB = new long[half];
    highA = new long[half]; highB = new long[half];

    // init arrays with 0
    for(int i=0; i<half; i++)
        lowA[i] = lowB[i] = midA[i] = midB[i] = highA[i] = highB[i] = 0;

    // init low coefficients to new arrays  - VECTORIZED
    for(int i=0; i<half; i++){
        lowA[i] = A[i];
        lowB[i] = B[i];
    }

    // init high coefficients
    // VECTORIZED
    for(int i=half; i<size; i++){
        highA[i - half] = A[i];
        highB[i - half] = B[i];
    }

    // init mid coefficients
    // VECTORIZED
    for(int i=0; i<half; i++){
        midA[i] = lowA[i] + highA[i];
        midB[i] = lowB[i] + highB[i];
    }

    // compute the parts of result
    long *z0 = karatsuba(lowA, lowB, half);
    long *z1 = karatsuba(midA, midB, half);
    long *z2 = karatsuba(highA, highB, half);

    // init the result array
    auto *result = new long[2*size-1];
    for (int i = 0; i < 2*size-1; i++)
        result[i] = 0;

    // compute the result
    // VECTORIZED
    for(int i=0; i<2*half-1; i++)
        result[i + 2 * half] += z2[i];

    // VECTORIZED
    for(int i=0; i<2*half-1; i++)
        result[i + half] += z1[i] - z2[i] - z0[i];

    // VECTORIZED
    for(int i=0; i<2*half-1; i++)
        result[i] += z0[i];

    return result;
}

int main(int argc, char* argv[]) {
    int size_A, size_B;
    long *A = readPolynomial(argv[1], &size_A);
    long *B = readPolynomial(argv[2], &size_B);

    #ifdef KARATSUBA
        std::ofstream outfile ("times/karatsuba_" + std::to_string(size_A) + "deg_times.txt");
    #else
        std::ofstream outfile ("times/naive_" + std::to_string(size_A) + "deg_times.txt");
    #endif

    for( int i = START_LIMIT; i <= END_LIMIT; i += STEP ) {

        ::naive_limit = i;

        // Record start time
        auto start = omp_get_wtime();

        // multiply the polynomials
        #ifdef KARATSUBA
            long *result = karatsuba(A, B, size_A);
        #else
            long *result = naive(A, B, size_A, size_B);
        #endif

        // stop timer
        auto finish = omp_get_wtime();

        // Record end time
        double elapsed = finish - start;

        // print the result if DEBUG is defined
        #ifdef DEBUG
                print(result, size_A + size_B - 1);
        #endif

        // print the elapsed time
        PRINT << "Elapsed time"<< ::naive_limit << ": " << elapsed << std::endl;
        outfile << elapsed << ::std::endl;
    }

    return 0;
}
