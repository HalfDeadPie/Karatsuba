#include <iostream>
#include <sstream>
#include <fstream>

#define PRINT std::cout
#define ENDLINE std::endl

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

// naive multiplication
long* naive(const long *A, const long *B, int size_A, int size_B){
    long *result = new long[size_A + size_B - 1];

    // initialization
    for(int iterator=0; iterator<size_A + size_B - 1; iterator++)
        result[iterator] = 0;

    // naive multiplication
    for (int i = 0; i < size_A; i++) {
        for (int j = 0; j < size_B; ++j) {
            result[i+j] += A[i] * B[j];
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

// recursive karatsuba algorithm
long* karatsuba(const long *A,const long *B, int size) {
    long *lowA, *highA, *lowB, *highB, *midA, *midB;

    // find the limit
    if (size == 1)
        return naive(A, B, size, size);

    // compute the half for splitting the polynomial
    int half = size / 2;

    // if size is odd number
    if(size % 2 == 1)
        half++;

    // prepare arrays for splitted parts
    lowA = new long[half];
    lowB = new long[half];
    midA = new long[half];
    midB = new long[half];
    highA = new long[half];
    highB = new long[half];

    // init arrays
    for(int i=0; i<half; i++) lowA[i] =  lowB[i] = midA[i] = midB[i] = highA[i] = highB[i] = 0;

    // init low coefficients
    for(int i=0; i<half; i++){
        lowA[i] = A[i];
        lowB[i] = B[i];
    }

    // init high coefficients
    for(int i=half; i<size; i++){
        highA[i-half] = A[i];
        highB[i-half] = B[i];
    }

    // init mid coefficients
    for(int i=0; i<half; i++){
        midA[i] = lowA[i] + highA[i];
        midB[i] = lowB[i] + highB[i];
    }

    long *z0 = karatsuba(lowA, lowB, half);
    long *z1 = karatsuba(midA, midB, half);
    long *z2 = karatsuba(highA, highB, half);

    // compute the result
    auto *result = new long[2*size-1];
    for (int i = 0; i < 2*size-1; i++) result[i] = 0;

    for(int i=0; i<2*half-1; i++) {
        result[i + 2 * half] += z2[i];
        result[i + half] += z1[i] - z2[i] - z0[i];
        result[i] += z0[i];
    }

    return result;
}

int main(int argc, char* argv[]) {
    int size_A, size_B;
    long* A = readPolynomial(argv[1], &size_A);
    long* B = readPolynomial(argv[2], &size_B);

    //multiply the polynomials
    long* result = karatsuba(A, B, size_A);

    print(result, size_A + size_B - 1);
    return 0;
}