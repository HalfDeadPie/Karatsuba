#include <iostream>
#include <sstream>
#include <utility>
#include <vector>
#include <iterator>
#include <typeinfo>
#include <fstream>
#include "Polynom.h"
#include <chrono>

#define DEBUG
//#define DEBUG_DEEP
#define LIMIT 3

// read polynomial from input file
std::vector<long> readPolynomial(const std::string &file){
    std::vector<long> result;
    std::ifstream inputFile(file);
    if (inputFile) {
        long value;

        // read the elements in the file into a vector
        while ( inputFile >> value ) {
            result.push_back(value);
        }
    }
    return result;
}

// iterative naive algorithm
Polynom multiply_naive(Polynom a, Polynom b){
    // get boundaries
    unsigned long m = a.getSize();
    unsigned long n = b.getSize();

    // create vector of double with proper size for new polynom
    std::vector<long> coef (m+n-1, 0);

    // all terms of first polynom
    for( unsigned long i = 0; i < m; i++ ){
        //all terms of second polynom
        for( unsigned long j = 0; j < n; j++ ){
            coef[i+j] += a.getAt(i) * b.getAt(j);
        }
    }

    // return the result polynom
    return Polynom (coef, coef.size() - 1);
}

// fill the vector Di
std::vector<long> fill_di(Polynom a, Polynom b){
    std::vector<long> result (a.getSize(), 0);
    for( unsigned long i = 0; i < a.getSize(); i++ ){
        result[i] = a.getAt(i) * b.getAt(i);
    }
    return result;
}

// start calculation
unsigned long calcStart(unsigned long position, unsigned long size){
    // if position is bigger than size of single polynom
    if(position>=size){
        // return proper position
        return position + 1 - size;
    }else{
        // return the start of polynom
        return 0;
    }
}

// end calculation
unsigned long calcEnd(unsigned long position){
    return (position+1)/2;
}

// non-recursive karatsuba algorithm
// https://eprint.iacr.org/2006/224.pdf
Polynom karatsuba(Polynom a, Polynom b){
    // get size of both polynoms (their number of coefficients is equal)
    unsigned long size = a.getSize();

    // create empty coefficient vector with proper size and fill it with 0
    std::vector<long> result (2 * size - 1, 0);

    // fill Di vector with Ai * Bi
    std::vector<long> D = fill_di(a, b);

    // set the first coefficient
    result[0] = D[0];

    //set the last coefficient
    result[2 * (size - 1)] = D[size - 1];

    // for all coefficients of result vector
    for (unsigned long position=1; position < 2*(size-1); position++){
        // for even coefficient add Di/2
        if ( position % 2 == 0)
            result[position] += D[position>>1];

        // calculate start position in polynom
        unsigned long start = calcStart(position, size);

        // calculate end position in polynom
        unsigned long end = calcEnd(position);

        // inner loop: sum (Dst) - sum (Ds + Dt) where s+t=i
        for(unsigned long inner = start; inner<end; inner++){
            result[position] += ( a.getAt(inner)+a.getAt(position-inner) ) * ( b.getAt(inner)+b.getAt(position-inner) );
            result[position] -= (D[inner] + D[position-inner]);
        }
    }

    // create polynom with computed coefficients
    return Polynom (result, result.size() - 1);
}

// decides which multiplication will be used
Polynom multiply(Polynom a, Polynom b){
    if( (a.getDegree() != b.getDegree()) || a.getDegree() < LIMIT || b.getDegree() < LIMIT){
#ifdef DEBUG
        std::cout << "Naive multiplication" << std::endl;
#endif
        return multiply_naive(a, b);
    }else{
#ifdef DEBUG
        std::cout << "Karatsuba multiplication" << std::endl;
#endif
        return karatsuba(a, b);
    }
}

int main(int argc, char* argv[]) {
    std::vector<long> poly_vec1 = readPolynomial(argv[1]);
    std::vector<long> poly_vec2 = readPolynomial(argv[2]);
    Polynom a (poly_vec1, poly_vec1.size() - 1), b (poly_vec2,poly_vec2.size() - 1);

    auto start = std::chrono::high_resolution_clock::now();
    Polynom result = multiply(a, b);
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";

#ifdef DEBUG_DEEP
    // print the polynom
    for(long actual : result.getCoefficients()){
        std::cout << actual << " ";
    }
    std::cout << std::endl;
#endif

    return 0;
}