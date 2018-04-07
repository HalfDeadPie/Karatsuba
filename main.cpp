#include "static.h"

// read polynomial from input file
TYPE* readPolynomial(const std::string &file, int *size){
    int degree;
    TYPE *result = nullptr;

    std::ifstream inputFile(file);
    if (inputFile) {
        inputFile >> degree;
        *size = degree;
        result = new TYPE[degree];
        for(int i=0; i<degree; i++){
            inputFile>>result[i];
        }
    }
    return result;
}

// print polynomial
void print(TYPE *A, int size){
    for (int i = 0; i < size ; i++) {
        PRINT<<A[i]<<" ";
    }
    PRINT<<std::endl;
}

// Compare two numbers and return the smaller one
TYPE min(int a, int b) {
    return a > b ? b : a;
}

// Print the time resuls to standard output and times-file output
void printResult(double elapsed, int size){
    std::string variable_parameter;
    std::string name;

    #ifdef KARATSUBA
    name = "karatsuba_";
    #else
    name = "naive_";
    #endif

    #ifdef VAR_NAIVE_LIMIT
    variable_parameter ="naive_limit_";
    #endif

    #ifdef VAR_TASK_LIMIT
    variable_parameter = "task_limit_";
    #endif

    #ifdef VAR_FOR_THREADS
    variable_parameter =  "for_threads_";
    #endif

    #ifdef VAR_FOR_LIMIT
    variable_parameter = "for_limit_";
    #endif

    #ifdef VAR_NAIVE_THREADS
    variable_parameter =  "naive_threads";
    #endif

    #ifdef VAR_TILING_FACTOR
    variable_parameter =  "tiling_factor";
    #endif

    #ifdef XEON
    variable_parameter = "xeon_for_threads+naive_threads_";
    #endif

    std::ofstream outfile ("times/" + name + variable_parameter + std::to_string(size) + "deg_times.txt");

    PRINT << "Elapsed time: " << elapsed << std::endl;
    outfile<<variable_parameter<<" "<<elapsed<<::std::endl;
}

// naive algorithm
TYPE* naive(const TYPE *A, const TYPE *B, int size_A, int size_B){
    auto *result = new TYPE[size_A + size_B - 1];

    // naive multiplication
    #ifdef NAIVE_PARALELISM
            #pragma omp parallel for collapse(2) shared(result) num_threads(::NAIVE_THREADS)
            for (int i = 0; i < size_A; i += ::TILING_FACTOR) {
                for (int j = 0; j < size_B; j += ::TILING_FACTOR) {
		            int boundA = min(size_A, ::TILING_FACTOR + i);
                    int boundB = min(size_B, ::TILING_FACTOR + j);
                    for (int m = i; m < boundA; m++)
                        for (int n = j; n < boundB; n++)
                            result[m + n] = A[m] * B[n];
                }
            }
    #else
        for (int i = 0; i < size_A; i++) {
                // VECTORIZED
                for (int j = 0; j < size_B; j++) {
                    result[i + j] = A[i] * B[j];
                }
            }
    #endif

    return result;
}

// recursive Karatsuba algorithm with OpenMP directives
TYPE* karatsuba(const TYPE *A,const TYPE *B, int size, int depth) {
    TYPE *__restrict__ lowA, *__restrict__ highA, *__restrict__ lowB, *__restrict__ highB, *__restrict__ midA, *__restrict__ midB;

    // when the size of polynomial is bellow the limit, naive algorithm is called
    if (size <= ::NAIVE_LIMIT)
        return naive(A, B, size, size);

    // compute the half for splitting the polynomial
    // if size is odd number
    int half = size / 2;

    if(size % 2 == 1)
        half++;


    // prepare arrays for splitted parts
    lowA = new TYPE[half];  lowB = new TYPE[half];
    midA = new TYPE[half];  midB = new TYPE[half];
    highA = new TYPE[half]; highB = new TYPE[half];

    // init arrays with 0
    /*for(int i=0; i<half; i++)
        lowA[i] = lowB[i] = midA[i] = midB[i] = highA[i] = highB[i] = 0;*/

    // init low coefficients to new arrays  - VECTORIZED
    for(int i=0; i<half; i++){
        lowA[i] = A[i];
        lowB[i] = B[i];
    }

    // init high coefficients
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

    TYPE *__restrict__ z0, *__restrict__ z1, *__restrict__ z2;

    // compute the parts of result
    #ifdef TASK_PARALLELISM

        #pragma omp task shared(z0) if(half > ::TASK_LIMIT)
        z0 = karatsuba(lowA, lowB, half, depth + 1);

        #pragma omp task shared(z1) if(half > ::TASK_LIMIT)
        z1 = karatsuba(midA, midB, half, depth + 1);

        // job for calling thread
        z2 = karatsuba(highA, highB, half, depth + 1);

        #pragma omp taskwait

    #else

        z0 = karatsuba(lowA, lowB, half, depth+1);

        z1 = karatsuba(midA, midB, half, depth+1);

        z2 = karatsuba(highA, highB, half, depth+1);

    #endif

    // init the result array
    TYPE *__restrict__ result = new TYPE[2*size-1];

    for (int i = 0; i < 2*size-1; i++)
        result[i] = 0;

    // compute the result

    #ifdef FOR_PARALLELISM

        int i,j,k;
        #pragma omp parallel  default(shared) private(i, j, k) num_threads(::FOR_THREADS) if(depth <= ::PARALLEL_FOR_LIMIT)
        {

            //VECTORIZED
            #pragma omp for
            for (i = 0; i < 2 * half - 1; i++)
                result[i + 2 * half] += z2[i];

            // VECTORIZED
            #pragma omp for
            for (j = 0; j < 2 * half - 1; j++)
                result[j + half] += z1[j] - z2[j] - z0[j];

            // VECTORIZED
            #pragma omp for
            for (k = 0; k < 2 * half - 1; k++)
                result[k] += z0[k];
        }

    #else

        // VECTORIZED
        for (int i = 0; i < 2 * half - 1; i++)
            result[i + 2 * half] += z2[i];

        // VECTORIZED
        for (int i = 0; i < 2 * half - 1; i++)
            result[i + half] += z1[i] - z2[i] - z0[i];

        // VECTORIZED
        for (int i = 0; i < 2 * half - 1; i++)
            result[i] += z0[i];

    #endif

    return result;
}



int main(int argc, char* argv[]) {
    int size_A, size_B;
    TYPE *A = readPolynomial(argv[1], &size_A);
    TYPE *B = readPolynomial(argv[2], &size_B);
    int threads [8] = { 1, 2, 4, 6, 8, 12, 24, 32 };

    for( int i = START_LIMIT; i <= END_LIMIT; i += STEP ) {

        #ifdef VAR_NAIVE_LIMIT
            ::NAIVE_LIMIT = static_cast<int>(pow(2, i));
        #endif

        #ifdef VAR_TASK_LIMIT
            ::TASK_LIMIT = i;
        #endif

        #ifdef VAR_FOR_THREADS
            ::FOR_THREADS =  threads[i];
        #endif

	    #ifdef VAR_FOR_LIMIT
	        ::PARALLEL_FOR_LIMIT = i;
        #endif

        #ifdef VAR_NAIVE_THREADS
            ::NAIVE_THREADS = threads[i];
        #endif

        #ifdef VAR_TILING_FACTOR
            ::TILING_FACTOR = static_cast<int>(pow(2, i));
        #endif

	    #ifdef XEON
		    ::NAIVE_THREADS = static_cast<int>(pow(2, i));
		    ::FOR_THREADS = static_cast<int>(pow(2, i));
	    #endif

        // Record start time
        auto start = omp_get_wtime();

        // multiply the polynomials
        TYPE *result;
        #ifdef KARATSUBA
            #ifdef TASK_PARALLELISM
                #pragma omp parallel
                #pragma omp single
            #endif
            result = karatsuba(A, B, size_A, 0);
        #else
            result = naive(A, B, size_A, size_B);
        #endif

        // stop timer
        auto finish = omp_get_wtime();

        // Record end time
        double elapsed = finish - start;

        // print the result if DEBUG is defined
        #ifdef DEBUG
            print(result, size_A + size_B - 1);
        #endif

        printResult(elapsed, size_A);
    }

    return 0;
}
