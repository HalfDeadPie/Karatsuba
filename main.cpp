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

// save the result to output file
void saveResult(std::string &file, const TYPE *polynomial, TYPE size){
    std::ofstream myfile;
    myfile.open (file);
    myfile << 2 * size - 1 << "\n";
    for(int i=0; i<2 * size - 1; i++){
        myfile << polynomial[i] << " ";
    }
    myfile.close();
}

// compare two numbers and return the smaller one
TYPE min(int a, int b) {
    return a > b ? b : a;
}

/* computes the checksum based on equality of coefficients on random positions
 * but the result must be saved first. The best way how to save result is to
 * define SAVE in static.h and run the naive sequent algorithm*/
void checksum(std::string &file, const TYPE *polynomial, TYPE size){
    int realSize;
    int mySize = 2 * size - 1;
    TYPE *realResult = readPolynomial(file, &realSize);

    if(realSize != mySize)
        PRINT<<"\n\n!!!WRONG CHECKSUM (SIZE)!!!\n\n";

    int position = 0;
    for(int i = 0; i < NUM_OF_CHECKS; i++){
        position = rand() % mySize;
        if(polynomial[position] != realResult[position])
            PRINT<<"\n\n!!!WRONG CHECKSUM (RESULT)!!!\n\n";
    }

}

// print the time resuls to standard output and times-file output
void printStats(double elapsed, const TYPE *polynomial, int size){
    std::string variable_parameter, name;
    int column=0;

    if  (::VERSION == 1)
        variable_parameter += "seq_";
    else if (::VERSION == 2)
        variable_parameter += "omp_";
    else if (::VERSION == 3)
        variable_parameter += "acc_";
    else if (::VERSION == 4)
        variable_parameter += "seq-iterative_";
    else if (::VERSION == 5)
        variable_parameter += "acc-iterative_";


    #ifdef KARATSUBA
    name = "karatsuba_";
    #else
    name = "naive_";
    #endif

    #ifdef VAR_NAIVE_LIMIT
    variable_parameter ="naive_limit_";
    column = ::NAIVE_LIMIT;
    #endif

    #ifdef VAR_TASK_LIMIT
    variable_parameter = "task_limit_";
    column = ::TASK_LIMIT;
    #endif

    #ifdef VAR_FOR_THREADS
    variable_parameter =  "for_threads_";
    column = ::FOR_THREADS;
    #endif

    #ifdef VAR_FOR_LIMIT
    variable_parameter = "for_limit_";
    column = ::PARALLEL_FOR_LIMIT;
    #endif

    #ifdef VAR_NAIVE_THREADS
    variable_parameter =  "naive_threads";
    column = ::NAIVE_THREADS;
    #endif

    #ifdef VAR_TILING_FACTOR
    variable_parameter =  "tiling_factor";
    column = ::TILING_FACTOR;
    #endif

    #ifdef XEON
    variable_parameter = "xeon_for_threads+naive_threads_";
    #endif

    #ifdef VAR_NUM_GANGS
    variable_parameter += "_num_gangs_";
    column = ::NUM_GANGS;
    #endif

    #ifdef VAR_NUM_WORKERS
    variable_parameter += "_num_workers_";
    column = ::NUM_WORKERS;
    #endif

    #ifdef VAR_VECTOR_LENGTH
    variable_parameter += "_vector_length_";
    column = ::VECTOR_LENGTH;
    #endif

    #ifdef SAVE
    if(::VERSION == 1) {
        std::string resultOutput = "results/result_" + std::to_string(size) + ".txt";
        saveResult(resultOutput, polynomial, size);
    }
    #endif

    #ifdef CHECKSUM
        std::string resultOutput = "results/result_" + std::to_string(size) + ".txt";
        checksum(resultOutput, polynomial, size);
    #endif

    if(!::outfile.is_open())
        ::outfile.open("times/" + name + variable_parameter + std::to_string(size) + "deg_times.txt");


    PRINT << "Elapsed time: " << elapsed << std::endl;
    ::outfile<<column<<" "<<elapsed<<::std::endl;
}

// set global parameter on the value
void setGlobalParameter(int i){
    #ifdef VAR_NAIVE_LIMIT
    ::NAIVE_LIMIT = static_cast<int>(pow(2, i));
    #endif

    #ifdef VAR_TASK_LIMIT
    ::TASK_LIMIT = i;
    #endif

    #ifdef VAR_FOR_THREADS
    ::FOR_THREADS =  ::threads[i];
    #endif

    #ifdef VAR_FOR_LIMIT
    ::PARALLEL_FOR_LIMIT = i;
    #endif

    #ifdef VAR_NAIVE_THREADS
    ::NAIVE_THREADS = ::threads[i];
    #endif

    #ifdef VAR_TILING_FACTOR
    ::TILING_FACTOR = static_cast<int>(pow(2, i));
    #endif

    #ifdef XEON
    ::NAIVE_THREADS = static_cast<int>(pow(2, i));
		    ::FOR_THREADS = static_cast<int>(pow(2, i));
    #endif

    #ifdef VAR_NUM_GANGS
    ::NUM_GANGS = static_cast<int>(pow(2, i));
    #endif

    #ifdef VAR_NUM_WORKERS
    ::NUM_WORKERS = static_cast<int>(pow(2, i));
    #endif

    #ifdef VAR_VECTOR_LENGTH
    ::VECTOR_LENGTH = static_cast<int>(pow(2, i));
    #endif

}

//---------------------------------------------------------

// sequence version of naive algorithm
TYPE* seq_naive(const TYPE *A, const TYPE *B, int size_A, int size_B){
    auto *result = new TYPE[size_A + size_B - 1];
    for(int i=0; i<size_A + size_B - 1; i++) result[i] = 0;

    for (int i = 0; i < size_A; i++) {
        // VECTORIZED
        for (int j = 0; j < size_B; j++) {
            result[i + j] += A[i] * B[j];
        }
    }

    return result;
}

// sequence recursive version of Karatsuba algorithm
TYPE* seq_karatsuba(const TYPE *A, const TYPE *B, int size, int depth) {
    TYPE *__restrict__ lowA, *__restrict__ highA, *__restrict__ lowB, *__restrict__ highB, *__restrict__ midA, *__restrict__ midB;

    // when the size of polynomial is bellow the limit, omp_naive algorithm is called
    if (size <= ::NAIVE_LIMIT)
        return seq_naive(A, B, size, size);

    // compute the half for splitting the polynomial
    // if size is odd number
    int half = size / 2;

    if(size % 2 == 1)
        half++;

    // prepare arrays for splitted parts
    lowA = new TYPE[half];  lowB = new TYPE[half];
    midA = new TYPE[half];  midB = new TYPE[half];
    highA = new TYPE[half]; highB = new TYPE[half];

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

    z0 = seq_karatsuba(lowA, lowB, half, depth + 1);
    z1 = seq_karatsuba(midA, midB, half, depth + 1);
    z2 = seq_karatsuba(highA, highB, half, depth + 1);

    // init the result array
    auto *__restrict__ result = new TYPE[2*size-1];
    for(int i=0; i<2*size-1; i++) result[i] = 0;

    // compute the result
    // VECTORIZED
    for (int i = 0; i < 2 * half - 1; i++)
        result[i + 2 * half] += z2[i];

    // VECTORIZED
    for (int i = 0; i < 2 * half - 1; i++)
        result[i + half] += z1[i] - z2[i] - z0[i];

    // VECTORIZED
    for (int i = 0; i < 2 * half - 1; i++)
        result[i] += z0[i];

    return result;
}

//---------------------------------------------------------

// parallel version of naive algorithm for OpenMP
TYPE* omp_naive(const TYPE *A, const TYPE *B, int size_A, int size_B){
    auto *__restrict__ result = new TYPE[size_A + size_B - 1];
    for(int i=0; i<size_A+size_B-1; i++) result[i] = 0;

    // omp_naive multiplication
    #ifdef NAIVE_PARALELISM
            #pragma omp parallel for collapse(2) shared(result) num_threads(::NAIVE_THREADS)
            for (int i = 0; i < size_A; i += ::TILING_FACTOR) {
                for (int j = 0; j < size_B; j += ::TILING_FACTOR) {
		            int boundA = min(size_A, ::TILING_FACTOR + i);
                    int boundB = min(size_B, ::TILING_FACTOR + j);
                    for (int m = i; m < boundA; m++)
                        for (int n = j; n < boundB; n++)
                            result[m + n] += A[m] * B[n];
                }
            }
    #else
        for (int i = 0; i < size_A; i++) {
                // VECTORIZED
                for (int j = 0; j < size_B; j++) {
                    result[i + j] += A[i] * B[j];
                }
            }
    #endif

    return result;
}

// parallel recursive version of Karatsuba algorithm for OpenMP
TYPE* omp_karatsuba(const TYPE *A, const TYPE *B, int size, int depth) {
    TYPE *__restrict__ lowA, *__restrict__ highA, *__restrict__ lowB, *__restrict__ highB, *__restrict__ midA, *__restrict__ midB;

    // when the size of polynomial is bellow the limit, omp_naive algorithm is called
    if (size <= ::NAIVE_LIMIT)
        return omp_naive(A, B, size, size);

    // compute the half for splitting the polynomial
    // if size is odd number
    int half = size / 2;

    if(size % 2 == 1)
        half++;

    // prepare arrays for splitted parts
    lowA = new TYPE[half];  lowB = new TYPE[half];
    midA = new TYPE[half];  midB = new TYPE[half];
    highA = new TYPE[half]; highB = new TYPE[half];

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
    #endif
    z0 = omp_karatsuba(lowA, lowB, half, depth + 1);

    #ifdef TASK_PARALLELISM
    #pragma omp task shared(z1) if(half > ::TASK_LIMIT)
    #endif
    z1 = omp_karatsuba(midA, midB, half, depth + 1);

    // job for calling thread
    z2 = omp_karatsuba(highA, highB, half, depth + 1);

    #ifdef TASK_PARALLELISM
    #pragma omp taskwait
    #endif

    // init the result array
    auto *__restrict__ result = new TYPE[2*size-1];
    for(int i=0; i<2*size-1; i++) result[i] = 0;

    // compute the result

    #ifdef FOR_PARALLELISM

        int i;
        #pragma omp parallel shared(z0,z1,z2,result) private(i) num_threads(::FOR_THREADS) if(depth <= ::PARALLEL_FOR_LIMIT)
        {

            //VECTORIZED
            #pragma omp for
            for (i = 0; i < 2 * half - 1; i++)
                result[i + 2 * half] += z2[i];

            // VECTORIZED
            #pragma omp for
            for (i = 0; i < 2 * half - 1; i++)
                result[i + half] += z1[i] - z2[i] - z0[i];

            // VECTORIZED
            #pragma omp for
            for (i = 0; i < 2 * half - 1; i++)
                result[i] += z0[i];
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

//---------------------------------------------------------

// parallel version of naive algorithm for OpenACC
TYPE* acc_naive(const TYPE *A, const TYPE *B, int size_A, int size_B){
    TYPE *__restrict__ result = new TYPE[size_A + size_B - 1];

	for(int i=0; i<size_A + size_B - 1; i++) result[i] = 0;

	#pragma acc kernels num_gangs(1024) num_workers(32) copy(result[0:2*size_A-1]) copyin(A[0:size_A], B[0:size_B])
	{
        #pragma acc loop gang
        for (int i = 0; i < size_A; i++) {
        	#pragma acc loop worker
            for (int j = 0; j < size_B; j++) {
                result[i + j] += A[i] * B[j];
            }
        }
    }

    return result;
}


// parallel version of Karatsuba algorithm for OpenACC
TYPE* acc_karatsuba(const TYPE *A, const TYPE *B, int size, int depth) {
    TYPE *__restrict__ lowA, *__restrict__ highA, *__restrict__ lowB, *__restrict__ highB, *__restrict__ midA, *__restrict__ midB;

    // when the size of polynomial is bellow the limit, omp_naive algorithm is called
    if (size <= ::NAIVE_LIMIT){
        return seq_naive(A, B, size, size);
		//return acc_naive(A, B, size, size);
	}

    // compute the half for splitting the polynomial
    // if size is odd number
    int half = size / 2;

    if(size % 2 == 1)
        half++;

    // prepare arrays for splitted parts
    lowA = new TYPE[half];  lowB = new TYPE[half];
    midA = new TYPE[half];  midB = new TYPE[half];
    highA = new TYPE[half]; highB = new TYPE[half];

    // init low coefficients to new arrays  - VECTORIZED
   	for (int i = 0; i < half; i++) {
    	lowA[i] = A[i];
        lowB[i] = B[i];
    }


   	// init high coefficients
   	for (int i = half; i < size; i++) {
   		highA[i - half] = A[i];
   		highB[i - half] = B[i];
   	}

	// init mid coefficients
    for (int i = 0; i < half; i++) {
    	midA[i] = lowA[i] + highA[i];
        midB[i] = lowB[i] + highB[i];
    }


    TYPE *__restrict__ z0, *__restrict__ z1, *__restrict__ z2;

    // compute the parts of result
    z0 = acc_karatsuba(lowA, lowB, half, depth + 1);
    z1 = acc_karatsuba(midA, midB, half, depth + 1);
    z2 = acc_karatsuba(highA, highB, half, depth + 1);

    // init the result array
    auto *__restrict__ result = new TYPE[2*size-1];
    for(int i=0; i<2*size-1; i++) result[i] = 0;

	if(depth < 3)
	{
    	// compute the result
    	#pragma acc data copy(result[0 : 2 * size - 1])
	    {
    		#pragma acc parallel vector_length(::VECTOR_LENGTH)
    		{
    			#pragma acc loop vector
    			for (int i = 0; i < 2 * half - 1; i++)
	        		result[i + 2 * half] += z2[i];

    			#pragma acc loop vector
    			for (int i = 0; i < 2 * half - 1; i++)
	        		result[i + half] += z1[i] - z2[i] - z0[i];

    			#pragma acc loop vector
    			for (int i = 0; i < 2 * half - 1; i++)
        			result[i] += z0[i];
        	}
		}
	}
	else
	{
    	for (int i = 0; i < 2 * half - 1; i++)
        	result[i + 2 * half] += z2[i];
    	for (int i = 0; i < 2 * half - 1; i++)
        	result[i + half] += z1[i] - z2[i] - z0[i];
       	for (int i = 0; i < 2 * half - 1; i++)
        	result[i] += z0[i];

	}
    return result;
}

//---------------------------------------------------------

//  iterative version of Karatsuba algorithm
TYPE* iter_karatsuba(const TYPE *A, const TYPE *B, TYPE size){
    // https://eprint.iacr.org/2006/224.pdf

    // create empty coefficient vector with proper size and fill it with 0
    auto *__restrict__ result = new TYPE[2 * (size - 1)];


    // fill Di vector with Ai * Bi
    TYPE *__restrict__ D = new TYPE[2*size-1];
    for( TYPE i = 0; i < size; i++ ) D[i] = A[i] * B[i];

    // set the first and last coefficients
    result[0] = D[0];
    result[2 * (size - 1)] = D[size - 1];

    // for all coefficients of result vector
    for (TYPE position=1; position < 2* (size-1); position++){
        // for even coefficient add Di/2
        if ( position % 2 == 0)
            result[position] += D[position >> 1];

        // calculate start position in polynom
        TYPE start = (position >= size) ? (position % size ) + 1 : 0;

        // calculate end position in polynom
        TYPE end = (position + 1) >> 1;

        // inner loop: sum (Dst) - sum (Ds + Dt) where s+t=i
        for(TYPE inner = start; inner < end; inner++){
            result[position] += ( A[inner] + A[position - inner] ) * ( B[inner] + B[position - inner] );
            result[position] -= ( D[inner] + D[position-inner] );
        }
    }

    return result;
}

// iterative version of Karatsuba algorithm using the OpenACC directives
TYPE* acc_iter_karatsuba(const TYPE  *A,  const TYPE *B, TYPE size) {
    // create empty coefficient vector with proper size and fill it with 0
    TYPE *__restrict__ result = new TYPE[2 * size - 1];

    // fill Di vector with Ai * Bi
    TYPE *__restrict__ D = new TYPE[2 * size - 1];
    //#pragma acc parallel loop
    for (TYPE i = 0; i < size; i++) D[i] = A[i] * B[i];

    // set the first and last coefficients
    result[0] = D[0];
    result[2 * (size - 1)] = D[size - 1];

    //#pragma acc kernels num_gangs(::NUM_GANGS) num_workers(::NUM_WORKERS) copy (result[0:2*size-1]) copyin(A[0:size], B[0:size], D[2*size-1])
    //#pragma acc parallel copy (result[0:2*size-1]) copyin(A[0:size], B[0:size], D[2*size-1])
    //#pragma acc kernels num_workers(::NUM_WORKERS) copy(result[0:2*size-1]) copyin(A[0:size], B[0:size], D[2*size-1])
    //#pragma acc parallel vector_length(1024) copy(result[0:2*size-1]) copyin(A[0:size], B[0:size], D[0:2*size-1])
    #pragma acc parallel copy(result[0:2*size-1]) copyin(A[0:size], B[0:size], D[0:2*size-1])
    {
        #pragma acc loop gang
        for(TYPE position = 1; position < 2 * (size - 1); position++) {
            // for even coefficient add Di/2
            if (position % 2 == 0)
                result[position] += D[position >> 1];

            TYPE start = (position >= size) ? (position % size) + 1 : 0;
            TYPE end = (position + 1) >> 1;

            // inner loop: sum (Dst) - sum (Ds + Dt) where s+t=i
           	//#pragma acc loop
            for (TYPE inner = start; inner < end; inner++) {
                result[position] += (A[inner] + A[position - inner]) * (B[inner] + B[position - inner]);
                result[position] -= (D[inner] + D[position - inner]);
            }
        }
    }

	return result;
}


//---------------------------------------------------------


int main(int argc, char* argv[]) {

 	int size_A, size_B;
    TYPE *A = readPolynomial(argv[1], &size_A);
    TYPE *B = readPolynomial(argv[2], &size_B);

    for( int i = START_LIMIT; i <= END_LIMIT; i += STEP ) {
        setGlobalParameter(i);

        // Record start time
        auto start = omp_get_wtime();

        // multiply the polynomials
        TYPE *result;

        if (::VERSION==1) {
            PRINT << "Sequence version:" << std::endl;
            #ifdef KARATSUBA
                result = seq_karatsuba(A, B, size_A, 0);
            #else
                result = seq_naive(A, B, size_A, size_B);
            #endif
        }

        else if (::VERSION==2) {
            PRINT << "OpenMP version:" << std::endl;
            #ifdef KARATSUBA
                #ifdef TASK_PARALLELISM
                    #pragma omp parallel
                    #pragma omp single
                #endif
                result = omp_karatsuba(A, B, size_A, 0);
            #else
                result = omp_naive(A, B, size_A, size_B);
            #endif
        }

		else if(::VERSION==3) {
            PRINT << "OpenACC version:" << std::endl;
            #ifdef KARATSUBA
                result = acc_karatsuba(A, B, size_A, 0);
            #else
                result = acc_naive(A, B, size_A, size_B);
            #endif
        }

        else if(::VERSION==4) {
            PRINT << "Sequence iterative version:" << std::endl;
            #ifdef KARATSUBA
                result = iter_karatsuba(A, B, size_A);
            #else
                result = seq_naive(A, B, size_A, size_B);
            #endif

        } else if(::VERSION==5) {
            PRINT << "OpenACC iterative version:" << std::endl;
            #ifdef KARATSUBA
            result = acc_iter_karatsuba(A, B, size_A);
            #else
            result = acc_naive(A, B, size_A, size_B);
            #endif
        }

        // stop timer
        auto finish = omp_get_wtime();

        // Record end time
        double elapsed = finish - start;

        // print the result if DEBUG is defined
        #ifdef DEBUG
            print(result, size_A + size_B - 1);
        #endif

        printStats(elapsed, result, size_A);
    }

    ::outfile.close();
    return 0;
}
