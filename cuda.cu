#include "static_cuda.h"

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

/* computes the checksum based on equality of coefficients on random positions
 * but the result must be saved first. The best way how to save result is to
 * define SAVE in static.h and run the naive sequent algorithm*/
void checksum(std::string &file, const TYPE *polynomial, TYPE size){
    // THIS IS NOT REAL CHECKSUM

	int realSize;
    int mySize = 2 * size - 1;
    TYPE *realResult = readPolynomial(file, &realSize);

    // check size
	if(realSize != mySize)
        PRINT<<"\n\n!!!WRONG CHECKSUM (SIZE)!!!\n\n";
	
	// check first coefficients
	if(polynomial[0] != realResult[0])
		PRINT<<PRINT<<"\n\n!!!WRONG CHECKSUM (FIRST COEFFICIENT)!!!\n\n";
	
	// check last coefficients
	if(polynomial[mySize - 1] != realResult[mySize - 1]){
		PRINT<<"\n\n!!!WRONG CHECKSUM (LAST COEFFICIENT)!!!\n";
		PRINT<<polynomial[mySize - 1]<<" "<<realResult[mySize - 1]<<std::endl<<std::endl;
	}            
	
	// check even coefficient
	int position = 0;
  	position = mySize / 2;
  	if(position % 2 == 1) position++;
  	if(polynomial[position] != realResult[position] || polynomial[mySize-3] != realResult[mySize-3]){
  		 PRINT<<"\n\n!!!WRONG CHECKSUM (EVEN COEFFICIENT)!!!\n\n";
  	}
        
	// check random coefficients
    for(int i = 0; i < NUM_OF_CHECKS; i++){
        position = rand() % mySize;
        if(polynomial[position] != realResult[position])
            PRINT<<"\n\n!!!WRONG CHECKSUM (RESULT)!!!\n\n";
    }
}

// print the time resuls to standard output and times-file output
void printStats(double elapsed, double computingTime, const TYPE *polynomial, int size){
    std::string variable_parameter, name;

   	variable_parameter += "cuda_";
    
    #ifdef KARATSUBA
    name = "karatsuba_";
    #else
    name = "naive_";
    #endif

    #ifdef CHECKSUM
        std::string resultOutput = "results/result_" + std::to_string(size) + ".txt";
        checksum(resultOutput, polynomial, size);
    #endif

    if(!::outfile.is_open())
        ::outfile.open("times/" 
        + name 
        + variable_parameter 
        + std::to_string(size) 
        + "deg_"
        + std::to_string(::BLOCK_X) + "_times.txt");


    PRINT << "Total/computing time: " << elapsed << " / " << computingTime << std::endl;
    ::outfile<<::BLOCK_X<<"\t"<<computingTime<<"\t"<<elapsed<<::std::endl;
}


// computes the number of needed blocks depending on the threads and size
int blockRatio(int size, int threads){
    return ( (size + (threads - 1)) / threads );
}

//---------------------------------------------------------

__global__ void kernel_naive(TYPE *A, TYPE *B, TYPE *result, TYPE size_A, TYPE size_B, TYPE resultSize){
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	int j = blockDim.y * blockIdx.y + threadIdx.y;
	
	if( i >= 0 && i < size_A) {
		if( j >= 0 && j < size_B ){
			atomicAdd( (result + i + j) , (A[i] * B[j]) ); 
		}
	}
}

// sequence version of naive algorithm
TYPE* cuda_naive(const TYPE *A, const TYPE *B, int size_A, int size_B, double *computingTime){
    int resultSize = size_A + size_B - 1;
    TYPE *__restrict__ hostResult = new TYPE[size_A + size_B - 1];
    
	int blockX = ::BLOCK_X;
	int blockY = ::BLOCK_Y;

	// blocks and grid initialisation
    dim3 block( blockX, blockY );
    dim3 grid( blockRatio(resultSize, blockX) , blockRatio(resultSize, blockY) );

    TYPE *__restrict__ devA, *__restrict__ devB, *__restrict__ devResult;
    cudaMalloc((void**)&devA, size_A * sizeof(TYPE));
    cudaMalloc((void**)&devB, size_B * sizeof(TYPE));
    cudaMalloc((void**)&devResult, resultSize * sizeof(TYPE)); 
    
   	// hostA -> devA
	cudaMemcpy(devA, A, size_A * sizeof(TYPE), 
		cudaMemcpyHostToDevice);
	
	// hostB -> devB
	cudaMemcpy(devB, B, size_B * sizeof(TYPE), 
		cudaMemcpyHostToDevice);

	// start the timer without copying data
	auto start = std::chrono::high_resolution_clock::now();
	
	// K E R N E L S
	kernel_naive <<< grid, block >>> (devA, devB, devResult, size_A, size_B, resultSize);
	
	cudaDeviceSynchronize();

	// Record end time without copying data
	auto finish = std::chrono::high_resolution_clock::now();
	        
	// Record end time
	std::chrono::duration<double> elapsed = finish - start;
	
	//PRINT<<"Elapsed time (no-copy): "<<elapsed.count()<<std::endl;
	*computingTime = elapsed.count();

	// devResult -> hostResult
	cudaMemcpy(hostResult, devResult, resultSize * sizeof(TYPE), 
		cudaMemcpyDeviceToHost);
	    
    return hostResult;
}

//---------------------------------------------------------

// compute the D according to iterative algorithm
__global__ void kernel_D(TYPE *A, TYPE *B, TYPE *D, TYPE size){
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < size)
    	D[i] = A[i] * B[i];
}

// set the first and last coefficients of the result
__global__ void kernel_res_bounds(TYPE *result, TYPE *D, TYPE size, TYPE resultSize){
	//result[resultSize - 1] = D[size - 1];
	result[0] = D[0];
	result[resultSize - 1] = D[size - 1];
}

// addition to even coefficients according to iterative algorithm
__global__ void kernel_res_even(TYPE *result, TYPE *D, TYPE resultSize){
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if(i % 2 == 0 && i > 0 && i < resultSize - 1)
		atomicAdd(result + i, D[ i / 2 ]);
}

// the main computing of the result using for-loop (this seems to be faster)
__global__ void kernel_res_main(TYPE *A, TYPE *B, TYPE *D, TYPE *result, TYPE size, TYPE resultSize){
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	if( i > 0 && i < resultSize - 1){
        TYPE start = (i >= size) ? (i % size ) + 1 : 0;
        TYPE end = (i + 1) / 2;
  
        for(TYPE inner = start; inner < end; inner++){
            atomicAdd(result + i, (  (A[inner] + A[i - inner]) * (B[inner] + B[i - inner]) ) );
            atomicSub(result + i, ( D[inner] + D[i-inner] ) );
        }		
	}
}


// the main computing of the result using grid (this seems to be slower)
__global__ void kernel_res_nested(TYPE *A, TYPE *B, TYPE *D, TYPE *result, TYPE size, TYPE resultSize){
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	int j = blockDim.y * blockIdx.y + threadIdx.y;    
    
	if(i>=resultSize-1 || j>size) return;

	TYPE start = (i >= size) ? (i % size ) + 1 : 0;
    TYPE end = (i + 1) >> 1;
    
	if( i<resultSize-1 && j>=start && j<end ){
		atomicAdd(result+i, ( (A[j] + A[i - j] ) * ( B[j] + B[i - j]) ) );
		atomicSub(result+i,  ( D[j] + D[i - j] ) );
	}
	
}

//  iterative version of Karatsuba algorithm
TYPE* cuda_karatsuba(const TYPE *A, const TYPE *B, TYPE size, double *computingTime){
 	// compute the size of the result
	TYPE resultSize = 2 * size - 1;
	
	// the number of threads for one dimensional computing
	int threadsOne = ::NUM_THREADS;

	// dimension of blocks for grid
	int blockX = ::BLOCK_X;
	int blockY = ::BLOCK_Y;

	// blocks and grid initialisation
    dim3 block( blockX, blockY );
    dim3 grid( blockRatio(resultSize, blockX) , blockRatio(resultSize, blockY) );

	    
    // create empty coefficient vector with proper size and fill it with 0
    TYPE *__restrict__ hostResult = new TYPE[resultSize];
	TYPE *__restrict__ hostD = new TYPE[size];
	
	// create and allocate space for polynomials and the result
	TYPE *__restrict__ devA, *__restrict__ devB, *__restrict__ devD, *__restrict__ devResult;
	cudaMalloc((void**)&devA, size * sizeof(TYPE));
	cudaMalloc((void**)&devB, size * sizeof(TYPE));
	cudaMalloc((void**)&devD, size * sizeof(TYPE));
	cudaMalloc((void**)&devResult, resultSize * sizeof(TYPE));
	
	// hostA -> devA
	cudaMemcpy(devA, A, size * sizeof(TYPE), 
		cudaMemcpyHostToDevice);
	// hostB -> devB
	cudaMemcpy(devB, B, size * sizeof(TYPE),
		cudaMemcpyHostToDevice);

	// start the timer without copying data
	auto start = std::chrono::high_resolution_clock::now();

	// K E R N E L S	
	
	// compute the D
	kernel_D <<<blockRatio(size, threadsOne) , threadsOne>>> (devA, devB, devD, size);
	
	// ini the first and last coefficent of the result
	kernel_res_bounds <<<1, 1>>> (devResult, devD, size, resultSize);
	
	// add to even coefficients of result
	kernel_res_even <<<blockRatio(resultSize, threadsOne), threadsOne>>> (devResult, devD, resultSize);
	
	// compute the main part of the result 
	kernel_res_main <<<blockRatio(resultSize, threadsOne), threadsOne>>> (devA, devB, devD, devResult, size, resultSize);
	//kernel_res_nested <<<grid, block>>> (devA, devB, devD, devResult, size, resultSize);

	// wait for the kernels
	cudaDeviceSynchronize();

	// Record end time without copying data
	auto finish = std::chrono::high_resolution_clock::now();
	        
	// Record end time
	std::chrono::duration<double> elapsed = finish - start;
	
	//PRINT<<"Elapsed time (no-copy): "<<elapsed.count()<<std::endl;
	*computingTime = elapsed.count();

	// devResult -> hostResult
	cudaMemcpy(hostResult, devResult, resultSize * sizeof(TYPE), 
		cudaMemcpyDeviceToHost);
	    
    return hostResult;
}


//---------------------------------------------------------


int main(int argc, char* argv[]) {

 	int size_A, size_B;
    TYPE *A = readPolynomial(argv[1], &size_A);
    TYPE *B = readPolynomial(argv[2], &size_B);
	double computingTime;

    // Record start time
    auto start = std::chrono::high_resolution_clock::now();

    // multiply the polynomials
    TYPE *result;

	// M U L T I P LI C A T I O N
	#ifdef KARATSUBA
    PRINT<<"Karatsuba version: CUDA"<<std::endl;
    result = cuda_karatsuba(A, B, size_A, &computingTime);
    #else
    PRINT<<"Naive version: CUDA"<<std::endl;
    result = cuda_naive(A, B, size_A, size_B, &computingTime);
    #endif

    // Record end time
    auto finish = std::chrono::high_resolution_clock::now();

    // Record end time
    std::chrono::duration<double> elapsed = finish - start;

    // print the result if DEBUG is defined
    #ifdef DEBUG
    print(result, 2 * size_A - 1);
    #endif

    printStats(elapsed.count(), computingTime, result, size_A);

    ::outfile.close();
    return 0;
}
