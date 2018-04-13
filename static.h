//
// Created by halfdeadpie on 7.4.18.
//
#ifndef UNTITLED_STATIC_H
#define UNTITLED_STATIC_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <chrono>
#include <omp.h>
#include <cmath>
#include <cstring>

//---------------------------------------------------------

/*
 * 1 - sequence version
 * 2 - OpenMP version
 * 3 - OpenACC version
 * 4 - sequence iterative version
 * 5 - OpenACC iterative version
 */
int VERSION = 5;

//---------------------------------------------------------

//#define XEON

//---------------------------------------------------------

// VARIABILE - only one!
//#define VAR_NAIVE_LIMIT
//#define VAR_TASK_LIMIT
//#define VAR_FOR_THREADS
//#define VAR_FOR_LIMIT
//#define VAR_NAIVE_THREADS
//#define VAR_TILING_FACTOR

//-----------------------------------------------------------

/* Main loop configuration - may be used for changing the selected parameter defining its "VAR" directives*/
#define START_LIMIT 0
#define END_LIMIT 0
#define STEP 1

//-----------------------------------------------------------

// Paralelisation modes
//#define NAIVE_PARALELISM
//#define TASK_PARALLELISM
//#define FOR_PARALLELISM

//-----------------------------------------------------------

/*If defined, Karatsuba algorithm is called. Otherwise, the omp_naive algorithm computes all multiplications*/
#define KARATSUBA

//-----------------------------------------------------------

/*Limit value that is compared to actual size of multiplicated polynomials. When their size is equal or bellow
 * this value, the polynomials are multiplicated using the omp_naive method*/
int NAIVE_LIMIT = 256;

//-----------------------------------------------------------

// Defines the limit for size of polynomial for creating new tasks
int TASK_LIMIT = 256;

//-----------------------------------------------------------

int FOR_THREADS = 12;
int PARALLEL_FOR_LIMIT = 3;

//-----------------------------------------------------------

int NAIVE_THREADS = 24;
int TILING_FACTOR = 1024;

//-----------------------------------------------------------

// Macro used for printing to standard output
#define PRINT std::cout

// type of coefficient
#define TYPE int

/*Debug the result of multiplication - if defined, the result is printed to standard output*/
#define DEBUG

std::ofstream outfile;

int threads [8] = { 1, 2, 4, 6, 8, 12, 24, 32 };

#endif //UNTITLED_STATIC_H
