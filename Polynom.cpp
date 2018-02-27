//
// Created by halfdeadpie on 23.2.18.
//

#include "Polynom.h"

#include <utility>

// constructor
Polynom::Polynom(std::vector<long> coef, unsigned long deg) {
    coefficients = std::move(coef);
    degree = deg;
}

// get------------------------------------
// return degree
unsigned long Polynom::getDegree() {
    return degree;
}

// return coefficients
std::vector<long> Polynom::getCoefficients(){
    return coefficients;
};

long Polynom::getAt(unsigned long position){
    return Polynom::getCoefficients()[position];
}

unsigned long Polynom::getSize(){
    return Polynom::getCoefficients().size();
}

// set------------------------------------
// set degree
void Polynom::setDegree(unsigned long deg){
    degree = deg;
}

// set coeffiecints
void Polynom::setCoefficients(std::vector<long> coef){
    coefficients = std::move(coef);
}
