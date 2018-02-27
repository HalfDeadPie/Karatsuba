//
// Created by halfdeadpie on 23.2.18.
//

#ifndef UNTITLED_POLYNOM_H
#define UNTITLED_POLYNOM_H

#include <vector>
#include <utility>

class Polynom {
private:
    unsigned long degree;
    std::vector<long> coefficients;

public:
    // constructor
    explicit Polynom(std::vector<long>, unsigned long);

    // degree
    unsigned long getDegree();
    void setDegree(unsigned long);

    // coefficients
    std::vector<long> getCoefficients();
    void setCoefficients(std::vector<long>);
    long getAt(unsigned long);
    unsigned long getSize();
};

#endif //UNTITLED_POLYNOM_H
