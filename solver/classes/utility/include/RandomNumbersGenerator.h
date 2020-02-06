//
// Created by matteo on 28/12/2019.
//

#ifndef LABS_TSP_RANDOMNUMBERSGENERATOR_H
#define LABS_TSP_RANDOMNUMBERSGENERATOR_H

#include <algorithm>
#include <iostream>

using namespace std;

class RandomNumbersGenerator {

public:

    static int getRandomInt(int lb = 0, int ub = RAND_MAX);

    static double getRandomDouble(int lb = 0, int ub = RAND_MAX);
};


#endif //LABS_TSP_RANDOMNUMBERSGENERATOR_H
