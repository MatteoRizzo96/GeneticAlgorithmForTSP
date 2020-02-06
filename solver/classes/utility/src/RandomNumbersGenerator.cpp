#include <random>
#include "../include/RandomNumbersGenerator.h"

using namespace std;

int RandomNumbersGenerator::getRandomInt(int lb, int ub) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(lb, ub);
    return dis(gen);
}

double RandomNumbersGenerator::getRandomDouble(int lb, int ub) {
    double base = getRandomInt(lb, ub);
    return base == ub ? base : base + double(getRandomInt(0, 9)) / 10;
}
