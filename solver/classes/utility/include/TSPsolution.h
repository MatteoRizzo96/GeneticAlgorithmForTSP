#ifndef LABS_TSP_TSPSOLUTION_H
#define LABS_TSP_TSPSOLUTION_H

#include <vector>

using namespace std;

// TSP Solution representation: ordered sequence of nodes (path representation)

class TSPsolution {
public:
    vector<int> sequence;

    TSPsolution() {};

    explicit TSPsolution(int n);

    TSPsolution(const TSPsolution &tspSol);

    TSPsolution(const vector<int> &seq);

    static TSPsolution getRandomSolution(int n);

public:

    void print() const;

    TSPsolution &operator=(const TSPsolution &right);
};


#endif //LABS_TSP_TSPSOLUTION_H
