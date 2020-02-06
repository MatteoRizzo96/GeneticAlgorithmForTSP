#ifndef LABS_TSP_TabuSearchSolver_H
#define LABS_TSP_TabuSearchSolver_H

#include <vector>
#include <iostream>

#include "../../utility/include/TSPsolution.h"
#include "TSPsolver.h"

// Structure representing substring reversal move

typedef struct move {
    int from;
    int to;
} TSPMove;

// Solves a TSP problem by neighbourhood search and 2-opt moves

class TabuSearchSolver : public TSPsolver {
public:
    void solve() final;

    double getObjFunValue() final;

    TSPsolution getSolution() final;

    void setParameters(const vector<double> &costs, int size);

    void logSolution() {};

private:
    double infinite = 1e10;
    int tabuLength = 10;
    vector<int> tabuList;
    TSPsolution bestSol = TSPsolution(N + 1);

    void initTabuList(int n);

    static TSPsolution &swap(TSPsolution &tspSol, const TSPMove &move);

    [[nodiscard]] double evaluate(const TSPsolution &sol) const;

    double findBestNeighbor(const TSPsolution &currSol, int currIter, double aspiration, TSPMove &move);
};

#endif //LABS_TSP_TabuSearchSolver_H
