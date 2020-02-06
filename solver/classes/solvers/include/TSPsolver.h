#ifndef LABS_TSP_TSPSOLVER_H
#define LABS_TSP_TSPSOLVER_H


#include <map>
#include <vector>
#include "../../utility/include/TSPsolution.h"

using namespace std;

class TSPsolver {

public:
    int N;
    vector<double> C;
    vector<vector<double>> M;
    bool log;

    virtual void solve() = 0;

    virtual void logSolution() = 0;

    virtual void setParameters(const vector<double> &costs, int size) = 0;

    virtual double getObjFunValue() = 0;

    virtual TSPsolution getSolution() = 0;

    virtual ~TSPsolver() = default;

protected:
    void costsToMatrix() {
        M = vector<vector<double>>(N, vector<double>(N));
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) M[i][j] = C[i * N + j];
    }
};


#endif //LABS_TSP_TSPSOLVER_H
