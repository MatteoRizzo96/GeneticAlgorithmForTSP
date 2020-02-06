#ifndef LAB01_TSP_TSPSOLVER_H
#define LAB01_TSP_TSPSOLVER_H

#include <map>
#include "../../utility/include/ProblemInstanceGenerator.h"
#include "TSPsolver.h"
#include "../../utility/include/cpxmacro.h"

class CPLEXsolver : public TSPsolver {

public:
    Env environment;
    Prob problem;

    explicit CPLEXsolver(const vector<double> &costs = {}, int size = 0, bool logSolutions = true);

    void setParameters(const vector<double> &costs, int size) final;

    void solve() final;

    void logSolution() final;

    double getObjFunValue() final;

    TSPsolution getSolution() final;

    ~CPLEXsolver() final;

private:

    void setupConstraints(const pair<map<pair<int, int>, int>, map<pair<int, int>, int>> &);

    pair<map<pair<int, int>, int>, map<pair<int, int>, int>> setupVariables();
};


#endif //LAB01_TSP_TSPSOLVER_H
