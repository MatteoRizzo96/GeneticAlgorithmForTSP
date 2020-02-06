#ifndef LABS_TSP_TSPSOLVERSCRIPT_H
#define LABS_TSP_TSPSOLVERSCRIPT_H


#include <utility>
#include "ProblemInstanceGenerator.h"
#include "../../solvers/include/TSPsolver.h"

class TSPsolverScript {
public:
    int numTrials;
    TSPsolver *solver;
    vector<int> size;
    vector<double> time;
    string problemName;
    string solverName;
    string fileName;
    string pathToSize = "../data_analysis/data/comparisons/size/";
    string pathToTrial = "../data_analysis/data/comparisons/trial/";
    string pathToSolutions = "../data_analysis/data/solutions/";
    ProblemInstanceGenerator problemGenerator;
    ProblemInstanceGenerator::ProblemType problem;

    pair<double, double> computeSolution(const vector<double> &costs, int size, TSPsolver *s = nullptr);

    pair<double, double> searchSolution(double t);

    void solveProblemInstances();

    void compareSolversByTrial(const vector<TSPsolver *> &solvers);

    void compareSolversBySize(const vector<TSPsolver *> &solvers);

    TSPsolverScript(TSPsolver *s = nullptr);
};


#endif //LABS_TSP_TSPSOLVERSCRIPT_H
