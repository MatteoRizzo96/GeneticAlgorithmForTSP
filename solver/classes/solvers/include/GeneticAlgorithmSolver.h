#ifndef LABS_TSP_GENETICALGORITHMSOLVER_H
#define LABS_TSP_GENETICALGORITHMSOLVER_H


#include "TSPsolver.h"

class GeneticAlgorithmSolver : public TSPsolver {
public:
    void solve() final;

    void logSolution() final {};

    void setParameters(const vector<double> &costs, int size) final;

    double getObjFunValue() final;

    TSPsolution getSolution() final;

private:
    int K;
    int maxIter;
    int tolerance;
    int eliteSize;
    int populationSize;
    double mutationPercentage;
    double mutationRate;
    vector<pair<int, double>> fitness;
    vector<TSPsolution> population;
    pair<int, double> bestSol;

    void nextGeneration();

    void breed();

    void mutate();

    void educate();

    void updateFitness();

    double calculateFitness(const TSPsolution &sol) const;

    static TSPsolution twoOptSwap(vector<int> s, int i, int j);

    TSPsolution tournamentSelection();

    TSPsolution mutateIndividual(TSPsolution individual);

    TSPsolution recombine(const TSPsolution &i1, const TSPsolution &i2);
};


#endif //LABS_TSP_GENETICALGORITHMSOLVER_H
