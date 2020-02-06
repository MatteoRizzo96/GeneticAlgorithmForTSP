#include "../include/GeneticAlgorithmSolver.h"
#include "../../utility/include/RandomNumbersGenerator.h"
#include "../../utility/include/Params.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <utility>

using namespace std;
using namespace rapidjson;

// -- Public --

void GeneticAlgorithmSolver::solve() {
    // Compute the initial generation
    nextGeneration();
    pair<int, double> initialSol = fitness[0];
    cout << "\nThe initial solution is: ";
    population[initialSol.first].print();
    cout << "\n The initial objective function value is: " << 1 / initialSol.second << endl;

    // Iterate the breeding of the population for the given number of generations and tolerance
    vector<pair<int, double>> iterations;
    int toleranceThreshold = tolerance;
    for (int i = 0; i < maxIter && tolerance; i++) {
        // Get the best solution in the new generation (i.e. a pair "sol id"-"fitness score")
        nextGeneration();
        pair<int, double> currentSol = fitness[0];
        iterations.push_back(currentSol);
        cout << "\n * Iteration (" << i + 1 << ") - Objective function has value: " << 1 / currentSol.second << endl;

        // If no improvement has been made, decrease the tolerance else update the best values
        if (currentSol.second <= bestSol.second) {
            tolerance--;
        } else {
            bestSol = currentSol;
            tolerance = toleranceThreshold;
        }
    }

    cout << "\nThe final solution is: ";
    population[bestSol.first].print();
    cout << "\nThe objective function value has decreased from: " << 1 / initialSol.second << " to "
         << 1 / bestSol.second << "\n" << endl;
}

void GeneticAlgorithmSolver::setParameters(const vector<double> &costs, int size) {
    // Initialize the main parameters
    C = costs;
    N = size;
    bestSol = make_pair(0, 0.0);
    costsToMatrix();

    // Initialize the specific parameters of the solver from json
    Document d = Params::getAllParams();
    K = max(2, int(N * d["genetic_algorithm"]["k"].GetDouble()));
    maxIter = d["genetic_algorithm"]["max_iter"].GetInt();
    tolerance = d["genetic_algorithm"]["tolerance"].GetInt();
    mutationRate = d["genetic_algorithm"]["mutation_rate"].GetDouble();
    populationSize = d["genetic_algorithm"]["population_size"].GetInt();
    eliteSize = int(populationSize * d["genetic_algorithm"]["elite_percentage"].GetDouble());
    mutationPercentage = d["genetic_algorithm"]["mutation_percentage"].GetDouble();
    fitness = vector<pair<int, double>>(populationSize);

    // Randomly initialize the population
    population = vector<TSPsolution>(populationSize);
    for (auto &i : population) i = TSPsolution::getRandomSolution(N);

    // Initialize the fitness scores
    updateFitness();
}

double GeneticAlgorithmSolver::getObjFunValue() { return 1 / fitness[0].second; }

TSPsolution GeneticAlgorithmSolver::getSolution() { return population[bestSol.first]; }

// -- Private --

/** Evaluate the fitness of a solution in O(N)
  * @param sol: solution to be evaluated
  * @return the evaluation of the solution
  */
double GeneticAlgorithmSolver::calculateFitness(const TSPsolution &sol) const {
    double total = 0.0;
    for (unsigned i = 0; i < sol.sequence.size() - 1; i++) total += M[sol.sequence[i]][sol.sequence[i + 1]];
    return total == 0 ? total : double(1 / total);
}

/** Update the fitness scores of all the elements in the population and sort them in descending order in O(populationSize * N) */
void GeneticAlgorithmSolver::updateFitness() {
    for (int i = 0; i < populationSize; i++) fitness[i] = pair<int, double>(i, calculateFitness(population[i]));
    sort(fitness.begin(), fitness.end(), [](const auto &a, const auto &b) { return a.second > b.second; });
}

/** Ordered crossover recombination in O(2N). Example:
  *
  *     idx  | 0 1 2 3 4 5
  *     ------------------
  *     ind1 | 0 1 3 2 4 0
  *     ind2 | 0 4 1 3 2 0
  *     ------------------
  *
  * idx: minIdx = 3, maxIdx = 4 => crossover seq 2 4
  * Output: 0 1 3 2 4 0
  *               ^ ^
  * @param i1: the first individual to recombine
  * @param i2: the second individual to recombine
  * @return a child solution that is the recombination of its parents
  */
TSPsolution GeneticAlgorithmSolver::recombine(const TSPsolution &i1, const TSPsolution &i2) {
    int idxMin = RandomNumbersGenerator::getRandomInt(1, N - 2);
    int idxMax = RandomNumbersGenerator::getRandomInt(idxMin + 1, N - 1);
    vector<int> childSequence = i1.sequence;
    vector<int> seq = vector<int>(childSequence.begin() + idxMin, childSequence.begin() + idxMax + 1);
    int childIdx = 1;
    for (int i = 1; i < N; i++) {
        if (childIdx >= idxMin && childIdx <= idxMax) childIdx += idxMax - idxMin + 1;
        if (std::find(seq.begin(), seq.end(), i2.sequence[i]) == seq.end()) {
            childSequence[childIdx] = i2.sequence[i];
            childIdx++;
        }
    }
    return TSPsolution(childSequence);
}

/** K-tournament selection in O(K)
  * @return: the member of the population that won the tournament
  * */
TSPsolution GeneticAlgorithmSolver::tournamentSelection() {
    pair<int, double> winner = fitness.back();
    for (int i = 0; i < K; i++) {
        auto participant = fitness[RandomNumbersGenerator::getRandomInt(0, populationSize - 1)];
        if (participant.second > winner.second) winner = participant;
    }
    return population[winner.first];
}

/** Recombine the individuals and update the population in O(populationSize * (2K + 2N)).
  * The fittest individuals of the previous generation are preserved via elitism
  */
void GeneticAlgorithmSolver::breed() {
    for (int i = 0; i < eliteSize; i++) population[i] = population[fitness[i].first];
    for (int i = eliteSize; i < populationSize; i++)
        population[i] = recombine(tournamentSelection(), tournamentSelection());
}

/** Swap mutates the given individual according to the given mutation rate in O(N).
  * @param individual: the individual to be mutated
  * @return the mutated individual
  * */
TSPsolution GeneticAlgorithmSolver::mutateIndividual(TSPsolution individual) {
    for (int i = 1; i < N; i++) {
        if (RandomNumbersGenerator::getRandomInt(0, 100) / 100.0 < mutationRate) {
            int j = RandomNumbersGenerator::getRandomInt(1, N - 1);
            int a = individual.sequence[i];
            int b = individual.sequence[j];
            individual.sequence[i] = b;
            individual.sequence[j] = a;
        }
    }
    return individual;
}

/** Swap mutate some individuals in the population according to the given mutation percentage in O(populationSize * N).
  * This is done to contrast genetic drift. Note that the mutated individuals will be at most
  * mutationPercentage of the population size, but could be less since some individual may mutate more than once
  * */
void GeneticAlgorithmSolver::mutate() {
    for (int i = 0; i < populationSize * mutationPercentage; i++) {
        int idx = RandomNumbersGenerator::getRandomInt(eliteSize, populationSize - 1);
        population[idx] = mutateIndividual(population[idx]);
    }
}

/** Perform a 2-opt swap on the given solution in O(N)
 * @param s: a sequence of nodes (TSPsolution.sequence)
 * @param i: the min index of the swap sequence
 * @param j: the max index of the swap sequence
 * @return a 2-opt swapped solution
 */
TSPsolution GeneticAlgorithmSolver::twoOptSwap(vector<int> s, int i, int j) {
    vector<int> seq = vector<int>(std::move(s));
    reverse(seq.begin() + i, seq.begin() + j);
    return TSPsolution(seq);
}

/** Educate the population performing a local search in O(populationSize * N^3)*/
void GeneticAlgorithmSolver::educate() {
    for (int k = 0; k < populationSize; k++)
        for (int i = 1; i < N - 1; i++)
            for (int j = i + 1; j < N - 1; j++) {
                TSPsolution sol = twoOptSwap(population[k].sequence, i, j);
                if (calculateFitness(sol) > calculateFitness(population[k])) population[k] = sol;
            }
}

/** Compute the new generation by breeding and mutating */
void GeneticAlgorithmSolver::nextGeneration() {
    breed();
    educate();
    mutate();
    updateFitness();
}
