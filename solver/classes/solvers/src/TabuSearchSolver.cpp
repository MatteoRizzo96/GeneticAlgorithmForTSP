#include "../include/TabuSearchSolver.h"
#include "../../utility/include/RandomNumbersGenerator.h"
#include <iostream>

// -- PUBLIC --

void TabuSearchSolver::solve() {
    try {
        tabuList.reserve(N);
        initTabuList(N);

        TSPMove move;
        TSPsolution currSol(TSPsolution::getRandomSolution(N));

        double bestValue, currValue;
        bestValue = currValue = evaluate(currSol);

        bool stop = false;
        int maxIter = 100;

        for (int iter = 1; iter < maxIter && !stop; iter++) {

            if (N < 20) currSol.print();

            cout << " (" << iter << ") value " << currValue << "\t(" << evaluate(currSol) << ")";

            // Aspired improvement (to improve over bestValue)
            double aspiration = bestValue - currValue;

            // Calculate the cost with respect to the best neighbor
            double bestNeighValue = currValue + findBestNeighbor(currSol, iter, aspiration, move);

            // Stop because all neighbors are tabu
            if (bestNeighValue >= infinite) {
                cout << "\tmove: NO legal neighbour" << endl;
                stop = true;
                continue;
            }

            cout << "\tmove: " << move.from << " , " << move.to;

            // Update the tabu list
            tabuList[currSol.sequence[move.from]] = iter;
            tabuList[currSol.sequence[move.to]] = iter;

            currSol = swap(currSol, move);
            currValue = bestNeighValue;

            // Update the incumbent solution if sensibly better than the current
            if (currValue < bestValue - 0.01) {
                bestValue = currValue;
                bestSol = currSol;
                cout << "\t***";
            }
            cout << endl;
        }
    }
    catch (exception &e) { cout << ">> EXCEPTION: " << e.what() << endl; }
}

void TabuSearchSolver::setParameters(const vector<double> &costs, int size) {
    C = costs;
    N = size;
    costsToMatrix();
}

TSPsolution TabuSearchSolver::getSolution() { return bestSol; }

double TabuSearchSolver::getObjFunValue() { return evaluate(bestSol); }

// -- PRIVATE --

/** Perform a swap move (corresponding to 2-opt)
  * @param tspSol: solution to be perturbed
  * @param move: move to perform
  * @return (into param tspSol) the perturbed solution
  */
TSPsolution &TabuSearchSolver::swap(TSPsolution &tspSol, const TSPMove &move) {
    TSPsolution tmpSol(tspSol);
    for (int i = move.from; i <= move.to; ++i) tspSol.sequence[i] = tmpSol.sequence[move.to - (i - move.from)];
    return tspSol;
}

/** Explore the neighborhood. Determine the NON-TABU (or satisfying aspiration) move yielding
  * the best 2-opt neighbor solution. Aspiration criteria: 'neighCostVariation' better than
  * 'aspiration' (notice that 'aspiration' has been set such that if 'neighCostVariation'
  * is better than 'aspiration' than we have a new incumbent solution)
  * @param tsp: TSP data
  * @param currSol: center solution
  * @return (as side effect to move) the selected move (steepest descent strategy)
  * @return the incremental cost with respect to currSol
  */
double TabuSearchSolver::findBestNeighbor(const TSPsolution &currSol,
                                                   int currIter,
                                                   double aspiration,
                                                   TSPMove &move) {
    double bestCostVariation = infinite;

    for (uint a = 1; a < currSol.sequence.size() - 2; a++) {
        int h = currSol.sequence[a - 1];
        int i = currSol.sequence[a];
        for (uint b = a + 1; b < currSol.sequence.size() - 1; b++) {
            int j = currSol.sequence[b];
            int l = currSol.sequence[b + 1];

            // Compute the 2-opt cost variation
            double neighCostVariation = -M[h][i] - M[j][l] + M[h][j] + M[i][l];

            // If tabu and not aspiration criteria
            if ((currIter - tabuList[i] <= tabuLength) && (currIter - tabuList[j] <= tabuLength) &&
                neighCostVariation >= aspiration - 0.01) {
                continue;
            }

            // Update the best cost variation
            if (neighCostVariation < bestCostVariation) {
                bestCostVariation = neighCostVariation;
                move.from = a;
                move.to = b;
            }
        }
    }

    if (bestCostVariation >= infinite) cout << "\nNo best neighbor has been found!" << endl;

    return bestCostVariation;
}

/** Evaluate a solution
  * @param sol: solution to be evaluated
  * @return the value of the solution
  */
double TabuSearchSolver::evaluate(const TSPsolution &sol) const {
    double total = 0.0;
    for (uint i = 0; i < sol.sequence.size() - 1; i++) total += M[sol.sequence[i]][sol.sequence[i + 1]];
    return total;
}

/** Initialize a tabu list with n random element
  * @param n: the number of elements in the tabu list
  */
void TabuSearchSolver::initTabuList(int n) { for (int i = 0; i < n; i++) tabuList.push_back(-tabuLength - 1); }

