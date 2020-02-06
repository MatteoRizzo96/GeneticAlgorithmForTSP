#include <iostream>
#include <vector>
#include "classes/utility/include/ProblemInstanceGenerator.h"
#include "classes/utility/include/TSPsolverScript.h"
#include "classes/solvers/include/CPLEXsolver.h"
#include "classes/solvers/include/TabuSearchSolver.h"
#include "classes/solvers/include/GeneticAlgorithmSolver.h"
#include "classes/utility/include/Params.h"

int main(int argc, char const *argv[]) {
    try {
        TSPsolver *solver = nullptr;
        map<string, int> comparisonType = {{"trial", 1},
                                           {"size",  2},
                                           {"none",  3}};
        map<string, int> solverType = {{"cplex",   1},
                                       {"tabu",    2},
                                       {"genetic", 3}};
        vector<TSPsolver *> solvers = vector<TSPsolver *>({new TabuSearchSolver(),
                                                           new GeneticAlgorithmSolver(),
                                                           new CPLEXsolver()});
        switch (comparisonType[Params::getAllParams()["comparison"].GetString()]) {
            case 1:
                TSPsolverScript().compareSolversByTrial(solvers);
                break;
            case 2:
                TSPsolverScript().compareSolversBySize(solvers);
                break;
            case 3:
                cout << "\nTSP SOLVER. Currently using ";
                switch (solverType[Params::getAllParams()["solver"].GetString()]) {
                    case 1:
                        cout << "CPLEX" << endl;
                        solver = new CPLEXsolver();
                        break;
                    case 2:
                        cout << "Tabu search" << endl;
                        solver = new TabuSearchSolver();
                        break;
                    case 3:
                        cout << "Genetic Algorithm" << endl;
                        solver = new GeneticAlgorithmSolver();
                        break;
                    default:
                        throw std::logic_error("unsupported solver. Supported solvers are: 'cplex', 'tabu', 'genetic'");
                }
                cout << "\n----------------------------------\n" << endl;
                TSPsolverScript(solver).solveProblemInstances();
                break;
            default:
                throw std::logic_error("unsupported comparison. Supported comparisons are: 'trial', 'size', 'none'");
        }
    }
    catch (exception &e) { cout << ">> EXCEPTION: " << e.what() << endl; }
    return 0;
}