#include "../include/TSPsolverScript.h"
#include "../include/Params.h"
#include "../../solvers/include/CPLEXsolver.h"

#include <iostream>
#include <utility>
#include <vector>
#include <list>
#include <chrono>
#include <fstream>

using namespace std;
using namespace rapidjson;


TSPsolverScript::TSPsolverScript(TSPsolver *s) {
    Document d = Params::getAllParams();

    // Initialize the main general parameters
    solver = s;
    solverName = d["solver"].GetString();
    problemName = d["problem"].GetString();
    fileName = string(d["filename"].GetString());
    problem = ProblemInstanceGenerator::problem[d["problem"].GetString()];
    problemGenerator = ProblemInstanceGenerator();
    numTrials = d["num_trials"].GetInt();

    cout << "Solving instance of type \"" << ProblemInstanceGenerator::problemName[problem] << "\"\n" << endl;
    cout << "Target size(s): ";

    // Initialize the target sizes
    if (problemName == "from_file") {
        ifstream data("data/" + fileName);
        size = vector<int>(1);
        data >> size[0];
        cout << to_string(size[0]) << endl;
    } else {
        double currentSize = d["size"]["start_size"].GetInt();
        double endSize = d["size"]["end_size"].GetInt();
        double step = d["size"]["step"].GetInt();
        do {
            size.push_back(currentSize);
            currentSize += step;
        } while (currentSize <= endSize);
        cout << "from " << to_string(size[0]) << " to " << to_string(size.back()) << " with step " << step << endl;
    }

    // Initialize the target times
    double currentTime = d["time"]["start_time"].GetDouble();
    double endTime = d["time"]["end_time"].GetDouble();
    double step = d["time"]["step"].GetDouble();
    do {
        time.push_back(currentTime);
        currentTime *= step;
    } while (currentTime < endTime);

    cout << "\nTarget time(s):\n" << endl;
    for (auto t : time) cout << "      * " << t << "s" << endl;
    cout << "\n----------------------------------\n" << endl;
}

/** Compute a solution to the TSP according to the given parameters and
  * records the time of execution
  * @param costs: the cost matrix of the instance of the problem to be solved
  * @param n: the size of the problem
  * @param s: optional solver. If not specified, the class solver is used
  * @return a pair (objective function value, time of execution)
  */
pair<double, double> TSPsolverScript::computeSolution(const vector<double> &costs, int n, TSPsolver *s) {
    if (s) solver = s;
    solver->setParameters(costs, n);
    auto start = chrono::high_resolution_clock::now();
    solver->solve();
    auto finish = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = finish - start;
    return make_pair(solver->getObjFunValue(), elapsed.count());
}

/** Search the greatest size of a problem instance solvable within the target time
  * @param t: a target time
  * @return a pair (objective function value, time of execution)
  */
pair<double, double> TSPsolverScript::searchSolution(double t) {
    ofstream file(pathToSolutions +
                  solverName + "__" +
                  to_string(size[0]) + "_" + to_string(size.back()) + "__" +
                  to_string(unsigned(t * 1000)) + "__" +
                  to_string(numTrials) + "_trials__" +
                  problemName + ".csv");
    file << "size,solutions" << endl;

    pair<double, double> prevSol = make_pair(-1.0, -1.0);
    bool stop = false;

    for (int i = 0; i < int(size.size()) && !stop; i++) {
        cout << "\tSolving the problem instance of size " << size[i] << "...\n" << endl;
        file << to_string(size[i]) + ",";
        for (int j = 0; j < numTrials; j++) {
            auto sol = computeSolution(problemGenerator.generateProblem(problem, size[i]), size[i]);
            file << to_string(sol.first) + ":" + to_string(sol.second) << " ";
            if (sol.second > t) {
                if (prevSol.first == -1)
                    cout << "\n\tNo solution found! The elapsed time for base size is: " << sol.second << "s" << endl;
                else
                    cout << "\n\tSOLUTION FOUND at size " << size[i] << "! Elapsed time: " << sol.second << "s" << endl;
                stop = true;
            } else prevSol = sol;
        }
        file << endl;
    }
    if (!stop) cout << "\n\tTested all possible sizes! " << endl;
    return prevSol;
}

/** Solve instances of the problem of increasing size until meeting the target times */
void TSPsolverScript::solveProblemInstances() {
    list<pair<double, double>> solutions;
    for (double t : time) {
        cout << "\n  --> Looking for the size of the greatest instance solvable within target time '" << t
             << "'...\n" << endl;
        solutions.emplace_back(searchSolution(t));
    }
}

/** Compare the given solvers solving increasing instances of the current problem
  * @param solvers: the vector of solvers to be compared
  */
void TSPsolverScript::compareSolversBySize(const vector<TSPsolver *> &solvers) {
    cout << "\nTSP SOLVER. Performing size comparison among solvers\n" << endl;
    cout << "----------------------------------\n" << endl;
    string pn = problemName == "from_file" ? fileName.substr(0, fileName.find('.')) : problemName;
    string filename = to_string(numTrials) + "_trials__" +
                      to_string(size[0]) + "_" + to_string(size.back()) + "_size__" +
                      pn + "__comparison.csv";
    std::ofstream file(pathToSize + filename);
    file << "size,tabu,genetic,cplex\n";
    for (int n : size) {
        cout << "\nSIZE (" << n << ")\n" << endl;
        file << to_string(n);
        vector<double> p = problemGenerator.generateProblem(problem, n);
        for (TSPsolver *s: solvers) {
            file << ",";
            for (int i = 0; i < numTrials; i++) {
                pair<double, double> sol = computeSolution(p, n, s);
                file << to_string(sol.first) << ":" << to_string(sol.second) << " ";
            }
        }
        file << endl;
    }
}

/** Compare the given solvers solving numTrials identical instances of the current problem
  * @param solvers: the vector of solvers to be compared. CPLEX is not present since it is
  * execute only once
  */
void TSPsolverScript::compareSolversByTrial(const vector<TSPsolver *> &solvers) {
    cout << "\nTSP SOLVER. Performing trials comparison among solvers\n" << endl;
    cout << "----------------------------------\n" << endl;
    string pn = problemName == "from_file" ? fileName.substr(0, fileName.find('.')) : problemName;
    std::ofstream file(pathToTrial + to_string(numTrials) + "_trials__" + pn + "__comparison.csv");
    file << "tabu_obj,tabu_time,genetic_obj,genetic_time,cplex_obj,cplex_time\n";
    int n = size[0];
    for (int i = 0; i < numTrials; i++) {
        cout << "\nTRIAL (" << i + 1 << ")\n" << endl;
        for (TSPsolver *s: solvers) {
            auto sol = computeSolution(problemGenerator.generateProblem(problem, n), n, s);
            file << to_string(sol.first) << "," << to_string(sol.second) << ",";
        }
        file << endl;
    }
}
