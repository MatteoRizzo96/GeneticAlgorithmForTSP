#include <cstdio>
#include <iostream>
#include <utility>
#include <vector>
#include "../../utility/include/cpxmacro.h"
#include "../include/CPLEXsolver.h"

using namespace std;

// Error status and message buffer
int status;
char errmsg[BUF_SIZE];

CPLEXsolver::CPLEXsolver(const vector<double> &costs, int size, bool logSolutions) {
    C = costs;
    N = size;
    log = logSolutions;
    try {
        DECL_ENV(env)
        environment = env;
        DECL_PROB(env, lp)
        problem = lp;
    }
    catch (exception &e) { cout << ">> EXCEPTION: " << e.what() << endl; }
}

/** Set the parameters for the current instance of the TSP
  * @param costs: the matrix of the costs
  * @param size: the number of nodes of the problem
  */
void CPLEXsolver::setParameters(const vector<double> &costs, int size) {
    C = costs;
    N = size;
    try {
        CPXfreeprob(environment, &problem);
        DECL_ENV(env)
        environment = env;
        DECL_PROB(env, lp)
        problem = lp;
    }
    catch (exception &e) { cout << ">> EXCEPTION: " << e.what() << endl; }
}


/** Set up the constraints for the TSP
  * @return a pair of maps containing the names and the values of the x and y variables
  */
void CPLEXsolver::setupConstraints(const pair<map<pair<int, int>, int>, map<pair<int, int>, int>> &varsIdx) {
    auto xVars = varsIdx.first;
    auto yVars = varsIdx.second;

    // (10) [sum_i (x_ik) - sum_j (x_kj) = 1], with (i, k), (k, j) in A and j != 0, for each k in N \ {0}
    int numXVars = 2 * N - 3;
    for (int k = 1; k < N; k++) {
        int index = 0;
        // Declaration of the variables and their coefficients (for the left-hand side)
        vector<int> idx(numXVars, 0);
        vector<double> coefficient(numXVars, 0);
        // sum_i (x_ik)
        for (int i = 0; i < N; i++) {
            if (k == i) continue;
            idx[index] = xVars.find(make_pair(i, k))->second;
            coefficient[index] = 1;
            index++;
        }
        // - sum_j (x_kj)
        for (int j = 1; j < N; j++) {
            if (k == j) continue;
            idx[index] = xVars.find(make_pair(k, j))->second;
            coefficient[index] = -1;
            index++;
        }
        double rhs = 1;
        char sense = 'E';
        int matbeg = 0;
        CHECKED_CPX_CALL(CPXaddrows, environment, problem, 0, 1, idx.size(), &rhs, &sense, &matbeg, &idx[0],
                         &coefficient[0],
                         nullptr, nullptr);
    }

    // (11) [sum_j (y_ij) = 1], for each i in N
    for (int i = 0; i < N; i++) {
        // Declaration of the variables and their coefficients (for the left-hand side)
        int numVars = N - 1;
        vector<int> idx(numVars);
        vector<double> coefficient(numVars, 1);
        int index = 0;
        for (int j = 0; j < N; j++) {
            if (i == j) continue;
            idx[index] = yVars.find(make_pair(i, j))->second;
            index++;
        }
        double rhs = 1;
        char sense = 'E';
        int matbeg = 0;
        CHECKED_CPX_CALL(CPXaddrows, environment, problem, 0, 1, idx.size(), &rhs, &sense, &matbeg, &idx[0],
                         &coefficient[0],
                         nullptr, nullptr);
    }

    // (12) [sum_i (y_ij) = 1], for each j in N
    for (int j = 0; j < N; j++) {
        // Declaration of the variables and their coefficients (for the left-hand side)
        int numVars = N - 1;
        vector<int> idx(numVars);
        vector<double> coefficient(numVars, 1);
        int index = 0;
        for (int i = 0; i < N; i++) {
            if (i == j) continue;
            idx[index] = yVars.find(make_pair(i, j))->second;
            index++;
        }
        double rhs = 1;
        char sense = 'E';
        int matbeg = 0;
        CHECKED_CPX_CALL(CPXaddrows, environment, problem, 0, 1, idx.size(), &rhs, &sense, &matbeg, &idx[0],
                         &coefficient[0],
                         nullptr, nullptr);
    }

    // (13) [x_ij <= (|N| - 1) * y_ij] for each (i, j) in A, j != 0
    for (int i = 0; i < N; i++) {
        for (int j = 1; j < N; j++) {
            if (i == j) continue;
            vector<int> idx(2);
            idx[0] = xVars.find(make_pair(i, j))->second;
            idx[1] = yVars.find(make_pair(i, j))->second;
            vector<double> coefficient(2);
            coefficient[0] = 1.0;
            coefficient[1] = -(N - 1);
            double rhs = 0;
            char sense = 'L';
            int matbeg = 0;
            CHECKED_CPX_CALL(CPXaddrows, environment, problem, 0, 1, idx.size(), &rhs, &sense, &matbeg, &idx[0],
                             &coefficient[0],
                             nullptr, nullptr);
        }
    }
}

/** Set up the variables for the TSP
  * @return a pair of maps containing the names and the values of the x and y variables
  */
pair<map<pair<int, int>, int>, map<pair<int, int>, int>> CPLEXsolver::setupVariables() {
    const int NAME_SIZE = 512;
    char name[NAME_SIZE];

    // Set up a var to index the x and y variables of the problem
    int currentVarPos = 0;

    // (14) x_ij = amount of the flow shipped from i to j, for each (i, j) in A, with j != 0
    map<pair<int, int>, int> xVars;
    for (int i = 0; i < N; i++) {
        for (int j = 1; j < N; j++) {
            if (i == j) continue;
            char xtype = 'C';
            double lb = 0.0;
            double ub = CPX_INFBOUND;
            double obj = 0.0;
            snprintf(name, NAME_SIZE, "x_%d_%d", i, j);
            char *xname = (char *) (&name[0]);
            CHECKED_CPX_CALL(CPXnewcols, environment, problem, 1, &obj, &lb, &ub, &xtype, &xname);
            xVars.insert({make_pair(i, j), currentVarPos});
            currentVarPos++;
        }
    }

    // (15) y_ij = 1 if arc (i, j) ships some flow, 0 otherwise, for each (i, j) in A
    map<pair<int, int>, int> yVars;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j) continue;
            char ytype = 'B';
            double lb = 0.0;
            double ub = 1.0;
            double obj = C[i * N + j];
            snprintf(name, NAME_SIZE, "y_%d_%d", i, j);
            char *yname = (char *) (&name[0]);
            CHECKED_CPX_CALL(CPXnewcols, environment, problem, 1, &obj, &lb, &ub, &ytype, &yname);
            yVars.insert({make_pair(i, j), currentVarPos});
            currentVarPos++;
        }
    }

    return make_pair(xVars, yVars);
}


void CPLEXsolver::solve() {
    try {
        auto varsIdx = setupVariables();
        setupConstraints(varsIdx);
        CHECKED_CPX_CALL(CPXmipopt, environment, problem);
        if (log) logSolution();
    }
    catch (exception &e) { cout << ">> EXCEPTION: " << e.what() << endl; }
}

double CPLEXsolver::getObjFunValue() {
    double objval;
    try { CHECKED_CPX_CALL(CPXgetobjval, environment, problem, &objval); }
    catch (exception &e) {
        cout << ">> EXCEPTION: " << e.what() << endl;
        objval = -1.0;
    }
    return objval;
}

TSPsolution CPLEXsolver::getSolution() { return TSPsolution(0); }

void CPLEXsolver::logSolution() {
    // Get the number of variables in the model
    int n = CPXgetnumcols(environment, problem);
    // Get the value of the objective function
    double objval = getObjFunValue();
    cout << "      Objval: " << objval << endl;
    // Get the values of the variables
    vector<double> varVals;
    varVals.resize(n);
    CHECKED_CPX_CALL(CPXgetx, environment, problem, &varVals[0], 0, n - 1);
    for (int i = 0; i < n; ++i) cout << "      var in position " << i << " : " << varVals[i] << endl;
    CHECKED_CPX_CALL(CPXsolwrite, environment, problem, "TSP.sol");
    CHECKED_CPX_CALL(CPXwriteprob, environment, problem, "TSP.lp", nullptr);
}

CPLEXsolver::~CPLEXsolver() {
    CPXfreeprob(environment, &problem);
    CPXcloseCPLEX(&environment);
}

