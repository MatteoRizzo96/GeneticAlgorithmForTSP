#include <iostream>
#include <random>
#include "../include/TSPsolution.h"
#include "../include/RandomNumbersGenerator.h"

/** Constructor: builds a standard solution as the sequence <0, 1, 2, 3 ... n-1, 0>
  * @param n: the size of the solution
  * */
TSPsolution::TSPsolution(const int n) {
    sequence.reserve(n + 1);
    for (int i = 0; i < n; ++i) sequence.push_back(i);
    sequence.push_back(0);
}

/** Copy constructor: builds a solution from another
  * @param tspSol: a TSP solution
  * */
TSPsolution::TSPsolution(const TSPsolution &tspSol) {
    sequence.reserve(tspSol.sequence.size());
    for (int i : tspSol.sequence) sequence.push_back(i);
}


/** Constructor: builds a solution from a sequence
  * @param tspSol: a TSP solution
  * */
TSPsolution::TSPsolution(const vector<int> &seq) {
    sequence.reserve(seq.size());
    for (int i : seq) sequence.push_back(i);
}

/** Initialize a solution as a random sequence by random swaps
  * @param n: the length of the solution to be initialized
  * @return the randomly initialized solution
  */
TSPsolution TSPsolution::getRandomSolution(int n) {
    vector<int> seq = vector<int>(n - 1);
    for (int i = 0; i < n - 1; i++) seq[i] = i + 1;
    shuffle(seq.begin(), seq.end(), mt19937(random_device()()));
    seq.insert(seq.begin(), 0);
    seq.push_back(0);
    return TSPsolution(seq);
}

/** Assignment overloading: copy a solution into another one
  * @param right TSP solution to get into
  * */
TSPsolution &TSPsolution::operator=(const TSPsolution &right) {
    if (this == &right) return *this;
    if (sequence.size() != right.sequence.size()) sequence = vector<int>(right.sequence.size(), 0);
    for (unsigned i = 0; i < sequence.size(); i++) sequence[i] = right.sequence[i];
    return *this;
}

void TSPsolution::print() const {
    cout << "### ";
    for (int i : sequence) cout << i << " ";
    cout << " ###" << endl;
}
