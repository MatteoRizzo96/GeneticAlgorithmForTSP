#ifndef LAB01_TSP_PROBLEMINSTANCEGENERATOR_H
#define LAB01_TSP_PROBLEMINSTANCEGENERATOR_H

#include <vector>
#include <string>
#include <map>

using namespace std;

class ProblemInstanceGenerator {
public:

    // Desired final number of holes in the board
    int numHoles = 0;

    // Dimensions of the board (width x height)
    pair<double, double> boardSize = make_pair(200.0, 200.0);

    // Possible types of problem instances
    enum ProblemType {
        fromFile, uniformGrid, randomGrid, randomPolygons, randomDistribution, randomSymmetricDistribution
    };

    // Mapping from enum to problem name as string
    static vector<string> problemName;

    // Mapping from string to enum
    static map<string, ProblemType> problem;

    // Flag to trigger the printing
    bool verbose;

    vector<double> generateProblem(ProblemType type, int nh, bool v = false);

private:

    // Utility

    void printCostMatrix(vector<double> costs);

    // From file

    vector<double> readFromFile(const string &filename);

    // Random distributions

    vector<double> generateRandomSymmetricDistribution();

    vector<double> generateRandomDistribution();

    // Random polygons

    void writePolygonsOnFile(const vector<vector<pair<double, double>>> &polygons);

    vector<double> generateRandomPolygonsCosts(const vector<pair<double, double>> &vertices);

    static bool isDisjoint(pair<pair<double, double>, pair<double, double>> box,
                           const vector<pair<pair<double, double>, pair<double, double>>> &boxes);

    static pair<pair<double, double>, pair<double, double>> getBoundingBox(double r, pair<double, double> c);

    static vector<pair<double, double>> getPolygon(int n = 4,
                                                   double r = 1,
                                                   double theta = 0,
                                                   pair<double, double> c = make_pair(0, 0));

    vector<vector<pair<double, double>>> generatePolygons(const vector<int> &numVerticesPolygons);

    static vector<int> generatePolygonsSplit(int k);

    vector<double> generateRandomPolygons();
};


#endif //LAB01_TSP_PROBLEMINSTANCEGENERATOR_H
