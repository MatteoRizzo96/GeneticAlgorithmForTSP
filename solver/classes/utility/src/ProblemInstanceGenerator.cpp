#include <iostream>
#include <algorithm>
#include <cmath>
#include <map>
#include <fstream>
#include "../include/ProblemInstanceGenerator.h"
#include "../include/RandomNumbersGenerator.h"
#include "../include/Params.h"

vector<string> ProblemInstanceGenerator::problemName = {
        "from file",
        "random polygons",
        "random distribution",
        "random symmetric distribution"
};

map<string, ProblemInstanceGenerator::ProblemType> ProblemInstanceGenerator::problem = {
        {"from_file",                     ProblemInstanceGenerator::ProblemType::fromFile},
        {"random_polygons",               ProblemInstanceGenerator::ProblemType::randomPolygons},
        {"random_distribution",           ProblemInstanceGenerator::ProblemType::randomDistribution},
        {"random_symmetric_distribution", ProblemInstanceGenerator::ProblemType::randomSymmetricDistribution}
};

vector<double>
ProblemInstanceGenerator::generateProblem(ProblemType type, int nh, bool v) {
    // Initialize the main parameters
    numHoles = nh;
    verbose = v;
    rapidjson::Document d = Params::getAllParams();
    // Generate the problem of the given type
    switch (type) {
        case fromFile:
            return readFromFile("data/" + string(d["filename"].GetString()));
        case randomPolygons:
            return generateRandomPolygons();
        case randomDistribution:
            return generateRandomDistribution();
        case randomSymmetricDistribution:
            return generateRandomSymmetricDistribution();
        default:
            cout << "Problem type not supported! Returning default unitary costs";
            return vector<double>(numHoles, 1);
    }
}

// Utility

/** Pretty print the cost matrix to the standard output
  * @param costs: the cost matrix
  */
void ProblemInstanceGenerator::printCostMatrix(vector<double> costs) {
    cout << "\n Printing the cost matrix for " << numHoles << "...\n" << endl;
    for (int i = 0; i < numHoles; i++) {
        string line = "|\t" + to_string(i) + "\t|\t";
        for (int j = 0; j < numHoles; j++) line += to_string(costs[numHoles * i + j]) + "\t";
        cout << line << " |\n";
    }
    cout << endl;
}

// -- From file --

/** Read a cost matrix from file
  * @param filename: the name of the file that stores the cost matrix
  * @return the cost matrix that has been read
  */
vector<double> ProblemInstanceGenerator::readFromFile(const string &filename) {
    ifstream in(filename);
    in >> numHoles;
    vector<double> costs = vector<double>(pow(numHoles, 2));
    for (int i = 0; i < pow(numHoles, 2); i++) in >> costs[i];
    printCostMatrix(costs);
    return costs;
}

// -- Random distribution --

/** Generate a matrix of costs as a random distribution of integers
  * @return: a vector of doubles representing the costs of the problem
  * */
vector<double> ProblemInstanceGenerator::generateRandomDistribution() {
    vector<double> costs = vector<double>(pow(numHoles, 2));
    for (double &cost : costs) cost = RandomNumbersGenerator::getRandomInt(0, 100);
    return costs;
}

/** Generate a symmetric matrix of costs as a random distribution of integers
  * @return: a vector of doubles representing the costs of the problem
  * */
vector<double> ProblemInstanceGenerator::generateRandomSymmetricDistribution() {
    vector<vector<double>> costsMatrix = vector<vector<double>>(numHoles, vector<double>(numHoles));
    // Initialize the matrix above the diagonal
    for (int i = 0; i < numHoles; i++)
        for (int j = i; j < numHoles; j++)
            costsMatrix[i][j] = RandomNumbersGenerator::getRandomInt(1, 100);
    // Initialize the matrix below the diagonal
    for (int i = 0; i < numHoles; i++) for (int j = 0; j < i; j++) costsMatrix[i][j] = costsMatrix[j][i];
    // Convert the matrix into a vector
    vector<double> costs = vector<double>(pow(numHoles, 2));
    for (int i = 0; i < numHoles; i++) for (int j = 0; j < numHoles; j++) costs[i * numHoles + j] = costsMatrix[i][j];
    return costs;
}

// -- Random polygons --

void ProblemInstanceGenerator::writePolygonsOnFile(const vector<vector<pair<double, double>>> &polygons) {
    std::ofstream file("../data_analysis/polygons/" + to_string(numHoles) + "holesPolygons.csv");
    file << "num_vertices,vertices\n";
    for (const auto &polygon: polygons) {
        file << polygon.size() << ",";
        for (auto point : polygon) file << point.first << ":" << point.second << " ";
        file << "\n";
    }
}

/** Generate the costs for the generated set of holes distributed as random polygons
  * as the euclidean distances between any vertex v1 and v2
  * @param vertices: the whole set of vertices of all the generated polygons
  * @return: a vector of doubles representing the costs of the problem
  * */
vector<double> ProblemInstanceGenerator::generateRandomPolygonsCosts(const vector<pair<double, double>> &vertices) {
    vector<double> costs = vector<double>(pow(numHoles, 2));
    int idx = 0;
    for (auto v1 : vertices) {
        string line = "|\t";
        for (auto v2 : vertices) {
            costs[idx] = sqrt(pow(v1.first - v2.first, 2) + pow(v1.second - v2.second, 2));
            line += to_string(costs[idx]) + "\t";
            idx++;
        }
        cout << line + "\t|\n";
    }
    return costs;
}

/** Check if a given square intersects one of the elements of a list of squares
  * @param box: the square to be checked against the list
  * @param boxes: the list of squared
  * @return: a boolean flag that is true if an intersection exists, false otherwise
  * */
bool ProblemInstanceGenerator::isDisjoint(pair<pair<double, double>, pair<double, double>> box,
                                          const vector<pair<pair<double, double>, pair<double, double>>> &boxes) {
    auto tl1 = box.first;
    auto br1 = box.second;
    for (const auto &b: boxes) {
        auto tl2 = b.first;
        auto br2 = b.second;
        if (!(tl1.first > br2.first || tl2.first > br1.first || tl1.second < br2.second || tl2.second < br1.second))
            return false;
    }
    return true;
}

/** Return the top-left and bottom-right coordinates of the square boxing the polygon
  * @param r: the radius of the polygon
  * @param c: the coordinates of the centre of the polygon
  * @return: a couple containing the top-left and bottom-right coordinates of the bounding box
  * */
pair<pair<double, double>, pair<double, double>>
ProblemInstanceGenerator::getBoundingBox(double r, pair<double, double> c) {
    return make_pair(make_pair(c.first - r, c.second + r), make_pair(c.first + r, c.second - r));
}

/** Generate the set of vertices of a polygon.
  * @param n: the number of sides of the polygon
  * @param r: the radius of the polygon
  * @param theta: the angle of orientation of the polygon
  * @param c: the coordinates of the centre of the polygon
  * @return: a vector of coordinates representing the vertices of the polygon
  * */
vector<pair<double, double>>
ProblemInstanceGenerator::getPolygon(int n, double r, double theta, pair<double, double> c) {
    vector<pair<double, double>> vertices = vector<pair<double, double>>(n);
    for (int i = 0; i < n; i++)
        vertices[i] = make_pair(r * cos(2 * M_PI * i / n + theta) + c.first,
                                r * sin(2 * M_PI * i / n + theta) + c.second);
    return vertices;

}

/** Generate a set of random polygons.
  * @param numVerticesPolygons: the number of vertices of each polygon to be generated
  * @param polygons: the vector of generated polygons
  * @param boxes: the boxes (i.e. squares) that contain the generated polygons
  * */
vector<vector<pair<double, double>>>
ProblemInstanceGenerator::generatePolygons(const vector<int> &numVerticesPolygons) {
    // Initialize a vector of polygons, and the boxes that contain them
    vector<vector<pair<double, double>>> polygons;
    vector<pair<pair<double, double>, pair<double, double>>> boxes;
    // Iterate over the number of vertices of each polygon
    for (int numVertices : numVerticesPolygons) {
        bool disjoint = false;
        while (!disjoint) {
            // Randomly generate the data for the polygon
            double angle = RandomNumbersGenerator::getRandomDouble(0, int(2 * M_PI));
            double radius = RandomNumbersGenerator::getRandomDouble(1, int(boardSize.first / 4));
            pair<double, double> center = make_pair(
                    RandomNumbersGenerator::getRandomDouble(ceil(radius), floor(boardSize.first - radius)),
                    RandomNumbersGenerator::getRandomDouble(ceil(radius), floor(boardSize.second - radius)));

            auto box = getBoundingBox(radius, center);
            if (isDisjoint(box, boxes)) {
                auto polygon = getPolygon(numVertices, radius, angle, center);
                polygons.push_back(polygon);
                boxes.push_back(box);
                disjoint = true;
            }
        }
    }
    return polygons;
}

/** Split the number of holes in chunks of random length corresponding
  * to the number of vertices of each polygon
  * @param k: the maximum overall number of vertices of the polygons
  */
vector<int> ProblemInstanceGenerator::generatePolygonsSplit(int k) {
    vector<int> numVerticesPolygons;
    while (k > 0) {
        int n = RandomNumbersGenerator::getRandomInt(3, k);
        if (k - n >= 3 || k - n == 0) {
            numVerticesPolygons.push_back(n);
            k -= n;
        }
    }
    return numVerticesPolygons;
}

/** Generate a random number of regular polygons with random number of sided and center.
  * The number of polygons generated is based on the desired number of holes. However,
  * if the number of holes is lower than 3, no polygon can be generated and just numHoles
  * random points are generated.
  * @return: a vector of doubles representing the costs of the problem
  * */
vector<double> ProblemInstanceGenerator::generateRandomPolygons() {
    // If the number of holes is lower than 3, no polygon can be generated
    if (numHoles < 3) return generateRandomDistribution();
    // Generate a random polygon for each chunk
    auto polygons = generatePolygons(generatePolygonsSplit(numHoles));
    writePolygonsOnFile(polygons);
    // Flatten all the vertices of the polygons into a single vector
    vector<pair<double, double>> vertices;
    for (auto p : polygons) vertices.insert(end(vertices), begin(p), end(p));
    // Generate the costs as the euclidean distance of the vertices
    cout << "Generated " << to_string(polygons.size()) << " random polygon(s)\n" << endl;
    return generateRandomPolygonsCosts(vertices);
}

