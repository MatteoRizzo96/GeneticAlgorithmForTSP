# A genetic algorithm for the TSP

This is the project for the [course of Methods and Models for Combinatorial Optimizazion](https://www.math.unipd.it/~luigi/courses/metmodoc/metmodoc.html) held at University of Padua by professors Luigi de Giovanni and Marco di Summa. You can find extensive documentation for the project on the `docs/Report.pdf` file.

The software features a C++ implementation of a Genetic Algorithm for solving the TSP. This has been extensively compared to an exact method (via IBM CPLEX) and to another heuristic method (Tabu Search).

### Installation

The project is built using Cmake. In order to run the software, make sure that the current working directory is set to `solver` and then run the following commands inside `solver`:

```
cmake .
make
cd build
./Labs-TSP
```

### Configuration

The user can tune the parameters of the execution via the `config/params.json` file.
