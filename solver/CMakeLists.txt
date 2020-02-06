cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_MODULE_PATH "cmake")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -DIL_STD -m64 -Wall -Werror")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Og -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG -O3 -flto")
set(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem ")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../solver/build")
set(LINKER_OPTIONS -flto -Wl,--no-as-needed)

if (CMAKE_SYSTEM_NAME STREQUAL Darwin)
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -save-temps=obj")
endif ()

# ======================
# Project specifications
# ======================

set(My_Project_Title "Labs-TSP")
set(My_Project_Description "TSP solver implementation using CPLEX APIs and Genetic Algorithm.")
set(My_Project_Author "Matteo Rizzo")
set(My_Project_Revised_Date "30-10-2019")

set(My_Project_Copyright_Date "2019-2020")
set(My_Project_Copyright "Copyright © ${My_Project_Copyright_Date} ${My_Project_Author}.")

string(TOUPPER ${My_Project_Title} My_Project_Title_Caps)
string(REGEX REPLACE " " "_" My_Project_Title_Caps ${My_Project_Title_Caps})

project(${My_Project_Title})

# =================
# Configure Project
# =================

find_package(Cplex)
include_directories(SYSTEM ${CPLEX_INCLUDE_DIRS})

set(SOURCE_FILES
		main.cpp
		classes/solvers/src/CPLEXsolver.cpp
		classes/solvers/include/CPLEXsolver.h
		classes/utility/src/ProblemInstanceGenerator.cpp
		classes/utility/include/ProblemInstanceGenerator.h
		classes/utility/src/TSPsolution.cpp
		classes/utility/include/TSPsolution.h
		classes/solvers/include/TSPsolver.h
		classes/utility/src/TSPsolverScript.cpp
		classes/utility/include/TSPsolverScript.h
		classes/solvers/src/TabuSearchSolver.cpp
		classes/solvers/include/TabuSearchSolver.h
		classes/utility/src/RandomNumbersGenerator.cpp
		classes/utility/include/RandomNumbersGenerator.h
		classes/solvers/src/GeneticAlgorithmSolver.cpp
		classes/solvers/include/GeneticAlgorithmSolver.h
		classes/utility/src/Params.cpp
		classes/utility/include/Params.h)

add_executable(${My_Project_Title} ${SOURCE_FILES})
target_link_libraries(${My_Project_Title} PRIVATE ${CPLEX_LIBRARIES})
target_link_libraries(${My_Project_Title} PRIVATE dl)


