
# RKO: Random-key Optimization

This is an implementation of the Random-key Optimization (RKO) to solve combinatorial optimization problems.

The C++ code of this algorithm has been designed to be easy of reuse. Users can only implement specific functions (read and decoder). 


## References

When using this algorithm in academic studies, please refer to the following works:

Title: A random-key GRASP for combinatorial optimization
Authors: Antonio A. Chaves, Mauricio G.C. Resende, Ricardo M.A. Silva
ArXiv: https://arxiv.org/submit/5626817

A random-key GRASP for combinatorial optimization. Antonio A. Chaves, Mauricio G.C. Resende, Ricardo M.A. Silva. In 15th Metaheuristics International Conference, June 4th to June 7th, 2024, Lorient, France.


## Scope

This code has been designed to solve the traveling salesman problem (TSP). Users need to configure only Problem.cpp file.


## Running the algorithm (linux or mac)

* Enter the Program directory: `cd Program`
* Run the make command: `make`
* Run the BRKGA_QL: `./runTest ../Instances/TSP-tests-debug.csv #MH 0`, where #MH is the metaheuristic chosen


## Code structure

The code structure is documented in [1] and organized in the following manner:

* **SPECIFIC_CODE:**
    * **Problem.cpp**: Contains data structure of the problem, the input function, and the decoder.

* **GENERAL_CODE:**
    * **MHs.cpp**: Contains all of the metaheuristics' mechanisms.
    * **Main.cpp**: Contains the main function to start the algorithm.
    * **Data.h**: Represents the data structures of MH and QL.
    * **Output.h**: Stores the outputs functions, including the best solution found and statistical analysis of the MH.

## File testScenario.csv is the input data problem and each line consists in:

- Instance Name
- Run mode (0 = debug, prints in the screen; 1 = run, prints in files)
- Number of implemented decoders
- Maximum rumtime (in seconds)
- Maximum number of runs
- Number of threads used in OpenMP
- Optimal or lower bound solution (if is known), 0 otherwise

Users need to create a folder named "Instances/ProblemName", where the instances must be; Users also need to create a folder named "Results" where the results files are writed.
