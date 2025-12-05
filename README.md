Monte Carlo Simulation of the Inspection Game
This repository contains the Fortran 90 source code used for the finite population analysis in the paper:"Finite Population Dynamics Resolve the Central Paradox of the Inspection Game" European Physical Journal B (EPJ B), 2025 Authors: Bianca Ishikawa and José Fontanari 
Overview:
This code simulates the stochastic evolutionary dynamics of the Inspection Game in a finite population of Citizens and Inspectors. It implements an imitation dynamics process (pairwise comparison) where agents update their strategies based on the success of others.
The simulation tracks the system until it reaches an absorbing state (fixation), where all individuals in a population adopt the same strategy. By running thousands of independent realizations, the code estimates the fixation probability of different outcomes (e.g., extinction of crime vs. dominance of crime) as a function of the game parameters.
Strategies
Citizens:
Crime (C): Commit a crime (gain g, penalty p if caught).
No Crime (NC): Law-abiding behavior (payoff 0).
Inspectors:
Inspect (I): Pay cost k to inspect (reward r if criminal caught).
No Inspect (NI): Do not inspect (payoff 0).
Files
inspection_game.f90: The main Fortran 90 source code.
README.md: This documentation file.
Requirements
To compile and run this code, you need a Fortran compiler. The code is standard Fortran 90 and should be compatible with most modern compilers, such as:
GNU Fortran (gfortran): Open source and widely available (part of GCC).
Intel Fortran Compiler (ifort): Commercial compiler (often free for academic/research use).
Compilation and Execution
Using gfortran (Linux/macOS/Windows MinGW)
Compile:Open your terminal and navigate to the directory containing inspection_game.f90. Run the following command:gfortran -O3 inspection_game.f90 -o inspection_sim
(The -O3 flag enables high-level optimization, which is recommended for Monte Carlo simulations).Run:Execute the compiled binary:./inspection_sim
(On Windows, run inspection_sim.exe)
Parameters
The physical parameters of the Inspection Game are set within the PROGRAM Inspection_Game block in inspection_game.f90. Key parameters include:
N: Population size of Citizens (default: 1000).
M: Population size of Inspectors (default: 100).
Nsample: Number of independent Monte Carlo realizations per parameter set (default: 10,000).
Ntime: Maximum time steps allowed before forcing a stop (default: 1,000,000).
rr: Reward for catching a criminal ($r$).kk: Cost of inspection ($k$).
p: Penalty for being caught ($p$).
The code automatically sweeps the parameter g (Gain from crime) relative to p using the loop variable gg.
Output Data
The program appends results to a data file named simg_N1000_M100_r4_p100.dat (filename may vary based on hardcoded strings in the source).Columns in the output file:
rr: Reward parameter ($r$).
ix0y0: Count of fixations at No Crime, No Inspect.
ix0ym: Count of fixations at No Crime, All Inspect (Crime Extinction).
ixny0: Count of fixations at All Crime, No Inspect (Crime Dominance).
ixnym: Count of fixations at All Crime, All Inspect.
tx0y0: Average time to fixation for (No Crime, No Inspect).
tx0ym: Average time to fixation for (No Crime, All Inspect).
txny0: Average time to fixation for (All Crime, No Inspect).
txnym: Average time to fixation for (All Crime, All Inspect).
X0: Initial number of criminals.
Y0: Initial number of inspecting inspectors.
N: Total citizens.
M: Total inspectors.
g: Gain from crime parameter.
p: Penalty parameter.
Nsample: Total samples run.
License
This project is open-source. 
Please cite the associated EPJ B paper if you use this code for academic research.
