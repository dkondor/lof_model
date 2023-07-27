# ODE model

This folder contains a simple analytical (ODE) model of landscape of fear interactions and related analysis.

## Compilation

The main code calculating the solutions is written in C++ and relies on the [GNU Scientific Library](https://www.gnu.org/software/gsl/). On a typical Linux system, it can be compiled as follows:
```
g++ -o pm3 pop_model3.cpp -lgsl -lgslcblas -O3 -march=native -lm
```

## Running the simulation

Running the resulting program will calculate one solution. Model parameters, initial conditions and noise levels can be given as command line parameters; results are written on standard output. Example:
```
./pm3 -t l -a 0.2 -b 0.03 -d 0.01 -c 1 -r 0.0135 -K 10000 -R 0.1 -x 0.1 -y 100 -z 10 -T 3000 -Nx 0.05 -Ny 0.05 -Nz 0.05 > lof1.out
```
Command line parameters:
 - `-a`, `-b`, `-c`, `-d`, `-r`, `-K`, `-R`: set the model parameters (notation according to Eqs. (1)-(3) in our paper)
 - `-x`, `-y`, `-z`: set initial values of variables, where `x` is the level of conflict, `y` is the population in the general landscape and `z` is population in the regufe
 - `-T`: total simulation tim
 - `-Nx`, `-Ny`, `-Nz`: (relative) level of noise to apply to each variable (in each year)
 - `-s`: seed to use for the PRNG (used when generating noise)

## Analysis

The script [lof_plots1.gp](lof_plots1.gp) contains commands for the [Gnuplot](http://gnuplot.info) plotting tool that run the simulation with the main parameters and create the main figures in our paper. The additional R script [lof_acf_simple.r](lof_acf_simple.r) calculates ACFs in the case of noisy simulation runs.

