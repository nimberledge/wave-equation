# Solving the Wave Equation using Concurrency and Domain Decomposition

This is a coursework project aimed at solving the Wave Equation in 2D. The Concurrent framework is to use MPI with domain decomposition in order to efficiently parallelize.

## Usage
```./make_run.sh <N>``` is a bash-script that compiles the code and runs it  on ```N``` MPI processes for a problem defined by ```src/conf_file.txt```. After the code finishes executing, a Python script is run to collate the output files of the program into a series of PNGs. This Python script also provides a means to convert the PNGs into an animated GIF, though due to compression issues, this is not the default behaviour. 

 In general, the code can be compiled for MPI by executing ```make```, and then run with ```mpiexec -n N bin/parwave <filename>``` (on ```N``` MPI processes). The filename is optional and allows one to specify the configuration file.

## Domain Decomposition / Concurrency
The details of the numerical model and the concurrent model are found in the ```Domain_Decomposition_for_the_Wave_Equation.pdf``` file, along with run-time analysis on an HPC system.

## Configuration settings
By editing a configuration file, one can set the bounds of the domain, the resolution of the final images (with dx and dy). One can also modify the output time step, wave speed constant, and the bound of the time domain. They also specify the type of boundary conditions, initial conditions, and whether or not to actually write outputs (no writes for benchmarks).

## Boundary conditions
The code implements 3 main kinds of boundary conditions - 
1. Dirichlet Boundaries where the value of displacement is 0 at domain edges
2. Neumann Boundaries where the gradient of displacement is 0 at domain edges
3. Periodic Boundaries where the domain edges are treates as continuous, i.e. the left edge wraps around to the right side of the domain.

In addition to this, users can specify a bool function to define internal boundaries on the domain, as seen in ```src/parwave.cpp```. Users can add more than 1 internal boundary by simply defining each obstacle independently and then combining them with a logical or operator. This is also demonstrated in ```src/parwave.cpp```, as a compound internal boundary. 

## Initial conditions
The initial conditions are pre-defined functions found in ```src/initialconditions.h``` that allow for a few initial conditions of interest. These are defined as mathematical functions $f: \mathbb{R}^2 \rightarrow \mathbb{R}$, and applied to the domain. To add an initial condition, one could add the function to ```src/initialconditions.h``` and update the array found at the end of the file. This can then be selected via its index in the configuration file.

## Results
Using a point disturbance as an initial condition, here are some results of the program. All inputs to the program can be found in the ```images/``` directory as text files.

### Dirichlet BCs, no internal BCs
![Link to image](/images/dirichlet_pd.gif "Dirichlet PD")

### Neumann BCs, no internal BCs
![Link to image](/images/neumann_pd.gif "Neumann PD")

### Periodic BCs, no internal BCs
![Link to image](/images/periodic_pd.gif "Periodic PD")

### Neumann BCs, with 2 internal obstacles
![Link to image](/images/internal_boundary.gif "Internal Boundaries PD")



