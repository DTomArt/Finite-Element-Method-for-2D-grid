# Finite-Element-Method-for-2D-grid

## Basic idea

Program was made for FEM university course - it's a console application written in C++.

The software is solving an undetermined problem heat exchange in the two-dimensional system. Simulation of boundary conditions (temperature t on the surface and perpendicular heat flux q according to the law of convection) will determine the approximate temperature field at the nodes of the mesh and print out the maximum and minimum temperatures in a non-stationery system, i.e. at any given moment of the time simulation.

## Compiling and running program

To run the program just run the file:  
`./a.out`  
You can change initial parameters for grid(mesh) in main.cpp. In this case you have to recompile.  
To do so, just compile all the files with:  
`g++ *.cpp`

## Program structure

Program consists of header files `.h` which contain structures made basing on FEM theory with thorough documentation and `.cpp` files, which contain implementation of functions and methods from header files.
