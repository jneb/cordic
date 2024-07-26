# cordic
Demonstration of the working of the CORDIC algorithms.

## Introduction
CORDIC is a set of hardware optimized implmentations of many calculator functions.
This implementation emphasizes on decimal calculations.

### The CORDIC loop
All functions use the same loop that consists of only the following operations:
* addition and subtraction
* shifting (per digit) of values
* table lookup
* check the sign of a value
This is why CORDIC is still used today in hardware optimized evaluation on FPGAs.

## Implemented functions
This program implements:
* multiplication and division
* sin, cos in degrees and radians
* sinh and cosh
* tan and tanh (using two CORDIC operations)
* 
