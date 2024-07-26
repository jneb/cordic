# cordic
Demonstration of the working of the CORDIC algorithms.

## Introduction
CORDIC is a set of hardware optimized implmentations of many calculator functions.
This implementation emphasizes on decimal calculations.

### The CORDIC loop
All functions use the same loop that consists of only the following operations:
* addition and subtraction
* shifting (per digit) of values
* table lookup (only four tables of 6 entries for the whole program!)
* check the sign of a value
This is why CORDIC is still used today in hardware optimized evaluation on FPGAs.

## Implemented functions
This program implements:
* multiplication and division
* sin, cos in degrees and radians
* atan in degrees and radians
* sinh and cosh
* atanh
* tan and tanh (using two CORDIC operations)
* hypot (sqrt(x**2+y**2) using two CORDIC operations
* exp, log, sqrt

## Goal of the file
The goal is to showcase the working of CORDIC.
That's why there is a very limited class Fixed that defines fixed point numbers with only a few operations:
- addition
- subtraction, negation
- shifting to the right (with proper rounding)
- comparing with a given value giving the difference as a number of digits (for accuracy measurements)
- comparing "greater than 0" and nothing else
- (left) multiplication with -1, 0 and 1

All the operations above are calculated using Fixed numbers.
Setup is using regular floats, of course.
