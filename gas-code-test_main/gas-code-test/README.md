`gas.cpp` contains a simple 1D Lagrangian gas dynamics code.

Please complete the following tasks:

1. Parallelise the loops - we recommend using OpenMP (https://www.openmp.org/) but you are free to choose any alternative method.
2. Explore performance with varying thread counts and varying problem
   sizes (within the bounds of compute available to you, if that is only e.g.
   two or four CPU cores that is not a problem).

There is a simple Makefile that assumes you are using gcc but any c++ compiler that supports
C++11 should work (https://en.cppreference.com/w/cpp/compiler_support/11) just as well.

The provided `plot()` function can be called to print CSV output to standard
out. These files can be compared directly (e.g. with `diff`), or plotted with
the included `plot.py` Python script. Two `.png` files are included showing the
density field for 100 cell and 100,000 cell runs.  Any changes should not
perceptibly alter the output for a fixed number of cells (very minor deviations
of the order of floating-point round-off may occur depending on compiler, flags
and system).

This exercise is expected to take ~2 hours.
