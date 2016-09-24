----------------------------------
      TORIC CODE SIMULATOR

	Tobias Hoelzer
	September 2013
	Institute for Quantum Information
	RWTH Aachen University
	Bachelor Thesis: Study of [[4,2,2]]-concatenated toric code for a scalable circuit-QED architecture


	based on the works of

      Martin Suchara
      UC Berkeley, 2012
----------------------------------



VERSION DESCRIPTION:

Toric code threshold simulator that assumes that syndrome measurements are perfect.

This version is only slightly different than the one from Martin Suchara.
Only Z and Y errors are removed.


1.1 COMPILATION

To compile type make clean; make.
To run the simulation type ./qecsim.


Some compilers may experience problems with the order of commands in the Makefile (especially the -ltr linker command)



1.2 SIMULATION PARAMETERS

The simulation parameters are specified in in/parameters.txt. The file has the following format:

XSize = 2 4 6
iterations = 100000
pMin = 0.001
pMax = 0.1
pStep = 0.001


The first line specifies the lattice sizes that the simulator evaluates. The second line specifies the number of repetitions of the simulation for each lattice size. For example, 100000 means that a random error is generated 100000 times and the simulator attempts to correct the error. The output is the percentage of the generated errors that were corrected successfully.

The following three lines specify how is the error in each iteration generated. For example, pMin = 0.001, pMax = 0.1, pStep = 0.001 means that a series of simulations is performed where the error is generated using the depolarizing error model starting with p=0.1% and following to p=10.0% by increasing the error in 0.1% increments.


1.3 SIMULATION OUTPUT

The simulation result is saved in the out directory in separate files for each lattice size. The directory also contains scripts that visualize the data. For visualization, copy the .txt files to the expt1 directory and type ./get_plot to produce a .ps and .pdf file with the visualizations.
