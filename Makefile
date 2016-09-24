all:
	g++ -Wall dijkstra.cc operator.cc lattice.cc simulation.cc qecsim.cc blossom5-v2.02.src/misc.o blossom5-v2.02.src/PMduals.o blossom5-v2.02.src/PMexpand.o blossom5-v2.02.src/PMinit.o blossom5-v2.02.src/PMinterface.o blossom5-v2.02.src/PMmain.o blossom5-v2.02.src/PMrepair.o blossom5-v2.02.src/PMshrink.o blossom5-v2.02.src/GEOM/GPMinit.o blossom5-v2.02.src/GEOM/GPMinterface.o blossom5-v2.02.src/GEOM/GPMkdtree.o blossom5-v2.02.src/GEOM/GPMmain.o blossom5-v2.02.src/MinCost/MinCost.o -O3 -o qecsim -lrt -fopenmp
clean:
	rm -f *.o qecsim *~
