#ifndef MYCONST_H
#define MYCONST_H

//Debug Martin
static const bool debug1 = false;
//First Debugger to see it running
static const bool debug2 = false;
//produce lattice, correction algorith lines, qubits
static const bool debug3 = false;
//correction algoritm working?
//be careful, dumps core if not!
static const bool debug4 = false;
//Random Number Debug
static const bool debug5 = false;
//Cell Initialisizer Debug
static const bool debug6 = false;
//Success() Debug wrong stuff
static const bool debug7 = false;
//Success() Debug Z string
static const bool debug8 = false;
//Success() Debug correct stuff
static const bool debug9 = false;
//Runtime?
static const bool RUNTIME = true;


//2DDEBUG makes the simulation 2Dimensional with one in time and one in space
//Simply a 2D cross section through spacetime
//Doesnt that sound like fun?
static const bool DEBUG2D = false;

//debug10 debugs dijkstra algorithm
static const bool debug10 = false;



//Want to see a progressbar?
//not shure, if recommended
static const bool PROGRESSBAR = false;


//Pmax is the maximum p determine_p_per_step is looking for
static const double PMAX = 0.5;
//PSTEP and DELTAPWUNSCH determine the precision of the
static const double PSTEP = 0.00001;
static const double DELTAPWUNSCH = 0.00005;

#endif
