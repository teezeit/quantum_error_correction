/*
qecsim.cc
Tobias Hoelzer
September 2013
Institute for Quantum Information
RWTH Aachen University
Bachelor Thesis: Study of [[4,2,2]]-concatenated toric code for a scalable circuit-QED architecture



Based on the work of:
qecsim.cc
Martin Suchara
November 2011
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include "simulation.h"
#include "constants.h"


using namespace std;


int main(int argc, char **argv)
{

    if (argc != 1)
    {
        cout << "Parameters will be ignored." << endl;
    }

    Simulation mySimulation;
    mySimulation.run();

    return 0;
}
