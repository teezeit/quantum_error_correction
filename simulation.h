#ifndef SIMULATION_H
#define SIMULATION_H

#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <assert.h>
#include <vector>
#include <set>
#include <map>
#include <fstream>
#include <climits>
#include <bitset>
#include <time.h>
#include <inttypes.h>
#include <iomanip>      // std::setprecision
#include "operator.h"
#include "lattice.h"




using namespace std;

struct Distance_and_ID
{
    float** Distance_Array;
    int number_of_IDs;
};


class Simulation
{
public:
    vector<int> XSize;
    int iterations;
    double pMin, pMax, pStep, q_sqr_factor, q_oct_factor;


    Simulation (void);
    ~Simulation (void) { }
    void run (void);

    struct Distance_and_ID Distance(int newXSize, int newT, double p_err, double q_SQR_err, double q_OCT_err);


    //for the error in time
    //is not really used in the final simulation
    //(only printed to screen)
    double faculty (double n);
    double ncr (double n, double k);
    double determine_p_per_step (double p_G, int T_in);
    double determine_p_effective (double p_G, int T_in);
};


#endif
