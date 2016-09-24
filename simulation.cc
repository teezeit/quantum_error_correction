#include "simulation.h"
#include "constants.h"
#include <omp.h>

// Constructor reads parameters from a config file.
Simulation::Simulation (void)
{

    string tmp;
    string item;
    int pos;
    ifstream f_conf("in/parameters.txt");

    getline(f_conf,tmp);
    pos = tmp.find(" = ");
    tmp = tmp.substr (pos + 3);
    stringstream ss(tmp);
    while(ss >> item)
    {
        XSize.push_back(atoi(item.c_str()));
    }

    getline(f_conf,tmp);
    pos = tmp.find(" = ");
    tmp = tmp.substr (pos + 3);
    iterations = atoi(tmp.c_str());

    getline(f_conf,tmp);
    pos = tmp.find(" = ");
    tmp = tmp.substr (pos + 3);
    pMin = atof(tmp.c_str());

    getline(f_conf,tmp);
    pos = tmp.find(" = ");
    tmp = tmp.substr (pos + 3);
    pMax = atof(tmp.c_str());

    getline(f_conf,tmp);
    pos = tmp.find(" = ");
    tmp = tmp.substr (pos + 3);
    pStep = atof(tmp.c_str());

    getline(f_conf,tmp);
    pos = tmp.find(" = ");
    tmp = tmp.substr (pos + 3);
    q_sqr_factor = atof(tmp.c_str());

    getline(f_conf,tmp);
    pos = tmp.find(" = ");
    tmp = tmp.substr (pos + 3);
    q_oct_factor = atof(tmp.c_str());

    f_conf.close();
}


// Simulate a simple decoding algorithm. Run simulation and save result.
void Simulation::run (void)
{

    // Initialize the random number generator seed
    srand (time(0));
    // Initialize constant random number seed
    //srand(8844);
    //Initialize Runtime TImer

    clock_t start, end;
    start = clock();


    //loop for different XSizes

    for (int i = 0; i < (int) XSize.size(); i++)
    {
        //data file
        ofstream f_out;
        stringstream ss;
        //add a 0 to the name if XSize is smaller than 0
        if(XSize[i] < 10)
        {
            ss << "out/results_0" << XSize[i] << ".txt";
        }
        else
        {
            ss << "out/results_" << XSize[i] << ".txt";

        }
        //open data file stream
        f_out.open (ss.str().c_str());


        //Print Information to screen

        cout   << "RUNNING SIMULATION:   XSize=" << XSize[i] << "   iterations=" << iterations << "   pMin=" << pMin << "   pMax=" << pMax << "   pStep=" << pStep << endl;


        cout << "ERROR [%]:\t\t\t\t\t\t\tSUCCESS [%]"<<endl<<"Level\tEffective\tper_Step\tSQR\tOCT\tSuccess"<< endl;


        // loop for different p's

        for (double p = pMin; p <= pMax + 0.00001; p += pStep)
        {

            //###################
            //configuration

            //number of succesful iterations
            int cntSucc = 0;


            //Error probability per time step
            double p_per_step = p;//determine_p_per_step(p,XSize[i]);

            //SquareSyndrome measurement error probability
            double q_SQR = p_per_step*q_sqr_factor;
            //OctogonSyndrome measurement error probability
            double q_OCT = p_per_step*q_oct_factor;

            //effective error probability
            double p_effective = determine_p_effective(p_per_step,XSize[i]);


            //Measuerement and Error Repetitions
            int Timesteps = XSize[i];


            //####################

            //Calculate all the distances in the history array

            struct Distance_and_ID returnDistandID = Distance(XSize[i], Timesteps, p_per_step,  q_SQR, q_OCT);
            float** D = returnDistandID.Distance_Array;
            int numberofids = returnDistandID.number_of_IDs;



            //loop the iterations with
            //enabled multithreading

            //omp_set_num_threads(2);
#pragma omp parallel for schedule(dynamic) default(shared) reduction(+:cntSucc)
            for (int iterate = 1; iterate <= iterations; iterate++ )
            {
                //There is a lot of debugging in this loop
                //important tasks are marked with ######

                //test the parallel computing
                //cout << omp_get_num_threads() << endl;

                if (debug4)
                {
                    cout<<"iteration: "<<iterate<<endl;

                }

                if (debug1) cout << "1. ERROR WAS GENERATED:" << endl;


                //############
                //Build a lattice
                Lattice myLattice (XSize[i],Timesteps);
                //############


                if (debug1) myLattice.printState();
                if (debug3)
                {
                    cout<<"Before Correction:"<<endl;
                    myLattice.printStateQB();
                }
                if (debug1) cout << "2. Z ERRORS CORRECTED USING SITE STABILIZERS:" << endl;


                //############
                //evolute and correct the lattice
                myLattice.evolute_and_correct(p_per_step,q_SQR,q_OCT,D);
                //############


                if (DEBUG2D)
                {
                    cout<<"Corrected"<<endl;
                    cout<<"   "<<": ";
                    myLattice.print2DState();
                }

                if (debug1) myLattice.printState();
                if (debug3)
                {
                    cout<<"After Correction:"<<endl;
                    myLattice.printStateQB();
                }

                //#########
                //do one last round of perfect correction
                myLattice.perfect_correction();
                //#########



                if (debug4)
                {
                    if ( !myLattice.is_corrected() )
                    {
                        cout<<"Not Corrected"<<endl;
                        assert(0);
                    }
                }


                //##########
                //Check success of the recovery
                if (myLattice.success ())
                {
                    cntSucc++;
                    //##########
                    if (debug1) cout << "SIMULATION REPORTS SUCCESS!!!" << endl;

                    if(debug9)
                    {
                        char in;
                        cout<< "Correction Succeded"<<endl;
                        cin>>in;
                    }
                }
                else
                {

                    if(debug7)
                    {
                        char in;
                        cout<< "Correction Failed"<<endl;
                        cin>>in;
                    }
                    if (debug1) cout << "SIMULATION REPORTS FAILURE!!!" << endl;
                }

                //Progressbar
                if (PROGRESSBAR && !(iterations%100))
                {
                    float progress = (float)iterate/(float)iterations;
                    cout<<"\rprogress[%]: "<<setprecision(3)<<100* progress;
                    //if(iterate == iterations-1){cout<<"\r";}

                }

            }// end of loop iterations


            //Print to Screen
            if(PROGRESSBAR)
            {
                cout<<"\r                 \r";
            }
            cout<< 100 * p <<"\t"<<setprecision(5)<<p_effective * 100<<"\t\t" <<p_per_step * 100<<  "\t\t"<<100*q_SQR<<"\t"<<100*q_OCT<<"\t\t"<< (double) (100 * cntSucc) / iterations  << endl;

            //Print to file
            f_out << 100 * p_per_step<< " " <<  (double) (100 * cntSucc) / iterations << endl;


            //clean memory (the distances)
            for (  int j = 0; j < numberofids; j++)
            {
                delete [] D[j];
            }
            delete [] D;

        }//end of loop p

        //close file
        f_out.close();

    }//end of loop XSize




    //Runtime Timer end
    //print Runtime information
    end = clock();
    if(RUNTIME)
    {
        cout<<endl<<"Computation Ended. Resume"<<endl<<"Starttime: "<<start<<endl;
        cout<<"Endtime:"<<end<<endl;
        cout<<"Clocks_per_Sec: "<<CLOCKS_PER_SEC<<endl;
        cout<<"Runtime ca.: " <<(double) (end/CLOCKS_PER_SEC)<<" s"<<endl;
    }

}//end of run()





//Determines the effective error probability for a given p and T per qubit
//If an error occurred 2,4,6,... times the qubit is error free
//So only if an error occures 1,3,5,.. times it really is an error
//This sums up to an effective error probability
double Simulation::determine_p_effective (double p_per_step, int T_in)
{



    //hier errechne wahrscheinlichkeit nach T schritten
    double p_eff_T = 0.0000f;

    for (int k = 0; k <= T_in; k++)
    {
        //only sum up odd k's
        if (k%2)
        {
            double k_loop =  k;

            double q = 1 - p_per_step;
            double l = T_in -k_loop;

            double prod1 = ncr(T_in,k_loop);
            double prod2 = pow(p_per_step,k_loop);
            double prod3 = pow(q,l);

            p_eff_T += prod1*prod3*prod2;

        }//endofif
    }//endoffor


    return p_eff_T;

}//






//for a given p and T we have an effective error probability
//Determines the probability p_per_step we need for an effective
//error probability after T steps
double Simulation::determine_p_per_step (double p_G, int T_in)
{

    double p_step = PSTEP;
    double p_max = PMAX;
    double T = (double) T_in;
    //runs


    //Probability nach T-steps
    vector <double> vec_p_eff_T;
    //Effektive Einzelwahrscheinlichkeit per step (Thats what we are looking for)
    vector <double> vec_p_per_step;



    for (double p_per_step = 0.0000f; p_per_step <= p_max; p_per_step  += p_step)
    {


        //hier calculate probability after T steps
        double p_eff_T = 0.0000f;

        for (int k = 0; k <= T; k++)
        {
            //only sum up odd k's
            if (k%2)
            {
                double k_loop =  k;

                double q = 1- p_per_step;
                double l = T-k_loop;

                double prod1 =ncr(T,k_loop);
                double prod2 = pow(p_per_step,k_loop);
                double prod3 = pow(q,l);

                p_eff_T += prod1*prod3*prod2;
                /*
                cout<<"kloop:"<<k_loop<<endl<<" ncr("<<T<<","<<k<<")="<<prod1<<endl<<"pow("<<p<<","<<k_loop<<")="<<prod2<<endl<<"pow("<<1-p<<","<<T-k_loop<<")="<<prod3<<endl;
                cout<<"p_eff="<<p_eff<<endl;
                */

            }//endofif
        }//endoffor

        vec_p_per_step.push_back(p_per_step);
        vec_p_eff_T.push_back(p_eff_T);

    }

    /*
    	//debug
    	cout<<"p\tp_eff"<<endl;
    	for (int i = 0; i < abs((int)prob.size()); i++) {

    		cout<<prob[i]<<"\t"<<prob_eff[i]<<endl;

    	}
    */

    //now check which p_real fits to which p_0

    double delta_p_wunsch = DELTAPWUNSCH;


    //nun suche ein p das das gewuenschte p gesamt ergibt
    while (true)
    {

        for (int i = 0; i < (int)vec_p_per_step.size(); i++)
        {

            double delta_p = p_G - vec_p_eff_T[i];
            //absolute value
            if (delta_p < 0)
            {
                delta_p *= -1;
            }

            if (delta_p < delta_p_wunsch )
            {
                return vec_p_per_step[i];
            }

        }

        //if there is no p found that fits

        delta_p_wunsch += 0.000005;


        if (delta_p_wunsch > 0.1)
        {
            cout<<"NO PROPER P FOUND!"<<endl;
            break;

        }


    }







    return p_G/ ( (double) T);




}//end of p_per_step











//Well, faculty
double Simulation::faculty (double n)
{
    if (n == 0)
    {
        return 1;
    }
    else
    {
        double n_minus_1 = n-1;
        n = n * faculty(n_minus_1);
    }

    return n;

}

//and n chose k
double Simulation::ncr (double n, double k)
{

    double n_min_k = n - k;

    double bin =  faculty(n) / (faculty(n_min_k) * faculty(k));

    return bin;


}




//This guy uses a Dijkstra algorithm to calculate all the distances
//from all vertices to all other vertices in the lattice and stores it
//in the struct. The second struct variable stores the dimensions of the
//Distancearray for the delete

//So it basically returns an array with all important distances

struct Distance_and_ID Simulation::Distance (int newXSize, int newT, double p_err, double q_SQR_err, double q_OCT_err)
{


    int xSize = 2*newXSize;
    int ySize = 2*newXSize;
    int T = newT;
    double p = p_err;
    double q_SQR = q_SQR_err;
    double q_OCT = q_OCT_err;




    //Vertices
    vector<Vertex*> vertices;
    //3d ID Tracker
    int ***IDTracker;
    IDTracker = new int** [xSize];

    int D_ID = 0;
    //fill history with false;
    //fill IDTracker with IDs
    if(debug10)
    {
        cout<<"Create history and ID Tracker:"<<endl;
    }

    for (int i = 0; i < xSize; i++)
    {
        IDTracker[i] = new int *[ySize];
        for (int j = 0; j < ySize; j++)
        {
            IDTracker[i][j] = new int [T+1];
            for (int t = 0; t < T+1; t++)
            {
                if(debug10)
                {
                    cout<<i<<j<<t<<" D_ID: "<<D_ID<<endl;
                }
                IDTracker[i][j][t] = D_ID;
                Vertex* v = new Vertex (i, j, t);
                vertices.push_back(v);
                D_ID++;
            }
        }
    }



    // Calculate relative edge weights of time and space edges
    double deltaXY = 1000;
    double deltaT_SQR = 1000;
    double deltaT_OCT = 1000;


    if (p!=0) deltaXY = log ((1 - p) / p);
    if (q_SQR!=0) deltaT_SQR = log ((1 - q_SQR) / q_SQR);
    if (q_OCT!=0) deltaT_OCT = log ((1 - q_OCT) / q_OCT);



//DIJKSTRA ALGORITHM
    //initialize size of adjacency list
    if(debug10)
    {
        cout<<"Create adjacency list lenght=D_ID= "<<D_ID<<endl<<endl;
    }
    adjacency_list_t adjacency_list(D_ID);

// Add all weighted edges
// remember to insert edges both ways for an undirected graph
    for (int t = 0; t < (T+1); t++)
    {
        for (int i = 0; i < xSize; i++)
        {
            for (int j = 0; j < ySize; j++)
            {
                //XOCTAGons do not have connections
                if(!(i%2) && !(j%2))
                {
                    if(debug10)
                    {
                        cout<<i<<j<<" XOctagon continue"<<endl;
                    }
                    continue;
                }

                //ZOctagon
                if (i%2 && j%2)
                {
                    if(debug10)
                    {
                        cout<<i<<j<<" ZOctagon: ";
                    }

                    //get the ID
                    int hereID = IDTracker[i][j][t];
                    //get the neigbor ids
                    int neighborID1 = IDTracker[(i+1)%xSize][j][t];
                    int neighborID2 = IDTracker[(i+xSize-1)%xSize][j][t];
                    int neighborID3 = IDTracker[i][(j+1)%ySize][t];
                    int neighborID4 = IDTracker[i][(j+ySize-1)%ySize][t];

                    if(debug10)
                    {
                        cout<<"hereID: "<<hereID<<" neighborIDS: "<<neighborID1<<" "<<neighborID2<<" "<<neighborID3<<" "<<neighborID4<<" ";
                    }

                    adjacency_list[hereID].push_back(neighbor(neighborID1, deltaXY));
                    adjacency_list[hereID].push_back(neighbor(neighborID2, deltaXY));
                    adjacency_list[hereID].push_back(neighbor(neighborID3, deltaXY));
                    adjacency_list[hereID].push_back(neighbor(neighborID4, deltaXY));

                    if(t!= 0)
                    {
                        int neighborID5 = IDTracker[i][j][t-1];
                        adjacency_list[hereID].push_back(neighbor(neighborID5, deltaT_OCT));
                        if(debug10)
                        {
                            cout<<neighborID5<<" ";
                        }

                    }
                    if(t!= T)
                    {
                        int neighborID6 = IDTracker[i][j][t+1];
                        adjacency_list[hereID].push_back(neighbor(neighborID6, deltaT_OCT));
                        if(debug10)
                        {
                            cout<<neighborID6<<" ";
                        }
                    }

                    if(debug10)
                    {
                        cout<<endl;
                    }


                    //horizontal connected Square
                }
                else if (i%2 && !(j%2))
                {
                    if(debug10)
                    {
                        cout<<i<<j<<"horizontal square: ";
                    }


                    //get the ID
                    int hereID = IDTracker[i][j][t];
                    //get the neigbor ids
                    int neighborID3 = IDTracker[i][(j+1)%ySize][t];
                    int neighborID4 = IDTracker[i][(j-1+ySize)%ySize][t];
                    //add the neighbors to the list
                    adjacency_list[hereID].push_back(neighbor(neighborID3, deltaXY));
                    adjacency_list[hereID].push_back(neighbor(neighborID4, deltaXY));

                    if(debug10)
                    {
                        cout<<"hereID: "<<hereID<<"neighborIDS: "<<neighborID3<<" "<<neighborID4<<" ";
                    }
                    if(t!= 0)
                    {

                        int neighborID5 = IDTracker[i][j][t-1];
                        adjacency_list[hereID].push_back(neighbor(neighborID5, deltaT_SQR));
                        if(debug10)
                        {
                            cout<<neighborID5<<" ";
                        }
                    }
                    if(t!= T)
                    {
                        int neighborID6 = IDTracker[i][j][t+1];
                        adjacency_list[hereID].push_back(neighbor(neighborID6, deltaT_SQR));
                        if(debug10)
                        {
                            cout<<neighborID6<<" ";
                        }
                    }
                    if(debug10)
                    {
                        cout<<endl;
                    }

                    //vertical connected square
                }
                else if (!(i%2) && j%2)
                {
                    if(debug10)
                    {
                        cout<<i<<j<<"vertical square: ";
                    }

                    //get the ID
                    int hereID = IDTracker[i][j][t];
                    //get the neigbor ids
                    int neighborID1 = IDTracker[(i+1)%xSize][j][t];
                    int neighborID2 = IDTracker[(i+xSize-1)%xSize][j][t];
                    //add the neigbors to the list
                    adjacency_list[hereID].push_back(neighbor(neighborID1, deltaXY));
                    adjacency_list[hereID].push_back(neighbor(neighborID2, deltaXY));

                    if(debug10)
                    {
                        cout<<"hereID: "<<hereID<<"neighborIDS: "<<neighborID1<<" "<<neighborID2<<" ";
                    }


                    if(t!= 0)
                    {
                        int neighborID5 = IDTracker[i][j][t-1];
                        adjacency_list[hereID].push_back(neighbor(neighborID5, deltaT_SQR));
                        if(debug10)
                        {
                            cout<<neighborID5<<" ";
                        }

                    }
                    if(t!= T)
                    {
                        int neighborID6 = IDTracker[i][j][t+1];
                        adjacency_list[hereID].push_back(neighbor(neighborID6, deltaT_SQR));
                        if(debug10)
                        {
                            cout<<neighborID6<<" ";
                        }
                    }
                    if(debug10)
                    {
                        cout<<endl;
                    }

                }
                else
                {
                    cout<<"Assert"<<endl;
                    assert(0);
                }


            }//end of i
        }//end of j
    }//end of t
//FINISHEd adjecency list


    //Now initialize the edgeLenghts and the IDStorage

    float** D_edgeLengths;

    D_edgeLengths = new float*[D_ID];
    for (int i = 0; i < D_ID; i++)
    {
        D_edgeLengths[i] = new float[D_ID];
        for (int j = 0; j < i; j++)
        {
            D_edgeLengths[i][j] = -1.0;
        }
    }

    /*
    double D_edgeLengths[D_ID][D_ID];

    for (int i = 0; i < D_ID; i++) {
        for (int j = 0; j < D_ID; j++) {
          D_edgeLengths[i][j] = -1.0;
        }
      }

    */

    // Add all weighted edges (u,v) with weight w in complete graph G
    //initialize parallel computing
    //every thread uses the same lattice
    //this is possible, because every thread operates on a different pointer
#pragma omp parallel for schedule(dynamic) default(shared)
    for (int i = 0; i < D_ID; i++)
    {
        //if XOctagon, kill dat thing
        if (!(vertices[i]->x%2) && !(vertices[i]->y%2))
        {
            continue;
        }

        //for every start vertex i get the ID
        int startID = i;
        //and calculate all the distances
        std::vector<weight_t> min_distance;
        std::vector<vertex_t> previous;
        DijkstraComputePaths(startID, adjacency_list, min_distance, previous);

        //Remember, the Matrix is Triangular

        for (int j = 0; j < i; j++)
        {
            //if i=j or XOctagon
            if (i == j || (!(vertices[j]->x%2) && !(vertices[j]->y%2)))
            {
                continue;
            }

            //for every end vertex j get the ID
            int endID = j;
            //And get the minimum distance
            float dist = (float) min_distance[endID];


            //then save it
            D_edgeLengths[i][j] = dist;



        }//end of j
    }//end of i



//Delete Vertices
    for (int i = 0; i < D_ID; i++)
    {
        delete vertices[i];
    }


// Delete IDTracker
    for (int i = 0; i < xSize; i++)
    {
        for (int j = 0; j < ySize; j++)
        {
            delete IDTracker[i][j];
        }
        delete IDTracker[i];
    }
    delete IDTracker;

    struct Distance_and_ID returnDistandID;

    returnDistandID.Distance_Array = D_edgeLengths;
    returnDistandID.number_of_IDs = D_ID;

    return returnDistandID;

}//end of Distance





//end of file simulation.cc
