#include "lattice.h"
#include "constants.h"


// Constructor generates a unit cell. In case of the toric code this is a single qubit.
//Cell::Cell (double & p) {
Cell::Cell (void)
{
    ErrInit    = new Operator (1);
    //ErrInit->generateRandom (p);
    // ErrCurrent is a copy of the random ErrInit
    //ErrCurrent = new Operator (*ErrInit);
    ErrCurrent = new Operator(1);
    ErrGuess   = new Operator (1);

}


// Destructor for a unit cell.
Cell::~Cell (void)
{
    delete ErrInit;
    delete ErrCurrent;
    delete ErrGuess;
}


// Constructor of a vertex representing an operator in the transformed lattice.
Vertex::Vertex (int newX, int newY, int newT)
{
    x  = newX;
    y  = newY;
    t  = newT;
}


// Constructor generates the lattice.
Lattice::Lattice(int newXSize, int newT)
{

    xSize = 2*newXSize;
    ySize = 2*newXSize;
    T	= newT;


// zSize: 0 is the syndrome level, 1,2,3,4 are the qubit indices
    zSize = 5;


    cells = new Cell***[xSize];
    for (int i = 0; i < xSize; i++)
    {
        cells[i] = new Cell**[ySize];
        for (int j = 0; j < ySize; j++)
        {
            cells[i][j] = new Cell *[zSize];
            for ( int k = 0; k < zSize; k++)
            {
                //if k= 0     OR  i,j both odd OR    both even DO NOT generate a new cell (qubit)
                if ( k == 0 || (i%2 && j%2)  || (!(i%2) && !(j%2)))
                {
                    continue;
                }

                //otherwise, if only one of i,j is odd  and k > 0 generate a new cell (qubit)
                if(debug6)
                {
                    cout<< i<<j<<k<<" ";
                }
                cells[i][j][k] = new Cell ();

            }//end of k
        }//end of j
    }// end of i


}


// Destructor of the lattice.
Lattice::~Lattice(void)
{

    for (int i = 0; i < xSize; i++)
    {
        for (int j = 0; j < ySize; j++)
        {
            for (int k = 0; k < zSize; k++)
            {
                if ((k == 0 )  || (i%2 && j%2) || (!(i%2) && !(j%2)))
                    continue;
                delete cells[i][j][k];
            }
            delete[] cells[i][j];
        }
        delete[] cells[i];
    }
    delete cells;

}




// Calculates Octaeder plaquette syndrome at (around) the specified coordinates.
bool Lattice::OCTsyndromeZZZZZZZZ (int xLoc, int yLoc, double q_err)
{
//both i,j odd


//cout<<" xLoc "<< xLoc <<"  xSize "<< xSize <<" yLoc  " <<yLoc <<" ySize  "<< ySize << " (xLoc % 2)  "<< (xLoc % 2)<<"   (yLoc % 2) "<< (yLoc % 2)<<endl;


    assert (xLoc < xSize && yLoc < ySize && (xLoc % 2) && (yLoc % 2));

    //q could be modeled here as well
    double q = q_err;


    bool result = false;

    Operator *P = new Operator(8);
    P->ops[0]  = Z;
    P->ops[1]  = Z;
    P->ops[2]  = Z;
    P->ops[3]  = Z;
    P->ops[4]  = Z;
    P->ops[5]  = Z;
    P->ops[6]  = Z;
    P->ops[7]  = Z;

    Operator *plaquette = new Operator(8);

//Remember QB indices are
//1--2
//|  |
//3--4

    /*/DEBUG
    //test where the plaquette tests the cells
    cout << "xLoc  " << xLoc << "    yLoc-1  " << yLoc-1 << endl;
    cout << "xLoc+1modxSiyze  " << (xLoc+1)%xSize << "    yLoc  " << yLoc << endl;
    cout << "xLoc  " << xLoc << "    yLoc+1modySize  " << (yLoc+1)%ySize << endl;
    cout << "xLoc-1  " << xLoc-1 <<"  yLoc  " << yLoc <<  endl;
    */



//One step left (if i is line, y is column) to the square, get QB 2,4
    plaquette->ops[0]  = cells[xLoc][(yLoc-1+ySize)%ySize][2]->ErrCurrent->ops[0];
    plaquette->ops[1]  = cells[xLoc][(yLoc-1+ySize)%ySize][4]->ErrCurrent->ops[0];
    //One step down get QB 1,2
    plaquette->ops[2]  = cells[(xLoc+1)%xSize][yLoc][1]->ErrCurrent->ops[0];
    plaquette->ops[3]  = cells[(xLoc+1)%xSize][yLoc][2]->ErrCurrent->ops[0];
// One step right, get QB 1,3
    plaquette->ops[4]  = cells[xLoc][(yLoc+1)%ySize][1]->ErrCurrent->ops[0];
    plaquette->ops[5]  = cells[xLoc][(yLoc+1)%ySize][3]->ErrCurrent->ops[0];
//One step up get QB 3,4
    plaquette->ops[6]  = cells[(xLoc-1+xSize)%xSize][yLoc][3]->ErrCurrent->ops[0];
    plaquette->ops[7]  = cells[(xLoc-1+xSize)%xSize][yLoc][4]->ErrCurrent->ops[0];



	//Does it commute?
    result = !P->commute(*plaquette);

    //generate noisy syndrome measurement with propapility q
    double r = (double)rand() / ((double)RAND_MAX + 1);
    //if it shall have noise, then flip it
    if (r < q)
    {
        result = !result;
    }


    delete P;
    delete plaquette;
    return result;
}

// Calculates SQUARE plaquette syndrome at (around) the specified coordinates.
bool Lattice::SQsyndromeZZZZ (int xLoc, int yLoc, double q_err)
{
//exactly one of i,j odd and one even
    assert (xLoc < xSize && yLoc < ySize && ( ((xLoc % 2) && !(yLoc % 2)) || ( !(xLoc % 2) && (yLoc % 2))) );


    //q could be modeled here as well
    double q = q_err;
    bool result = false;

    Operator *P = new Operator(4);
    P->ops[0]  = Z;
    P->ops[1]  = Z;
    P->ops[2]  = Z;
    P->ops[3]  = Z;


    Operator *plaquette = new Operator(4);

//Remember QB indices are
//1--2
//|ij|
//3--4
    plaquette->ops[0]  = cells[xLoc][yLoc][1]->ErrCurrent->ops[0];
    plaquette->ops[1]  = cells[xLoc][yLoc][2]->ErrCurrent->ops[0];
    plaquette->ops[2]  = cells[xLoc][yLoc][3]->ErrCurrent->ops[0];
    plaquette->ops[3]  = cells[xLoc][yLoc][4]->ErrCurrent->ops[0];



    result = !P->commute(*plaquette);




    //generate noisy syndrome measurement with propapility q
    double r = (double)rand() / ((double)RAND_MAX + 1);
    //if it shall have noise, then flip it
    if (r < q)
    {
        result = !result;
    }


    delete P;
    delete plaquette;
    return result;
}


// Calculates a signature which is a concatentation of all site and plaquette syndromes in the lattice.
vector<bool> Lattice::getSignature (void)
{

    // The signature starts out empty
    vector<bool> signature;



    // Add ZZZZ syndromes
    for (int i = 1; i < xSize; i+=2)
    {
        for (int j = 1; j < ySize; j+=2)
        {
            // signature.push_back(SQsyndromeZZZZ(i,j));
        }
    }

    return signature;
}


//  Transform X or Z errors at the specified location.
void Lattice::correct ( int xLoc, int yLoc, int QubitNr)
{
    if (debug4)
    {
        cout<< "xsize,ysixze,ysize "<<xSize<<ySize<<zSize<<endl;
        cout<< "xLoc,yLoc,QBNr "<<xLoc<<yLoc<<QubitNr<<endl;
    }
    assert (xLoc < xSize && yLoc < ySize && QubitNr < zSize);

//X*X=I
    pauli p = X;

    // Correct the error
    cells[xLoc][yLoc][QubitNr]->ErrCurrent->ops[0] = cells[xLoc][yLoc][QubitNr]->ErrCurrent->ops[0] * p;
    cells[xLoc][yLoc][QubitNr]->ErrGuess->ops[0] = cells[xLoc][yLoc][QubitNr]->ErrGuess->ops[0] * p;

    if (debug1)
    {
        cout << " -> correcting " << p << " at x=" << xLoc << " y=" << yLoc << " QB_Nr="<<QubitNr<<endl;
    }
    if (debug3)
    {
        cout << " -> correcting " << p << " at x=" << xLoc << " y=" << yLoc << " QB_Nr="<<QubitNr<<endl;
    }

}//end of correct()



//evolutes the lattice in time
void Lattice::evolute(double p_err)
{
    //qubit error
    double p = p_err;

    for (int i = 0; i < xSize; i++)
    {
        for (int j = 0; j < ySize; j++)
        {
            for (int k = 0; k < zSize; k++)
            {


				//produces 2D cross section through the spacetime
				//very nice debugging feature
                if (DEBUG2D)
                {

                    //only second line nr 1 gets errors
                    if (i==1 && !(j%2))  //j=1,3,5is syndrome,0,2,4 is qubit
                    {
                        cells[i][j][1]->ErrInit->generate_random_evolution(p);
                        //and copy into Current
                        cells[i][j][1]->ErrCurrent->ops[0] = cells[i][j][1]->ErrInit->ops[0];
                        cells[i][j][2]->ErrInit->generate_random_evolution(p);
                        //and copy into Current
                        cells[i][j][2]->ErrCurrent->ops[0] = cells[i][j][2]->ErrInit->ops[0];
                    }




                }
                else
                {
                    //if k= 0     OR  i,j both odd OR    both even DO NOTHING
                    if ( k == 0 || (i%2 && j%2)  || (!(i%2) && !(j%2)))
                    {
                        continue;
                    }
                    //otherwise, throw error at qubit
                    cells[i][j][k]->ErrInit->generate_random_evolution(p);
                    //and copy into Current
                    cells[i][j][k]->ErrCurrent->ops[0] = cells[i][j][k]->ErrInit->ops[0];


                }

                //if that doesnt work, we need a copy method for operators (seems to work:-)
            }//endk
        }//endj
    }//endi









}//end of evolute










// Calls the mathcing algorithm to find pairs of squares with nontrivial syndromes,
// and calls the correction subroutine.
void Lattice::evolute_and_correct (double p_err, double q_err_SQR, double q_err_OCT, float **Dist)
{

    //proper algorithm for q
    double p = p_err;
    double q_SQR = q_err_SQR;
    double q_OCT = q_err_OCT;


    /*
     General Structure:
     
    find proper algorithm for q, (is maybe done in simulation.cc)

    create history

    loop:
    evolute,measure,save

    create vertices

    match vertices

    correct


    */

    //Create 3d	history AND 3d ID Tracker
    bool*** history;
    history = new bool** [xSize];

    int*** IDTracker;
    IDTracker = new int** [xSize];

    int ID = 0;
    //fill history with false;
    //fill IDTracker with IDs
    if(debug10)
    {
        cout<<"Create history and ID Tracker:"<<endl;
    }

    for (int i = 0; i < xSize; i++)
    {
        history[i] = new bool *[ySize];
        IDTracker[i] = new int *[ySize];
        for (int j = 0; j < ySize; j++)
        {
            history[i][j] = new bool [T+1];
            IDTracker[i][j] = new int [T+1];
            for (int t = 0; t < T+1; t++)
            {
                if(debug10)
                {
                    cout<<i<<j<<t<<" ID: "<<ID<<endl;
                }
                history[i][j][t] = false;
                IDTracker[i][j][t] = ID;
                ID++;

            }
        }
    }


    //now evolute the lattice, measure after each time step and save the syndromes in the history.

    if (DEBUG2D)
    {
        cout<<"Evolution"<<endl;
    }


    //last T slice must be error free!
    for (int t = 0; t < T; t++)
    {

        //evolute
        this->evolute(p);


        if (DEBUG2D)
        {
            cout<<"T="<<t<<": ";
            print2DState();
        }

        //measure and save

        for (int i = 0; i < xSize; i++)
        {
            for (int j = 0; j < ySize; j++)
            {

                // if i,j both odd =square do sqr Zsyndrome, if that is true
                // OR
                // if on Zoctagon (NOT Xoctagon) do oct Zsyndrome, if that is true
                if ( (which_cell_type (i,j) == sqr && SQsyndromeZZZZ(i,j,q_SQR)) || ( (i%2) && (j%2) && OCTsyndromeZZZZZZZZ(i,j,q_OCT)))
                {
                    //save in history
                    history[i][j][t] = true;
                    if (debug3)
                    {
                        cout<<"history i,j,t = "<<i<<" "<<j<<" "<<t<<" = "<<history[i][j][t]<<endl;
                    }
                }


            }//end of j
        }//end of i

    }//end of for t evolute and measure


//Maybe here, debug print history

    //Print history
    if (DEBUG2D)
    {
        cout<<"History"<<endl;
        for (int t = 0; t < T+1; t++)
        {

            cout<<"T="<<t<<": ";
            for (int i = 0; i < xSize; i++)
            {
                for (int j = 0; j < ySize; j++)
                {


                    if (i==1)  //j=1,3,5is syndrome,0,2,4 is qubit
                    {
                        cout<< " "<<history[i][j][t]<<" ";


                    }

                }
                if (i==1)
                {
                    cout<<endl;
                }
            }//end of i

        }//endof t
    }






    if (DEBUG2D)
    {
        cout<<"Vertices"<<endl;
    }

    // Identify vertices that represent squares with change in syndrome
    vector<Vertex*> vertices;
    for (int t = 0; t < (T+1); t++)
    {
        if (DEBUG2D)
        {
            cout<<"T="<<t<<": ";
        }

        for (int i = 0; i < xSize; i++)
        {
            for (int j = 0; j < ySize; j++)
            {

                bool is_vertex = false;
                //is on Syndrome ij?
                // if i,j both odd(Zocto) OR exactly one odd and one even (Zsquare) ->check if vertex
                if ( (i%2 && j%2 ) ||  which_cell_type (i,j) == sqr )
                {
                    /*
                    in the first slice, every syndrome is a vertex.
                    otherwise, check the current t and the one before
                    if the syndrome changed, if->then it is considered a vertex
                    */
                    if((t==0 && history[i][j][t]) || (t!=0 && history[i][j][t]!=history[i][j][t-1]))
                    {
                        is_vertex = true;
                    }//end if is vertex

                }//end of if Syndrome ij

                //Print vertex history
                if (DEBUG2D)
                {
                    if (i==1)
                    {
                        if(is_vertex)
                        {
                            cout<<" "<<vertices.size();
                            //if(vertices.size()<10){cout<<" ";}
                            cout<<" ";
                        }
                        else
                        {
                            cout<<" . ";
                        }

                    }
                }



                //if it is a vertex, save it
                if (is_vertex)
                {
                    if (debug1) cout << "Is_Vertex at: "  << " i=" << i << " j= " << j << " t= " << t<< endl;
                    Vertex* v = new Vertex (i, j, t);
                    vertices.push_back(v);
                }//end if is vertex


            } //end of loop j
        }//end of loop i

        if (DEBUG2D)
        {
            cout<<endl;
        }

    }//end of loop t



//now prepare for the minimum weight matching

    int numVert = vertices.size();
    int edges = 0;

    float** edgeLengths;
    edgeLengths = new float*[numVert];
    for (int i = 0; i < numVert; i++)
    {
        edgeLengths[i] = new float[numVert];
        for (int j = 0; j < numVert; j++)
        {
            edgeLengths[i][j] = -1.0;
        }
    }



    // Add all weighted edges (u,v) with weight w in complete graph G
    for (int i = 0; i < numVert; i++)
    {
        //for every start vertex i get the ID
        int startID = IDTracker[vertices[i]->x][vertices[i]->y][vertices[i]->t];
        for (int j = 0; j < numVert; j++)
        {
            if (i == j)
            {
                continue;
            }

            //for every end vertex j get the ID
            int endID = IDTracker[vertices[j]->x][vertices[j]->y][vertices[j]->t];
            //And get the minimum distance
            float dist = 0.000f;
            if(startID>endID)
            {
                dist = Dist[startID][endID];
            }
            else
            {
                dist = Dist[endID][startID];

            }


            //then save it
            edgeLengths[i][j] = dist;
            edges++;



        }//end of j
    }//end of i




// Call minimum weight perfect matching
    struct PerfectMatching::Options options;
    struct GeomPerfectMatching::GPMOptions gpm_options;
    options.verbose = false;
    PerfectMatching *pm = new PerfectMatching(numVert,edges);

    for (int i = 0; i < numVert; i++)
    {
        for (int j = 0; j < numVert; j++)
        {
            double dist = edgeLengths[i][j];
            if (dist >= 0)
                pm->AddEdge(i,j,dist);
        }
    }



    pm->options = options;
    pm->Solve();



    if (debug4 || DEBUG2D)
    {
        for (int i = 0; i < numVert; i++ )
        {
            int j = pm->GetMatch(i);
            if (i<j)
            {
                cout<<"Matched: "<<i<<" and "<<j<<endl;
            }
        }
    }

    // Call the correction subroutine for all matched vertices
    for (int i = 0; i < numVert; i++ )
    {
        // i and j are matched
        int j = pm->GetMatch(i);
        //only correct:
        //-once
        // -those vertices whose connections has space components (not only time) <-> d_x OR d_y >0
        //-the last history slice is special:
        //a path is in space and time, then it is corrected BUT:
        // the last slice is without
        //error by default SO: only correct which matching partner  vertex is not in the last slice == correct in space AND time


        if ( !(vertices[i]->t == T && vertices[j]->t == T) && (i < j) && (j < numVert) && ( (vertices[i]->x != vertices[j]->x) || (vertices[i]->y != vertices[j]->y) )   )
        {
            if (debug4 ||DEBUG2D)
            {
                cout<<"Corrected: "<<i<<" to "<<j<<endl;
            }
            correctLine (vertices[i]->x, vertices[i]->y, vertices[j]->x, vertices[j]->y);


        }//end of if correct
    }//end of numVertices

    for (int i = 0; i < numVert; i++)
    {
        delete vertices[i];
    }

    delete pm;

    for (int i = 0; i < numVert; i++)
    {
        delete edgeLengths[i];
    }
    delete edgeLengths;


// Delete syndrome history AND IDTracker
    for (int i = 0; i < xSize; i++)
    {
        for (int j = 0; j < ySize; j++)
        {
            delete history[i][j];
            delete IDTracker[i][j];
        }
        delete history[i];
        delete IDTracker[i];
    }
    delete history;
    delete IDTracker;


}// end of evolute and correct





// Does one round of perfect syndrome measurement and corrects the lattice
// Without it, it might not be part of the stabilizer
//(errors could still be in it)
void Lattice::perfect_correction(void)
{

    if (debug3)
    {
        cout<<"Evolute and Correct "<<endl;
    }

    double q_SQR = 0;
    double q_OCT = 0;



    //Create 1d
    bool*** history;


    history = new bool** [xSize];

    //fill history with false;
    for (int i = 0; i < xSize; i++)
    {
        history[i] = new bool *[ySize];
        for (int j = 0; j < ySize; j++)
        {
            history[i][j] = new bool [T+1];
            for (int t = 0; t < 1; t++)
            {
                history[i][j][t] = false;
            }
        }
    }


    //measure and save
    for (int t = 0; t < 1; t++)
    {
        for (int i = 0; i < xSize; i++)
        {
            for (int j = 0; j < ySize; j++)
            {

                // if i,j both odd =square do sqr Zsyndrome, if that is true
                // OR
                // if on Zoctagon (NOT Xoctagon) do oct Zsyndrome, if that is true
                if ( (which_cell_type (i,j) == sqr && SQsyndromeZZZZ(i,j,q_SQR)) || ( (i%2) && (j%2) && OCTsyndromeZZZZZZZZ(i,j,q_OCT)))
                {
                    //save in history
                    history[i][j][t] = true;
                    if (debug3)
                    {
                        cout<<"history i,j,t = "<<i<<" "<<j<<" "<<t<<" = "<<history[i][j][t]<<endl;
                    }
                }


            }//end of j
        }//end of i
    }//end of t




    // Identify vertices that represent squares with change in syndrome
    vector<Vertex*> vertices;
    if (debug4)
    {
        cout<<"Vertices"<<endl;
    }

    for (int t = 0; t < 1; t++)
    {
        for (int i = 0; i < xSize; i++)
        {
            for (int j = 0; j < ySize; j++)
            {

                bool is_vertex = false;
                // if i,j both odd check if vertex
                if (( (i%2 && j%2 ) || ( (!(i%2) && j%2)  || (i%2 && !(j%2)) ) ) && history[i][j][t])
                {
                    is_vertex = true;

                }//end of is vertex






                //if it is a vertex, save it
                if (is_vertex)
                {
                    if (debug3) cout << "Is_Vertex at: "  << " i=" << i << " j= " << j << " t= " << t<< endl;
                    Vertex* v = new Vertex (i, j, t);
                    vertices.push_back(v);

                }

            }//end of j

        }//end of i
    }//end of t







    int numVert = vertices.size();
    int edges = 0;

    double** edgeLengths;
    edgeLengths = new double*[numVert];
    for (int i = 0; i < numVert; i++)
    {
        edgeLengths[i] = new double[numVert];
        for (int j = 0; j < numVert; j++)
        {
            edgeLengths[i][j] = -1.0;
        }
    }



    // Add all weighted edges (u,v) with weight w in complete graph G
    for (int i = 0; i < numVert; i++)
    {
        for (int j = 0; j < numVert; j++)
        {
            if (i == j)
            {
                continue;
            }



			//this is the distance in x y t
            double deltaX = min (abs (vertices[i]->x - vertices[j]->x), xSize - abs (vertices[i]->x - vertices[j]->x));
            //if both syndromes are both on the same odd column the distance in x increases by2
            //just a fancy specialty of this metrik (taxi-cab with holes)
            if ( (vertices[i]->y == vertices[j]->y) && (!(vertices[i]->y % 2)) )
            {
                deltaX = deltaX + 2;
            }

            double deltaY = min (abs (vertices[i]->y - vertices[j]->y), ySize - abs (vertices[i]->y - vertices[j]->y));
            //if the syndromes are both on the same odd line the distance in y increases by 2
            //just a fancy specialty of this metrik (taxi-cab with holes)
            if ( (vertices[i]->x == vertices[j]->x) && (!(vertices[i]->x % 2)) )
            {
                deltaY = deltaY + 2;
            }

            double deltaT = abs (vertices[i]->t - vertices[j]->t);

            double dist   = deltaX + deltaY + deltaT;
            edgeLengths[i][j] = dist;
            edges++;
            if (debug1) cout << "Matching: i=" << i << " j= " << j << " dist=" << dist << endl;
        }
    }



// Call minimum weight perfect matching
    struct PerfectMatching::Options options;
    struct GeomPerfectMatching::GPMOptions gpm_options;
    options.verbose = false;
    PerfectMatching *pm = new PerfectMatching(numVert,edges);

//cout<<"numVert: "<<numVert<<endl;

    for (int i = 0; i < numVert; i++)
    {
        for (int j = 0; j < numVert; j++)
        {
            double dist = edgeLengths[i][j];
            if (dist >= 0)
                pm->AddEdge(i,j,dist);
        }
    }



    pm->options = options;
    pm->Solve();

    if (debug4)
    {
        for (int i = 0; i < numVert; i++ )
        {
            int j = pm->GetMatch(i);
            if (i<j)
            {
                cout<<"Matched: "<<i<<" and "<<j<<endl;
            }
        }
    }

    // Call the correction subroutine for all matched vertices
    for (int i = 0; i < numVert; i++ )
    {
        // i and j are matched
        int j = pm->GetMatch(i);

        if (  (i < j) && (j < numVert) && ( (vertices[i]->x != vertices[j]->x) || (vertices[i]->y != vertices[j]->y) )   )
        {
            if (debug4)
            {
                cout<<"Corrected: "<<i<<" to "<<j<<endl;
            }
            correctLine (vertices[i]->x, vertices[i]->y, vertices[j]->x, vertices[j]->y);

        }
    }

    for (int i = 0; i < numVert; i++)
    {
        delete vertices[i];
    }

    delete pm;

    for (int i = 0; i < numVert; i++)
    {
        delete edgeLengths[i];
    }
    delete edgeLengths;


// Delete syndrome history
    for (int i = 0; i < xSize; i++)
    {
        for (int j = 0; j < ySize; j++)
        {
            delete history[i][j];
        }
        delete history[i];
    }
    delete history;

    if (debug3)
    {
        cout<<"Finished Correction "<<endl;
    }
}








// Given the locations of two syndromes with a non-trivial syndrome, correct errors on a line
// connecting the two syndromes.
void Lattice::correctLine ( int x1Loc, int y1Loc, int x2Loc, int y2Loc)
{

    // Initialize x1, y1, x2, y2 so that we can correct errors on a line in
    // the SE or SW direction coming from coordinates 1 to coordinates 2

//x y-->
//I
//v    EAST
//  00 01 02
//N 10 11 13 South
//  20 21 22
//  West

    int x1, y1, x2, y2;
    if ( (y1Loc < y2Loc && y2Loc - y1Loc < ySize - y2Loc + y1Loc) ||
            (y2Loc < y1Loc && y1Loc - y2Loc > ySize - y1Loc + y2Loc) )
    {
        x1 = x1Loc;
        y1 = y1Loc;
        x2 = x2Loc;
        y2 = y2Loc;
    }
    else
    {
        x1 = x2Loc;
        y1 = y2Loc;
        x2 = x1Loc;
        y2 = y1Loc;
    }


//A LONG WALK TO FREEDOM ALGORITHM

// go from x1|y1 -> x2|y2
// y2 = y1+(...)%ySize

//DIRECTION
    // Determine if we need to go up, down or neither
    vertical_direction v_dir = v_where_to_go(x1, x2);
//up,down, neither, wrong_direction

//HORIZONTAL Y-ACHSIS j
    // Determine if we need to start from sqr or oct
    celltype h_start_type = which_cell_type ( x1, y1);
//sqr,oct,wrong_type
    // Determine if we need to start from sqr or oct
    celltype h_end_type = which_cell_type ( x1, y2);

//VERTICAL X-ACHSIS i
// Determine if we need to start from sqr or oct
    celltype v_start_type = h_end_type;
    // Determine if we need to start from sqr or oct
    celltype v_end_type = which_cell_type ( x2, y2);

//Determine distance between syndromes

    int delta_x = delta_xy(x2,x1,xSize);
    int delta_y = delta_xy(y2,y1,ySize);








    /*old debugging
      if (debug1) cout << "Correcting line x1=" << x1 << " y1=" << y1 << " x2=" << x2 << " y2=" << y2 << " east=" << east << ":" << endl;
    */

    if (debug3)
    {
        cout <<endl<< "Correcting line x1=" << x1 << " y1=" << y1 << " x2=" << x2 << " y2=" << y2
             << " direction= ";
        if (v_dir==up)
        {
            cout<< " up "<<endl;
        }
        else if (v_dir==down)
        {
            cout<< " down "<<endl;
        }
        else if (v_dir==neither)
        {
            cout<< " neither "<<endl;
        }
        cout<< "Start Horizontal= ";
        if (h_start_type==sqr)
        {
            cout<< " square "<<endl;
        }
        else if (h_start_type==octo)
        {
            cout<< " 0ctogon "<<endl;
        }
        cout<< "End Horizontal= ";
        if (h_end_type==sqr)
        {
            cout<< " square "<<endl;
        }
        else if (h_end_type==octo)
        {
            cout<< " 0ctogon "<<endl;
        }
        cout<< "Start Vertical= ";
        if (v_start_type==sqr)
        {
            cout<< " square "<<endl;
        }
        else if (v_start_type==octo)
        {
            cout<< " 0ctogon "<<endl;
        }
        cout<< "End Vertical= ";
        if (v_end_type==sqr)
        {
            cout<< " square "<<endl;
        }
        else if (v_end_type==octo)
        {
            cout<< " 0ctogon "<<endl;
        }
        cout<< "Delta_y= "<<delta_xy(y2,y1,ySize)<<endl;
        cout<< "Delta_x= "<<delta_xy(x2,x1,xSize)<<endl;
    }
/////////////////////////////////////////////HORIZONTAL
    // correct in the horizontal direction first
//Qubits to correct in horizontal y Achsis
//y1 ++ to y2


// use correct( x1, y1, qubitcorr)

    //do we need to correct in horizontal direction?
    if (delta_y != 0)
    {

        if (debug4)
        {
            cout<<"delta_y != 0"<<endl;
        }
        //then start correction

        //exitvariable
        bool exit = false;
        //exit next time variable
        bool exitnextloop = false;


        //laufvariable
        int Y = y1;
        //countingvariable counting the amount of corrections
        int l = 1;

        //while exit is false
        while (!exit)
        {

            if (debug4)
            {
                cout<<"Y, y2: "<<Y<<" "<<y2<<endl;
            }


            //only if on an octagon
            if ( which_cell_type(x1,Y) == octo)
            {

                if (debug4)
                {
                    cout<<" is octagon "<<endl;
                }

                //initialize the stuff!
                bool first = false;
                bool last = false;
                correctiontype whatcorrection = wrong_correctiontype;
                correctiontypeQB whatcorrectionqubit = wrong_correctiontypeQB;


                //check if first
                if ( l == 1)
                {
                    first = true;
                }

                //check if last
                if ( delta_xy(Y,y2,ySize) <= 1 )
                {
                    last = true;
                }

                //Determine correctiontype
                if ( last && h_end_type == octo )
                {
                    whatcorrection = before;
                }
                else if (first && h_start_type)
                {
                    whatcorrection = next;
                }
                else
                {
                    whatcorrection = four;
                }

                //Determine correctiontypequbit

                if ( (whatcorrection == before) && ( v_dir == up) )
                {
                    whatcorrectionqubit = before_up;
                }
                else if ( (whatcorrection == before) && ( v_dir == down  || v_dir == neither) )
                {
                    whatcorrectionqubit = before_down;
                }
                else if ( (whatcorrection == next) && ( v_dir == up) )
                {
                    whatcorrectionqubit = next_up;
                }
                else if ( (whatcorrection == next) && ( v_dir == down  || v_dir == neither) )
                {
                    whatcorrectionqubit = next_down;
                }
                else if ( (whatcorrection == four) && ( v_dir == up) )
                {
                    whatcorrectionqubit = four_up;
                }
                else if ( (whatcorrection == four) && ( v_dir == down  || v_dir == neither) )
                {
                    whatcorrectionqubit = four_down;
                }
                else
                {
                    assert(0);
                }

                if(debug4)
                {
                    if (first)
                    {
                        cout<< " first ";
                    }
                    else
                    {
                        cout<<" not first ";
                    }
                    if (last)
                    {
                        cout<< " last ";
                    }
                    else
                    {
                        cout<<" not last ";
                    }

                    if (whatcorrectionqubit == before_up)
                    {
                        cout<< " before_up ";
                    }
                    else if (whatcorrectionqubit == before_down)
                    {
                        cout<< " before_down ";
                    }
                    else if(whatcorrectionqubit == next_up)
                    {
                        cout<< " next_up ";
                    }
                    else if(whatcorrectionqubit == next_down)
                    {
                        cout<< " next_down ";
                    }
                    else if(whatcorrectionqubit == four_up)
                    {
                        cout<< " four_up ";
                    }
                    else if(whatcorrectionqubit == four_down)
                    {
                        cout<< " four_down ";
                    }
                    else
                    {
                        assert(0);
                    }
                }




                switch (whatcorrectionqubit)
                {
                case before_up:
                    correct(x1,(Y+ySize-1)%ySize, 2);
                    Y = (Y+1)%ySize;
                    break;


                case before_down:
                    correct(x1,(Y+ySize-1)%ySize, 4);
                    Y = (Y+1)%ySize;
                    break;
                case next_up:
                    correct(x1,(Y+1)%ySize, 1);
                    Y = (Y+1)%ySize;
                    break;
                case next_down:
                    correct(x1,(Y+1)%ySize, 3);
                    Y = (Y+1)%ySize;
                    break;
                case four_down:
                    correct(x1,(Y+ySize-1)%ySize, 4);
                    correct((x1+1)%xSize,Y, 1);
                    correct((x1+1)%xSize,Y, 2);
                    correct(x1,(Y+1)%ySize, 3);
                    Y = (Y+1)%ySize;
                    break;
                case four_up:
                    correct(x1,(Y+ySize-1)%ySize, 2);
                    correct((x1+xSize-1)%ySize,Y, 3);
                    correct((x1+xSize-1)%ySize,Y, 4);
                    correct(x1,(Y+1)%ySize, 1);
                    Y = (Y+1)%ySize;
                    break;
                default:
                    assert(0);
                }

                //if we were at an octagon it shure must have corrected something
                // (otherwise assert(0);)-->so counter up!
                l++;

            }
            else if ( which_cell_type(x1,Y) == sqr )
            {

                Y = (Y+1)%ySize;

            }
            else
            {
                assert(0);
            }// end if octo or square

            // this makes the loop run a last time after Y hits y2
            // No, I do not want to talk about it!
            if (exitnextloop)
            {
                exit = true;
            }
            if (Y == y2)
            {
                exitnextloop = true;
            }




        } //endofwhileY
    }//endofifdeltay


/////////////////////////////////////////////VERTICAL
    // correct in the vertical direction
//Qubits to correct in vertical x Achsis
// go from x1|y2 to x2|y2
// x1 --> x2


// use correct( x1, y1, qubitcorr)

    //do we even need to correct in vertical direction?
    if (delta_x != 0)
    {
        //then start correction

        //exitvariable
        bool exit = false;
        //exit next time variable
        bool exitnextloop = false;

        //laufvariable
        int X = x1;
        //countingvariable counting the amount of corrections
        int l = 1;

        // while exit is false
        while (!exit)
        {

            //only if on an octagon
            if ( which_cell_type(X,y2) == octo)
            {


                //initialize that stuff
                bool first = false;
                bool last = false;
                correctiontype whatcorrection = wrong_correctiontype;
                correctiontypeQB whatcorrectionqubit = wrong_correctiontypeQB;



                //check if first
                if ( l == 1)
                {
                    first = true;
                }
                //check if last
                if ( delta_xy(X,x2,xSize) <= 1 )
                {
                    last = true;
                }

                //Determine correctiontype
                if ( last && v_end_type == octo )
                {
                    whatcorrection = before;
                }
                else if (first && v_start_type)
                {
                    whatcorrection = next;
                }
                else
                {
                    whatcorrection = four;
                }

                //Determine correctiontypequbit

                if ( (whatcorrection == before) && ( v_dir == up) )
                {
                    whatcorrectionqubit = before_up;
                }
                else if ( (whatcorrection == before) && ( v_dir == down ) )
                {
                    whatcorrectionqubit = before_down;
                }
                else if ( (whatcorrection == next) && ( v_dir == up) )
                {
                    whatcorrectionqubit = next_up;
                }
                else if ( (whatcorrection == next) && ( v_dir == down ) )
                {
                    whatcorrectionqubit = next_down;
                }
                else if ( (whatcorrection == four) && ( v_dir == up) )
                {
                    whatcorrectionqubit = four_up;
                }
                else if ( (whatcorrection == four) && ( v_dir == down  ) )
                {
                    whatcorrectionqubit = four_down;
                }
                else
                {
                    assert(0);
                }





                switch (whatcorrectionqubit)
                {
                case before_up:
                    correct((X+xSize+1)%xSize,y2, 1);
                    X = (X+xSize-1)%xSize;
                    break;
                case before_down:
                    correct((X+xSize-1)%xSize,y2, 3);
                    X = (X+1)%xSize;
                    break;
                case next_up:
                    correct((X+xSize-1)%xSize,y2, 3);
                    X = (X+xSize-1)%xSize;
                    break;
                case next_down:
                    correct((X+xSize+1)%xSize,y2, 1);
                    X = (X+1)%xSize;
                    break;
                case four_down:
                    correct((X+xSize-1)%xSize,y2, 3);
                    correct(X,(y2+ySize-1)%ySize, 2);
                    correct(X,(y2+ySize-1)%ySize, 4);
                    correct((X+xSize+1)%xSize,y2, 1);
                    X = (X+1)%xSize;
                    break;
                case four_up:
                    correct((X+xSize-1)%xSize,y2, 3);
                    correct(X,(y2+ySize-1)%ySize, 2);
                    correct(X,(y2+ySize-1)%ySize, 4);
                    correct((X+xSize+1)%xSize,y2, 1);
                    X = (X+xSize-1)%xSize;
                    break;
                default:
                    assert(0);
                }

                //if we were at an octagon it shure must have corrected something
                // (otherwise assert(0);)-->so counter up!
                l++;

            }
            else if ( which_cell_type(X,y2) == sqr )
            {

                switch (v_dir)
                {
                case up:
                    X = (X+xSize-1)%xSize;
                    break;
                case down:
                    X = (X+1)%xSize;
                    break;
                case neither:
                    assert(0);
                default:
                    assert(0);
                }


            }
            else
            {
                assert(0);
            }// end if octo or square


            // this makes the loop run a last time after X hits x2
            // No, I do not want to talk about it!
            if (exitnextloop)
            {
                exit = true;
            }
            if (X == x2)
            {
                exitnextloop = true;
            }




        } //endofwhileY
    }//endofifdeltax


}// end of function I guess...




// Determines if we need to go up, down or neither in the lattice
vertical_direction Lattice::v_where_to_go (int x1, int x2)
{

    vertical_direction v_dir = wrong_direction;

    if ( (x1 < x2 && x2 - x1 <  xSize - x2 + x1) ||
            (x2 < x1 && x1 - x2 >  xSize - x1 + x2) )
        v_dir = down;
    else if (x1 == x2)
        v_dir = neither;
    else
        v_dir = up;

    return v_dir;
}




// Determines if the celltype is square or an octagon
celltype Lattice::which_cell_type ( int x, int y)
{

    int i = x;
    int j = y;

    celltype type = wrong_type;

    if ( (i%2 && j%2)  || (!(i%2) && !(j%2)) )
    {
        type = octo;
    }
    else if ( (!(i%2) && j%2)  || (i%2 && !(j%2)) )
    {
        type = sqr;
    }
    else
    {
        type = wrong_type;
    }

    return type;

}//end of which_cell_type



//determines distance between coordinates in x direction
int Lattice::delta_xy(int xy1, int xy2, int Size)
{

    int d = min (abs(xy2-xy1), Size-abs(xy2-xy1) );

    return d;

}//endof delta_x



//returns true if there is no error syndrome anywhere on the lattice

bool Lattice::is_line()
{


    bool is_line_var = false;
    int syndromecounter = 0;


    // Identify vertices that represent squares/octagons with non-trivial syndrome


    for (int i = 0; i < xSize; i++)
    {
        for (int j = 0; j < ySize; j++)
        {
            bool syndrome = false;
// if i,j both odd do OCT-Z-syndrome (detects X Errors)

            if ( ((i%2) && (j%2)) )
            {
                // syndrome = OCTsyndromeZZZZZZZZ(i,j);
            }
//if exactly  one of i and j is odd and one is even do SQ-Z-Syndrome

            if (  (i%2 && !(j%2))   ||   (!(i%2) && j%2) )
            {
                //syndrome = SQsyndromeZZZZ(i,j);
            }


// if one of the above has a syndrome
            if (syndrome)
            {
                syndromecounter++;
            }
//otherwise go further
        }
    }

    if (syndromecounter == 2)
    {
        is_line_var = true;
    }


    return is_line_var;


}




//generates exactly 1 line of errors on the lattice

void Lattice::generateline()
{
    int x1 =  -1;
    int y1 = -1;

    int x2 =  -1;
    int y2 = -1;




    while(true)
    {

        x1 =  rand()%xSize;	//element [0,xSize]
        y1 = rand()%ySize;	//element [0,ySize]

        x2 =  rand()%xSize;	//element [0,xSize]
        y2 = rand()%ySize;	//element [0,ySize]

        //no syndrome at X Octagons erlaubt
        if ( (!(x1%2) && !(y1%2)) || (!(x2%2) && !(y2%2)))
        {
            continue;
        }



        //if we have generated two different points stop the generator
        if ((x1 != x2 || y1 != y2)    )
        {
            break;
        }

    }




    correctLine( x1, y1, x2, y2);


}//end of generateline



//returns true if there is no error syndrome anywhere on the lattice

bool Lattice::is_corrected()
{



    // Identify vertices that represent squares/octagons with non-trivial syndrome


    for (int i = 0; i < xSize; i++)
    {
        for (int j = 0; j < ySize; j++)
        {
            bool syndrome = false;
// if i,j both odd do OCT-Z-syndrome (detects X Errors)

            if ( ((i%2) && (j%2)) )
            {
                syndrome = OCTsyndromeZZZZZZZZ(i,j,0);
            }
//if exactly  one of i and j is odd and one is even do SQ-Z-Syndrome

            if (  (i%2 && !(j%2))   ||   (!(i%2) && j%2) )
            {
                syndrome = SQsyndromeZZZZ(i,j,0);
            }


// if one of the above has a syndrome
            if (syndrome)
            {
                return false;
            }
//otherwise go further
        }
    }



//if there was no syndrome anywhere
    return true;


}







// Generates logical operator  Z1, or Z2 on this lattice which do anticommute with nontrivil Xloops
Operator* Lattice::getLogical(string whichOp)
{

    //number of qubits
    Operator *myOp = new Operator (xSize*ySize/2*4);

    for (int i = 0, l = 0; i < xSize; i++)
    {
        for (int j = 0; j < ySize; j++)
        {
            for (int k = 0; k < zSize; k++)
            {

                if ( which_cell_type(i,j) == octo || (k == 0) )
                {
                    if(debug8)
                    {
                        cout<<"continue"<<endl;
                    }
                    continue;

                }

                if (whichOp == "Z1")
                {
                    //first line squares
                    if (((i == 0) && which_cell_type(i,j) == sqr && (k == 1)) || ((i == 0) && which_cell_type(i,j) == sqr && (k == 2)) )
                    {
                        if(debug8)
                        {
                            cout<<"Z"<<endl;
                        }
                        myOp->ops[l] = Z;
                    }
                    else
                    {
                        if(debug8)
                        {
                            cout<<"I"<<endl;
                        }
                        myOp->ops[l] = I;
                    }
                }
                else if (whichOp == "Z2")
                {
                    //first column squares
                    if ((j == 0 && which_cell_type(i,j) == sqr && k == 1) || (j == 0 && which_cell_type(i,j) == sqr && k == 3))
                    {
                        if(debug8)
                        {
                            cout<<"Z"<<endl;
                        }
                        myOp->ops[l] = Z;
                    }
                    else
                    {
                        if(debug8)
                        {
                            cout<<"I"<<endl;
                        }
                        myOp->ops[l] = I;
                    }
                }
                else
                {
                    assert (0);
                }
                l++;

            }//endofk
        }//endofj
    }//endofi

    return myOp;
}


// Examimes ErrInit and ErrGuess to determine if error correction succeeds.
// Multiplies ErrGuess with operator p.
bool Lattice::success(void)
{

    // Calculate O = E * E_guessed
    Operator E(0), EGuessP(0);
    for (int i = 0; i < xSize; i++)
    {
        for (int j = 0; j < ySize; j++)
        {
            for (int k = 0; k < zSize; k++)
            {
                if ((k == 0) || (i%2 && j%2) || (!(i%2) && !(j%2)) )
                {
                    continue;
                }

                E.pushBack (cells[i][j][k]->ErrInit->ops[0]);
                EGuessP.pushBack (cells[i][j][k]->ErrGuess->ops[0]);
                //ECurrent.pushBack (cells[i][j]->ErrCurrent->ops[0]);
            }
        }
    }
    Operator O = E * EGuessP;

    /*
    	Debug

      Operator  ECurrent(0);
      //O and ECurrent should be the same
    	cout<<endl<<"E: ";
    	E.printState();
    	cout<<endl<<"G: ";
    	EGuessP.printState();
    	cout<<endl<<"O: ";
    	O.printState();
    	cout<<endl<<"C: ";
    	ECurrent.printState();
    	cout<<endl<<endl;
    */

    // Generate logical operators  Z1, Z2

    Operator *Z1 = getLogical("Z1");
    Operator *Z2 = getLogical("Z2");

    assert ( Z1->commute(*Z2));
    if (debug7)
    {
        cout<<"E0: ";
        E.printState();
        cout<<endl;
        cout<<"EG: ";
        EGuessP.printState();
        cout<<endl;
        cout<<"Ob: ";
        O.printState();
        cout<<endl;
        cout<<"Z1: ";
        Z1->printState();
        cout<<endl;

    }

    // Determine result of error correction
    bool result = true;
//if ( !O.commute(*Z1) || !O.commute(*Z2) ) {
    if ( !O.commute(*Z2)  )
    {
        result = false;
    }

    delete Z1;
    delete Z2;

    return result;


}


// Print the current 2Dstate of the lattice. The content of ErrCurrent and ZSyndrome for all qubits in i==1 is printed.
void Lattice::print2DState(void)
{

    for (int i = 0; i < xSize; i++)
    {
        for (int j = 0; j < ySize; j++)
        {




            if (i==1 && !(j%2))  //j=1,3,5is syndrome,0,2,4 is qubit
            {
                cells[i][j][1]->ErrCurrent->printState();
                if (SQsyndromeZZZZ (i, j, 0))
                {
                    cout<<"-";
                }
                else
                {
                    cout<<"+";
                }
                cells[i][j][2]->ErrCurrent->printState();

            }
            else if (i==1 && j%2)
            {
                if (OCTsyndromeZZZZZZZZ (i, j, 0))
                {
                    cout<<" - ";
                }
                else
                {
                    cout<<" + ";
                }
            }

        }//endofj
    }//endofi
    cout<<endl;
}






// Print the current state of the lattice. The content of ErrCurrent for all qubits is printed.
void Lattice::printState(void)
{

    for (int i = 0; i < xSize; i++)
    {
        for (int j = 0; j < ySize; j++)
        {
            for (int k = 0; k < zSize; k++)
            {
                if ((k == 0) || (i%2 && j%2) || (!(i%2) && !(j%2)))
                    continue;
                cout << "Error at location X=" << i << " Y=" << j << " Z=" << k <<": ";
                cells[i][j][k]->ErrCurrent->printState();
                cout<< endl;
            }
        }
    }
}

// Print the current state of the QB lattice. The content of ErrCurrent for all qubits is printed.
void Lattice::printStateQB(void)
{
    cout<<"Print Lattice State"<<endl;

    for (int i = 0; i < xSize; i++)
    {
        for (int l = 0; l < 3; l++)
        {
            for (int j = 0; j < ySize; j++)
            {

                if ( (i%2 && j%2) || (!(i%2) && !(j%2)))
                {

                    if (l==1)
                    {
                        if(i%2 && j%2)
                        {
                            cout<< "    ";
                            // if there is syndrome print -
                            //if(OCTsyndromeZZZZZZZZ(i,j)){cout<<"-";} else {cout<<"+";}
                            cout<<"    ";
                            continue;
                        }
                        else
                        {
                            cout<< "         ";
                            continue;
                        }

                    }
                    else
                    {
                        cout<< "         ";
                        continue;
                    }
                }



                if (l==0)
                {
                    cout<<" ";
                    cells[i][j][1]->ErrCurrent->printState();
                    cout<<" --- ";
                    cells[i][j][2]->ErrCurrent->printState();
                    cout<<" ";


                }
                else if (l==1)
                {
                    cout<<"  | ";
                    //if there is syndrome print -
                    //if(SQsyndromeZZZZ(i,j)){cout<<"-";} else {cout<<"+";}
                    cout <<" |  ";
                }
                else if (l==2)
                {
                    cout<<" ";
                    cells[i][j][3]->ErrCurrent->printState();
                    cout<<" --- ";
                    cells[i][j][4]->ErrCurrent->printState();
                    cout<<" ";
                }










            }//end of j
            cout<<endl;
        }// end of l
        cout<<endl;



    }//end of i
}//end of print






//end of file lattice.cc
