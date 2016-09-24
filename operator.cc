#include "operator.h"
#include "constants.h"


// Defines products of elementary Pauli operators.
pauli operator*(pauli p, pauli q)
{

    pauli res = I;

    if (p == I)
        res = q;
    else if (q == I)
        res = p;
    else if ((p == Y && q == Z) || ((p == Z && q == Y)))
        res = X;
    else if ((p == X && q == Z) || ((p == Z && q == X)))
        res = Y;
    else if ((p == X && q == Y) || ((p == Y && q == X)))
        res = Z;

    return res;
}


// Constructor generates identity operator of size newSize.
Operator::Operator(int newSize)
{

    size = newSize;

    for (int i = 0; i < size; i ++)
    {
        ops.push_back (I);
    }
}


// Copy constructor. Constructs a new operator that is an exact copy of p.
Operator::Operator (const Operator & p)
{

    size = p.size;

    for (int i = 0; i < size; i++)
    {
        ops.push_back (p.ops[i]);
    }
}


// Push back a Pauli operator to ops.
void Operator::pushBack (pauli & op)
{

    ops.push_back (op);
    size++;
}


// Retunrs true if two operators commute. Otherwise false.
bool Operator::commute (Operator & op)
{

    assert (size == op.size);
    int result = 1;

    for (int i = 0; i < size; i++)
    {
        pauli p1 = ops[i];
        pauli p2 = op.ops[i];

        if (p1 != I && p2 != I && (p1 != p2))
            result *= -1;
    }

    if (result == 1)
        return true;
    else
        return false;
}



// Randomize this operator. With probability p, each operation is X, Y, Z. Otherwise I.
void Operator::generateRandom (double & p)
{

    assert (p < 1.000001);

//gives always the same number with srand()
//srand(time(0));

    for (int i = 0; i < size; i++)
    {

        double randomnumber = (double) rand();
        double r = randomnumber / ((double)RAND_MAX + 1);
        if(debug5)
        {
            cout<<" r: "<<r<<" p: "<<p;
        }
        if (r < p)
        {
            ops[i] = X;
        }
        else
        {
            ops[i] = I;
        }

        if(debug5)
        {
            cout<<" Err: "<<ops[0]<<endl;
        }
    }
}

void Operator::generate_random_evolution(double & p)
{
    assert (p < 1.000001);

    for (int i = 0; i < size; i++)
    {

        double r = (double)rand() / ((double)RAND_MAX + 1);
        if (r < p)
        {
            ops[i] = ops[i]*X;
        }/*else
	//stupid to really multiply it...
	//ops[i] = ops[i]*I;
	continue;
	*/
    }//end of for


}






// Operator * overloaded.
Operator operator*(const Operator & p, const Operator & q)
{

    assert (p.size == q.size);
    Operator res (p.size);

    for (int i = 0; i < p.size; i++)
    {
        res.ops[i] = p.ops[i] * q.ops[i];
    }

    return res;
}


// Operator == overloaded.
bool operator==(const Operator & p, const Operator & q)
{

    if (p.size != q.size)
        return false;

    bool res = true;

    for (int i = 0; i < p.size; i++)
    {
        if (p.ops[i] != q.ops[i])
            res = false;
    }

    return res;
}


// Operator < overloaded. Needed so that operators can be inserted into a set.
bool operator<(const Operator & p, const Operator & q)
{

    if (p.size < q.size)
        return true;

    if (p.size > q.size)
        return false;

    for (int i = 0; i < p.size; i++)
    {
        if (p.ops[i] < q.ops[i])
            return true;
        if (p.ops[i] > q.ops[i])
            return false;
    }

    return false;
}


// Is this operator the identity?
bool Operator::identity (void)
{

    bool res = true;

    for (int i = 0; i < size; i++)
    {
        if (ops[i] != I)
            res = false;
    }

    return res;
}


// Print this operator on screen.
void Operator::printState (void)
{

    for (int i = 0; i < size; i++)
    {
        if (i != 0)
            cout << " ";
        switch(ops[i])
        {
        case I:
            cout << "I";
            break;
        case X:
            cout << "X";
            break;
        case Y:
            cout << "Y";
            break;
        case Z:
            cout << "Z";
            break;
        default:
            assert (0);
        }
    }


}


//end of file operator.cc
