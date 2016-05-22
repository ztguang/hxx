#include "hxx.hpp"

hxx_matrices::hxx_matrices (int N_in, int M_in)
: N_ (N_in), M_ (M_in), A (N_in * N_in), B (N_in * M_in), Pi (N_in)
{
}

hxx_matrices::hxx_matrices (initializer_list<double> A_in,
                            initializer_list<double> B_in,
                            initializer_list<double> Pi_in)
: A (A_in), B (B_in), Pi (Pi_in), N_ (Pi_in.size ()), M_ (B_in.size ()/Pi_in.size ())
{
}

double&
hxx_matrices::a (int i, int j)
{
    return A[N_*i + j];
}

double&
hxx_matrices::b (int i, int k)
{
    return B[N_*i + k];
}

double&
hxx_matrices::p (int i)
{
    return Pi[i];
}

int
hxx_matrices::N () const
{
    return N_;
}

int hxx_matrices::M () const
{
    return M_;
}

