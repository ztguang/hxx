#include "hxx.hpp"

#include <random>
using std::random_device;
using std::ranlux48;
using std::uniform_real_distribution;

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

double
hxx_matrices::a (int i, int j) const
{
    return A[N_*i + j];
}

double
hxx_matrices::b (int i, int k) const
{
    return B[M_*i + k];
}

double
hxx_matrices::p (int i) const
{
    return Pi[i];
}

double&
hxx_matrices::a (int i, int j)
{
    return A[N_*i + j];
}

double&
hxx_matrices::b (int i, int k)
{
    return B[M_*i + k];
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

void
hxx_matrices::randomize ()
{
    random_device rd;
    ranlux48 eng (rd());
    uniform_real_distribution<> urdst (0., 1.);

    for (int i = 0; i < N_; ++i) {
        double Arow_denom = 0.;

        for (int j = 0; j < N_; ++j) {
            Arow_denom += A[N_*i + j] = urdst (eng);
        }

        for (int j = 0; j < N_; ++j) {
            A[N_*i + j] /= Arow_denom;
        }
    }

    for (int i = 0; i < N_; ++i) {
        double Brow_denom = 0.;
        
        for (int k = 0; k < M_; ++k) {
            Brow_denom += B[M_*i + k] = urdst (eng);
        }
        
        for (int k = 0; k < M_; ++k) {
            B[M_*i + k] /= Brow_denom;
        }
    }

    double Pi_denom = 0.;

    for (int i = 0; i < N_; ++i) {
        Pi_denom += Pi[i] = urdst (eng);
    }
    
    for (int i = 0; i < N_; ++i) {
        Pi[i] /= Pi_denom;
    }
}
