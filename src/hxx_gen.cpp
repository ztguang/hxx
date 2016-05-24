#include "hxx.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
using std::abs;
using std::cerr;
using std::partial_sum;
using std::random_device;
using std::lower_bound;

hxx_gen::hxx_gen (initializer_list<double> A,
                  initializer_list<double> B,
                  initializer_list<double> Pi)
: reng (rd()), urd (0., 1.)
 {
    a.assign (A);
    b.assign (B);
    p.assign (Pi);

    N = Pi.size ();
    M = B.size () / N;

    partial_sum (p.begin (), p.end (), p.begin ());
    if (abs (p.back () - 1.) > 1.e-6) {
        cerr << "Error: hxx_gen.cpp: hxx_gen: cumulative probability does not sum to 1.: Pi\n";
    }
    p.back () = 1.;

    for (int i = 0; i < N; ++i) {
        partial_sum (a.begin () + N*i, a.begin () + N*(i+1), a.begin () + N*i);
        if (abs (a[N*(i+1) - 1] - 1.) > 1.e-6) {
            cerr << "Error: hxx_gen.cpp: hxx_gen: cumulative probability does not sum to 1.: A\n";
        }
        a[N*(i+1) - 1] = 1.;
    }
    for (int i = 0; i < N; ++i) {
        partial_sum (b.begin () + M*i, b.begin () + M*(i+1), b.begin () + M*i);
        if (abs (b[M*(i+1) - 1] - 1.) > 1.e-6) {
            cerr << "Error: hxx_gen.cpp: hxx_gen: cumulative probability does not sum to 1.: B\n";
        }
        b[M*(i+1) - 1] = 1.;
    }

    t = -1;
    q = -1;
}

pair<int, int>
hxx_gen::operator ()()
{
    int Ot;

    if (-1 != t) {
        q = lower_bound (a.begin () + N*q, a.begin () + N*(q+1), urd (reng)) - a.begin () - N*q;

        ++t;

        Ot = lower_bound (b.begin () + M*q, b.begin () + M*(q+1), urd (reng)) - b.begin () - M*q;
    } else {
        q = lower_bound (p.begin (), p.end (), urd (reng)) - p.begin ();
        t = 0;

        Ot = lower_bound (b.begin () + M*q, b.begin () + M*(q+1), urd (reng)) - b.begin () - M*q;
    }

    return pair<int, int>(q, Ot);
}
