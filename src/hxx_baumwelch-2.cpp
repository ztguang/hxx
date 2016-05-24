#include "hxx.hpp"

#include <algorithm>
#include <iostream>
#include <random>
#include <utility>
#include <vector>
using std::cin;
using std::cout;
using std::nth_element;
using std::ranlux48;
using std::swap;
using std::uniform_int_distribution;
using std::vector;

void output_lambda (const hxx_matrices&);

double median (double* p, size_t n)
{
    nth_element (p, p + n/2, p + n);

    if (1 & n) {
        return p[n/2];
    } else {
        double b = p[n/2];
        nth_element (p, p + n/2 - 1, p + n);
        double a = p[n/2 - 1];
        return (a + b)/2.;
    }
}

int main (int argc, char *argv[])
{
    vector<double> raw;

    double d_in;

    while (cin >> d_in) {
        raw.push_back (d_in);
    }

    for (auto i: raw) {
        cout << i << " ";
    }
    cout << "\n";

    size_t T = raw.size ();

    {
        vector<double> raw_copy (raw);

        double m = median (&raw[0], T);

        cout << "Median = " << m << "\n";
    }

    vector<int> O (T);
    ranlux48 eng;
    uniform_int_distribution<size_t> uidst (0, T - 1);

//    vector<int> Q (T);
    vector<double> Gamma, Xi;

    hxx_matrices lambda1 (2, 2), lambda2 (2, 2);
    hxx_matrices* lambda_in = &lambda1;
    hxx_matrices* lambda_out = &lambda2;

    for (int s = 0; s < 10; ++s) {
        double split = raw[uidst (eng)];

        for (int t = 0; t < T; ++t) {
            if (raw[t] < split) {
                O[t] = 0;
            } else {
                O[t] = 1;
            }
        }

        lambda_in->randomize ();
        output_lambda (*lambda_in);
        cout << " split=" << split;
        cout << "\n";

        for (int i = 0; i < 10; ++i) {
            hxx_forwardbackward (O, *lambda_in, Xi, Gamma);

            hxx_baumwelch (O, Xi, Gamma, *lambda_out);

            output_lambda (*lambda_out);

            cout << "\n";

            swap (lambda_in, lambda_out);
        }
    }
}


void output_lambda (const hxx_matrices& lambda) {
    int N = lambda.N ();
    int M = lambda.M ();

    cout << "A = {";
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            cout << " " << lambda.a (i, j);
        }
    }
    cout << " } ";

    cout << "B = {";
    for (int i = 0; i < N; ++i) {
        for (int k = 0; k < M; ++k) {
            cout << " " << lambda.b (i, k);
        }
    }
    cout << " } ";

    cout << "Pi = {";
    for (int i = 0; i < N; ++i) {
        cout << " " << lambda.p (i);
    }
    cout << " } ";
}

