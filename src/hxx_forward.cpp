#include "hxx.hpp"

#include <algorithm>
#include <cmath>
#include <vector>
using std::accumulate;
using std::nan;
using std::transform;
using std::vector;

double hxx_forward (vector<int>::const_iterator first,
                    vector<int>::const_iterator last,
                    initializer_list<double> A_in,
                    initializer_list<double> B_in,
                    initializer_list<double> Pi_in)
{
    const vector<double> A(A_in);
    const vector<double> B(B_in);
    const vector<double> Pi(Pi_in);

    const int N = Pi.size ();
    const int M = B.size () / N;

    auto a = [N, &A](int i, int j) { return A[N*i + j]; };
    auto b = [M, &B](int i, int k) { return B[M*i + k]; };
    auto p = [&Pi](int i) { return Pi[i]; };

    vector<double> alpha_bar (N);
    vector<double> alpha_hat (N);
    double alpha_bar_sum;
    double c;
    double bigCT = 1.;

    if (first == last) {
        return nan("");
    }

    for (int i = 0; i < N; ++i) {
        alpha_bar[i] = p(i) * b(i, *first);
    }

    alpha_bar_sum = accumulate (alpha_bar.begin (), alpha_bar.end (), 0.);
    if (alpha_bar_sum > 0.) {
        c = 1. / alpha_bar_sum;
        bigCT *= c;
        transform (alpha_bar.begin (), alpha_bar.end (), alpha_hat.begin (), [c](double v) { return v * c; });
    } else {
        alpha_hat = alpha_bar;
    }
    ++first;

    for (; first != last; ++first) {
        for (int j = 0; j < N; ++j) {
            alpha_bar[j] = 0.;

            for (int i = 0; i < N; ++i) {
                alpha_bar[j] += alpha_hat[i] * a(i, j) * b(j, *first);
            }
        }

        alpha_bar_sum = accumulate (alpha_bar.begin (), alpha_bar.end (), 0.);
        if (alpha_bar_sum > 0.) {
            c = 1. / alpha_bar_sum;
            bigCT *= c;
            transform (alpha_bar.begin (), alpha_bar.end (), alpha_hat.begin (), [c](double v) { return v * c; });
        } else {
            alpha_hat = alpha_bar;
        }
    }

    return 1. / bigCT;
    //return accumulate (alpha_hat.begin (), alpha_hat.end (), 0.) / bigCT;
}

