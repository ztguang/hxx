#include "hxx.hpp"

#include <cmath>
#include <utility>
using std::nan;
using std::pair;

void hxx_viterbi (vector<int>::const_iterator it,
                  vector<int>::const_iterator last,
                  initializer_list<double> A_in,
                  initializer_list<double> B_in,
                  initializer_list<double> Pi_in,
                  vector<int>& q_star)
{
    vector<double> A (A_in);
    vector<double> B (B_in);
    vector<double> Pi (Pi_in);

    hxx_viterbi (it, last, A, B, Pi, q_star);
}

void hxx_viterbi (vector<int>::const_iterator it,
                  vector<int>::const_iterator last,
                  const vector<double> A,
                  const vector<double> B,
                  const vector<double> Pi,
                  vector<int>& q_star)
{
    const int T = last - it;
    const int N = Pi.size ();
    const int M = B.size () / N;

    auto a = [N, &A](int i, int j) { return A[N*i + j]; };
    auto b = [M, &B](int i, int k) { return B[M*i + k]; };
    auto p = [&Pi](int i) { return Pi[i]; };

    vector<double> Delta (T * N);
    vector<int> Psi (T * N);

    auto delta = [N, &Delta](int t, int i) ->double& { return Delta[N*t + i]; };
    auto psi = [N, &Psi](int t, int i) ->int& { return Psi[N*t + i]; };

    if (it == last) {
        q_star.clear ();
        return;
    }

    for (int i = 0; i < N; ++i) {
        delta(0, i) = p(i) * b(i, *it);
        psi(0, i) = -1; // Unused;
    }
    ++it;

    for (int t = 1; t < T; ++t, ++it) {
        for (int j = 0; j < N; ++j) {
            double max_val = -1.;
            int max_ind = -1;
            for (int i = 0; i < N; ++i) {
                double tmp = delta (t - 1, i) * a (i, j);
                if (tmp > max_val) {
                    max_val = tmp;
                    max_ind = i;
                }
            }
            delta (t, j) = max_val * b (j, *it);
            psi (t, j) = max_ind;
        }
    }

    double P_star = -1.;
    int q_star_T = -1;
    for (int i = 0; i < N; ++i) {
        if (delta (T - 1, i) > P_star) {
            P_star = delta (T - 1, i);
            q_star_T = i;
        }
    }

    q_star.resize (T);
    q_star[T - 1] = q_star_T;

    for (int t = T - 2; t > -1; --t) {
        q_star[t] = psi (t + 1, q_star[t + 1]);
    }
}
