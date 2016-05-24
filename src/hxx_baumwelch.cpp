#include "hxx.hpp"

#include <vector>
using std::vector;

void
hxx_baumwelch (const vector<int>& O,
               const vector<double>& Eta,
               const vector<double>& Gamma,
               hxx_matrices& biglambda)
{
    auto a = [&biglambda](int i, int j) ->double& { return biglambda.a (i, j); };
    auto b = [&biglambda](int i, int j) ->double& { return biglambda.b (i, j); };
    auto p = [&biglambda](int i) ->double& { return biglambda.p (i); };

    int N = biglambda.N ();
    
    auto eta = [N, &Eta](int t, int i, int j) { return Eta[N*(N*t + i) + j]; };
    auto gamma = [N, &Gamma](int t, int i) { return Gamma[N*t + i]; };

    int T = O.size ();

    if (0 == T) {
        return;
    }

    for (int i = 0; i < N; ++i) {
        p (i) = gamma (0, i);
    }

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            double aij_numer = 0.;
            for (int t = 0; t < T - 1; ++t) {
                aij_numer += eta (t, i, j);
            }

            double aij_denom = 0.;
            for (int t = 0; t < T - 1; ++t) {
                aij_denom += gamma (t, i);
            }

            a (i, j) = aij_numer / aij_denom;
        }
    }

    for (int j = 0; j < N; ++j) {
        for (int k = 0; k < M; ++k) {
            double bik_numer = 0.;
            for (int t = 0; t < T; ++t) {
                if (O[t] == k) {
                    bik_numer += gamma (t, j);
                }
            }

            double bik_denom = 0.;
            for (int t = 0; t < T; ++t) {
                bik_denom += gamma (t, j);
            }

            b (j, k) = bik_numer / bik_denom;
        }
    }
}
