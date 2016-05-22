#include "hxx.hpp"

#include <cassert>
#include <cmath>
using std::nan;

double
hxx_forwardbackward (const vector<int>& O,
                     hxx_matrices& biglambda,
                     vector<double>& Eta,
                     vector<double>& Gamma)
{
    auto a = [&biglambda](int i, int j) { return biglambda.a (i, j); };
    auto b = [&biglambda](int i, int j) { return biglambda.b (i, j); };
    auto p = [&biglambda](int i) { return biglambda.p (i); };
    int N = biglambda.N ();
    int M = biglambda.M ();

    int T = O.size ();

    if (0 == T) {
        return nan ("");
    }

    // Forward algorithm.
    vector<double> alpha_bar (N);
    vector<double> Alpha_Hat (T * N);
    auto alpha_hat = [N, &Alpha_Hat](int t, int i) ->double& { return Alpha_Hat[N*t + i]; };
    double alpha_bar_sum = 0.;
    vector<double> c (T);
    double bigCT = 1.;

    for (int i = 0; i < N; ++i) {
        alpha_bar[i] = p (i) * b (i, O[0]);
        alpha_bar_sum += alpha_bar[i];
    }

    assert (alpha_bar_sum > 0.);

    c[0] = 1. / alpha_bar_sum;
    for (int i = 0; i < N; ++i) {
        alpha_hat(0, i) = c[0] * alpha_bar[i];
    }
    bigCT = c[0];

    for (int t = 1; t < T; ++t) {
        for (int j = 0; j < N; ++j) {
            alpha_bar[j] = 0.;
            
            for (int i = 0; i < N; ++i) {
                alpha_bar[j] += alpha_hat (t - 1, i) * a (i, j) * b (j, O[t]);
            }
        }

        alpha_bar_sum = 0.;

        for (int i = 0; i < N; ++i) {
            alpha_bar_sum += alpha_bar[i];
        }

        assert (alpha_bar_sum > 0.);
        c[t] = 1. / alpha_bar_sum;
        bigCT *= c[t];

        for (int i = 0; i < N; ++i) {
            alpha_hat (t, i) = c[t] * alpha_bar[i];
        }
    }

    // Backward algorithm.
    vector<double> beta_bar (N);
    vector<double> Beta_Hat (T * N, 0.);
    auto beta_hat = [N, &Beta_Hat](int t, int i) ->double& { return Beta_Hat[N*t + i]; };

    for (int i = 0; i < N; ++i) {
        beta_hat (T - 1, i) = c[T - 1];
    }

    for (int t = T - 2; t > -1; --t) {
        for (int j = 0; j < N; ++j) {
            beta_bar[j] = 0.;
            for (int i = 0; i < N; ++i) {
                beta_bar[j] += a (i, j) * b (j, O[t + 1]) * beta_hat (t + 1, i);
            }
        }

        for (int i = 0; i < N; ++i) {
            beta_hat (t, i) = c[t] * beta_bar[i];
        }
    }

    Eta.resize (T * N * N);
    Gamma.resize (T * N);

    auto eta = [N, &Eta](int t, int i, int j) ->double& { return Eta[N*(N*t + i) + j]; };
    auto gamma = [N, &Gamma](int t, int i) ->double& { return Gamma[N*t + i]; };

    for (int t = 0; t < T - 1; ++t) {
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                eta (t, i, j) = alpha_hat (t, i) * a (i, j) * b (j, O[t + 1]) * beta_hat (t + 1, j);
            }
        }
    }
    
    for (int t = 0; t < T; ++t) {
        for (int i = 0; i < N; ++i) {
            gamma (t, i) = alpha_hat (t, i) * beta_hat (t, i) / c[t];
        }
    }

    return bigCT;
}
