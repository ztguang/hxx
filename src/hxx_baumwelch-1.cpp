#include "hxx.hpp"

#include <algorithm>
#include <iostream>
#include <vector>
using std::cout;
using std::random_shuffle;
using std::vector;

void output_biglambda (const hxx_matrices&);

int main (int argc, char* argv[])
{
    auto A = {.5, .5,
              .5, .5};

    auto B = {.75, .25,
              .25, .75};

    //auto B = {.99, .01,
    //          1., .0};

    auto Pi = {.1, .9};

    hxx_matrices biglambda (A, B, Pi);

    hxx_gen hg (A, B, Pi);

    pair<int, int> tmp;
    int T = 100000;
    vector<int> Q (T);
    vector<int> O (T);

    for (int t = 0; t < T; ++t) {
        tmp = hg ();
        Q[t] = tmp.first;
        O[t] = tmp.second;
    }

    vector<double> Xi, Gamma;

    //random_shuffle (O.begin (), O.end ());

    hxx_forwardbackward (O, biglambda, Xi, Gamma);

    hxx_matrices newbiglambda (A, B, Pi);

    hxx_baumwelch (O, Xi, Gamma, newbiglambda);

    for (int t = 4; t < T; ++t) {
        double likelihood = hxx_forward (O.cbegin () + t - 4, O.cbegin () + t + 1, A, B, Pi);
        cout << Q[t] << " " << O[t] << " " << likelihood << "\n";
    }

    cout << "Original biglambda\n";
    output_biglambda (biglambda);
    cout << "New biglambda\n";
    output_biglambda (newbiglambda);
}

void output_biglambda (const hxx_matrices& biglambda) {
    int N = biglambda.N ();
    int M = biglambda.M ();

    cout << "A = {";
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            cout << " " << biglambda.a (i, j);
        }
    }
    cout << " }\n";

    cout << "B = {";
    for (int i = 0; i < N; ++i) {
        for (int k = 0; k < M; ++k) {
            cout << " " << biglambda.b (i, k);
        }
    }
    cout << " }\n";

    cout << "Pi = {";
    for (int i = 0; i < N; ++i) {
        cout << " " << biglambda.p (i);
    }
    cout << " }\n";
}
