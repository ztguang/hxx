/*
[root@bupt src]# make
[root@bupt src]# ./hxx_baumwelch-1 0.5 0.5 0.5 0.5 0.75 0.35 0.15 0.75 0.1 0.9
*/

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
//ztg add
    double a1 = atof(argv[1]), a2 = atof(argv[2]), a3 = atof(argv[3]), a4 = atof(argv[4]);
    double b1 = atof(argv[5]), b2 = atof(argv[6]), b3 = atof(argv[7]), b4 = atof(argv[8]);
    double p1 = atof(argv[9]), p2 = atof(argv[10]);

//ztg alter
//-----------------------------------
/*
    auto A = {.5, .5,
              .5, .5};

    auto B = {.75, .25,
              .25, .75};

    //auto B = {.99, .01,
    //          1., .0};

    auto Pi = {.1, .9};

//*/
    auto A = {a1, a2,
              a3, a4};

    auto B = {b1, b2,
              b3, b4};

    auto Pi = {p1, p2};
//-----------------------------------


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
