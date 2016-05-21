#include "hxx.hpp"

#include <iostream>
#include <vector>
using std::cout;
using std::vector;

int main (int argc, char* argv[])
{
    auto A = {.5, .5,
              .5, .5};

    //auto B = {.75, .25,
    //          .25, .75};

    auto B = {.99, .01,
              .5, .5};

    auto Pi = {.1, .9};

    hxx_gen hg (A, B, Pi);

    pair<int, int> tmp;
    int T = 100;
    vector<int> Q (T);
    vector<int> O (T);

    for (int t = 0; t < T; ++t) {
        tmp = hg ();
        Q[t] = tmp.first;
        O[t] = tmp.second;
    }

    vector<int> q_star;

    hxx_viterbi (O.cbegin (), O.cend (), A, B, Pi, q_star);

    for (int t = 4; t < T; ++t) {
        double likelihood = hxx_forward (O.cbegin () + t - 4, O.cbegin () + t + 1, A, B, Pi);
        cout << Q[t] << " " << O[t] << " " << likelihood << " " << q_star[t] << "\n";
    }

    int num_diff = 0;

    for (int t = 0; t < T; ++t) {
        if (Q[t] != q_star[t]) {
            ++num_diff;
        }
    }
    cout << "Generalized error: " << static_cast<double>(num_diff)/T << "\n";
}
