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
              1., .0};

    auto Pi = {.5, .5};

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

    for (int t = 0; t < T; ++t) {
        double likelihood = hxx_forward (O.cbegin (), O.cbegin () + t + 1, A, B, Pi);
        cout << Q[t] << " " << O[t] << " " << likelihood << "\n";
    }
}
