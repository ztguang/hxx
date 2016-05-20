#include "hxx.hpp"

#include <iostream>
using std::cout;

int main (int argc, char* argv[])
{
    auto A = {.5, .5,
              .5, .5};

    auto B = {.75, .25,
              .25, .75};

    auto Pi = {.5, .5};

    hxx_gen hg (A, B, Pi);

    pair<int, int> tmp;

    for (int i = 0; i < 100; ++i) {
        tmp = hg ();
        cout << tmp.first << " " << tmp.second << "\n";
    }
}
