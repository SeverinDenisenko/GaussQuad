#include <iostream>
#include <algorithm>
#include <cstdint>
#include "gauss.h"

int main()
{
    auto res = solve_polynomial_bernoulli({1, -12, 12, 80});

    std::for_each(res.begin(), res.end(), [&](const auto &item)
    {
        std::cout << item << " ";
    });

    std::cout << std::endl;

    return 0;
}
