#include <iostream>
#include <algorithm>
#include <cstdint>
#include "gauss.h"

int main()
{
    auto res = calculate_roots_gauss(11);

    std::for_each(res.begin(), res.end(), [&](const auto &item)
    {
        std::cout << item << " ";
    });

    std::cout << std::endl;

    return 0;
}
