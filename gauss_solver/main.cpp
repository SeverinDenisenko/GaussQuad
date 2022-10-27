#include <iostream>
#include "gauss.h"

int main()
{

    auto res = get_legendre_polynomial(5);

    std::for_each(res.begin(), res.end(), [&](const auto &item)
    {
        std::cout << item << " ";
    });

    std::cout << std::endl;

    return 0;
}
