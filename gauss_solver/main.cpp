#include <iostream>
#include "gauss.h"

int main()
{
    std::vector<double> test = {10, 0, 3, -1, 0};

    auto [res, residual] = divide_polynomial_by(test, 3);

    std::for_each(res.begin(), res.end(), [&](const auto &item)
    {
        std::cout << item << " ";
    });

    std::cout << std::endl << residual << std::endl;

    return 0;
}
