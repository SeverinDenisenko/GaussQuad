#include <iostream>
#include <algorithm>
#include <fstream>
#include <numeric>

#include "gauss.h"

int main(int argc, char **argv)
{
    if (argc != 2){
        std::cerr << "Expected one input argument!" << std::endl;
        return 1;
    }

    uint32_t n = std::stoul(argv[1]);

    auto roots = calculate_roots_gauss(n);
    auto coefficients = calculate_coefficients_gauss(roots);

    std::ofstream out("quad" + std::to_string(n) + ".dat");

    for (uint32_t i = 0; i < n; ++i)
    {
        out << coefficients.at(i) << " " << roots.at(i) << std::endl;
    }

    out.close();

    return 0;
}
