//
// Created by Severin on 27.10.2022.
//

#include "gauss.h"

std::tuple<std::vector<double>, double> divide_polynomial_by(const std::vector<double>& p, double c)
{
    uint32_t n = p.size();
    std::vector<double> result(n - 1);
    
    result[0] = p[0];

    for (uint32_t i = 1; i < n - 1; ++i)
    {
        result[i] = p[i] + c * result[i - 1];
    }

    double residual = p[n - 1] + c * result[n - 2];

    return {result, residual};
}

