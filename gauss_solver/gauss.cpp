//
// Created by Severin on 27.10.2022.
//

#include "gauss.h"
#include <iostream>

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

std::vector<double> get_legendre_polynomial(uint32_t n)
{
    std::vector<double> res(n + 1);

    switch (n)
    {
        case 0:
            res[0] = 1;
            break;
        case 1:
            res[0] = 1;
            res[1] = 0;
            break;
        default:
            std::vector<double> p_n_1 = get_legendre_polynomial(n - 1);
            std::vector<double> p_n_2 = get_legendre_polynomial(n - 2);

            // effectively multiples by x
            p_n_1.insert(p_n_1.end(), 0);

            p_n_2.insert(p_n_2.begin(), 0);
            p_n_2.insert(p_n_2.begin(), 0);

            for (uint32_t i = 0; i < n + 1; ++i)
            {
                res[i] = (2 * (double) n - 1) / n * p_n_1[i] - (double) (n - 1) / n * p_n_2[i];
            }
            break;
    }

    return res;
}

