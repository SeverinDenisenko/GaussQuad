//
// Created by Severin on 27.10.2022.
//

#include "gauss.h"

#include <random>
#include <algorithm>
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

std::vector<double> solve_polynomial_bernoulli(const std::vector<double> &p)
{
    std::vector<double> roots;
    std::vector<double> _p = p;

    // iterate throw polynomial, dividing by biggest root
    while (_p.size() != 2){
        std::vector<double> difference_equation(_p.size());

        // fill with random values
        std::random_device rd;
        std::default_random_engine eng(rd());
        std::uniform_real_distribution<double> randomizer(-1.0, 1.0);

        for (auto& item : difference_equation)
        {
            item = randomizer(eng);
        }

        double root = difference_equation.at(difference_equation.size() - 1) /
                      difference_equation.at(difference_equation.size() - 2);
        double root1 = 0;

        // iterate throw the difference equation
        while (fabs(root - root1) > 10e-10){
            // find next y
            double y = 0;
            for (uint32_t i = 1; i < difference_equation.size(); ++i)
            {
                y -= _p.at(i) * difference_equation.at(difference_equation.size() - i);
            }

            // shift difference equation
            std::shift_left(difference_equation.begin(), difference_equation.end(), 1);
            difference_equation.at(difference_equation.size() - 1) = y;

            root1 = root;
            root = difference_equation.at(difference_equation.size() - 1) /
                   difference_equation.at(difference_equation.size() - 2);

        }

        double residual = 0;
        roots.push_back(root);

        std::tie(_p, residual) = divide_polynomial_by(_p, root);

        difference_equation.pop_back();
    }

    // Find the last root
    roots.push_back(- _p.at(1) / _p.at(0));

    return roots;
}
