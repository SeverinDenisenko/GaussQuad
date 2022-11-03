//
// Created by Severin on 27.10.2022.
//

#include "gauss.h"

#include <random>
#include <algorithm>
#include <iostream>

std::tuple<std::vector<double>, double> divide_polynomial_by(const std::vector<double> &p, double c)
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
    while (_p.size() != 2)
    {
        std::vector<double> difference_equation(_p.size());

        // fill with random values
        std::random_device rd;
        std::default_random_engine eng(rd());
        std::uniform_real_distribution<double> randomizer(-1.0, 1.0);

        for (auto &item: difference_equation)
        {
            item = randomizer(eng);
        }

        double root = difference_equation.at(difference_equation.size() - 1) /
                      difference_equation.at(difference_equation.size() - 2);
        double root1 = 0;

        // iterate throw the difference equation
        while (fabs(root - root1) > 10e-10)
        {
            // find next y
            double y = 0;
            for (uint32_t i = 1; i < difference_equation.size(); ++i)
            {
                y -= _p.at(i) * difference_equation.at(difference_equation.size() - i);
            }
            y /= _p.at(0);

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
    roots.push_back(-_p.at(1) / _p.at(0));

    return roots;
}

std::vector<double> calculate_roots_gauss(uint32_t n)
{
    std::vector<double> roots;

    std::vector<double> legendre = get_legendre_polynomial(n);

    if (n % 2 == 1)
    {
        // For odd order Legendre polynomial (a_n = 0) can be divided by x and 0.0 root excluded
        legendre.pop_back(); // divide by x

        roots.push_back(0.0);
    }

    // For even order Legendre polynomial can be made substitution y = x^2
    std::vector<double> mod_legendre;

    for (uint32_t i = 0; i < legendre.size(); i += 2)
    {
        mod_legendre.push_back(legendre.at(i));
    }

    std::vector<double> square_roots = solve_polynomial_bernoulli(mod_legendre);

    // Back substitution
    for (auto item: square_roots)
    {
        roots.push_back(-sqrt(item));
        roots.push_back(sqrt(item));
    }

    return roots;
}

void SolveLeadElement(double **A, std::vector<double> &X);
void SolveLeadElement(double **A, std::vector<double> &X)
{
    int n = X.size();

    int *x_order = (int *) malloc(sizeof(int) * n);

    for (int i = 0; i < n; i++)
    {
        x_order[i] = i;
    }

    for (int j = 0; j < n; j++)
    {
        // Serch for biggest element and change order //

        int max_i = 0;
        int max_j = 0;
        double max = 0;

        for (int z = j; z < n; z++)
        {
            for (int t = j; t < n; t++)
            {
                if (fabs(A[z][t]) > max)
                {
                    max = fabs(A[z][t]);
                    max_i = z;
                    max_j = t;
                }
            }
        }

        //for (int i = 0; i < n; i++)
        //{
        //    std::swap(A[i][max_j], A[i][j]);
        //}
        //std::swap(A[max_i], A[j]);
        //std::swap(x_order[j], x_order[max_j]);

        for (int i = 0; i < n; i++)
        {
            double tmp = A[i][max_j];
            A[i][max_j] = A[i][j];
            A[i][j] = tmp;
        }

        double *tmp = A[j];
        A[j] = A[max_i];
        A[max_i] = tmp;

        int tmp1 = x_order[j];
        x_order[j] = x_order[max_j];
        x_order[max_j] = tmp1;


        ////////////////////////////////////////////////

        for (int i = j + 1; i < n; i++)
        {
            double div = A[i][j] / A[j][j];

            for (int k = 0; k < n + 1; k++)
            {
                A[i][k] = A[i][k] - div * A[j][k];
            }
        }
    }

    X[n - 1] = A[n - 1][n] / A[n - 1][n - 1];

    for (int i = n - 2; i >= 0; i--)
    {
        double sum = 0;
        for (int j = i + 1; j < n; j++)
        {
            sum = sum + A[i][j] * X[j];
        }
        X[i] = (A[i][n] - sum) / A[i][i];
    }

    // Restore the order of X //

    for (int i = 0; i < n; i++) {
        int next = i;

        while (x_order[next] >= 0) {

            //std::swap(X[i], X[x_order[next]]);
            double tmp = X[i];
            X[i] = X[x_order[next]];
            X[x_order[next]] = tmp;

            int temp = x_order[next];

            x_order[next] -= n;
            next = temp;
        }
    }

    free(x_order);

    ////////////////////////////
}

std::vector<double> calculate_coefficients_gauss(const std::vector<double>& roots)
{
    std::vector<double> res(roots.size());

    // Solve linear equation with Full-piv
    uint32_t n = roots.size();

    auto **A = (double **)malloc(sizeof(double *) * roots.size());
    for (uint32_t i = 0; i < n; ++i)
    {
        A[i] = (double *)malloc(sizeof(double) * (n + 1));
    }

    for (uint32_t k = 0; k < n; ++k)
    {
        for (uint32_t i = 0; i < n; ++i)
        {
            A[k][i] = pow(roots.at(i), k);
        }
    }

    for (uint32_t k = 0; k < n; ++k)
    {
        if (k % 2 == 1){
            A[k][n] = 0;
        } else {
            A[k][n] = 2 / ((double)k + 1);
        }
    }

    SolveLeadElement(A, res);

    return res;
}
