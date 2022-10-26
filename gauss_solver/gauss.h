//
// Created by Severin on 27.10.2022.
//

#ifndef NUMERICAL_TASK_8_GAUSS_H
#define NUMERICAL_TASK_8_GAUSS_H

#include <vector>
#include <tuple>

/**
 * Divides polynomial p(x) by x - c
 * @param p dividend
 * @param c divider
 * @return quotient and residual
 */
std::tuple<std::vector<double>, double> divide_polynomial_by(const std::vector<double>& p, double c);

/**
 * Solves p(x) = 0
 * @param p
 * @return vector of roots
 */
std::vector<double> solve_polynomial_bernoulli(const std::vector<double>& p);

/**
 * Calculates Legendre polynomial coefficients for power of n
 * @param n power of a Legendre polynomial
 * @return Legendre polynomial
 */
std::vector<double> get_legendre_polynomial(uint32_t n);

/**
 * Calculate optimal points in Gauss integration method
 * @param n number of the points
 * @return vector of the points
 */
std::vector<double> calculate_roots_gauss(uint32_t n);

/**
 * Calculate optimal coefficients in Gauss integration method
 * @param n number of the points
 * @return vector of the coefficients
 */
std::vector<double> calculate_coefficients_gauss(uint32_t n);

#endif //NUMERICAL_TASK_8_GAUSS_H
