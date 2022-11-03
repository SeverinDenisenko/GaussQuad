#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

double func(double x);
double func(double x){
    return exp(x);
}

double integrate(double (*func)(double), uint32_t n, double a, double b);
double integrate(double (*func)(double), uint32_t n, double a, double b){

    std::ifstream method("../gauss_solver/quad" + std::to_string(n) + ".dat");

    if (!method.good()){
        std::cerr << "Can,t find file!" << std::endl;
        exit(1);
    }

    std::vector<double> roots(n);
    std::vector<double> coefficients(n);

    for (uint32_t i = 0; i < n; ++i)
    {
        method >> coefficients[i] >> roots[i];
    }

    method.close();

    double ans = 0;

    for (uint32_t i = 0; i < n; ++i)
    {
        ans += func((roots.at(i) + 1) * (b - a) / 2 + a) * coefficients.at(i);
    }

    ans *= (b - a) / 2;

    return ans;
}

int main(int argc, char **argv)
{
    if (argc != 4){
        std::cerr << "Expected three input arguments!" << std::endl;
        return 1;
    }

    uint32_t n = std::stoul(argv[1]);
    double a = std::stod(argv[2]);
    double b = std::stod(argv[3]);

    double res = integrate(func, n, a, b);

    std::cout << res << std::endl;

    return 0;
}
