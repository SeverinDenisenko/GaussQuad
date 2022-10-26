# Gaussian quadrature evaluator

To find Gaussian quadrature parameters:

```
make gauss_solver/Makefile build
./gauss_solver/solver {n}
```

To integrate:
```
make gauss_integrator/Makefile build
./gauss_integrator/integrator {n} {a} {b}
```

Where:

* n - number of nodes
* a - lower limit
* b - higher limit
