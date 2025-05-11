/**
# Dimensions for the Poisson solver 

We solve
$$
\nabla^2\psi = \omega
$$
with $\psi$ the streamfunction (dimensions [2,-1]) and $\omega$ the
vorticity (dimensions [0,-1]). */

#include "poisson.h"

/**
 * @brief Tests the Poisson solver with scalar fields of varying physical dimensions.
 *
 * Initializes a computational grid and sets up streamfunction (`psi`) and vorticity (`omega`) fields with specified dimensions. The Poisson equation is solved twice: first with initial field values and dimensions, then with updated values and dimensions, verifying the solver's handling of physical units.
 *
 * @return int Exit status code.
 */
int main()
{
  init_grid (1);
  size (1. [1]);
  scalar psi[], omega[];
  foreach() {
    psi[] = 0. [2,-1];
    omega[] = 1. [0,-1];
  }
  poisson (psi, omega);

  foreach() {
    psi[] = 0. [2,-2];
    omega[] = 2. [0,-2];
  }
  poisson (psi, omega);
}
