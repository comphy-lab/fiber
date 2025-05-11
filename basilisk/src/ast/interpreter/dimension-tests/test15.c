#include "layered/hydro.h"
double G0 = 9.81 [1,-2], H0 = 10. [1];
double c = G*H0;

/**
 * @brief Initializes simulation parameters and starts the layered hydrodynamics run.
 *
 * Sets the simulation domain size and gravitational acceleration, then executes the main simulation loop.
 * @return int Exit status code.
 */
int main()
{
  L0 = 1. [1];
  G = G0;
  run();
}

event init (i = 0)
{
  foreach()
    u.x[] = L0; // fixme: should fail here
}

event logfile (t = 10 [0,1]);
