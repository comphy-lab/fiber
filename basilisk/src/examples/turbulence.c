/**
# Decaying two-dimensional turbulence

We solve the two-dimensional incompressible Euler equations using a
vorticity--streamfunction formulation. */

#include "navier-stokes/stream.h"

/**
 * @brief Initializes the simulation grid and starts the two-dimensional turbulence simulation.
 *
 * The domain is a unit square. The grid resolution is set by the first command-line argument, or defaults to 256 if not specified.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line argument strings.
 * @return int Exit status code.
 */

int main (int argc, char * argv[]) {
  init_grid (argc > 1 ? atoi(argv[1]) : 256);
  run();
}

/**
The initial condition for vorticity is just a white noise in the range
$[-1:1]$ .*/

event init (i = 0) {
  double a = 1. [0,-1];
  foreach (cpu)
    omega[] = a*noise();
}

/**
We generate images of the vorticity field every 4 timesteps up to
$t=1000$. We fix the colorscale to $[-0.3:0.3]$.

![Evolution of the vorticity](turbulence/omega.mp4)(autoplay loop) */

event output (t += 1; t <= 1000) {
#if !BENCHMARK
  output_ppm (omega, min = -0.3, max = 0.3, file = "omega.mp4");
#endif
}

/**
On GPUs we have the option to display the vorticity field in real
time. */

#if _GPU && SHOW
event display (i++)
  output_ppm (omega, min = -0.3, max = 0.3, fps = 30, fp = NULL);
#endif

/**
## See also

[Benchmark on GPUs](/src/grid/gpu/Benchmarks.md#two-dimensional-turbulence)
*/
