/**
Checks that minmod2() is always undefined. */

#include "utils.h"

/**
 * @brief Initializes a grid, sets up a scalar field, and tests calling an undefined function.
 *
 * Initializes a computational grid with one level, declares a scalar field, and sets the interpreter verbosity. Iterates over all grid cells, and if the previous cell's scalar value is nonzero, attempts to call the undefined function `minmod2` with neighboring scalar values. This is intended to verify that `minmod2()` is always undefined in the interpreter context.
 *
 * @return int Exit status code.
 */
int main()
{
  init_grid (1);
  scalar s[];
  {
    interpreter_verbosity (2);
    foreach() {
      if (s[-1]) {
	double a = minmod2 (s[-1], s[], s[1]);
	if (a);
      }
    }
  }
}
