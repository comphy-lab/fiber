/**
 * @brief Tests interpreter behavior for reporting unset variable usage in constraints.
 *
 * Declares unset variables and assigns them to another variable, verifying that unset value usage is only reported when the value is subsequently used in a constraint. Initializes a grid, assigns values to a scalar field based on cell coordinates, and applies constraints to test detection logic.
 *
 * @return int Exit status code.
 */

int main()
{
  { // interpreter_verbosity (3);
    double s; // unset
    double val = s; // should not report this
    double s1; // unset
    val = s1; // should report this ...
    val == 1 [0]; // ... since it is used here
    init_grid(1);
    scalar s[];
    foreach() {
      if (x > 12)
	s[] = 33.;
      else
	s[] = 10;
      val = s[];
    }
    val == 2;
  }
}
