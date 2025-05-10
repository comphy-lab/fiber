/**
 * @brief Tests reporting of unset variable usage in constraint expressions.
 *
 * Demonstrates that unset variables assigned but not used in constraints are not reported, while unset variables used in constraint expressions trigger a report. Also includes conditional assignments within a loop and checks constraint evaluation after assignments.
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
