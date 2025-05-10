/**
 * @brief Tests dimension consistency in branched expressions.
 *
 * Initializes a grid and evaluates conditional assignments to verify that both branches of a conditional expression have matching dimensions. Demonstrates valid and invalid cases, including an intentional dimension mismatch to trigger an error.
 *
 * @return int Exit status code.
 */

int main()
{
  init_grid (1);
  {
    double a, b;

    if (a)
      b = 1. [1];
    else
      b = 2.;
    display_value (b);

    scalar s[];
    foreach() {
      if (a)
	s[] = 1. [1];
      else
      	s[] = 0.;
      display_value (s[]);
    }

    if (a)
      b = 1. [1];
    else
      b = 1. [2]; // should give an error since [2] != [1]

  }
}
