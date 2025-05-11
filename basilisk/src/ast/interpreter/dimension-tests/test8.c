/**
 * @brief Tests dimension propagation in conditional expressions within a computational grid.
 *
 * Initializes a grid and assigns values to a scalar field using both if-else and ternary conditional expressions, verifying that dimensions are handled consistently. Also checks that multiplicative conditional expressions with constants are treated as dimensionless and that conditional branches with dimensioned values enforce consistent dimensions.
 *
 * Outputs values for manual or automated verification of dimension correctness.
 *
 * @return int Exit status code.
 */

int main()
{
  init_grid (1);
  scalar f[];
  foreach() {
    double a;
    if (a)
      f[] = 1. [1];
    else
      f[] = 2.;
  }
  foreach()
    f[] = f[] < 3. ? 4. : f[];

  /**
  Multiplicative conditional expressions with a constant must be
  dimensionless. So [d] should be [c] since [(b > 0 ? 1 : -2)] = [0]
  since it is a multiplicative constant. */

  double b, c = 7, d;
  d = c*(b > 0 ? 1 : -2);
  display_value (d);

  /**
  Here, there is only one constant (3) in the conditional expression,
  whose dimension must [e] == [-1]. */
  
  double e = 4 [-1];
  d = c*(b > 0 ? 3 : e);
  display_value (d);

  /**
  Values of undefined conditional expressions must have the same
  dimensions (equivalent to [test1.c]() for undefined branching). */

  double d = b ? 1 [2] : 2;
  3. + d; // [3] should be [2]
}
