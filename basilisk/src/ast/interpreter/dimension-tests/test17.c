/**
 * @brief Demonstrates arithmetic operations with dimensioned constant expressions.
 *
 * Initializes double precision variables and performs addition with a dimensioned constant to test handling of such expressions.
 *
 * @return int Always returns 0.
 */

int main()
{
  double a = 1. + 2.;
  double c = 1. + 2.*sqrt(2.);
  double b = 4 [1];
  a + b;
  c + b;
}
