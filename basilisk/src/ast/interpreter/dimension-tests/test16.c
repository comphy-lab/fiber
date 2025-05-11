/**
 * @brief Demonstrates dimensioning of constant expressions and constants in C.
 *
 * Initializes variables using dimensioned constant expressions and constants, performs arithmetic operations, and displays results to illustrate the behavior of dimensioning in constant expressions.
 *
 * @return int Exit status code.
 */

int main()
{

  /**
  Constant expressions can be dimensioned. */
  
  double a = (1. + 2.)[1];
  display_value (a);
  double b = 3. + a;

  /**
  Constants in constant expressions can also be dimensioned (but this
  is silly and confusing). */
  
  double c = 4.[1] + 5.;
  double d = 6. + c;
}
