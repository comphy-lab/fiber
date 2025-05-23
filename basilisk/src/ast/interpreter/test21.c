/**
 * @brief Demonstrates evaluation and display of various expressions and variables.
 *
 * Exercises arithmetic operations, type conversions, variable assignments, constant usage, and conditional logic by displaying the results of multiple expressions using the display_value function.
 *
 * @return int Exit status code.
 */

int main()
{
  int a = 1;
  display_value (a);
  display_value (1);
  display_value (1*2);
  display_value (1*2 + 4);
  display_value (1 << 2 + 4);
  display_value (1*a + 4);
  display_value ((1*2 + 4)/32.);
  display_value (sqrt(2.));
  const double b = 3;
  display_value (2*b);
  a = b;
  display_value (a);

  {
    double a, b;
    if (a)
      b = 1.;
    else
      b = 2.;
    display_value (b);
  }

  display_value (sq(4));
}
