/**
 * @brief Returns the constant value 1.0 regardless of input.
 *
 * @param x Input value (unused).
 * @return double The constant value 1.0.
 */

double func (double x)
{
  return 1.;
}

/**
 * @brief Entry point of the program.
 *
 * Tests the behavior of combining an invalid array access expression with a function returning a constant double.
 *
 * @return int Exit status code.
 */
int main()
{
  {
    double a = 1 [1] + func (1);
  }
}
