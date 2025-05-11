/**
 * @brief Returns the sum of the input and a dimension-matching constant.
 *
 * Adds 1 to the input value, where the constant 1 inherits the dimension of the argument `x`.
 * This allows the function to operate correctly with inputs of varying dimensions.
 *
 * @param x Input value, whose dimension determines the dimension of the constant added.
 * @return double The result of `x + 1`, with the same dimension as `x`.
 */

double sum (double x)
{
  return x + 1.; // 1. has the dimension of x, so will change between
		 // the calls below
}

/**
 * @brief Executes test cases to verify dimension propagation in function calls.
 *
 * Initializes variables and calls the sum function with arguments that have attached dimensions, demonstrating that constants within functions can inherit and adapt to the dimension of their arguments.
 *
 * @return int Exit status code.
 */
int main()
{
  double a = 2;
  a += sum (1. [1]);
  sum (2. [2]);
#if 0
  double b = 1.;
  sum (sq(b));
#endif
}
