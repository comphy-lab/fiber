/**
 * @brief Returns the sum of a local constant and the input, preserving the dimension of the input.
 *
 * The local constant is initialized to 1.0 and inherits the dimension of the parameter `b`.
 *
 * @param b Input value whose dimension is used for the local constant.
 * @return double The sum of the local constant and `b`, with the same dimension as `b`.
 */

double func (double b) {
  double a = 1.; // this local constant has the dimension of b
  return a + b;
}

/**
 * @brief Runs dimension tests for the interpreter.
 *
 * Sets the interpreter verbosity level and calls the test function `func` with arguments of different dimensions to verify correct handling of local constants with varying dimensions.
 *
 * @return int Exit status code.
 */
int main()
{
  {
    interpreter_verbosity (4);
    func (1. [1]);
    func (1. [2]);
  }
}
