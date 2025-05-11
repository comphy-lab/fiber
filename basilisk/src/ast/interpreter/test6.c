/**
 * @brief Demonstrates access to elements of a fixed-size float array parameter.
 *
 * This function is used to test or illustrate the handling of named parameters and array argument passing in function calls.
 *
 * @param fc Array of three floats passed as a named parameter.
 */

void bidule (float fc[3])
{
  fc[0];
  fc[1];
  fc[2];
}

/**
 * @brief Tests named parameter passing with array initialization in function calls.
 *
 * Sets the interpreter verbosity level and calls the `bidule` function using named parameter syntax with an inline array initializer to verify interpreter support for these features.
 *
 * @return int Exit status code.
 */
int main()
{
  {
    interpreter_verbosity (4);
    bidule (fc = {1,2,3});
  }
}
