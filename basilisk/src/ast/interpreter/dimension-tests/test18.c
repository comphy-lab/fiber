/**
 * @brief Tests conditional assignment and modification of a double return value.
 *
 * Returns a double whose value depends on the condition `a` and subsequent checks, used to verify dimension consistency in return values.
 *
 * @return double The resulting value after conditional assignments.
 */

double func1()
{
  double ret;
  if (a)
    ret = 1.;
  if (ret == 1.)
    ret = 2.;
  return ret;
}

/**
 * @brief Returns 2.0 if the condition 'a' is true, otherwise returns an uninitialized value.
 *
 * If the global variable 'a' is true, sets a local variable to 1.0 and returns 2.0.
 * If 'a' is false, returns the value of an uninitialized local variable, which may lead to undefined behavior.
 *
 * @return double The result based on the condition 'a'; may be undefined if 'a' is false.
 */
double func2()
{
  double ret;
  if (a)
    ret = 1.;
  if (ret == 1.)
    return 2.;
  return ret;
}

/**
 * @brief Returns a coord with its x component conditionally set based on the value of a.
 *
 * Iterates over each dimension, assigning 1.0 to ret.x if a is true, then updates ret.x to 2.0 if it equals 1.0. Returns the resulting coord.
 *
 * @return coord The resulting coord after conditional assignments in each dimension.
 */
coord func3()
{
  coord ret;
  foreach_dimension() {
    if (a)
      ret.x = 1.;
    if (ret.x == 1.)
      ret.x = 2.;
  }
  return ret;
}

/**
 * @brief Runs tests to verify dimension consistency in function return values.
 *
 * Sets interpreter verbosity and asserts that the return values of `func1`, `func2`, and the `x` component of `func3` match expected values across dimensions.
 *
 * @return int Exit status code.
 */
int main()
{
  {
    interpreter_verbosity (3);
    func1() == 3. [1];
    func2() == 3. [2];
    foreach_dimension()
      func3().x == 3. [3];
  }
}
