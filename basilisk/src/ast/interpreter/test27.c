/**
 * @brief Returns the first value in a loop or -1 if the loop does not execute.
 *
 * Iterates from 0 to 2 and immediately returns the current index on the first iteration.
 * If the loop does not execute, returns -1.
 *
 * @return int The first loop index (0), or -1 if the loop is skipped.
 */

int func()
{
  for (int i = 0; i < 3; i++)
    return i;
  return -1;
}

/**
 * @brief Entry point of the program.
 *
 * Calls func() and passes its return value to display_value() for output.
 *
 * @return int Exit status code.
 */
int main()
{
  display_value (func());
}
