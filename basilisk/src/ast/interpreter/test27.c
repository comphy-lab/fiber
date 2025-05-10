/**
 * @brief Returns the first value of a loop or -1 if the loop does not execute.
 *
 * Iterates from 0 to 2 and returns the current index on the first iteration. If the loop does not execute, returns -1.
 *
 * @return int The value 0 if the loop runs; otherwise, -1.
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
 * Calls func() and passes its return value to display_value().
 *
 * @return int Exit status of the program.
 */
int main()
{
  display_value (func());
}
