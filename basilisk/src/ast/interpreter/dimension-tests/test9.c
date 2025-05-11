/**
 * @brief Tests conversion of floating-point values to boolean in conditional statements.
 *
 * Verifies that nonzero double values are treated as true by the interpreter in `if` conditions,
 * both for variables and direct floating-point literals.
 *
 * @return int Exit status code.
 */

int main()
{
  double a = 0.1;
  if (a)
    display_value (1); // must display 1
  if (0.1 || 0.2)
    display_value (1); // must display 1
}
