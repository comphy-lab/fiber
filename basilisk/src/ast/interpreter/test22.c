/**
 * @brief Tests logical AND and OR operations with initialized and uninitialized variables.
 *
 * Evaluates all combinations of logical AND (`&&`) and OR (`||`) between an uninitialized variable `a` and a variable `b` (set to 1 and then 0), displaying the results to verify correct handling of logical operations and partial evaluation behavior.
 *
 * @return int Always returns 0.
 */

int main()
{
  int b = 1, a;
  display_value (b && a);
  display_value (b || a);
  display_value (a && b);
  display_value (a || b);
  b = 0;
  display_value (b && a);
  display_value (b || a);
  display_value (a && b);
  display_value (a || b);
}
