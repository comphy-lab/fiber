/**
# Variable scopes */

double c;

/**
 * @brief Displays the values of the given argument and the global variable c.
 *
 * Calls display_value on the provided double x and on the global variable c.
 *
 * @param x The double value to display.
 */
double func (double x)
{
  display_value (x);
  display_value (c);
}

/**
 * @brief Entry point for testing variable scope and function calls with nested blocks.
 *
 * Demonstrates the visibility and lifetime of local and global variables by declaring variables in nested blocks, invoking `display_value` on them, and calling `func` with different arguments.
 * No value is returned.
 */
int main()
{ // 2
  func (c);
  { // 3
    double a;
    display_value (a);
    func (a);
    { // 4
      double b;
      display_value (b);
      display_value (c);
      func (c);
    }
    display_value (a);
  }
}
