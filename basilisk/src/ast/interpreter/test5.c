/**
(Probably) checks for potential scope troubles */

int b = 1;

/**
 * @brief Assigns the value of the global variable `b` to the parameter `a`.
 *
 * This function demonstrates assignment from a global variable, potentially illustrating variable shadowing or scope behavior.
 */
void func (int a)
{
  a = b;
}

/**
 * @brief Entry point of the program demonstrating variable shadowing and scope.
 *
 * Sets interpreter verbosity, declares a local struct variable that shadows a global variable, and calls a function to illustrate scope resolution.
 *
 * @return int Exit status code.
 */
int main()
{
  {
    interpreter_verbosity (4);
    struct {
      double x;
    } b;
    func (1);
  }
}
