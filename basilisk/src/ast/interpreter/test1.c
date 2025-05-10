/**
 * @brief Placeholder function accepting an integer argument.
 *
 * Intended for use in interpreter loop tests; performs no operations.
 *
 * @param i Integer value passed to the function.
 */

void func (int i) {}

/**
 * @brief Tests various loop constructs and structure usage for the interpreter.
 *
 * Executes while, do-while, and for loops, each calling a test function with loop variables. Also tests structure initialization and the sizeof operator by passing results to the test function. Intended as a functional test for interpreter loop and structure handling.
 *
 * @return int Exit status code.
 */
int main()
{
{
  interpreter_verbosity (4);
  
  int i = 0;
  while (i < 3) {
    func (i);
    i++;
  }

  do {
    i--;
    func (i);
  } while (i > 0);
  
  for (int j = 0; j < 3; j++)
    func (j);

  /**
  And some structure tests which have nothing to do with loops. */
  
  func (sizeof (scalar));

  scalar s;
  s = (scalar){1};
  func (s.i);
 }
}
