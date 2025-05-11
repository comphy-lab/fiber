/**
 * @brief Tests correct resolution of external variable declarations.
 *
 * Verifies that an external declaration within a nested block refers to the global variable rather than a local variable with the same name.
 *
 * @return int Exit status code.
 */

int main()
{
  {
    interpreter_verbosity (5);
    int uf = 3;
    {
      extern int uf;
      fprintf (uf); // must be: fprintf (4)
    }
  }
}

int uf = 4;
