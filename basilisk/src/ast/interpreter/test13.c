/**
 * @brief Tests initialization and accessibility of structure members.
 *
 * Verifies that structure members are properly initialized to zero and remain accessible, even when conditional assignments are involved.
 *
 * @return int Always returns 0.
 */

int main()
{
  {
    interpreter_verbosity (4);
    coord a = {0}, c = {1};
    double b = 0;
    if (b)
      a = c;
    if (a.y); // must be defined and zero
  }
}
