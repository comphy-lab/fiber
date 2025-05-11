/**
 * @brief Tests variable definedness after modification under undefined conditions.
 *
 * Demonstrates that if a variable modified within a conditional branch depending on an uninitialized value is restored to its original state, its use is not undefined; but if the variable is modified in only one branch, subsequent use may be undefined.
 *
 * @return int Always returns 0.
 */
  
int main()
{
  {
    interpreter_verbosity (2);
    int a = 0, b;
    if (b) {
      ++a;
      --a;
    }
    else
      a = 0;
    if (a); // this is not undefined since a went back to its original value

    a = 0;
    if (b)
      ++a;
    else
      a = 0;
    if (a); // this is undefined since one of the branches modifies a
  }
}
