/**
 * @brief Tests that all elements of a dimensioned array, including those beyond the interpreter's default initialization range, have consistent dimension metadata.
 *
 * Initializes a double array of size 100 with a dimensioned value. Uses a conditional expression to compare the dimensions of the first and last elements, triggering an error if they differ. This verifies that uninitialized elements beyond the default iteration limit retain the same dimension as initialized ones.
 *
 * @return int Always returns 0.
 */

int main()
{

  /**
  The default maximum number of iterations of the interpreter is 32,
  so the last elements of the array will not be initialised. */
  
  double a[100];
  for (int i = 0; i < 100; i++)
    a[i] = 2. [1];
  
  /**
  The undefined condition below will return an error if a[0] and a[99]
  don't have the same dimensions. */
  
  double cond, val;
  if (cond)
    val = a[0];
  else
    val = a[99];
}
