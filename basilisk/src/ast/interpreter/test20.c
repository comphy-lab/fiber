/**
 * @brief Tests memory allocation and freeing under an undefined conditional.
 *
 * Declares an uninitialized variable and uses it as a condition to allocate and free memory, intended for use with Valgrind to detect issues related to undefined behavior and memory management.
 *
 * @return int Exit status code.
 */

int main()
{
  double a;
  if (a) {
    double * b = malloc (2*sizeof (double));
    b[0] = 1.;
    b[1] = 2.;
    free (b);
  }
}
