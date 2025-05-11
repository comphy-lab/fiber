/**
 * @brief Displays the values of a one-dimensional array of two doubles.
 *
 * Iterates over the provided array and calls display_value on each element.
 *
 * @param d Array of two double values to be displayed.
 */

void myfunc1 (double d[2])
{
  for (int i = 0; i < 2; i++)
    display_value (d[i]);
}

/**
 * @brief Iterates over a 2x2 array of doubles and displays each element.
 *
 * Calls display_value on each element of the provided two-dimensional array.
 *
 * @param a A 2x2 array of doubles to be displayed.
 */
void myfunc (double a[2][2])
{
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      display_value (a[i][j]);
}

/**
 * @brief Demonstrates initialization and usage of one- and two-dimensional arrays.
 *
 * Initializes variables and arrays, then tests array passing and iteration by displaying their elements using helper functions.
 * Exercises both postfix (compound literal) and declarator-based array initialization in C.
 *
 * @return int Exit status code.
 */
int main()
{{  interpreter_verbosity (2);
    double x1 = 1, y1 = 2, x2 = 3, y2 = 4;
    myfunc1 ((double[2]){x1, y1}); // this is a postfix unidimensional array
    myfunc ((double[2][2]) {{x1, y1},{x2, y2}}); // this is a postfix multidimensional array
    double b[2][2] = {{x1, y1},{x2, y2}}; // this is an init_declarator multidimensional array
    for (int i = 0; i < 2; i++)
      for (int j = 0; j < 2; j++)
	display_value (b[i][j]);
}}
