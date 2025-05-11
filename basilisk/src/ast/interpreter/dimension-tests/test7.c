/**
Checks minmod2() */

#include "utils.h"

/**
 * @brief Tests the minmod2 function with unset double values and displays the results.
 *
 * Initializes three double variables with specific values, unsets them, and passes them to minmod2.
 * The result is displayed using display_value. This process is repeated with a second set of values.
 *
 * @return int Exit status code.
 */
int main()
{
  {
    double a = 1.[1], b = 2.[1], c = 3.[1];
    unset_double (&a);
    unset_double (&b);
    unset_double (&c);
    display_value (minmod2 (a, b, c));
    
    a = 1.[2], b = 2.[2], c = 3.[2];
    unset_double (&a);
    unset_double (&b);
    unset_double (&c);
    display_value (minmod2 (a, b, c));
  }
}
