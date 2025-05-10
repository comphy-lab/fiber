/**
 * @brief Demonstrates and verifies various mathematical and grid operations.
 *
 * Performs a sequence of floating-point mathematical computations, including multiplication, square root, power, absolute value, a custom square function, and sine. Initializes a computational grid, assigns values to scalar fields based on trigonometric calculations, and iterates over grid cells to perform additional arithmetic operations. Interpreter verbosity is set to a high level throughout for detailed output.
 *
 * @return int Returns 0 on successful execution.
 */

int main()
{
  double c;
  {
    interpreter_verbosity (4);
    double a = 1, b = 2;

    double c = a*b;
    double d = sqrt (c);
    double e = pow (d, 2);
    double f = fabs (e);
    double g = sq(f);
    double a = 1.;
    a *= g;
    double b = sin (a/g);
  }
  
  init_grid (1);

  scalar s[], f[];
  foreach() {
    interpreter_verbosity (4);
    double a = L0*sin (x/L0);
    s[] = a*a;
    f[] = 2;

  }

  foreach() {
    interpreter_verbosity (4);
    f[] + L0*L0;
    s[] + L0*L0;
  }
}
