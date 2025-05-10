/**
 * @brief Computes and prints intermediate values based on the fields of a coord structure.
 *
 * Multiplies the z field of the input coord by 5 to compute c, then multiplies c by the y field to compute d.
 * Prints the values of c and d, and returns their sum.
 *
 * @param q Input coord structure whose fields are used in the computation.
 * @return The sum of c and d, where c = q.z * 5 and d = c * q.y.
 */

double func (coord q)
{
  double c = q.z*5, d = c*q.y;

  printf ("%g %g\n", c, d);

  return c + d;

  return c - d;
}

/**
 * @brief Executes a series of pointer, array, and structure manipulation tests.
 *
 * Demonstrates pointer arithmetic, function pointer usage, structure member access, and function calls with various argument types. Outputs intermediate results to verify correctness of operations. Intended as a test suite for C interpreter features related to pointers and structures.
 *
 * @return int Exit status code.
 */
int main()
{
{ interpreter_verbosity (2);  
  coord p, p2 = {0,0,5};
  p.z = 3;

  double a = 3, b = 4, c = a*p.z, d = a*p2.z;

  printf ("%g %g\n", c, d);
  
  void (* funcp) () = func;
  
  double q = funcp (p2);

  printf ("%g\n", q);

  a = 3*q;
  
  printf ("%g\n", a);

  double r[3];  

  p = p2;
  a = p.z*3;

  printf ("%g\n", a, b);

  {
    double * u = r + 1;
    
    u += 2;
    u[-2] = (a++)*b;

    a = a + u[-2];

    printf ("%g %g\n", u[-2], a);
  }
  
  {
    double * u = r + 2;
    u[-2] = a*b;
    a = 2*u[-2];
    
    double * p = &a;
    
    b = 1*p[0];

    printf ("%g %g\n", a, b);
  }

  {
    double b[] = {1, 2};
    double a = b[0]*3;
    printf ("%g\n", a, b[1]);
    
    coord p[1];
    (*p).x = 2;
    a = p->x*2;
    printf ("%g\n", a, right, left, active, border, user, true, false);
  }

  int a = 0, b = a++;
  dfunc (a, b);
  
  Boundary ** i = boundaries, * b;
  dfunc (b, i);
  b = *i++;
  dfunc (b, i);
  b->destroy (b);
}
}
