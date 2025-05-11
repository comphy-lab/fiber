/**
 * @brief Computes a value based on the fields of a coord structure.
 *
 * Calculates two intermediate values using the `z` and `y` members of the input `coord`, then returns their sum.
 *
 * @param q The coord structure containing the input values.
 * @return double The computed sum of the intermediate values.
 */

double func (coord q)
{
  double c = q.z*5, d = c*q.y;
  return c + d;

  return c - d;
}

double func (coord p);

/**
 * @brief Executes tests involving function pointers, struct manipulation, and pointer arithmetic.
 *
 * Initializes and manipulates `coord` structures and doubles, demonstrates function pointer usage with structs, performs arithmetic operations, and exercises pointer arithmetic on arrays. Prints a computed value to standard output.
 *
 * @return int Exit status code.
 */
int main()
{{ interpreter_verbosity (4);
    
  coord p, p2 = {0,0,5};
  p.z = 3;
  
  double a = 3, b = 4, c = a*p.z, d = a*p2.z;

  void (* pfunc) (coord) = func;
  
  double q = pfunc (p2);

  a = 3*q;

  double r[3];  

  p = p2;
  a = p.z*3;

  printf ("%g\n", a);

  {
    double * u = r + 1;
    
    u += 2;
    u[-2] = (a++)*b;
    
    a = a + u[-2];
    
    c = 2*a;
  }
  
  {
    double * u = r + 2;
    u[-2] = a*b;
    a = 2*u[-2];
    
    double * p = &a;
    
    b = 1*p[0];
  }
}}
