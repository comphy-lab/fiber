/**
 * @brief Tests array and struct initialization and basic operations.
 *
 * Initializes a double array and a struct array, performs assignments and arithmetic operations to verify correct handling of arrays, struct field access, and pointer dereferencing.
 *
 * @return int Always returns 0.
 */

int main()
{{ interpreter_verbosity (4);
  double b[2] = {1, 2};
  double a = b[0]*3;
  a;
  coord p[1];
  (*p).x = b[1];
  a = p->x*2;
}}
