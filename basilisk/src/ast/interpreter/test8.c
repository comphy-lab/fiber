/**
 * @brief Tests interpreter behavior with undefined and conditional expressions.
 *
 * Declares an uninitialized variable and evaluates it in conditional statements to observe how the interpreter handles undefined conditions and control flow.
 *
 * @return int Exit status code.
 */

int main()
{{ interpreter_verbosity (4);
  double a, b = 2, c = 1;
  if (a)
    c++;
  else
    c--;
  c;
  if (c)
    b;
  else
    b*b;
}}
