/**
 * @brief Tests behavior when accessing an uninitialized structure member.
 *
 * Declares a structure variable, assigns a value to one member, and reads both the initialized and uninitialized members to observe the resulting values.
 *
 * @return int Exit status code.
 */

int main()
{{ interpreter_verbosity (4);
  coord p;
  p.x = 2.;
  double a, b;
  a = p.x;
  b = p.y;
}}
