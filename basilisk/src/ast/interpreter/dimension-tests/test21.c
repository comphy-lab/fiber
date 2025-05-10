/**
 * @brief Tests that array element dimensions are handled identically.
 *
 * Initializes a grid and verifies that arithmetic operations involving elements of double and float arrays with `L0` are processed without dimension mismatches.
 *
 * @return int Exit status code.
 */

int main()
{
  init_grid (1);
  double a[2] = {1,2};
  a[0] + L0;
  float b[2] = {3,4};
  b[0] + L0;
}
