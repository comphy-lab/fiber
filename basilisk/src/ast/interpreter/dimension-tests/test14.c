/**
 * @brief Tests the correct definition of point variables during face iteration.
 *
 * Initializes a grid of size 1 and iterates over each face, displaying the value of the point variable to verify its proper definition and accessibility.
 *
 * @return int Exit status code.
 */

int main()
{
  init_grid (1);
  foreach_face()
    display_value (point);
}
