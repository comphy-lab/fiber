/**
 * @brief Tests that field properties remain defined after conditional allocation of temporary fields.
 *
 * Initializes a computational grid and declares a scalar field. Conditionally allocates a temporary scalar field inside an uninitialized conditional block, then verifies that properties such as `datasize` and `_attribute[1].freed` are properly set and not left unset, regardless of the conditional allocation path.
 *
 * @return int Exit status code.
 */

int main()
{
  init_grid (1);
  scalar s[];
  int cond;
  if (cond) {
    scalar b[];
  }
  display_value (datasize); // must not be 'unset'
  display_value (_attribute[1].freed); // must not be 'unset'
}
