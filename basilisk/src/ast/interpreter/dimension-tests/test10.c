/**
 * @brief Tests that the expression 1./Delta has the correct dimension of [-1].
 *
 * Initializes a variable with dimension [1], iterates over a domain, and displays the value of 1./Delta to verify dimensional correctness.
 *
 * @return int Exit status code.
 */

int main()
{
  L0 = 1. [1];
  foreach()
    display_value (1./Delta);
}
