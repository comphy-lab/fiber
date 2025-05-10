/**
 * @brief Tests that scalar, vector, and tensor fields have default dimension L0.
 *
 * Initializes a single-cell grid, declares scalar, vector, and tensor fields, resets them, and verifies that all components have the default dimension L0 using dimensional assertions.
 *
 * @return int Returns 0 on successful completion.
 */

int main()
{
  init_grid (1);
  scalar s[];
  vector u[];
  tensor t[];
  { interpreter_verbosity (2);
    reset ({u, s, t}, 0.);
    foreach() {
      dimensional (s[] == L0);
      foreach_dimension() {
	dimensional (u.x[] == L0);
	dimensional (t.x.x[] == L0);
	dimensional (t.x.y[] == L0);
      }
    }
  }
}
