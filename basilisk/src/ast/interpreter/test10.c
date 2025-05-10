/**
 * @brief Tests interpreter handling of undefined conditions and branch evaluations in a grid loop.
 *
 * Initializes a scalar grid and iterates over its cells, evaluating conditional branches and variable assignments to verify correct interpreter behavior with potentially undefined values and conditions.
 */

int main()
{
  init_grid (1);
  scalar s[];
  foreach() {
    interpreter_verbosity (4);
    coord a = {33};
    if (max(0,s[]) > 0) {
      double q = 1;
      coord u = {0};
      if (q)
	a = u;
      q = a.x;
      if (q);	  
    }
    else {
      a.x; // this must be 33
      a.x = 2;	
    }
    if (a.x > 0) // this must be (unset)
      ;
  }
}
