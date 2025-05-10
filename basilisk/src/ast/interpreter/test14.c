/**
Checks Arrays. */

typedef int Elem;

/**
 * @brief Tests dynamic array creation, element appending, and pointer-based iteration.
 *
 * Creates a dynamic array of integers, appends elements, and verifies access to the internal data buffer using pointer arithmetic and iteration.
 *
 * @return int Exit status code.
 */
int main()
{
  {
    interpreter_verbosity (4);

    Array * a = array_new();
    Elem b = 1;
    array_append (a, &b, sizeof(Elem));

    Elem * j = (Elem *) a->p;
    for (int i = 0; i < a->len/sizeof(Elem); i++, j++)
      *j;

    b = 2;
    array_append (a, &b, sizeof(Elem));

    a->p;
    ((char *)a->p) + 5;
    
    j = (Elem *) a->p;
    for (int i = 0; i < a->len/sizeof(Elem); i++, j++) {
      j;
      *j;
    }
  }
}
