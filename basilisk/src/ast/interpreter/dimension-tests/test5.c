/**
Check the dimensions of implicit zeros in structures. */ 

struct Func {
  double a; // explicit (1)
  double b; // implicit (0)
};

/**
 * @brief Tests implicit zero initialization of struct members.
 *
 * Sets interpreter verbosity, creates a Func struct with only the first member initialized, and displays both members to verify that the second is zero-initialized.
 *
 * @return int Exit status code.
 */
int main()
{
  { interpreter_verbosity (2);
    struct Func p = {1};
    display_value (p.a);
    display_value (p.b);
  }
}
