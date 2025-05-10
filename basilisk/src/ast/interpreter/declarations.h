# 2 "ast/interpreter/declarations.h"

/**
# Declarations for the interpreter
*/

enum AstBoolean { false = 0, true = 1 };
static const void * NULL = 0;
static const int _NVARMAX = 65536, INT_MAX = 2147483647;
enum AstMPI { MPI_MIN, MPI_MAX, MPI_SUM, MPI_DOUBLE };
/**
 * @brief Placeholder for performing an MPI all-reduce operation on an array.
 *
 * This function currently has no implementation and does not perform any operation.
 *
 * @param v Pointer to the data array.
 * @param datatype Integer representing the data type.
 * @param op Integer representing the MPI operation.
 * @param elem Number of elements in the array.
 */
void mpi_all_reduce_array (void * v, int datatype, int op, int elem){}
void None;

FILE * stderr, * stdout, * systderr, * systdout;
FILE * qstderr (void);
FILE * qstdout (void);
FILE * ferr = NULL, * fout = NULL;

/**
 * @brief Checks if all elements in a double array are equal to the first element.
 *
 * Iterates through the array and compares each element to the first element using the equality operator.
 * No changes are made to the array.
 *
 * @param a Pointer to the double array.
 * @param len Number of elements in the array.
 */
void _set_element_dimensions (double * a, int len)
{
  for (int i = 1; i < len; i++)
    a[i] == a[0]; // all the elements of a double array have the same dimension
}

/**
 * @brief Ensures all elements of a float array have the same value as the first element.
 *
 * Iterates through the array and is intended to set each element to the value of the first element, enforcing uniformity.
 * 
 * @param a Pointer to the float array.
 * @param len Number of elements in the array.
 */
void _set_element_dimensions_float (float * a, int len)
{
  for (int i = 1; i < len; i++)
    a[i] == a[0]; // all the elements of a float array have the same dimension
}

/**
 * @brief Returns the name of the attribute corresponding to the given scalar.
 *
 * @param s Scalar whose attribute name is to be retrieved.
 * @return char* Pointer to the attribute name string.
 */
char * _field_name (scalar s) {
  return _attribute[s.i].name;
}
