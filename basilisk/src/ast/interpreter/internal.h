# 2 "ast/interpreter/internal.h"

/**
 * @brief Allocates a block of memory of the specified size.
 *
 * @param size Number of bytes to allocate.
 * @return Pointer to the beginning of the allocated memory block, or NULL if allocation fails.
 */

void * malloc (long size) {}
/**
 * @brief Allocates memory for an array of elements and initializes it to zero.
 *
 * Allocates enough memory for an array of nmemb elements of size bytes each, and sets all bytes in the allocated memory to zero.
 *
 * @param nmemb Number of elements to allocate.
 * @param size Size in bytes of each element.
 * @return Pointer to the allocated zero-initialized memory, or NULL if allocation fails.
 */
void * calloc (long nmemb, long size) {}
/**
 * @brief Changes the size of a previously allocated memory block.
 *
 * Resizes the memory block pointed to by ptr to the specified size. If ptr is NULL, behaves like malloc. If size is zero and ptr is not NULL, the memory is freed and NULL is returned.
 *
 * @param ptr Pointer to the memory block to resize, or NULL to allocate a new block.
 * @param size New size for the memory block in bytes.
 * @return void* Pointer to the resized memory block, or NULL if allocation fails.
 */
void * realloc (void * ptr, long size) {}
/**
 * @brief Deallocates memory previously allocated by memory allocation functions.
 *
 * Releases the memory block pointed to by ptr, making it available for future allocations.
 *
 * @param ptr Pointer to the memory block to be freed. If ptr is NULL, no action is taken.
 */
void free (void * ptr) {}
/**
 * @brief Fills a block of memory with a specified byte value.
 *
 * Sets the first n bytes of the memory area pointed to by s to the byte value c.
 *
 * @param s Pointer to the memory area to fill.
 * @param c Byte value to set.
 * @param n Number of bytes to set.
 * @return Pointer to the memory area s.
 */
void * memset (void * s, int c, long n) {}
/**
 * @brief Copies a specified number of bytes from the source to the destination buffer.
 *
 * Copies exactly @p n bytes from the memory area pointed to by @p src to the memory area pointed to by @p dest. The memory areas must not overlap.
 *
 * @param dest Pointer to the destination buffer.
 * @param src Pointer to the source buffer.
 * @param n Number of bytes to copy.
 * @return Pointer to the destination buffer.
 */
void * memcpy (void * dest, const void * src, long n) {}
/**
 * @brief Creates a duplicate of a null-terminated string.
 *
 * Allocates memory and copies the contents of the given string, returning a pointer to the new copy.
 *
 * @param s The null-terminated string to duplicate.
 * @return char* Pointer to the newly allocated duplicate string, or NULL if allocation fails.
 */
char * strdup (const char * s) {}
/**
 * @brief Copies a null-terminated string to a destination buffer.
 *
 * Copies the string pointed to by src, including the terminating null byte, to the buffer pointed to by dest.
 * The destination buffer must be large enough to receive the copy.
 *
 * @return Pointer to the destination buffer dest.
 */
char * strcpy (char *dest, const char *src) {}
/**
 * @brief Appends the source string to the destination string.
 *
 * Concatenates the null-terminated string pointed to by src to the end of the null-terminated string pointed to by dest. The destination buffer must be large enough to hold the resulting string.
 *
 * @param dest Pointer to the destination string buffer.
 * @param src Pointer to the source string to append.
 * @return Pointer to the destination string.
 */
char * strcat (char * dest, const char * src) {}
/**
 * @brief Returns the length of a null-terminated string.
 *
 * @param s Pointer to the null-terminated string.
 * @return long Number of characters in the string, excluding the null terminator.
 */
long strlen (const char * s) {}
/**
 * @brief Compares two strings lexicographically.
 *
 * @param s1 Pointer to the first null-terminated string.
 * @param s2 Pointer to the second null-terminated string.
 * @return An integer less than, equal to, or greater than zero if s1 is found, respectively, to be less than, to match, or be greater than s2.
 */
int strcmp (const char * s1, const char * s2) {}
/**
 * @brief Compares up to a specified number of characters of two strings.
 *
 * Compares at most `size` characters of the strings `s1` and `s2` and returns an integer indicating their lexical relationship.
 *
 * @param s1 First string to compare.
 * @param s2 Second string to compare.
 * @param size Maximum number of characters to compare.
 * @return int Less than, equal to, or greater than zero if `s1` is found to be less than, to match, or be greater than `s2`, respectively, within the first `size` characters.
 */
int strncmp (const char * s1, const char * s2, long size) {}

/**
 * @brief Returns the absolute value of a floating-point number.
 *
 * @param x The input value.
 * @return The non-negative value of x.
 */
double fabs (double x) {}
/**
 * @brief Computes the square root of a non-negative double-precision value.
 *
 * @param x The value for which to compute the square root.
 * @return The non-negative square root of x.
 */
double sqrt (double x) {}
/**
 * @brief Computes the exponential of a floating-point value.
 *
 * @param a The exponent.
 * @return double The value of e raised to the power of a.
 */
double exp (double a) {}
/**
 * @brief Computes the natural logarithm of a floating-point value.
 *
 * @param x The value for which to compute the natural logarithm.
 * @return The natural logarithm (base e) of x.
 */
double log (double x) {}
/**
 * @brief Computes the base-10 logarithm of a floating-point value.
 *
 * @param a The input value.
 * @return The base-10 logarithm of a.
 */
double log10 (double a) {}
/**
 * @brief Computes the sine of a given angle in radians.
 *
 * @param x Angle in radians.
 * @return double Sine of the angle.
 */
double sin (double x) {}
/**
 * @brief Computes the cosine of a given angle in radians.
 *
 * @param a Angle in radians.
 * @return double Cosine of the angle.
 */
double cos (double a) {}
/**
 * @brief Computes the tangent of a given angle in radians.
 *
 * @param x Angle in radians.
 * @return double Tangent of the angle.
 */
double tan (double x) {}
/**
 * @brief Computes the arc sine of a value.
 *
 * @param x Value whose arc sine is to be calculated, in the range [-1, 1].
 * @return double Arc sine of x, in radians.
 */
double asin (double x) {}
/**
 * @brief Computes the arc cosine (inverse cosine) of a value.
 *
 * @param a Value whose arc cosine is to be calculated, in the range [-1, 1].
 * @return double Arc cosine of a in radians, in the range [0, π].
 */
double acos (double a) {}
/**
 * @brief Computes the arc tangent of a value.
 *
 * @param x The value whose arc tangent is to be calculated.
 * @return double The principal value of the arc tangent of x, in radians.
 */
double atan (double x) {}
/**
 * @brief Computes the hyperbolic sine of a double-precision value.
 *
 * @param a The value for which to compute the hyperbolic sine.
 * @return The hyperbolic sine of the input value.
 */
double sinh (double a) {}
/**
 * @brief Computes the hyperbolic cosine of a value.
 *
 * @param a The value for which to compute the hyperbolic cosine.
 * @return double The hyperbolic cosine of the input.
 */
double cosh (double a) {}
/**
 * @brief Computes the hyperbolic tangent of a value.
 *
 * @param a The input value.
 * @return The hyperbolic tangent of a.
 */
double tanh (double a) {}
/**
 * @brief Computes the inverse hyperbolic sine of a value.
 *
 * @param a The value for which to compute the inverse hyperbolic sine.
 * @return The inverse hyperbolic sine of a.
 */
double asinh (double a) {}
/**
 * @brief Computes the inverse hyperbolic cosine of a value.
 *
 * @param a The value for which to compute the inverse hyperbolic cosine.
 * @return double The inverse hyperbolic cosine of the input.
 */
double acosh (double a) {}
/**
 * @brief Computes the inverse hyperbolic tangent of a value.
 *
 * @param a The value for which to compute the inverse hyperbolic tangent.
 * @return double The inverse hyperbolic tangent of the input.
 */
double atanh (double a) {}

/**
 * @brief Computes the arc tangent of y/x using the signs of both arguments to determine the correct quadrant.
 *
 * @param y The ordinate (vertical coordinate).
 * @param x The abscissa (horizontal coordinate).
 * @return double The angle in radians between the positive x-axis and the point (x, y), in the range [-π, π].
 */
double atan2 (double y, double x) {}
/**
 * @brief Raises a number to the power of another.
 *
 * Calculates x raised to the power of y (x^y) and returns the result.
 *
 * @param x The base value.
 * @param y The exponent value.
 * @return double The result of x raised to the power y.
 */
double pow (double x, double y) {}

/**
 * @brief Sets the verbosity level of the interpreter.
 *
 * Adjusts the amount of diagnostic or informational output produced by the interpreter based on the specified verbosity level.
 *
 * @param verbosity The desired verbosity level.
 */
void interpreter_verbosity (int verbosity) {}
/**
 * @brief Sets the maximum number of iterations allowed in the interpreter.
 *
 * Limits the number of iterations that interpreter-controlled loops or processes may execute.
 */
void interpreter_maximum_iterations (int maximum_iterations) {}
/**
 * @brief Displays the given value in a human-readable format.
 *
 * Intended for outputting or inspecting interpreter values during execution.
 */
void display_value (void * value) {}
/**
 * @brief Resets the value of a specified field element.
 *
 * Sets the value of the field element identified by name and block index to the given value.
 *
 * @param field Pointer to the field array to modify.
 * @param name Name identifying the field element.
 * @param val Value to assign to the field element.
 * @param block Index of the block containing the field element.
 */
void reset_field_value (real * field, const char * name, real val, int block) {}
