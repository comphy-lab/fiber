# 2 "ast/interpreter/internal.h"

/**
# Internally defined functions */

void * malloc (long size) {}
void * calloc (long nmemb, long size) {}
void * realloc (void * ptr, long size) {}
void free (void * ptr) {}
void * memset (void * s, int c, long n) {}
void * memcpy (void * dest, const void * src, long n) {}
char * strdup (const char * s) {}
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
 * @param s1 First null-terminated string to compare.
 * @param s2 Second null-terminated string to compare.
 * @return int An integer less than, equal to, or greater than zero if s1 is found, respectively, to be less than, to match, or be greater than s2.
 */
int strcmp (const char * s1, const char * s2) {}
/**
 * @brief Compares two strings up to a specified number of characters.
 *
 * Compares at most `size` characters of the strings `s1` and `s2`.
 * Returns an integer less than, equal to, or greater than zero if `s1` is found,
 * respectively, to be less than, to match, or be greater than `s2`.
 *
 * @param s1 First string to compare.
 * @param s2 Second string to compare.
 * @param size Maximum number of characters to compare.
 * @return int Negative if `s1` < `s2`, zero if equal, positive if `s1` > `s2`.
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
 * @brief Computes the square root of a non-negative number.
 *
 * @param x The value for which to compute the square root.
 * @return double The non-negative square root of x.
 *
 * If x is negative, the behavior is undefined.
 */
double sqrt (double x) {}
double exp (double a) {}
double log (double x) {}
double log10 (double a) {}
double sin (double x) {}
double cos (double a) {}
double tan (double x) {}
double asin (double x) {}
double acos (double a) {}
double atan (double x) {}
double sinh (double a) {}
double cosh (double a) {}
double tanh (double a) {}
double asinh (double a) {}
double acosh (double a) {}
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
 * Calculates x raised to the power of y (x^y).
 *
 * @param x The base value.
 * @param y The exponent value.
 * @return double The result of x raised to the power y.
 */
double pow (double x, double y) {}

/**
 * @brief Sets the verbosity level for the interpreter.
 *
 * Adjusts the amount of diagnostic or informational output produced during interpreter execution.
 *
 * @param verbosity The desired verbosity level.
 */
void interpreter_verbosity (int verbosity) {}
/**
 * @brief Sets the maximum number of iterations for the interpreter.
 *
 * Limits the number of iterations the interpreter will perform during execution.
 *
 * @param maximum_iterations The maximum allowed iteration count.
 */
void interpreter_maximum_iterations (int maximum_iterations) {}
/**
 * @brief Displays the provided value in a human-readable format.
 *
 * Intended for outputting or inspecting the contents of a value during interpretation or debugging.
 *
 * @param value Pointer to the value to be displayed.
 */
void display_value (void * value) {}
/**
 * @brief Resets the value of a specified field element.
 *
 * Sets the value of the field element identified by name and block index to the given value.
 *
 * @param field Pointer to the field array to modify.
 * @param name Name identifying the field element to reset.
 * @param val Value to assign to the field element.
 * @param block Index of the block containing the field element.
 */
void reset_field_value (real * field, const char * name, real val, int block) {}
