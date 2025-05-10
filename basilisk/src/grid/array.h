// Arrays

typedef struct {
  void * p;
  long max, len;
} Array;

/**
 * @brief Allocates and initializes a new dynamic array structure.
 *
 * @return Pointer to a newly created Array with zero length and capacity.
 */
Array * array_new()
{
  Array * a = qmalloc (1, Array);
  a->p = NULL;
  a->max = a->len = 0;
  return a;
}

/**
 * @brief Frees the memory used by an Array and its internal buffer.
 *
 * Releases both the data buffer and the Array structure itself.
 *
 * @param a Pointer to the Array to be freed.
 */
void array_free (Array * a)
{
  free (a->p);
  free (a);
}

/**
 * @brief Appends a block of data to the dynamic array, resizing if necessary.
 *
 * Copies `size` bytes from `elem` to the end of the array's buffer, expanding the buffer if there is insufficient space. Updates the array's length and returns a pointer to the newly appended data within the buffer.
 *
 * @param elem Pointer to the data to append.
 * @param size Number of bytes to append.
 * @return void* Pointer to the location in the array's buffer where the new data was added.
 */
void * array_append (Array * a, void * elem, size_t size)
{
  if (a->len + size >= a->max) {
    a->max += max (size, 4096);
    a->p = realloc (a->p, a->max);
  }
  memcpy (((char *)a->p) + a->len, elem, size);
  a->len += size;
  return (void *)(((char *)a->p) + a->len - size);
}

/**
 * @brief Shrinks the array's buffer to its current length and frees the Array structure.
 *
 * Reallocates the internal buffer to exactly the number of bytes currently used (`len`),
 * frees the Array structure, and returns a pointer to the resized buffer.
 *
 * @return void* Pointer to the resized buffer containing the array's data.
 */
void * array_shrink (Array * a)
{
  void * p = realloc (a->p, a->len);
  free (a);
  return p;
}
