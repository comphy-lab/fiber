#include <stdlib.h>
#include <string.h>
#include "allocator.h"

struct _Allocator {
  void ** m;
  int n;
  long len, maxlen;
};

/**
 * @brief Creates and initializes a new memory allocator.
 *
 * Allocates and returns a new Allocator instance with a default memory block size of 65,536 bytes.
 *
 * @return Pointer to the newly created Allocator.
 */
Allocator * new_allocator()
{
  Allocator * a = calloc (1, sizeof (Allocator));
  a->len = a->maxlen = 1 << 16;
  return a;
}

/**
 * @brief Allocates a zero-initialized memory block of the specified size from the allocator.
 *
 * If the current memory block does not have enough space, a new block is allocated. The returned pointer is always aligned within the allocator's managed memory and is zero-initialized.
 *
 * @param a Pointer to the allocator instance.
 * @param size Number of bytes to allocate.
 * @return void* Pointer to the allocated memory block.
 */
void * allocate (Allocator * a, long size)
{
  if (a->len + size >= a->maxlen) {
    a->n++;
    a->m = realloc (a->m, a->n*sizeof (void *));
    a->m[a->n - 1] = malloc (size > a->maxlen ? size : a->maxlen);
    a->len = 0;
  }
  void * p = (void *)(((char *)a->m[a->n - 1]) + a->len);
  memset (p, 0, size);
  a->len += size;
  return p;
}

/**
 * @brief Frees all memory managed by the allocator.
 *
 * Releases all memory blocks allocated by the allocator, as well as the allocator's internal structures.
 * After calling this function, the allocator pointer is no longer valid.
 */
void free_allocator (Allocator * a)
{
  for (int i = 0; i < a->n; i++)
    free (a->m[i]);
  free (a->m);
  free (a);
}
