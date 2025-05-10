#include <stdlib.h>
#include <string.h>
#include "stack.h"

struct _Stack {
  void * p;
  int n, size;
  void * (* push) (Stack *, void *);
  void * data;
};

Stack * stack_new (int size)
{
  Stack * s = calloc (1, sizeof (Stack));
  s->size = size;
  return s;
}

void stack_push (Stack * s, void * p)
{
  if (s->push)
    p = s->push (s, p);
  s->n++;
  s->p = realloc (s->p, s->n*s->size);
  char * dest = ((char *)s->p) + (s->n - 1)*s->size;
  memcpy (dest, p, s->size);
}

void * stack_pop (Stack * s)
{
  if (!s->n)
    return NULL;
  return ((char *)s->p) + --s->n*s->size;
}

/**
 * @brief Returns a pointer to the element at the specified index from the top of the stack.
 *
 * The index is zero-based, where 0 refers to the top element. Returns NULL if the stack is empty or if the index is out of bounds (negative or greater than the highest valid index).
 *
 * @param s Pointer to the stack.
 * @param i Index from the top of the stack (0-based).
 * @return void* Pointer to the element at the specified index, or NULL if out of bounds.
 */
void * stack_index (Stack * s, int i)
{
  if (!s->n || i > s->n - 1 || i < 0)
    return NULL;
  return ((char *)s->p) + (s->n - i - 1)*s->size;
}

void stack_destroy (Stack * s)
{
  free (s->p);
  free (s);
}

void * stack_set_push (Stack * s, void * push (Stack *, void *))
{
  void * old = s->push;
  s->push = push;
  return old;
}

void * stack_set_data (Stack * s, void * data)
{
  void * old = s->data;
  s->data = data;
  return old;
}

void * stack_get_data (const Stack * s)
{
  return s->data;
}
