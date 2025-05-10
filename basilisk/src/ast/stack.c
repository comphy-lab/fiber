#include <stdlib.h>
#include <string.h>
#include "stack.h"

struct _Stack {
  void * p;
  int n, size;
  void * (* push) (Stack *, void *);
  void * data;
};

/**
 * @brief Allocates and initializes a new stack for elements of a specified size.
 *
 * @param size Size in bytes of each stack element.
 * @return Pointer to the newly created stack, or NULL if allocation fails.
 */
Stack * stack_new (int size)
{
  Stack * s = calloc (1, sizeof (Stack));
  s->size = size;
  return s;
}

/**
 * @brief Pushes an element onto the stack.
 *
 * If a custom push function is set, it is applied to the element before insertion. The stack's internal buffer is resized as needed to accommodate the new element.
 */
void stack_push (Stack * s, void * p)
{
  if (s->push)
    p = s->push (s, p);
  s->n++;
  s->p = realloc (s->p, s->n*s->size);
  char * dest = ((char *)s->p) + (s->n - 1)*s->size;
  memcpy (dest, p, s->size);
}

/**
 * @brief Removes and returns the top element from the stack.
 *
 * @param s Pointer to the stack.
 * @return Pointer to the popped element's data, or NULL if the stack is empty.
 */
void * stack_pop (Stack * s)
{
  if (!s->n)
    return NULL;
  return ((char *)s->p) + --s->n*s->size;
}

/**
 * @brief Returns a pointer to the element at the specified index from the top of the stack.
 *
 * Index 0 refers to the top element. Returns NULL if the stack is empty or if the index is out of bounds.
 *
 * @param i Index of the element to access, counting from the top (0 is the top element).
 * @return void* Pointer to the requested element, or NULL if not found.
 */
void * stack_index (Stack * s, int i)
{
  if (!s->n || i > s->n - 1 || i < 0)
    return NULL;
  return ((char *)s->p) + (s->n - i - 1)*s->size;
}

/**
 * @brief Frees all memory associated with the stack.
 *
 * Releases the internal buffer and the stack object itself.
 */
void stack_destroy (Stack * s)
{
  free (s->p);
  free (s);
}

/**
 * @brief Sets a custom push function for the stack.
 *
 * Replaces the current push callback with the provided function, allowing custom behavior when elements are pushed onto the stack.
 *
 * @param push Pointer to a function that takes a Stack pointer and an element pointer, returning a possibly transformed element pointer.
 * @return The previous push function pointer.
 */
void * stack_set_push (Stack * s, void * push (Stack *, void *))
{
  void * old = s->push;
  s->push = push;
  return old;
}

/**
 * @brief Sets the user-defined data pointer for the stack.
 *
 * Replaces the stack's auxiliary data pointer with the provided value and returns the previous pointer.
 *
 * @param data New user-defined data pointer to associate with the stack.
 * @return void* The previous data pointer associated with the stack.
 */
void * stack_set_data (Stack * s, void * data)
{
  void * old = s->data;
  s->data = data;
  return old;
}

/**
 * @brief Retrieves the user-defined data pointer associated with the stack.
 *
 * @param s Pointer to the stack.
 * @return The auxiliary data pointer set for the stack, or NULL if not set.
 */
void * stack_get_data (const Stack * s)
{
  return s->data;
}
