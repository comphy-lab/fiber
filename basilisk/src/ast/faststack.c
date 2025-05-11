#include <stdlib.h>
#include <string.h>
#include "stack.h"
#include "khash.h"
#include "ast.h"
#include "symbols.h"

typedef struct {
  int n, ext;
} Index;

KHASH_MAP_INIT_STR(STR, Index *)

struct _Stack {
  void * p;
  int n, size;
  void * (* push) (Stack *, void *);
  void * data;
  khash_t(STR) * h;
};

/**
 * @brief Creates and initializes a new stack for elements of a specified size.
 *
 * Allocates a new Stack structure, sets the element size, and initializes the internal hash map for name-to-index tracking.
 *
 * @param size Size in bytes of each element to be stored in the stack.
 * @return Pointer to the newly created Stack, or NULL on allocation failure.
 */
Stack * stack_new (int size)
{
  Stack * s = calloc (1, sizeof (Stack));
  s->size = size;
  s->h = kh_init (STR);
  return s;
}

/**
 * @brief Reallocates memory for an array in fixed-size blocks.
 *
 * If the pointer is NULL or the number of elements is a multiple of the block size, reallocates the memory buffer to accommodate one additional block. Otherwise, returns the original pointer unchanged.
 *
 * @param ptr Pointer to the current memory block.
 * @param nmemb Current number of elements.
 * @param size Size of each element in bytes.
 * @param block Block size for allocation granularity.
 * @return void* Pointer to the reallocated memory, or the original pointer if no reallocation was performed.
 */
static void * realloc_block (void * ptr, size_t nmemb, size_t size, int block)
{
  if (!ptr || nmemb % block == 0)
    return realloc (ptr, (nmemb/block + 1)*block*size);
  return ptr;
}

/**
 * @brief Pushes an element onto the stack, tracking identifier or typedef names.
 *
 * If the element represents an AST node for an identifier or typedef name, extracts its name and updates the stack's internal hash map to track its position and linkage type. The element and its associated name pointer are stored together in the stack buffer. If a custom push function is set, it is invoked before processing.
 */
void stack_push (Stack * s, void * p)
{
  if (s->push)
    p = s->push (s, p);
  char * name = NULL;
  Ast * n = *((Ast **)p);
  if (n->sym == sym_IDENTIFIER || n->sym == sym_TYPEDEF_NAME) {
    AstTerminal * t = ast_terminal (n);
    if (!t->after)
      name = strdup (t->start);
    else {
      int len = t->after - t->start + 1;
      name = malloc (len + 1);
      name[len] = '\0';
      strncpy (name, t->start, len);
    }
    Index * index;
    int m = 0, ret;
    khiter_t k = kh_put (STR, s->h, name, &ret);
    assert (ret >= 0);
    if (ret == 0) {
      free (name);
      name = (char *) kh_key (s->h, k);
      index = kh_value (s->h, k);
      for (Index * i = index; i->n >= 0; i++) m++;
    }
    else
      index = NULL;
    index = realloc_block (index, m + 2, sizeof(Index), 10);
    index[m + 1].n = -1;
    index[m].n = s->n;
    index[m].ext = 0;
    if (n->sym == sym_IDENTIFIER) {
      Ast * type = ast_declaration_from_type (n);
      if (ast_schema (type, sym_declaration,
		      0, sym_declaration_specifiers,
		      0, sym_storage_class_specifier,
		      0, sym_EXTERN))
	index[m].ext = 1; // extern
      else if (ast_schema (type->parent, sym_external_declaration,
			   0, sym_declaration))
	index[m].ext = 2; // global
    }
    kh_value (s->h, k) = index;
  }
  s->n++;  
  s->p = realloc_block (s->p, s->n, s->size + sizeof(char *), 100);
  char * dest = ((char *)s->p) + (s->n - 1)*(s->size + sizeof(char *));
  memcpy (dest, p, s->size);
  memcpy (dest + s->size, &name, sizeof (char *));
}

/**
 * @brief Removes and returns the top element from the stack.
 *
 * If the popped element represents an identifier or typedef name, updates or removes its associated entry in the internal hash map and frees related memory as needed.
 *
 * @return Pointer to the popped element, or NULL if the stack is empty.
 */
void * stack_pop (Stack * s)
{
  if (!s->n)
    return NULL;
  void * p = ((char *)s->p) + --s->n*(s->size + sizeof(char *));
  Ast * identifier = *((Ast **)p);
  if (identifier->sym == sym_IDENTIFIER || identifier->sym == sym_TYPEDEF_NAME) {
    char * name = *((char **)(((char *)p) + s->size));
    khiter_t k = kh_get (STR, s->h, name);
    assert (k != kh_end (s->h));
    Index * index = kh_value (s->h, k);
    int n = 0;
    for (Index * i = index; i->n >= 0; i++) {
      assert (i->n <= s->n);
      n++;
    }
    assert (index[n - 1].n == s->n);
    if (n == 1) {
      free (index);
      free (name);
      kh_del (STR, s->h, k);
    }
    else
      index[n - 1].n = -1;
  }
  return p;
}

/**
 * @brief Returns a pointer to the stack element at the given index from the top.
 *
 * Retrieves the element at position `i` counting from the top of the stack (0 is the topmost element).
 * Returns NULL if the index is out of bounds.
 *
 * @param i Index from the top of the stack (0 = top element).
 * @return Pointer to the element at the specified position, or NULL if out of bounds.
 */
void * stack_index (Stack * s, int i)
{
  if (!s->n || i > s->n - 1 || i < 0)
    return NULL;
  return ((char *)s->p) + (s->n - i - 1)*(s->size + sizeof(char *));
}

/**
 * @brief Returns a pointer to the element at the specified index from the bottom of the stack.
 *
 * @param i Zero-based index from the bottom (0 = bottom element).
 * @return Pointer to the element at index i, or NULL if out of bounds.
 */
void * stack_indexi (Stack * s, int i)
{
  if (!s->n || i > s->n - 1 || i < 0)
    return NULL;
  return ((char *)s->p) + i*(s->size + sizeof(char *));
}

/**
 * @brief Finds the most relevant stack element associated with a given name.
 *
 * Looks up the stack element for the specified key using the internal hash map. If the latest entry is an extern declaration, searches backwards for a global declaration and returns that instead. Returns NULL if the key is not found.
 *
 * @param key The name to search for in the stack.
 * @return Pointer to the stack element associated with the key, or NULL if not found.
 */
void * fast_stack_find (Stack * s, const char * key)
{
  khiter_t k = kh_get (STR, s->h, key);
  if (k == kh_end (s->h))
    return NULL;
  Index * index = kh_value (s->h, k);
  int n = 0;
  for (Index * i = index; i->n >= 0; i++) n++;
  if (index[n - 1].ext == 1) // extern declaration
    for (int i = n - 2; i >= 0; i--)
      if (index[i].ext == 2) // return first global declaration
	return ((char *)s->p) + index[i].n*(s->size + sizeof(char *));
  return ((char *)s->p) + index[n - 1].n*(s->size + sizeof(char *));
}

/**
 * @brief Frees all resources associated with the stack.
 *
 * Releases memory for all stored name strings, index arrays, the internal hash map, the stack buffer, and the stack structure itself.
 */
void stack_destroy (Stack * s)
{
  const char * name;
  Index * index;
  kh_foreach (s->h, name, index,
	      (free ((void *) name), free (index)));
  kh_destroy (STR, s->h);
  free (s->p);
  free (s);
}

/**
 * @brief Sets a custom push function for the stack.
 *
 * Replaces the current push function with the provided one and returns the previous push function pointer.
 *
 * @param push Pointer to the new push function to use for stack operations.
 * @return Pointer to the previous push function.
 */
void * stack_set_push (Stack * s, void * push (Stack *, void *))
{
  void * old = s->push;
  s->push = push;
  return old;
}

/**
 * @brief Sets the auxiliary data pointer for the stack.
 *
 * Replaces the stack's current auxiliary data pointer with the provided value and returns the previous pointer.
 *
 * @param data The new auxiliary data pointer to associate with the stack.
 * @return The previous auxiliary data pointer.
 */
void * stack_set_data (Stack * s, void * data)
{
  void * old = s->data;
  s->data = data;
  return old;
}

/**
 * @brief Retrieves the auxiliary data pointer associated with the stack.
 *
 * @return The opaque data pointer stored in the stack, or NULL if none is set.
 */
void * stack_get_data (const Stack * s)
{
  return s->data;
}
