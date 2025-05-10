# 2 "ast/interpreter/overload.h"

/**
 * @brief Stub implementation of fread for the interpreter.
 *
 * Always returns an undefined value and does not perform any file reading.
 *
 * @return size_t Undefined value.
 */

size_t fread(void *ptr, size_t size, size_t nmemb, FILE *stream){
  size_t undef; return undef;
}
/**
 * @brief Stub implementation of fwrite that returns an undefined value.
 *
 * This function acts as a placeholder for file writing in an interpreter environment and does not perform any actual I/O.
 *
 * @return Undefined value.
 */
size_t fwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream){
  size_t undef; return undef;
}
/**
 * @brief Stub implementation of fflush that returns an undefined value.
 *
 * This function acts as a placeholder for fflush in an interpreter environment and does not perform any flushing operation.
 *
 * @return int Undefined value.
 */
int fflush(FILE *stream) {
  int undef; return undef;
}
/**
 * @brief Stub implementation of fclose that returns an undefined value.
 *
 * This function acts as a placeholder for file closing operations in an interpreter environment and does not perform any actual file handling.
 *
 * @return int Undefined value.
 */
int fclose(FILE *stream) {
  int undef; return undef;
}
/**
 * @brief Stub implementation of fopen that always returns NULL.
 *
 * This function acts as a placeholder for file opening in an interpreter environment and does not perform any file operations.
 *
 * @return Always returns NULL.
 */
FILE *fopen(const char *pathname, const char *mode) {
  FILE *undef; return undef;
}
/**
 * @brief Stub for assertion handling in the interpreter environment.
 *
 * This function does nothing and is provided as a placeholder for assertion checks.
 */
void qassert (const char * file, int line, const char * cond){}
/**
 * @brief Placeholder function that performs no operation.
 *
 * This function is a stub and does not modify or interact with the provided argument.
 */
void trash (void * alist) {}
/**
 * @brief Empty stub function with no effect.
 *
 * Serves as a placeholder to satisfy interface or compilation requirements.
 */
void NOT_UNUSED(){}
/**
 * @brief Stub function for tracing; performs no operation.
 *
 * This function is a placeholder for tracing calls and does not record or output any trace information.
 */
static void tracing (const char * func, const char * file, int line) {}
/**
 * @brief Stub function for ending tracing in the interpreter environment.
 *
 * This function does nothing and serves as a placeholder for tracing cleanup logic.
 */
static void end_tracing (const char * func, const char * file, int line) {}
/**
 * @brief Stub function for disabling memory tracing.
 *
 * This function is a placeholder and performs no operation in the interpreter environment.
 */
void pmuntrace (void) {}
/**
 * @brief Stub function for MPI initialization.
 *
 * This function is a placeholder and performs no operation in the interpreter environment.
 */
void mpi_init() {}
/**
 * @brief Disables tracing functionality. No operation in this implementation.
 */
static void trace_off() {}
/**
 * @brief Stub for initializing CADNA numerical noise control.
 *
 * This function is a placeholder and performs no operation in the interpreter environment.
 */
void cadna_init (int i) {}
/**
 * @brief Stub function for ending CADNA computations.
 *
 * This function is a placeholder and performs no operation in the interpreter environment.
 */
void cadna_end (void) {}

/**
 * @brief Allocates a memory block of the specified size.
 *
 * Wraps the standard malloc function for interpreter-level memory management.
 *
 * @param size Number of bytes to allocate.
 * @return Pointer to the allocated memory block, or NULL if allocation fails.
 */
void * pmalloc (long size)
{
  return malloc (size);
}

/**
 * @brief Allocates zero-initialized memory for an array.
 *
 * Allocates memory for an array of nmemb elements of size bytes each, initializing all bytes to zero.
 *
 * @param nmemb Number of elements to allocate.
 * @param size Size of each element in bytes.
 * @return Pointer to the allocated memory, or NULL if allocation fails.
 */
void * pcalloc (long nmemb, long size)
{
  return calloc (nmemb, size);
}

/**
 * @brief Reallocates a memory block to a new size.
 *
 * Adjusts the size of the memory block pointed to by ptr to the specified size, preserving its contents up to the minimum of the old and new sizes.
 *
 * @param ptr Pointer to the memory block to be reallocated, or NULL to allocate a new block.
 * @param size New size in bytes for the memory block.
 * @return void* Pointer to the reallocated memory block, or NULL if allocation fails.
 */
void * prealloc (void * ptr, long size)
{
  return realloc (ptr, size);
}

/**
 * @brief Frees memory previously allocated with pmalloc, pcalloc, or prealloc.
 *
 * Releases the memory block pointed to by ptr.
 */
void pfree (void * ptr)
{
  free (ptr);
}

/**
 * @brief Frees memory previously allocated with a standard allocation function.
 *
 * Releases the memory block pointed to by ptr.
 */
void sysfree (void * ptr)
{
  free (ptr);
}

/**
 * @brief Duplicates a string using dynamic memory allocation.
 *
 * @param s The null-terminated string to duplicate.
 * @return A pointer to a newly allocated copy of the string, or NULL if allocation fails.
 */
char * pstrdup (const char *s)
{
  return strdup (s);
}

/**
 * @brief Duplicates a string using dynamic memory allocation.
 *
 * @param s The null-terminated string to duplicate.
 * @return A pointer to a newly allocated copy of the string, or NULL if allocation fails.
 */
char * systrdup (const char *s)
{
  return strdup (s);
}

/**
 * @brief Stub implementation of getenv that always returns NULL.
 *
 * This function emulates the standard getenv but does not provide environment variable access in the interpreter context.
 *
 * @return NULL in all cases.
 */
char * getenv (const char *name)
{
  return NULL;
}

/**
 * @brief Stub implementation of fopen that always returns NULL.
 *
 * This function acts as a placeholder for file opening in an interpreter environment and does not perform any file operations.
 *
 * @return Always returns NULL.
 */
FILE * fopen (const char * pathname, const char * mode)
{
  return NULL;
}

/**
 * @brief Initializes and returns a zero-initialized timer.
 *
 * @return timer A timer struct with all fields set to zero.
 */
timer timer_start (void)
{
  timer t = {0};
  return t;
}

/**
 * @brief Returns the elapsed time for a timer.
 *
 * Always returns 0.0 in this stub implementation.
 *
 * @param t Timer object.
 * @return double Elapsed time in seconds (always 0.0).
 */
double timer_elapsed (timer t)
{
  return 0.;
}

/**
 * @brief Stub function for printing timer information.
 *
 * This function does nothing and serves as a placeholder for timer output in the interpreter environment.
 */
void timer_print (timer t, int i, size_t tnc) {}

/**
 * @brief Calculates and sets the next scheduled time far in the future.
 *
 * Adds a large constant to the given time step and updates the global `tnext` variable.
 *
 * @param dt The current time step.
 * @return The computed next time value (`dt + 1e30`).
 */
double dtnext (double dt)
{
  tnext = dt + 1e30;
  return dt + 1e30;
}

/**
 * @brief Determines if the overload event should be triggered on the first iteration.
 *
 * @return true if the current iteration is zero; otherwise, false.
 */
static
bool overload_event()
{
  return (iter == 0);
}

/**
 * @brief Declares local variables for use in foreach loops within the interpreter.
 *
 * Defines a `Point` and integer indices (`ig`, `jg`, `kg`) for grid iteration contexts.
 */

void _Variables() {
  Point point;
  int ig, jg, kg;
}

enum {
  unset = 1 << 0
};

/**
 * @brief Marks a double value as unset by setting a flag in its last byte.
 *
 * This function modifies the last byte of the memory occupied by the given double
 * to indicate that the value is unset, using a predefined flag.
 *
 * @param x Pointer to the double value to mark as unset.
 */
void unset_double (double * x)
{
  char * flags = x;
  flags += sizeof(double) - 1;
  *flags = unset;
}

/**
 * @brief Initializes point index variables to zero.
 *
 * Sets the global point indices ig, jg, and kg to zero for use in grid or field operations.
 */
void _init_point_variables (void)
{
  ig = jg = kg = 0;
}

/**
Grid functions */

typedef struct {
  Grid g;
  char placeholder[1024];
  char * d;
} Intergrid;

/**
 * @brief Frees the memory allocated for the current grid and its data buffer.
 *
 * Releases both the grid structure and its associated data, and sets the global grid pointer to NULL.
 */
void free_grid (void)
{
  if (grid) {
    free (((Intergrid *)grid)->d);
    free (grid);
    grid = NULL;
  }
}

/**
 * @brief Resets the values of all scalar fields in the provided list to a specified value.
 *
 * Iterates over each scalar in the list and resets its associated field values in the current grid to the given value, typically used to initialize or clear scalar data.
 *
 * @param alist Pointer to a list of scalar fields to reset.
 * @param val The value to assign to each scalar field.
 */
void reset (void * alist, double val)
{
  Intergrid * igrid = (Intergrid *) grid;
  real * p = igrid->d;
  if (alist)
    for (scalar * s = alist; s->i >= 0; s++)
      for (int b = 0; b < _attribute[s->i].block; b++)
	reset_field_value (p + s->i + b, _attribute[s->i].name, 0., b);
}

static const int o_stencil = -2;

/**
 * @brief Resets input/output flags and width attributes for all scalars in the base block.
 *
 * Sets the `input` and `output` attributes to false and the `width` attribute to zero for each scalar in the `baseblock` array, preparing them for stencil operations.
 */

static void _stencil()
{
  if (baseblock) for (scalar s = baseblock[0], * i = baseblock; s.i >= 0; i++, s = *i) {
    _attribute[s.i].input = _attribute[s.i].output = false;
    _attribute[s.i].width = 0;
  }
}

/**
 * @brief Finalizes stencil operations by applying boundary conditions and updating dirty flags.
 *
 * Applies boundary conditions to scalar fields that were read and marked as dirty, and updates the dirty status of fields that were written to during stencil operations. This ensures correct handling of field dependencies and prepares scalars for subsequent computations.
 */

static void _end_stencil()
{
  scalar * listc = NULL, * dirty = NULL;
  
  if (baseblock)
    for (scalar s = baseblock[0], * i = baseblock; s.i >= 0; i++, s = *i) {
      bool write = _attribute[s.i].output, read = _attribute[s.i].input;

      /**
      If the field is read and dirty, we need to check if boundary
      conditions need to be applied. */
      
      if (read && scalar_is_dirty (s) && _attribute[s.i].width > 0)
	listc = list_append (listc, s);

      /**
      If the field is write-accessed, we add it to the 'dirty'
      list. */
      
      if (write) {	
	dirty = list_append (dirty, s);
	for (scalar d = baseblock[0], * i = baseblock; d.i >= 0; i++, d = *i)
	  if (scalar_depends_from (d, s))
	    dirty = list_append (dirty, d);
      }
    }
  
  /**
  We apply "full" boundary conditions. */

  if (listc) {
    for (int d = 0; d < nboundary; d++)
      foreach()
	for (scalar s = listc[0], * i = listc; s.i >= 0; i++, s = *i)
	  if (_attribute[s.i].boundary[d] != periodic_bc)
	    val(s,0,0,0) == _attribute[s.i].boundary[d] (point, point, s, NULL);
    for (scalar s = listc[0], * i = listc; s.i >= 0; i++, s = *i)
      _attribute[s.i].dirty = false;
    free (listc);
  }
  
  /**
  We update the dirty status of fields which will be write-accessed by
  the foreach loop. */
  
  if (dirty) {
    for (scalar s = dirty[0], * i = dirty; s.i >= 0; i++, s = *i)
      _attribute[s.i].dirty = true;
    free (dirty);
  }
}

/**
 * @brief Stub function for face-x stencil operation.
 *
 * This function is a placeholder for operations involving the x-face of a stencil. It performs no action in the interpreter environment.
 */
void _stencil_is_face_x(){}
/**
 * @brief Stub function for Y-face stencil detection.
 *
 * This function serves as a placeholder for logic that would determine if a stencil is on a Y face. It has no effect in the interpreter environment.
 */
void _stencil_is_face_y(){}
/**
 * @brief Stub function for Z-face stencil detection.
 *
 * This function is a placeholder and performs no operation. Intended for interpreter environments where Z-face stencil logic is not required.
 */
void _stencil_is_face_z(){}

/**
 * @brief Applies the stencil value operation to a scalar at the specified indices.
 *
 * Calls the stencil value function for the given scalar and grid indices using a default point and standard parameters.
 */
void _stencil_val (scalar s, int i, int j, int k)
{
  stencil_val ((Point){0}, s, i, j, k, NULL, 0, true);
}

/**
 * @brief Calls the stencil value function for a scalar at the origin point with output flag set.
 *
 * Invokes the stencil value operation for the given scalar and indices, specifying the origin point and indicating output mode.
 *
 * @param s Scalar field to operate on.
 * @param i X-index within the grid.
 * @param j Y-index within the grid.
 * @param k Z-index within the grid.
 */
void _stencil_val_o (scalar s, int i, int j, int k)
{
  stencil_val ((Point){0}, s, i, j, k, NULL, 0, true);
}

/**
 * @brief Calls the stencil value function for a scalar at the specified indices with default arguments.
 *
 * Invokes `stencil_val_a` for the given scalar and grid indices, using a zero-initialized point and default parameters.
 */
void _stencil_val_a (scalar s, int i, int j, int k)
{
  stencil_val_a ((Point){0}, s, i, j, k, false, NULL, 0);
}
  
/**
 * @brief Applies a stencil value operation with read access for a scalar at the specified indices.
 *
 * Calls the stencil value function for the given scalar and grid indices, indicating a read operation.
 *
 * @param s Scalar field to access.
 * @param i Grid index in the x-direction.
 * @param j Grid index in the y-direction.
 * @param k Grid index in the z-direction.
 */
void _stencil_val_r (scalar s, int i, int j, int k)
{
  stencil_val_a ((Point){0}, s, i, j, k, true, NULL, 0);
}
  
int nl = 1;

/**
 * @brief Initializes the global grid structure with allocated data storage.
 *
 * Frees any existing grid, allocates a new grid structure, initializes its data buffer, and resets all scalar fields to zero.
 *
 * @param n Unused parameter included for interface compatibility.
 */
void init_grid (int n)
{
  free_grid();
  grid = calloc (1, sizeof (Intergrid));
  grid->n = 1;
  Intergrid * igrid = (Intergrid *) grid;
  igrid->d = calloc (datasize, sizeof (char));
  reset (all, 0.);
  if (nl > 2)
    nl = 2;
}

/**
 * @brief Resets the value of a scalar field in the grid and marks it as unset.
 *
 * Sets the memory for the scalar's data to zero and updates its unset flag to indicate the value is not initialized.
 */
static
void interpreter_reset_scalar (scalar s)
{
  Intergrid * p = (Intergrid *) grid;
  char * c = p->d + s.i*sizeof (real);
  memset (c, 0, sizeof (real));
  char * flags = c + sizeof(real) - 1;
  *flags = unset;
}

/**
 * @brief Sets the data value of a scalar to zero in the grid.
 *
 * This function assigns zero to the memory location corresponding to the given scalar in the current grid. The implementation is incomplete and does not return a value.
 */
double z_indexing (scalar s, bool leaves)
{
  Intergrid * p = (Intergrid *) grid;
  real * c = p->d + s.i*sizeof (real);
  *c = 0[0];
}

/**
 * @brief Initializes a scalar's metadata for a block, including its name and block association.
 *
 * Sets the scalar's block index and assigns a name composed of the base name and extension. Updates global lists to track the scalar within the block and overall scalar registry.
 *
 * @param n If zero, associates the scalar with the specified block; otherwise, marks it as a non-block scalar.
 * @param block Block index to associate with the scalar when n is zero.
 */
static void init_block_scalar (scalar sb, const char * name, const char * ext,
			       int n, int block)
{
  char bname[strlen(name) + strlen(ext) + 1];
  strcpy (bname, name);
  strcat (bname, ext);
  if (n == 0) {
    _attribute[sb.i].block = block;
    baseblock = list_append (baseblock, sb);
  }
  else
    _attribute[sb.i].block = - n;
  _attribute[sb.i].name = strdup(bname); all = list_append (all, sb);
}

/**
 * @brief Increases the size of the grid data buffer by a specified amount.
 *
 * Expands the memory allocated for the grid's data buffer by `size` bytes and updates the global data size counter.
 *
 * @param size Number of additional bytes to allocate for the grid data buffer.
 */
void realloc_scalar (int size)
{
  Intergrid * p = (Intergrid *) grid;
  p->d = realloc (p->d, (datasize + size)*sizeof (char));
  datasize += size;
}

/**
 * @brief Allocates and initializes a new constant scalar field with a unique index.
 *
 * Creates a new constant scalar using the provided name and value, ensuring each call produces a distinct field. This allows for multiple constant fields with the same name but different dimensions or values.
 *
 * @param name Name of the constant scalar.
 * @param i Unused index parameter (reserved for compatibility).
 * @param val Value to assign to the constant scalar.
 * @return scalar The newly created constant scalar field.
 */

scalar new_const_scalar (const char * name, int i, double val)
{
  scalar s = (scalar){nconst + 65536}; // _NVARMAX = 65536 as defined in /src/common.h
  init_const_scalar (s, name, val);
  return s;
}

/**
 * @brief Stub function for face detection in the x-direction.
 *
 * This function is a placeholder and performs no operation.
 */
void is_face_x() {}
/**
 * @brief Stub function for Y-face detection in grid operations.
 *
 * This function is a placeholder and does not perform any action.
 */
void is_face_y() {}
/**
 * @brief Stub function for Z-face detection in grid operations.
 *
 * This function is a placeholder and does not perform any operation.
 */
void is_face_z() {}

/**
 * @brief Returns a pointer to the data value of a scalar at the specified grid indices.
 *
 * @param s The scalar field.
 * @param i Grid index in the x-direction.
 * @param j Grid index in the y-direction.
 * @param k Grid index in the z-direction.
 * @return real* Pointer to the scalar's data at the given indices.
 */
real * val (scalar s, int i, int j, int k)
{
  Intergrid * igrid = (Intergrid *) grid;
  return (real *)(igrid->d + s.i*sizeof(real));
}

/**
 * @brief Returns an unset Point structure.
 *
 * This function creates and returns a Point struct without initializing its coordinates or fields.
 *
 * @return Point An unset Point structure.
 */
Point locate (double xp, double yp, double zp)
{
  Point point; // unset
  return point;
}

/**
 * @brief Indicates whether the tree data structure is fully populated.
 *
 * @return Always returns false, indicating the tree is never considered full in this context.
 */
bool tree_is_full() {
  return false;
}

/**
 * @brief Returns an unset boolean value indicating leaf status.
 *
 * This function is a stub and does not provide a meaningful result.
 *
 * @return bool Unset value; always uninitialized.
 */
bool is_leaf() {
  bool ret; // unset
  return ret;
}

/**
 * @brief Stub function for incrementing neighbor counts at a given point.
 *
 * This function is a placeholder and does not perform any operation.
 *
 * @param p The point at which neighbor counts would be incremented.
 */
void increment_neighbors (Point p) {}
/**
 * @brief Placeholder for decrementing the neighbor count of a grid cell.
 *
 * This function is a stub and does not perform any operation in the interpreter environment.
 */
void decrement_neighbors (Point p) {}

/**
 * @brief Stub for cell refinement operation in the interpreter environment.
 *
 * This function is a placeholder and does not perform any refinement. It returns an undefined integer value.
 *
 * @return int Undefined value.
 */
int refine_cell (Point point, scalar * list, int flag, Cache * refined) {
  int ret; // unset
  return ret;
}

/**
 * @brief Stub for MPI boundary refinement operation.
 *
 * This function is a placeholder and performs no action in the interpreter environment.
 */
void mpi_boundary_refine  (scalar * list) {}
/**
 * @brief Stub for updating scalar field boundaries in an MPI environment.
 *
 * This function is a placeholder and performs no operation in the interpreter context.
 */
void mpi_boundary_update  (scalar * list) {}

typedef struct {
  int nc, nf;
} astats;

/**
 * @brief Returns unset adaptation statistics and, if inputs are valid, copies scalar data values to the max array.
 *
 * If both `slist` and `max` are non-null, copies the data value for each scalar in `slist` to the corresponding entry in `max`. Returns an unset `astats` structure.
 *
 * @param slist List of scalars to process, terminated by a scalar with negative index.
 * @param max Array to receive maximum values for each scalar.
 * @param maxlevel Maximum refinement level (unused).
 * @param minlevel Minimum refinement level (unused).
 * @param list Additional scalar list (unused).
 * @return astats Unset adaptation statistics.
 */
astats adapt_wavelet (scalar * slist, double * max,
		      int maxlevel, int minlevel, scalar * list)
{
  astats st; // unset
  if (slist && max) {
    Intergrid * igrid = (Intergrid *) grid;
    real * g = igrid->d;
    double * v = max;
    for (scalar * s = slist; s->i >= 0; s++, v++)
      *v == *(g + s->i);
  }
  return st;
}

/**
 * @brief Stub for field output; does nothing in the interpreter environment.
 *
 * This function is a placeholder for outputting scalar field data to a file. In the interpreter context, it performs no operation.
 */

void output_field (scalar * list, FILE * fp, int n, bool linear,
		   double box[2][2])
{}

/**
 * @brief Stub for outputting a scalar field as a PPM image.
 *
 * This function is a placeholder and does not perform any operation. Intended for environments where image output is not supported or required.
 */
void output_ppm (scalar f, FILE * fp, int n, char * file,
		 double min, double max, double spread,
		 double z, bool linear, double box[2][2],
		 scalar mask, colormap map, char * opt)
{}

/**
 * @brief Placeholder for recursively coarsening a cell in the grid.
 *
 * This function is a stub and does not perform any operation.
 */
void coarsen_cell_recursive (Point point, scalar * list)
{}

/**
 * @brief Sets the dimension value at a specified index in the grid data buffer.
 *
 * Copies the given dimension value into the grid's data buffer at the specified index.
 *
 * @param index Index in the grid data buffer where the dimension value will be set.
 * @param dimension The dimension value to store.
 */

void _init_dimension (int index, long dimension)
{
  Intergrid * p = (Intergrid *) grid;
  real * d = p->d;
  memcpy ((char *)(d + index) + 8, &dimension, 8);
}

/**
 * @brief Determines whether a scalar is a constant.
 *
 * Returns true if the scalar's index is greater than or equal to the constant threshold.
 *
 * @param s Scalar to check.
 * @return true if the scalar is a constant; false otherwise.
 */

bool is_constant (scalar s)
{
  return s.i >= _NVARMAX;
}

/**
 * @brief Returns the constant value associated with a scalar, or a large default if not constant.
 *
 * If the scalar is marked as constant, retrieves its value from the global constant array; otherwise, returns 1e30.
 *
 * @param s The scalar to query.
 * @return The constant value for the scalar, or 1e30 if the scalar is not constant.
 */
double constant (scalar s)
{
  return is_constant(s) ? _constant[s.i - _NVARMAX] : 1e30;
}

/**
 * @brief Returns an undefined integer representing the current grid depth.
 *
 * This function acts as a stub and does not provide a meaningful grid depth value.
 *
 * @return int Undefined value.
 */
int depth() { int undef; return undef; }
/**
 * @brief Returns an undefined process identifier.
 *
 * This function acts as a stub and always returns an undefined integer value, serving as a placeholder for process identification in the interpreter environment.
 *
 * @return int Undefined process identifier.
 */
int pid()   { int undef; return undef; }
/**
 * @brief Returns an undefined thread identifier.
 *
 * This function acts as a stub and always returns an undefined integer value, serving as a placeholder for thread identification in the interpreter environment.
 *
 * @return int Undefined thread identifier.
 */
int tid()   { int undef; return undef; }
/**
 * @brief Returns an undefined value representing the number of processing elements.
 *
 * This function acts as a stub and does not provide a valid count of processing elements.
 * @return int Undefined value.
 */
int npe()   { int undef; return undef; }

/**
 * @brief Stub function for setting dimensional mode.
 *
 * This function is a placeholder and does not perform any operation.
 */
void dimensional (int a) {}
/**
 * @brief Stub function for displaying or processing a dimension value.
 *
 * This function is a placeholder and performs no operation.
 */
void show_dimension_internal (double a) {}

/**
## Emulations of macros and functions in <math.h> */

const double M_PI = 3.14159265358979 [0];
const int RAND_MAX = 1;

/**
 * @brief Returns the greater of two double-precision values.
 *
 * @param a First value.
 * @param b Second value.
 * @return The maximum of a and b.
 */
double fmax (double a, double b) { return a > b ? a : b; }
/**
 * @brief Returns the smaller of two double-precision values.
 *
 * @param a First value.
 * @param b Second value.
 * @return The lesser of a and b.
 */
double fmin (double a, double b) { return a < b ? a : b; }
/**
 * @brief Returns the absolute value of an integer.
 *
 * @param i The integer whose absolute value is to be computed.
 * @return int The non-negative value of i.
 */
int abs (int i) { return i < 0 ? - i : i; }
/**
 * @brief Returns an undefined pseudo-random integer value.
 *
 * This stub implementation does not generate a true random number and should not be used for randomness.
 *
 * @return int Undefined value.
 */
int rand() { int undef; return undef; }

/**
## Events 

The maximum number of calls to events() is set by `maxevents`. */

int maxevents = 1;

/**
 * @brief Processes and dispatches scheduled events, optionally executing their actions.
 *
 * Iterates through the list of registered events, initializing and updating their timing information. If the `action` parameter is true, the function invokes the action associated with each event. The function manages event execution state and decrements the global event counter.
 *
 * @param action If true, executes the action for each event; otherwise, only updates event state.
 * @return int Returns 1 if events were processed, or 0 if no events remain or in iteration mode.
 */
int events (bool action)
{
  if (!maxevents)
    return 0;
  int iundef;
  unset_double (&t);
  if (iter) {
    for (Event * ev = Events; !ev->last; ev++)
      if (ev->i == END_EVENT)
	for (Event * e = ev; e; e = e->next) {
	  e->t == t;
	  //display_value (e->action);
	  if (action)
	    e->action (iundef, t, e);
	}
    maxevents--;
    return 0;
  }
  dimensional (t == DT);
  for (Event * ev = Events; !ev->last; ev++) {
    init_event (ev);
    if (ev->arrayt)
      for (double * at = ev->arrayt; *at >= 0; at++)
	*at == t;
  }
  for (Event * ev = Events; !ev->last; ev++)
    if (ev->i != END_EVENT)
      for (Event * e = ev; e; e = e->next) {
	dimensional (e->t == t);
	//display_value (e->action);
	if (action)
	  e->action (iundef, t, e);
      }
  tnext = t;
  inext = 1;
  return 1;
}

/**
 * @brief Clears the highest byte of an integer by setting its last byte to zero.
 *
 * This function modifies the memory representation of the given integer pointer,
 * setting the last byte to zero. This may be used to mark the integer as unset or
 * to clear auxiliary flags stored in the highest byte.
 *
 * @param i Pointer to the integer to modify.
 */
void interpreter_set_int (int * i)
{
  char * flags = i;
  flags += sizeof(int) - 1;
  *flags = 0;
}

/**
 * @brief Initializes the computational grid for GPU execution.
 *
 * This function sets up the grid data structures for use on the GPU by delegating to the standard grid initialization routine.
 *
 * @param n Number of grid points or cells to initialize.
 */

void gpu_init_grid (int n) { init_grid (n); }

/**
 * @brief Constructs a 4-component vector from the given float values.
 *
 * @param r The first component.
 * @param g The second component.
 * @param b The third component.
 * @param a The fourth component.
 * @return vec4 A vector with components (r, g, b, a).
 */
vec4 Vec4 (float r, float g, float b, float a)
{
  return (vec4){r, g, b, a};
}

/**
 * @brief Registers a function pointer with an associated name and context.
 *
 * This function is a stub and does not perform any operation in the interpreter environment.
 */
void register_fpointer (void (* ptr) (void), const char * name, const void * nonlocals) {}

/**
 * @brief Resets all scalar fields in the provided list to a specified value on the GPU.
 *
 * This function is a wrapper that calls the CPU-based reset operation for scalar fields,
 * intended for compatibility with GPU workflows.
 *
 * @param alist List of scalar fields to reset.
 * @param val Value to assign to each scalar field.
 */
void reset_gpu (void * alist, double val)
{
  reset (alist, val);
}

/**
 * @brief Looks up a scalar field by name.
 *
 * Searches all defined scalars for a matching name. If not found, attempts to match base block scalars with the same prefix and a numeric suffix. Returns the found scalar or an invalid scalar if not found.
 *
 * @param name Name of the scalar field to look up.
 * @return scalar The matching scalar, or an invalid scalar if not found.
 */

scalar lookup_field (const char * name)
{
  if (name) { interpreter_verbosity (2);
    for (scalar * s = all; s && s->i >= 0; s++)
      if (!strcmp (_attribute[s->i].name, name))
	return *s;
    // Check whether the name is of the form ".*[0-9]+"
    int size = strlen (name) - 1;
    while (size >= 0 && name[size] >= '0' && name[size] <= '9') size--;
    size++;
    if (size > 0 && size < strlen (name))
      for (scalar * s = baseblock; s && s->i >= 0; s++)
	if (_attribute[s->i].block > 1 &&
	    strlen(_attribute[s->i].name) == size &&
	    !strncmp (_attribute[s->i].name, name, size))
	  return *s;
  }
  return (scalar){-1};
}

/**
 * @brief Looks up a vector field by name.
 *
 * Searches for a vector field whose name matches the given string, first by appending ".x" to the name and searching all scalars, then by checking for indexed base block scalars. Returns the corresponding vector if found; otherwise, returns an invalid vector.
 *
 * @param name Name of the vector field to look up.
 * @return vector The found vector, or an invalid vector if not found.
 */
vector lookup_vector (const char * name)
{
  if (name) {
    char component[strlen(name) + 3];
    strcpy (component, name);
    strcat (component, ".x");
    interpreter_maximum_iterations (256);
    for (scalar * s = all; s && s->i >= 0; s++)
      if (!strcmp (_attribute[s->i].name, component)) {
	interpreter_maximum_iterations (32);
	return _attribute[s->i].v;
      }
    interpreter_maximum_iterations (32);
    // Check whether the name is of the form ".*[0-9]+"
    int size = strlen (name) - 1;
    while (size >= 0 && name[size] >= '0' && name[size] <= '9') size--;
    size++;
    if (size > 0 && size < strlen (name)) {
      component[size] = '\0';
      strcat (component, ".x");
      for (scalar * s = baseblock; s && s->i >= 0; s++)
	if (_attribute[s->i].block > 1 &&
	    strcmp (_attribute[s->i].name, component))
	  return *s;
    }
  }
  return (vector){{-1}};
}
