# 2 "ast/interpreter/overload.h"

/**
# Function overloading for the interpreter
*/

size_t fread(void *ptr, size_t size, size_t nmemb, FILE *stream){
  size_t undef; return undef;
}
size_t fwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream){
  size_t undef; return undef;
}
int fflush(FILE *stream) {
  int undef; return undef;
}
int fclose(FILE *stream) {
  int undef; return undef;
}
FILE *fopen(const char *pathname, const char *mode) {
  FILE *undef; return undef;
}
void qassert (const char * file, int line, const char * cond){}
void trash (void * alist) {}
void NOT_UNUSED(){}
static void tracing (const char * func, const char * file, int line) {}
static void end_tracing (const char * func, const char * file, int line) {}
void pmuntrace (void) {}
void mpi_init() {}
static void trace_off() {}
void cadna_init (int i) {}
void cadna_end (void) {}

void * pmalloc (long size)
{
  return malloc (size);
}

void * pcalloc (long nmemb, long size)
{
  return calloc (nmemb, size);
}

void * prealloc (void * ptr, long size)
{
  return realloc (ptr, size);
}

void pfree (void * ptr)
{
  free (ptr);
}

void sysfree (void * ptr)
{
  free (ptr);
}

char * pstrdup (const char *s)
{
  return strdup (s);
}

char * systrdup (const char *s)
{
  return strdup (s);
}

char * getenv (const char *name)
{
  return NULL;
}

FILE * fopen (const char * pathname, const char * mode)
{
  return NULL;
}

timer timer_start (void)
{
  timer t = {0};
  return t;
}

double timer_elapsed (timer t)
{
  return 0.;
}

void timer_print (timer t, int i, size_t tnc) {}

double dtnext (double dt)
{
  tnext = dt + 1e30;
  return dt + 1e30;
}

static
bool overload_event()
{
  return (iter == 0);
}

/**
 * @brief Declares variables for use in foreach loops over grid points.
 *
 * Defines a `Point` structure and integer indices for grid traversal.
 */

void _Variables() {
  Point point;
  int ig, jg, kg;
}

enum {
  unset = 1 << 0
};

void unset_double (double * x)
{
  char * flags = x;
  flags += sizeof(double) - 1;
  *flags = unset;
}

/**
 * @brief Initializes point index variables to zero.
 *
 * Sets the global point indices `ig`, `jg`, and `kg` to zero for use in grid-based computations.
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

void free_grid (void)
{
  if (grid) {
    free (((Intergrid *)grid)->d);
    free (grid);
    grid = NULL;
  }
}

/**
 * @brief Resets all scalar fields in the provided list to zero within the grid.
 *
 * Iterates over each scalar in the list and resets its field values for all associated blocks using the grid's data buffer.
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
See the definition of foreach_stencil() in [/src/grid/stencils.h](). */

static void _stencil()
{
  if (baseblock) for (scalar s = baseblock[0], * i = baseblock; s.i >= 0; i++, s = *i) {
    _attribute[s.i].input = _attribute[s.i].output = false;
    _attribute[s.i].width = 0;
  }
}

/**
This is a simplified version of end_stencil() in [/src/grid/stencils.h](). */

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

void _stencil_is_face_x(){}
void _stencil_is_face_y(){}
void _stencil_is_face_z(){}

void _stencil_val (scalar s, int i, int j, int k)
{
  stencil_val ((Point){0}, s, i, j, k, NULL, 0, true);
}

void _stencil_val_o (scalar s, int i, int j, int k)
{
  stencil_val ((Point){0}, s, i, j, k, NULL, 0, true);
}

void _stencil_val_a (scalar s, int i, int j, int k)
{
  stencil_val_a ((Point){0}, s, i, j, k, false, NULL, 0);
}
  
void _stencil_val_r (scalar s, int i, int j, int k)
{
  stencil_val_a ((Point){0}, s, i, j, k, true, NULL, 0);
}
  
int nl = 1;

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
 * @brief Resets the data and status of a scalar field in the grid.
 *
 * Sets the scalar's data to zero and marks it as unset within the grid's data buffer.
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
 * @brief Sets the data for a scalar field to zero in the grid.
 *
 * This function accesses the data buffer for the given scalar field and sets its value to zero. The `leaves` parameter is currently unused.
 *
 * @param s The scalar field whose data is to be reset.
 * @param leaves Unused parameter.
 * @return Always returns zero.
 */
double z_indexing (scalar s, bool leaves)
{
  Intergrid * p = (Intergrid *) grid;
  real * c = p->d + s.i*sizeof (real);
  *c = 0[0];
}

/**
 * @brief Initializes metadata for a scalar block and appends it to relevant lists.
 *
 * Constructs the scalar's name from the provided base and extension, sets its block index, and adds it to the baseblock or all list as appropriate.
 *
 * @param sb Scalar to initialize.
 * @param name Base name for the scalar.
 * @param ext Extension to append to the base name.
 * @param n Index of the scalar within the block; if zero, marks as baseblock.
 * @param block Block index to assign if n is zero.
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
 * Reallocates the data buffer for the current grid to accommodate additional data and updates the total data size.
 *
 * @param size Number of bytes to add to the existing data buffer.
 */
void realloc_scalar (int size)
{
  Intergrid * p = (Intergrid *) grid;
  p->d = realloc (p->d, (datasize + size)*sizeof (char));
  datasize += size;
}

/**
We overload new_const_scalar() so that a new constant field is
allocated for each call. This is necessary when the same constant
field is used for values with different dimensions. */

scalar new_const_scalar (const char * name, int i, double val)
{
  scalar s = (scalar){nconst + 65536}; // _NVARMAX = 65536 as defined in /src/common.h
  init_const_scalar (s, name, val);
  return s;
}

void is_face_x() {}
/**
 * @brief Placeholder function for Y-face operations.
 *
 * This function is a stub and does not perform any action. It may be used as a placeholder for operations related to Y-oriented faces in grid or field management.
 */
void is_face_y() {}
/**
 * @brief Placeholder function for Z-face operations.
 *
 * This function is a stub and does not perform any actions. It may be used as a placeholder for future Z-face related logic in grid or field operations.
 */
void is_face_z() {}

/**
 * @brief Returns a pointer to the data value of a scalar field at the specified grid indices.
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
 * This function creates and returns a Point structure that is not initialized to any specific coordinates or values.
 *
 * @return Point An unset Point.
 */
Point locate (double xp, double yp, double zp)
{
  Point point; // unset
  return point;
}

bool tree_is_full() {
  return false;
}

bool is_leaf() {
  bool ret; // unset
  return ret;
}

void increment_neighbors (Point p) {}
void decrement_neighbors (Point p) {}

int refine_cell (Point point, scalar * list, int flag, Cache * refined) {
  int ret; // unset
  return ret;
}

void mpi_boundary_refine  (scalar * list) {}
void mpi_boundary_update  (scalar * list) {}

typedef struct {
  int nc, nf;
} astats;

/**
 * @brief Copies scalar field values to a max array for wavelet adaptation.
 *
 * If both the scalar list and max array are provided, assigns each scalar's value from the grid to the corresponding element in the max array. Returns unset adaptation statistics.
 *
 * @param slist List of scalar fields to process, terminated by a scalar with negative index.
 * @param max Array to receive the scalar values.
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
These functions are too expensive to interpret. */

void output_field (scalar * list, FILE * fp, int n, bool linear,
		   double box[2][2])
{}

/**
 * @brief Stub function for outputting a scalar field as a PPM image.
 *
 * This function is a placeholder and does not perform any output or processing.
 */
void output_ppm (scalar f, FILE * fp, int n, char * file,
		 double min, double max, double spread,
		 double z, bool linear, double box[2][2],
		 scalar mask, colormap map, char * opt)
{}

/**
 * @brief Placeholder for recursively coarsening a cell and its children.
 *
 * This function is a stub and does not perform any operations.
 */
void coarsen_cell_recursive (Point point, scalar * list)
{}

/**
 * @brief Initializes a grid dimension value at the specified index.
 *
 * Copies the given dimension value into the grid's data buffer at the provided index.
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
Emulations of macros in [/src/common.h](). */

bool is_constant (scalar s)
{
  return s.i >= _NVARMAX;
}

/**
 * @brief Returns the constant value associated with a scalar, or a large default if not constant.
 *
 * If the scalar is marked as constant, retrieves its value from the constant array; otherwise, returns 1e30.
 *
 * @param s Scalar to query.
 * @return double The constant value for the scalar, or 1e30 if the scalar is not constant.
 */
double constant (scalar s)
{
  return is_constant(s) ? _constant[s.i - _NVARMAX] : 1e30;
}

/**
 * @brief Returns an undefined or uninitialized depth value.
 *
 * This function serves as a stub and does not provide a meaningful depth value.
 * @return int Undefined depth.
 */
int depth() { int undef; return undef; }
/**
 * @brief Returns an undefined process identifier.
 *
 * This function is a stub and does not provide a valid process ID.
 *
 * @return int Undefined value representing the process ID.
 */
int pid()   { int undef; return undef; }
/**
 * @brief Returns an undefined thread identifier.
 *
 * This function is a stub and does not provide a valid thread ID.
 * @return int Undefined value.
 */
int tid()   { int undef; return undef; }
/**
 * @brief Returns an undefined value for the number of processing elements.
 *
 * This function serves as a stub and does not provide a valid count of processing elements.
 *
 * @return int Undefined value.
 */
int npe()   { int undef; return undef; }

/**
 * @brief Placeholder for dimension-related operations.
 *
 * This function currently performs no action and serves as a stub for future dimension handling logic.
 */
void dimensional (int a) {}
/**
 * @brief Placeholder function for displaying or processing a dimension value.
 *
 * This function currently has no implementation and serves as a stub for future dimension-related operations.
 */
void show_dimension_internal (double a) {}

/**
## Emulations of macros and functions in <math.h> */

const double M_PI = 3.14159265358979 [0];
const int RAND_MAX = 1;

/**
 * @brief Returns the greater of two double-precision values.
 *
 * @param a First value to compare.
 * @param b Second value to compare.
 * @return The maximum of a and b.
 */
double fmax (double a, double b) { return a > b ? a : b; }
/**
 * @brief Returns the smaller of two double-precision values.
 *
 * @param a First value.
 * @param b Second value.
 * @return The minimum of a and b.
 */
double fmin (double a, double b) { return a < b ? a : b; }
/**
 * @brief Returns the absolute value of an integer.
 *
 * @param i The integer whose absolute value is to be computed.
 * @return int The non-negative value of the input integer.
 */
int abs (int i) { return i < 0 ? - i : i; }
/**
 * @brief Returns an undefined pseudo-random integer value.
 *
 * This stub implementation does not generate a true random number and always returns an undefined value.
 * @return int Undefined integer value.
 */
int rand() { int undef; return undef; }

/**
## Events 

The maximum number of calls to events() is set by `maxevents`. */

int maxevents = 1;

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
 * @brief Clears the highest byte of an integer to zero.
 *
 * Sets the most significant byte of the integer pointed to by @p i to zero, which may be used to clear flag bits or metadata stored in that byte.
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
 * @brief Initializes the grid for GPU operations with the specified size.
 *
 * Calls the standard grid initialization routine for use in GPU contexts.
 *
 * @param n The number of grid points or cells to initialize.
 */

void gpu_init_grid (int n) { init_grid (n); }

/**
 * @brief Constructs a 4-component vector from four float values.
 *
 * @param r Red component.
 * @param g Green component.
 * @param b Blue component.
 * @param a Alpha component.
 * @return vec4 A vector containing the specified components.
 */
vec4 Vec4 (float r, float g, float b, float a)
{
  return (vec4){r, g, b, a};
}

/**
 * @brief Registers a function pointer for use in the interpreter environment.
 *
 * This function allows associating a function pointer with a name and optional non-local context.
 * The implementation is a stub and does not perform any registration.
 *
 * @param ptr Function pointer to register.
 * @param name Name to associate with the function pointer.
 * @param nonlocals Optional context for the function pointer.
 */
void register_fpointer (void (* ptr) (void), const char * name, const void * nonlocals) {}

/**
 * @brief Resets all scalar fields in the provided list to a specified value on the GPU.
 *
 * Calls the reset function for the given list of scalars, setting each field to the provided value.
 *
 * @param alist List of scalar fields to reset.
 * @param val Value to assign to each scalar field.
 */
void reset_gpu (void * alist, double val)
{
  reset (alist, val);
}

/**
 * @brief Searches for a scalar field by name.
 *
 * Looks up a scalar field with the specified name in the global scalar lists. If not found, also checks for block scalars whose names match the prefix of the input name followed by a numeric suffix. Returns an invalid scalar if no match is found.
 *
 * @param name Name of the scalar field to search for.
 * @return scalar The matching scalar field, or an invalid scalar if not found.
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
 * @brief Searches for a vector field by name and returns its vector handle.
 *
 * Attempts to locate a vector field whose base name matches the provided string, appending ".x" to search for the x-component. If not found, also checks for names with numeric suffixes. Returns an invalid vector if no match is found.
 *
 * @param name The base name of the vector field to search for.
 * @return vector The found vector handle, or an invalid vector if not found.
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
